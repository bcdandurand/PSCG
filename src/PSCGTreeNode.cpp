/*===========================================================================*
 * This file is part of the Bcps Linear Solver (PSCG).                       *
 *                                                                           *
 * PSCG is distributed under the Eclipse Public License as part of the       *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2017, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include <cassert>
#include <iostream>
#include <utility>
#include <cmath>
#include <vector>

#include "CoinUtility.hpp"
//#include "OsiRowCut.hpp"
//#include "OsiColCut.hpp"
//#include "OsiRowCutDebugger.hpp"
//#include "OsiCuts.hpp"

#include "AlpsKnowledge.h"
#include "AlpsEnumProcessT.h"
#include "AlpsKnowledgeBroker.h"
#include "AlpsTreeNode.h"

#include "BcpsBranchStrategy.h"

//#include "PSCGBranchObjectInt.h"
//#include "PSCGConstraint.h"
//#include "PSCGHelp.h"
//#include "PSCGObjectInt.h"
#include "PSCGParams.h"
//#include "PSCGSolution.h"
//#include "PSCGVariable.h"

#define REMOVE_SLACK 1

//#############################################################################

AlpsTreeNode*
PSCGTreeNode::createNewTreeNode(AlpsNodeDesc *&desc) const
{
    // Create a new tree node
    PSCGTreeNode *node = new PSCGTreeNode(desc);

    // Set solution estimate for this nodes.
    // solEstimate = quality_ + sum_i{min{up_i, down_i}}
    
    node->setSolEstimate(solEstimate_);

#ifdef PSCG_DEBUG_MORE
    printf("PSCG:createNewTreeNode: quality=%g, solEstimate=%g\n",
           quality_, solEstimate_);
#endif
    
    desc = NULL;

    return node;
}

//#############################################################################

// NOTE: if rampup,
// - parent must be explicit if not NULL,
// - this node is explicit.

int
PSCGTreeNode::process(bool isRoot, bool rampUp)
{
    int status = PSCG_OK;

    PSCGModel* model = dynamic_cast<PSCGModel*>(desc_->getModel());

    AlpsPhase phase = knowledgeBroker_->getPhase();
    
    //------------------------------------------------------
    // Check if this can be fathomed by objective cutoff.
    //------------------------------------------------------
    
    model->setActiveNode(this);
    model->addNumNodes();
    
    //------------------------------------------------------
    // Get model information and parameters.
    //------------------------------------------------------


    // Mark if this node is root or not.
    model->isRoot_ = isRoot;

    //======================================================
    // Restore, load and solve the subproblem.
    // (1) LP infeasible
    //     a. set status to be fathom.
    // (2) LP feasible
    //     a. MILP feasible. Check whether need update incumbent.
    //     b. LP feasible but not MIP feasible. Check whether can be 
    //        fathomed, if not, choose a branch variable.
    //======================================================

    //------------------------------------------------------
    // Extract info from this node and load subproblem into lp solver.
    //------------------------------------------------------
    
    installSubProblem(model);
    status = bound(model);
    setStatus(AlpsNodeStatusFathomed);

#if 0
    while (keepOn && (pass < maxPass)) {
        ++pass;
        keepOn = false;
        
        //--------------------------------------------------
        // Bounding to get the quality of this node.
        //--------------------------------------------------
        
        status = bound(model);
	if (pass == 1) {
	    int iter = model->solver()->getIterationCount();
	    model->addNumIterations(iter);
	}
        
        switch(status) {
        case PSCG_LP_OPTIMAL:
            // Check if IP feasible 
            feasibleIP = false;//model->feasibleSolution(numIntInfs, numObjInfs);
            
            if (feasibleIP) {         
                // IP feasible 
		
		if (quality_ < cutoff) {  
                    // Better than incumbent
                    // Update cutoff
                    cutoff = getKnowledgeBroker()->getIncumbentValue();
                }
                setStatus(AlpsNodeStatusFathomed);
		//break;
            }
            else {
                cutoff = getKnowledgeBroker()->getIncumbentValue();
                if (quality_ > cutoff) {
                    setStatus(AlpsNodeStatusFathomed);
		    //break;
                }
                needBranch = true;
                reducedCostFix(model);
                
                //------------------------------------------
                // Check if tailoff
                //------------------------------------------

                if (pass > 1) {
                    improvement = quality_ - preObjValue;
                    if (improvement > tailOffTol) {
                        // NOTE: still need remove slacks, although
                        //       tailoff.
                        keepOn = true;
                    }
                    
#ifdef PSCG_DEBUG_MORE
                    std::cout << "PROCESS: pass["<< pass << "], improvement=" 
                              << improvement << ", tailOffTol=" << tailOffTol
                              << std::endl;
#endif
                }
                else {
                    keepOn = true;
                }
                // Update previous objective value.
                preObjValue = quality_;

                if ( genConsHere &&
                     //(improvement > tailOffTol) && 
                     //(numRows > numStartRows) ) {
                     (numRows > numCoreRows) ) {   
                 
		    
                    
                }
            }
            
            break;
        case PSCG_LP_ABANDONED:
#ifdef PSCG_DEBUG
            assert(0);
#endif
            status = PSCG_ERR_LP;
            goto TERM_PROCESS;
        case PSCG_LP_DUAL_INF:
            // FIXME: maybe also primal infeasible
#ifdef PSCG_DEBUG
	    assert(0);
#endif
            status = PSCG_UNBOUND;
            goto TERM_PROCESS;
        case PSCG_LP_PRIMAL_INF:
            setStatus(AlpsNodeStatusFathomed);
            quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
            goto TERM_PROCESS;
        case PSCG_LP_DUAL_LIM:
            setStatus(AlpsNodeStatusFathomed);
            quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
            goto TERM_PROCESS;
        case PSCG_LP_PRIMAL_LIM:
        case PSCG_LP_ITER_LIM:
            /* Can say much, need branch */
            needBranch = true;
#ifdef PSCG_DEBUG
            assert(0);
#endif
            goto TERM_BRANCH;
            break;
        default:
#ifdef PSCG_DEBUG
            std::cout << "PROCESS: unknown status "  <<  status << std::endl;
            assert(0);
#endif
            break;
        }

        //--------------------------------------------------
        // Apply heuristics.
        //--------------------------------------------------
        
        //--------------------------------------------------
        // Generate constraints.
        //--------------------------------------------------
        
    }
#endif 
    
    //------------------------------------------------------
    // End of process()
    //------------------------------------------------------
    
    model->isRoot_ = false;
    return status;
}

//#############################################################################

std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
PSCGTreeNode::branch()
{

    //------------------------------------------------------
    // Change one var hard bound and record the change in nodedesc:
    // THINK: how about constraint bounds? When to update?
    // TODO: how about other SOS object, etc.?
    //------------------------------------------------------

    AlpsPhase phase = knowledgeBroker_->getPhase();

    double objVal = getQuality();

    PSCGNodeDesc* childDesc = NULL;  
      
    std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
	childNodeDescs;
    
    PSCGModel* model = dynamic_cast<PSCGModel*>(desc_->getModel());    
#if 0
    int numCols = model->getNumCols();


#ifdef PSCG_DEBUG_MORE
    // Debug survived old constraints.
    int currNumOldCons = model->getNumOldConstraints();
    for (int k = 0; k < currNumOldCons; ++k) {
	PSCGConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	std::cout << "BRANCH: DBG: "
		  << "k=" << k << ", len=" << aCon->getSize()
		  << ", node=" << index_ << std::endl;
    }
#endif

    //------------------------------------------------------
    // Get branching object. TODO: Assume integer branching object. 
    //------------------------------------------------------
    
    PSCGBranchObjectInt *branchObject =
	dynamic_cast<PSCGBranchObjectInt *>(branchObject_);
    int objInd = branchObject->getObjectIndex();

    double bValue = branchObject->getValue();

    PSCGObjectInt *obj = dynamic_cast<PSCGObjectInt *>(model->objects(objInd));
    int branchVar = obj->columnIndex();
    
#ifdef PSCG_DEBUG
    if ( (branchVar < 0) || (branchVar >= numCols) ) {
	std::cout << "ERROR: BRANCH(): branchVar = " << branchVar 
		  << "; numCols = " << numCols  << std::endl;
	throw CoinError("branch index is out of range", 
			"branch", "PSCGTreeNode");
    }
#endif

#ifdef PSCG_DEBUG
    printf("BRANCH(): on %d, phase %d\n", branchVar, phase);
    printf("DOWN: lb %g, up %g\n",
	   branchObject->getDown()[0], branchObject->getDown()[1]);
    printf("UP  : lb %g, up %g\n",
	   branchObject->getUp()[0], branchObject->getUp()[1]);
#endif

    PSCGNodeDesc* thisDesc = dynamic_cast<PSCGNodeDesc*>(desc_);

        
    //======================================================
    //------------------------------------------------------
    // Create down-branch node description.
    //------------------------------------------------------
    //======================================================
    
    childDesc = new PSCGNodeDesc(model);
    
    if (phase == AlpsPhaseRampup) {
	
	//--------------------------------------------------
	// Store a full description since each node will be the root of
	// a subtree.
	// NOTE: this desc must be explicit during rampup.
	//--------------------------------------------------	

	int index, k;
	int numModify = -1;
	double value;
	
	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];
		
	double *fVarSoftLB = NULL;
	double *fVarSoftUB = NULL;
	int *fVarSoftLBInd = NULL;
	int *fVarSoftUBInd = NULL;

	//--------------------------------------------------
	// Full hard variable bounds.
	//--------------------------------------------------

	numModify = thisDesc->getVars()->lbHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	    index = thisDesc->getVars()->lbHard.posModify[k];
	    assert(index == k);
	    value = thisDesc->getVars()->lbHard.entries[k];
	    fVarHardLB[k] = value;
	    fVarHardLBInd[k] = index;
	}
	
	numModify = thisDesc->getVars()->ubHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	    index = thisDesc->getVars()->ubHard.posModify[k];
	    assert(index == k);
	    value = thisDesc->getVars()->ubHard.entries[k];
	    fVarHardUB[k] = value;
	    fVarHardUBInd[k] = index;
	}
	
	// Branching bounds.
	fVarHardLB[branchVar] = branchObject->getDown()[0];
	fVarHardUB[branchVar] = branchObject->getDown()[1];


	childDesc->assignVarHardBound(numCols,
				      fVarHardLBInd,
				      fVarHardLB,
				      numCols,
				      fVarHardUBInd,
				      fVarHardUB);

	//--------------------------------------------------
	// Soft variable bounds.
	//--------------------------------------------------

	int numSoftVarLowers = thisDesc->getVars()->lbSoft.numModify;
	assert(numSoftVarLowers >= 0 && numSoftVarLowers <= numCols);
	if (numSoftVarLowers > 0) {
	    fVarSoftLB = new double [numSoftVarLowers];
	    fVarSoftLBInd = new int [numSoftVarLowers];
	    for (k = 0; k < numSoftVarLowers; ++k) {
		index = thisDesc->getVars()->lbSoft.posModify[k];
		value = thisDesc->getVars()->lbSoft.entries[k];
		fVarSoftLB[k] = value;
		fVarSoftLBInd[k] = index;
	    }
	}
		    
	int numSoftVarUppers = thisDesc->getVars()->ubSoft.numModify;
	assert(numSoftVarUppers >= 0 && numSoftVarUppers <= numCols);
	if (numSoftVarUppers > 0) {
	    fVarSoftUB = new double [numSoftVarUppers];
	    fVarSoftUBInd = new int [numSoftVarUppers];
	    for (k = 0; k < numSoftVarUppers; ++k) {
		index = thisDesc->getVars()->ubSoft.posModify[k];
		value = thisDesc->getVars()->ubSoft.entries[k];
		fVarSoftUB[k] = value;
		fVarSoftUBInd[k] = index;
	    }
	}

#ifdef PSCG_DEBUG_MORE
	// Print soft bounds.
	std::cout << "\nBRANCH: numSoftVarLowers=" << numSoftVarLowers
		  << ", numSoftVarUppers=" << numSoftVarUppers
		  << std::endl;
	for (k = 0; k < numSoftVarLowers; ++k) {
	    std::cout << "Col[" << fVarSoftLBInd[k] << "]: soft lb="
		      << fVarSoftLB[k] << std::endl;		    
	}
	std::cout << "------------------" << std::endl;
	for (k = 0; k < numSoftVarUppers; ++k) {
	    std::cout << "Col[" << fVarSoftUBInd[k] << "]: soft ub="
		      << fVarSoftUB[k] << std::endl;
	}
	std::cout << "------------------" << std::endl << std::endl;
#endif

	// Assign it anyway so to transfer ownership of memory(fVarSoftLBInd,etc.)
	childDesc->assignVarSoftBound(numSoftVarLowers, 
				      fVarSoftLBInd,
				      fVarSoftLB,
				      numSoftVarUppers, 
				      fVarSoftUBInd,
				      fVarSoftUB);

	//--------------------------------------------------
	// Full set of non-core constraints.
	// NOTE: non-core constraints have been saved in description
	//       when process() during ramp-up.
	//--------------------------------------------------

	BcpsObject **tempCons = NULL;
	int tempInt = 0;
	
	tempInt = thisDesc->getCons()->numAdd;
	if (tempInt > 0) {
	    tempCons = new BcpsObject* [tempInt];
	    for (k = 0; k < tempInt; ++k) {
		PSCGConstraint *aCon = dynamic_cast<PSCGConstraint *>
		    (thisDesc->getCons()->objects[k]);
                
		assert(aCon);
		assert(aCon->getSize() > 0);
		assert(aCon->getSize() < 100000);
		PSCGConstraint *newCon = new PSCGConstraint(*aCon);
		tempCons[k] = newCon;
	    }
	}
	    

#if 0
	else {
	    // No cons or only root cons.
	    tempInt = model->getNumOldConstraints();

	    if (tempInt > 0) {
		tempCons = new BcpsObject* [tempInt];
	    }
	    for (k = 0; k < tempInt; ++k) {
		PSCGConstraint *aCon = model->oldConstraints()[k];                
		assert(aCon);
		assert(aCon->getSize() > 0);
		assert(aCon->getSize() < 100000);
		PSCGConstraint *newCon = new PSCGConstraint(*aCon);
		tempCons[k] = newCon;
	    }
	}
#endif

#ifdef PSCG_DEBUG_MORE
	std::cout << "BRANCH: down: tempInt=" << tempInt <<std::endl;
#endif
	// Fresh desc, safely add.
	childDesc->setAddedConstraints(tempInt, tempCons);
    }
    else {
	
	//--------------------------------------------------
	// Relative: Only need to record hard var bound change. 
	// NOTE: soft var bound changes are Record after selectBranchObject.
	//--------------------------------------------------
	
	childDesc->setVarHardBound(1,
				   &branchVar,
				   &(branchObject->getDown()[0]),
				   1,
				   &branchVar,
				   &(branchObject->getDown()[1]));
    }

    childDesc->setBranchedDir(-1);
    childDesc->setBranchedInd(objInd);
    childDesc->setBranchedVal(bValue);

    // Copy warm start.
    CoinWarmStartBasis *ws = thisDesc->getBasis();
    CoinWarmStartBasis *newWs = new CoinWarmStartBasis(*ws);
    childDesc->setBasis(newWs);
    
    childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					    (childDesc),
					    AlpsNodeStatusCandidate,
					    objVal));

    //======================================================
    //------------------------------------------------------
    // Create up-branch node description.
    //------------------------------------------------------
    //======================================================
    
    childDesc = new PSCGNodeDesc(model);

    if (phase == AlpsPhaseRampup) {

	//--------------------------------------------------
	// Store a full description since each node will be the root of
	// a subtree.
	// NOTE: parent must be explicit during rampup.
	//--------------------------------------------------	

	int index, k;
	int numModify = -1;
	double value;
	
	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];
		
	double *fVarSoftLB = NULL;
	double *fVarSoftUB = NULL;
	int *fVarSoftLBInd = NULL;
	int *fVarSoftUBInd = NULL;

	//--------------------------------------------------
	// Full hard variable bounds.
	//--------------------------------------------------
	
	numModify = thisDesc->getVars()->lbHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	    index = thisDesc->getVars()->lbHard.posModify[k];
	    assert(index == k);
	    value = thisDesc->getVars()->lbHard.entries[k];
	    fVarHardLB[k] = value;
	    fVarHardLBInd[k] = index;
	}
	
	numModify = thisDesc->getVars()->ubHard.numModify;
	assert(numModify == numCols);
	for (k = 0; k < numModify; ++k) {
	    index = thisDesc->getVars()->ubHard.posModify[k];
	    assert(index == k);
	    value = thisDesc->getVars()->ubHard.entries[k];
	    fVarHardUB[k] = value;
	    fVarHardUBInd[k] = index;
	}
	
	// Branching bounds.
	fVarHardLB[branchVar] = branchObject->getUp()[0];
	fVarHardUB[branchVar] = branchObject->getUp()[1];

	childDesc->assignVarHardBound(numCols,
				      fVarHardLBInd,
				      fVarHardLB,
				      numCols,
				      fVarHardUBInd,
				      fVarHardUB);

	//--------------------------------------------------
	// Soft variable bounds.
	//--------------------------------------------------

	int numSoftVarLowers = thisDesc->getVars()->lbSoft.numModify;
	assert(numSoftVarLowers >= 0 && numSoftVarLowers <= numCols);
	if (numSoftVarLowers > 0) {
	    fVarSoftLB = new double [numSoftVarLowers];
	    fVarSoftLBInd = new int [numSoftVarLowers];
	    for (k = 0; k < numSoftVarLowers; ++k) {
		index = thisDesc->getVars()->lbSoft.posModify[k];
		value = thisDesc->getVars()->lbSoft.entries[k];
		fVarSoftLB[k] = value;
		fVarSoftLBInd[k] = index;
	    }
	}
		    
	int numSoftVarUppers = thisDesc->getVars()->ubSoft.numModify;
	assert(numSoftVarUppers >= 0 && numSoftVarUppers <= numCols);
	if (numSoftVarUppers > 0) {
	    fVarSoftUB = new double [numSoftVarUppers];
	    fVarSoftUBInd = new int [numSoftVarUppers];
	    for (k = 0; k < numSoftVarUppers; ++k) {
		index = thisDesc->getVars()->ubSoft.posModify[k];
		value = thisDesc->getVars()->ubSoft.entries[k];
		fVarSoftUB[k] = value;
		fVarSoftUBInd[k] = index;
	    }
	}

	// Assign it anyway so to transfer ownership of memory(fVarSoftLBInd,etc.)
	childDesc->assignVarSoftBound(numSoftVarLowers, 
				      fVarSoftLBInd,
				      fVarSoftLB,
				      numSoftVarUppers, 
				      fVarSoftUBInd,
				      fVarSoftUB);

	//--------------------------------------------------
	// Full set of non-core constraints.
	// NOTE: non-core constraints have been saved in description
	//       when process() during ramp-up.
	//--------------------------------------------------
	
	BcpsObject **tempCons = NULL;
	int tempInt = 0;
	
	tempInt = thisDesc->getCons()->numAdd;
	if (tempInt > 0) {
	    tempCons = new BcpsObject* [tempInt];
	
	    for (k = 0; k < tempInt; ++k) {
		PSCGConstraint *aCon = dynamic_cast<PSCGConstraint *>
		    (thisDesc->getCons()->objects[k]);
	    
		assert(aCon);
		assert(aCon->getSize() > 0);
		assert(aCon->getSize() < 1000);
		PSCGConstraint *newCon = new PSCGConstraint(*aCon);
		tempCons[k] = newCon;
	    }
	}
	
#if 0
	else {
	    // No cons or only root cons.
	    tempInt = model->getNumOldConstraints();
	    if (tempInt > 0) {
		tempCons = new BcpsObject* [tempInt];
	    }
	    for (k = 0; k < tempInt; ++k) {
		PSCGConstraint *aCon = model->oldConstraints()[k];                
		assert(aCon);
		assert(aCon->getSize() > 0);
		assert(aCon->getSize() < 1000);
		PSCGConstraint *newCon = new PSCGConstraint(*aCon);
		tempCons[k] = newCon;
	    }
	}
#endif

#ifdef PSCG_DEBUG_MORE
	std::cout << "BRANCH: up: tempInt=" << tempInt <<std::endl;
#endif
	// Fresh desc, safely add.
	childDesc->setAddedConstraints(tempInt, tempCons);
    }
    else {

	//--------------------------------------------------
	// Relative: Only need to record hard var bound change. 
	// NOTE: soft var bound changes are Record after selectBranchObject.
	//--------------------------------------------------

	childDesc->setVarHardBound(1,
				   &branchVar,
				   &(branchObject->getUp()[0]),
				   1,
				   &branchVar,
				   &(branchObject->getUp()[1]));
    }

    childDesc->setBranchedDir(1);
    childDesc->setBranchedInd(objInd);
    childDesc->setBranchedVal(bValue);
    
    // Copy warm start.
    CoinWarmStartBasis *newWs2 = new CoinWarmStartBasis(*ws);
    childDesc->setBasis(newWs2);
    
    childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					    (childDesc),
					    AlpsNodeStatusCandidate,
					    objVal));  

    // Change node status to branched.
#endif
    status_ = AlpsNodeStatusBranched;
    
    return childNodeDescs;
}

//#############################################################################

/* FIXME: need rewrite from scratch */
/* 0: find a branch var, -1 no branch var (should not happen) */

int PSCGTreeNode::selectBranchObject(PSCGModel *model, 
                                     bool& foundSol, 
                                     int numPassesLeft) 
{
    int bStatus = 0;
#if 0
    if(branchObject_) {
        delete branchObject_;
        branchObject_ = NULL;
    }
    
    //------------------------------------------------------
    // Get branching strategy.
    //------------------------------------------------------

    BcpsBranchStrategy *strategy = model->branchStrategy();
    if (!strategy) {
        throw CoinError("No branch strategy.", "process()","PSCGTreeNode");
    }

    //------------------------------------------------------
    // Create branching object candidates.
    //-----------------------------------------------------
    
    bStatus = strategy->createCandBranchObjects(numPassesLeft,
                                                model->getCutoff());
    
    //------------------------------------------------------
    // Select the best branching objects.
    //-----------------------------------------------------
    
    if (bStatus >= 0) {
        
       branchObject_ = strategy->bestBranchObject();

       if (branchObject_) {
	  // Move best branching object to node.

#ifdef PSCG_DEBUG_MORE
	  std::cout << "SELECTBEST: Set branching obj" << std::endl;
#endif
       }
       else {
#ifdef PSCG_DEBUG
	  std::cout << "ERROR: Can't find branching object" << std::endl;
#endif
	  assert(0);
       }

       // Set guessed solution value
       // solEstimate_ = quality_ + sumDeg;
    }
    
    if (!model->branchStrategy()) {
        delete strategy;
    }
#endif
    return bStatus;
}

//#############################################################################

int PSCGTreeNode::bound(BcpsModel *model) 
{
    int status = PSCG_OK;
    PSCGModel *m = dynamic_cast<PSCGModel *>(model);
    

    m->computeBound(5);
#if 0
    if (m->solver()->isAbandoned()) {
	status = PSCG_LP_ABANDONED;
    }
    else if (m->solver()->isProvenOptimal()) {
	status = PSCG_LP_OPTIMAL;
        PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);

        double objValue = m->solver()->getObjValue() *
            m->solver()->getObjSense();
        
        int dir = desc->getBranchedDir();
        if (dir != 0) {
            double objDeg = objValue - quality_;
            int objInd = desc->getBranchedInd();
            double lpX = desc->getBranchedVal();
            PSCGObjectInt *intObject = 
                dynamic_cast<PSCGObjectInt *>(m->objects(objInd));            

            intObject->pseudocost().update(dir, objDeg, lpX);
        }

        // Update quality of this nodes.
        quality_ = objValue;
    }
    else if (m->solver()->isProvenPrimalInfeasible()) {
	status = PSCG_LP_PRIMAL_INF;
    }
    else if (m->solver()->isProvenDualInfeasible()) {
	status = PSCG_LP_DUAL_INF;
    }
    else if (m->solver()->isPrimalObjectiveLimitReached()) {
	status = PSCG_LP_PRIMAL_LIM;
    }
    else if (m->solver()->isDualObjectiveLimitReached()) {
	status = PSCG_LP_DUAL_LIM;
    }
    else if (m->solver()->isIterationLimitReached()) {
	status = PSCG_LP_ITER_LIM;
    }
    else {
	std::cout << "UNKNOWN LP STATUS" << std::endl;
	assert(0);
    }
    
#endif
    return status;
}

//#############################################################################

int PSCGTreeNode::installSubProblem(BcpsModel *m)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    PSCGModel *model = dynamic_cast<PSCGModel *>(m);
    assert(model);
    
    PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);
    
    desc->installSubproblemFromNodeDesc();
#if 0
    int numModify = 0;

    int numCoreVars = model->getNumCoreVariables();
    int numCoreCons = model->getNumCoreConstraints();

    int numCols = model->solver()->getNumCols();
    int numRows = model->solver()->getNumRows();

    //double *varSoftLB = NULL;
    //double *varSoftUB = NULL;
    double *varHardLB = NULL;
    double *varHardUB = NULL;

    //double *conSoftLB = NULL;
    //double *conSoftUB = NULL;
    //double *conHardLB = NULL;
    //double *conHardUB = NULL;
    
    double *startColLB = model->startVarLB();
    double *startColUB = model->startVarUB();
    double *startRowLB = model->startConLB();
    double *startRowUB = model->startConUB();
    
    CoinFillN(startColLB, numCoreVars, -ALPS_DBL_MAX);
    CoinFillN(startColUB, numCoreVars, ALPS_DBL_MAX);
    CoinFillN(startRowLB, numCoreCons, -ALPS_DBL_MAX);
    CoinFillN(startRowUB, numCoreCons, ALPS_DBL_MAX);

#if 0
    //memcpy(startColLB, model->origVarLB(), sizeof(double) * numCoreVars);
    //memcpy(startColUB, model->origVarUB(), sizeof(double) * numCoreVars);
    //memcpy(startRowLB, model->origConLB(), sizeof(double) * numCoreCons);
    //memcpy(startRowUB, model->origConUB(), sizeof(double) * numCoreCons);
#endif

    int numOldCons = 0;
    int tempInt = 0;
    PSCGConstraint *aCon = NULL;

    int nodeID = -1;
    nodeID = getIndex();
    //std::cout << "nodeID=" << nodeID << std::endl;
#endif
    AlpsPhase phase = knowledgeBroker_->getPhase();

    //======================================================
    // Restore subproblem: 
    //  1. Remove noncore var/con
    //  2. Travel back to root and correct differencing to
    //     full var/con bounds into model->startXXX
    //  3. Set col bounds
    //  4. Set row bounds (is this necessary?)
    //  5. Add contraints except cores
    //  6. Add variables except cores
    //  7. Set basis (should not need modify)
    //======================================================


    //------------------------------------------------------
    // Remove old constraints from lp solver.
    //------------------------------------------------------
#if 0
    int numDelCons = numRows - numCoreCons;
    
#ifdef PSCG_DEBUG
    std::cout << "INSTALL: numDelCons = " << numDelCons << std::endl;
#endif
	
    if (numDelCons > 0) {
	int *indices = new int [numDelCons];
	if (indices == NULL) {
	    throw CoinError("Out of memory", "installSubProblem", "PSCGTreeNode");
	}
	
	for (i = 0; i < numDelCons; ++i) {
	    indices[i] = numCoreCons + i;
	}
	
	model->solver()->deleteRows(numDelCons, indices);
	delete [] indices;
	indices = NULL;
    }
    
    //--------------------------------------------------------
    // Travel back to a full node, then collect diff (add/rem col/row,
    // hard/soft col/row bounds) from the node full to this node.
    //----------------------------
    // Note: if we store full set of logic/agorithm col/row, then
    //       no diff are needed for col/row
    //--------------------------------------------------------
    
    //--------------------------------------------------------
    // Collect differencing bounds. Branching bounds of this node
    // are ALSO collected.
    //--------------------------------------------------------

    PSCGNodeDesc* pathDesc = NULL;
    AlpsTreeNode *parent = parent_;    
    
    /* First push this node since it has branching hard bounds. 
       NOTE: during rampup, this desc has full description when branch(). */
    model->leafToRootPath.push_back(this);
    
    if (phase != AlpsPhaseRampup) {
	while(parent) {
#ifdef PSCG_DEBUG_MORE
	    std::cout << "Parent id = " << parent->getIndex() << std::endl;
#endif     
	    model->leafToRootPath.push_back(parent);
	    if (parent->getExplicit()) {
		// Reach an explicit node, then stop.
		break;
	    }
	    else {
		parent = parent->getParent();
	    }
	}
    }
    
#ifdef PSCG_DEBUG_MORE
    std::cout << "INSTALL: path len = " << model->leafToRootPath.size()
	      << std::endl;
#endif
    
    //------------------------------------------------------
    // Travel back from this node to the explicit node to
    // collect full description.
    //------------------------------------------------------
    
    for(i = static_cast<int> (model->leafToRootPath.size() - 1); i > -1; --i) {

#ifdef PSCG_DEBUG_MORE
        if (index_ == 3487) {
            std::cout << "\n----------- NODE ------------" 
                      << model->leafToRootPath.at(i)->getIndex() << std::endl;
        }
#endif
	
	//--------------------------------------------------
	// NOTE: As away from explicit node, bounds become 
	//       tighter and tighter.
	//--------------------------------------------------
        
        pathDesc = dynamic_cast<PSCGNodeDesc*>((model->leafToRootPath.at(i))->
                                               getDesc());
        
        varHardLB = pathDesc->getVars()->lbHard.entries;
        varHardUB = pathDesc->getVars()->ubHard.entries;
        
	//--------------------------------------------------      
        // Adjust bounds according to hard var lb/ub.
	// If rampup or explicit, collect hard bounds so far.
	//--------------------------------------------------
	
        numModify = pathDesc->getVars()->lbHard.numModify;
	
#ifdef PSCG_DEBUG_MORE
	std::cout << "INSTALL: numModify lb hard = " << numModify << std::endl;
#endif
	
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->lbHard.posModify[k];
            value = pathDesc->getVars()->lbHard.entries[k];
      
#ifdef PSCG_DEBUG_MORE
	    printf("INSTALL: 1, col %d, value %g, startColLB %x\n", 
		   index, value, startColLB);
#endif      
	    // Hard bounds do NOT change according to soft bounds, so
	    // here need std::max.
            startColLB[index] = std::max(startColLB[index], value);
            
#ifdef PSCG_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 1, col %d, hard lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
#endif
        }

#ifdef PSCG_DEBUG_MORE
	std::cout << "INSTALL: numModify ub hard = " << numModify<<std::endl;
#endif

        numModify = pathDesc->getVars()->ubHard.numModify;
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->ubHard.posModify[k];
            value = pathDesc->getVars()->ubHard.entries[k];
            startColUB[index] = std::min(startColUB[index], value);
	    
#ifdef PSCG_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 2, col %d, hard lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            if (startColLB[index] > startColUB[index]) {
                    //assert(0);
            }
#endif
        }
        
        //--------------------------------------------------
        // Adjust bounds according to soft var lb/ub.
	// If rampup or explicit, collect soft bounds so far.
        //--------------------------------------------------

        numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef PSCG_DEBUG_MORE
	std::cout << "INSTALL: i=" << i << ", numModify soft lb="
		  << numModify << std::endl;
#endif
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->lbSoft.posModify[k];
            value = pathDesc->getVars()->lbSoft.entries[k];
            startColLB[index] = std::max(startColLB[index], value);
            
#ifdef PSCG_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 3, col %d, soft lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            
            if (startColLB[index] > startColUB[index]) {
		//assert(0);
            }
#endif
        }
        numModify = pathDesc->getVars()->ubSoft.numModify;
        
#ifdef PSCG_DEBUG_MORE
	std::cout << "INSTALL: i=" << i << ", numModify soft ub="
		  << numModify << std::endl;
#endif

        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->ubSoft.posModify[k];
            value = pathDesc->getVars()->ubSoft.entries[k];
            startColUB[index] = std::min(startColUB[index], value);
            
#ifdef PSCG_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 4, col %d, soft lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            
            if (startColLB[index] > startColUB[index]) {
		//assert(0);
            }
#endif
        }
        
        //--------------------------------------------------
        // TODO: Modify hard/soft row lb/ub.
        //--------------------------------------------------


        //--------------------------------------------------
        // Collect active non-core constraints at parent.
        //--------------------------------------------------

	//----------------------------------------------
	// First collect all generated cuts, then remove 
	// deleted.
	//----------------------------------------------
	
	tempInt = pathDesc->getCons()->numAdd;
	    
#ifdef PSCG_DEBUG_MORE
	std::cout << "\nINSTALL: numAdd = " << tempInt << std::endl;
#endif

	int maxOld = model->getOldConstraintsSize();
	
	for (k = 0; k < tempInt; ++k) {
	    aCon = dynamic_cast<PSCGConstraint *>
		(pathDesc->getCons()->objects[k]);
	    
	    assert(aCon);
	    assert(aCon->getSize() > 0);
	    assert(aCon->getSize() < 100000);
	    
#ifdef PSCG_DEBUG_MORE
	    std::cout << "INSTALL: cut  k=" << k 
		      << ", len=" <<aCon->getSize() 
		      << ", node="<< index_ << std::endl;
#endif
	    (model->oldConstraints())[numOldCons++] = aCon;
	    
	    if (numOldCons >= maxOld) {
		// Need resize
#ifdef PSCG_DEBUG_MORE
		std::cout << "INSTALL: resize, maxOld = " 
			  << maxOld << std::endl;
#endif
		maxOld *= 2;
		PSCGConstraint **tempCons = new PSCGConstraint* [maxOld];
		
		memcpy(tempCons, 
		       model->oldConstraints(), 
		       numOldCons * sizeof(PSCGConstraint *));
		
		model->delOldConstraints();
		model->setOldConstraints(tempCons);
		model->setOldConstraintsSize(maxOld);
	    }
	}
	    
	//----------------------------------------------
	// Remove those deleted. 
	// NOTE: model->oldConstraints_ stores all previously 
	// generated active constraints at parent.
	//----------------------------------------------
	
	tempInt = pathDesc->getCons()->numRemove;
            
	if (tempInt > 0) {
	    int tempPos;
	    int *tempMark = new int [numOldCons];
	    CoinZeroN(tempMark, numOldCons);
	    for (k = 0; k < tempInt; ++k) {
		tempPos = pathDesc->getCons()->posRemove[k];
#ifdef PSCG_DEBUG_MORE
		std::cout << "tempPos=" << tempPos 
			  << ", tempInt=" << tempInt 
			  << ", numOldCons=" << numOldCons << std::endl;
#endif
		tempMark[tempPos] = 1;
		
	    }
	    
	    tempInt = 0;
	    for (k = 0; k < numOldCons; ++k) {
		if (tempMark[k] != 1) {
		    // Survived.
		    (model->oldConstraints())[tempInt++]=
			(model->oldConstraints())[k];
		}
	    }
	    if (tempInt + pathDesc->getCons()->numRemove != numOldCons) {
		std::cout << "INSTALL: tempInt=" << tempInt
			  <<", numRemove="<<pathDesc->getCons()->numRemove
			  << ", numOldCons=" << numOldCons << std::endl;
		
		assert(0);
	    }
	    
	    // Update number of old non-core constraints.
	    numOldCons = tempInt;
	    delete [] tempMark;
	}
    } // EOF leafToRootPath.


    //--------------------------------------------------------
    // Debug variable bounds to be installed.
    //--------------------------------------------------------

#ifdef PSCG_DEBUG_MORE
    for (k = 0; k < numCols; ++k) {
        //if (index_ == -1) {
            printf("INSTALL: Col %d, \tlb %g,  \tub %g\n",
                   k, startColLB[k], startColUB[k]);
	    //}
        
        if (startColLB[k] > startColUB[k] + ALPS_GEN_TOL) {
            printf("INSTALL: Col %d, \tlb %g,  \tub %g\n",
                   k, startColLB[k], startColUB[k]);
            assert(0);
        }
    }
#endif   

    //--------------------------------------------------------
    // Clear path vector.
    //--------------------------------------------------------
    
    model->leafToRootPath.clear();
    assert(model->leafToRootPath.size() == 0);
    
    //--------------------------------------------------------
    // Adjust column bounds in lp solver  
    //--------------------------------------------------------
    
    for(i = 0; i < numCols; ++i) {
	model->solver()->setColBounds(i, startColLB[i], startColUB[i]); 
    }
    
    //--------------------------------------------------------
    // TODO: Set row bounds 
    //--------------------------------------------------------
    

    //--------------------------------------------------------
    // Add old constraints, which are collect from differencing.
    //--------------------------------------------------------
    
    // If removed cuts due to local cuts.

    model->setNumOldConstraints(numOldCons);
    
#ifdef PSCG_DEBUG
    std::cout << "INSTALL: after collecting, numOldCons = " << numOldCons 
	      << std::endl;
#endif

    if (numOldCons > 0) {
	const OsiRowCut ** oldOsiCuts = new const OsiRowCut * [numOldCons];
	for (k = 0; k < numOldCons; ++k) {
	    OsiRowCut * acut = 
		PSCGConstraintToOsiCut(model->oldConstraints()[k]);
	    oldOsiCuts[k] = acut;
	}
	model->solver()->applyRowCuts(numOldCons, oldOsiCuts);
	for (k = 0; k < numOldCons; ++k) {
	    delete oldOsiCuts[k];
	}
	delete [] oldOsiCuts;
	oldOsiCuts = NULL;
    }
    
    //--------------------------------------------------------
    // Add parent variables, which are collect from differencing.
    //--------------------------------------------------------
    

    //--------------------------------------------------------
    // Set basis
    //--------------------------------------------------------

    CoinWarmStartBasis *pws = desc->getBasis();

    if (pws != NULL) {
	model->solver()->setWarmStart(pws);

#ifdef PSCG_DEBUG
	printf("NODE %d: set warm start\n", getIndex());
	
	numCols = model->solver()->getNumCols();
	numRows = model->solver()->getNumRows();
	int nStr = pws->getNumStructural();
	int nArt = pws->getNumArtificial();
	
	if (numCols != nStr) {
	    std::cout << "nStr=" << nStr << ", numCols=" << numCols 
		      << std::endl;
	    assert(0);
	}
	std::cout << "nArt=" << nArt << ", numRows=" << numRows 
		  << std::endl;
	if (numRows != nArt) {
	    std::cout << "nArt=" << nArt << ", numRows=" << numRows 
		      << std::endl;
	    assert(0);
	}
#endif

    }  
#endif
    return status;
}

//#############################################################################

int 
PSCGTreeNode::generateConstraints(PSCGModel *model, OsiCuts & cutPool) 
{
    int status = PSCG_LP_OPTIMAL;
#if 0
    int i, numCGs;
    int preNumRowCons = 0;
    int preNumColCons = 0;
    int newCons = 0;
    int strategy = -2;
    int maxStrategy = -2;
    
    bool mustResolve = false;
    bool fullScan = true;
    
    double useTime;
    
    numCGs = model->numCutGenerators();
  
    
    for (i = 0 ; i < numCGs; ++i) {
	
	//----------------------------------------------------
	// Check if call this generator.
	//----------------------------------------------------
	
	strategy =  model->cutGenerators(i)->strategy();
	maxStrategy = ALPS_MAX(strategy, maxStrategy);
      
	bool useThis = false;
	if (strategy == -2) {
	    useThis = false;
	}
	else if (strategy == -1) {
	    if (model->isRoot_) useThis = true;
	}
	else if (strategy == 0) {
	    useThis = true;
	}
	else if (strategy > 0) {
	    // Num of nodes is set at the beginning of process().
	    int numNodes = model->getNumNodes();
	    if ((numNodes-1) % strategy == 0) {
		useThis = true;
	    }
	}
	
#ifdef PSCG_DEBUG_MORE
	std::cout<<"CUTGEN: " << model->cutGenerators(i)->name() 
		 <<": useThis="<<useThis
		 << ", strategy=" << strategy 
		 << ", num of nodes=" << model->getNumNodes()
		 <<std::endl;
#endif

	//----------------------------------------------------
	// Generator constraints.
	//----------------------------------------------------

	if (useThis) {
	    newCons = 0;
	    preNumRowCons = cutPool.sizeRowCuts();
	    preNumColCons = cutPool.sizeColCuts();
          
	    useTime = CoinCpuTime();
	    mustResolve = 
		model->cutGenerators(i)->generateCons(cutPool, fullScan);
	    useTime = CoinCpuTime() - useTime;

	    if (mustResolve) {
		// TODO: Only probing will return ture.
		status = bound(model);
		if (status == PSCG_LP_OPTIMAL) {
#ifdef PSCG_DEBUG
		    std::cout << "CUTGEN: after probing, this node survived."
			      << std::endl;
#endif             
		}
		else {
#ifdef PSCG_DEBUG
		    std::cout<<"CUTGEN: after probing, this node can fathomed."
			     << std::endl;
#endif
		    break;
		}
	    }
	    
	    //------------------------------------------------
	    // Modify control. 
	    // NOTE: only modify if user choose automatic.
	    //------------------------------------------------
	    
	    if ( (model->useCons() == 0) && 
		 (model->cutGenerators(i)->noConsCalls() > 30) ) {
		// disable.
		model->cutGenerators(i)->setStrategy(-2);
		maxStrategy = ALPS_MAX(strategy, maxStrategy);
	    }
	}
    }
    
    if ( (model->useCons() == 0) && (maxStrategy == -2) ){
	// Previously automatic, now all cut generators have been disabled.
	model->setUseCons(-2);
    }
#endif
    
    return status;
}

//#############################################################################

int 
PSCGTreeNode::applyConstraints(PSCGModel *model, 
                               OsiCuts & osiCutSet,
                               const double *solution)
{
    int status = PSCG_OK;
#if 0
    int i, k;
    
    //int numColumnCuts = osiCutSet.sizeColCuts() ;
    int numRowCuts = osiCutSet.sizeRowCuts();
    int numToAdd = numRowCuts;
    int numAdded = 0;
    
    if (numRowCuts > 0) {
	PSCGParams * PSCGPar = model->PSCGPar();
	double scaleConFactor = PSCGPar->entry(PSCGParams::scaleConFactor);
        
        if (numToAdd > 0) { 
	    
            OsiRowCut *rowCut = NULL;
	    
#ifdef PSCG_DEBUG
            printf("\nAPPLYCUT: before selecting, num of new cuts = %d\n",
                   numRowCuts);
#endif
            int numRowsNow = model->solver()->getNumRows();
            int numCols = model->solver()->getNumCols();
            CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis*>
                (model->solver()->getWarmStart());
            
            const OsiRowCut ** addCuts = new const OsiRowCut * [numToAdd];

            for (i = 0 ; i < numToAdd ; i++) {
		bool keep = true;
		
                rowCut = &(osiCutSet.rowCut(i));
                
		//------------------------------------------
		// Remove:
		//  - empty cuts
		//  - dense cuts
		//  - bad scaled cuts
                //  - weak cuts
		//  - parallel cuts
		//------------------------------------------

		const CoinPackedVector & rowVector = rowCut->row();
		int length = rowVector.getNumElements();
                bool check = true;
                
                while (check) {
                    //--------------------------------------                   
                    // Empty.
                    //--------------------------------------

                    if (length <= 0) {
                        keep = false;

#ifdef BLIS_DEBUG_MORE
                        std::cout << "APPLYCUT: A empty cut." << std::endl;
#endif
                        break;
                    }

                    //--------------------------------------
                    // Dense.
                    //--------------------------------------

                    if(length > model->getDenseConCutoff()){
                        keep = false;
#ifdef BLIS_DEBUG
                        std::cout << "APPLYCUT: A dense cut. length = " 
                                  << length << ", cutoff = " 
                                  << model->getDenseConCutoff() << std::endl;
#endif  
                        break;
                    }

                    //--------------------------------------
                    // Compuate scale factor.
                    //--------------------------------------

                    const double *elements = rowVector.getElements();
                    const int *indices = rowVector.getIndices();
                    
                    int index;
                    double activity = 0.0;
                    
                    double maxElem = 0.0;
                    double minElem = ALPS_DBL_MAX;
                    double scaleFactor;
                    
                    for (k = 0; k < length; ++k) {
                        if (fabs(elements[k]) > maxElem) {
                            maxElem = fabs(elements[k]);
                        }
                        if (fabs(elements[k]) < minElem) {
                            minElem = fabs(elements[k]);
                        }
                        index = indices[k];
                        activity += elements[k] * solution[k];
                    }
                    if(minElem != 0.0) {
                        scaleFactor = maxElem/minElem;
                    }
                    else {
                        assert(0);
                        scaleFactor = ALPS_DBL_MAX;
                    }
                    
#ifdef BLIS_DEBUG
                    std::cout << "APPLYCUT: scaleFactor=" << scaleFactor
                              << ", maxElem=" << maxElem 
                              << ", minElem=" << minElem << std::endl;
#endif
                    if (scaleFactor > scaleConFactor) {
#ifdef BLIS_DEBUG
                        std::cout<< "APPLYCUT: remove a bad scaled cut"
                                 << std::endl;
#endif
                        keep = false;
                        break;
                    }
                    
                    //--------------------------------------
                    // Weak.
                    //--------------------------------------

                    char cutSense = rowCut->sense();
                    double cutRhs = rowCut->rhs();
                    double violation = -1.0;
                    double range;
                    
                    switch(cutSense) {
                    case 'E':
                        violation = fabs(activity - cutRhs);
                        break;
                    case 'G':
                        violation = cutRhs - activity;
                        break;
                    case 'L':
                        violation = activity - cutRhs;
                    case 'R':
                        range = rowCut->range();
                        violation = ALPS_MAX(violation, activity - cutRhs);
                        violation = ALPS_MAX(violation, cutRhs-range-activity);
                        break;
                    case 'N':
                    default:
                        throw CoinError("Unknown cut sense", 
                                        "applyConstraint", "PSCGTreeNode");
                        break;
                    }
                    
                    if (violation < 1.0e-6) {
                        // Found a weak cuts.
#ifdef PSCG_DEBUG
                        std::cout<< "APPLYCUT: remove a weak cut, violation="
                                 << violation << std::endl;
#endif
                        keep = false;
                        break;
                    }
                    
                    //--------------------------------------
                    // Parallel cuts.
                    //--------------------------------------
                    
		    bool paral = parallel(model, 
					  &osiCutSet,
					  i,
					  rowCut);
		    if (paral) {
#ifdef PSCG_DEBUG
                        std::cout<< "APPLYCUT: remove a parallel"<< std::endl;
#endif
			keep = false;
			break;
		    }
                    
                    //--------------------------------------
                    // Check once and stop.
                    //--------------------------------------

                    check = false;
                }//while
                
		if (keep) {
                    addCuts[numAdded++] = &(osiCutSet.rowCut(i));
                }
                else {
                    osiCutSet.eraseRowCut(i);
                    --i;
                    --numToAdd;
                }


            }
	    

#ifdef PSCG_DEBUG
            printf("APPLYCUT: after selecting, num of new cuts = %d\n\n",
                   numAdded);
#endif
            
            //----------------------------------------------
            // Add cuts to lp and adjust basis.
            //----------------------------------------------

            model->solver()->applyRowCuts(numAdded, addCuts);
            delete [] addCuts;
            
            ws->resize(numRowsNow + numToAdd, numCols);
            for (i = 0 ; i < numToAdd; ++i) { 
                ws->setArtifStatus(numRowsNow + i,
                                   CoinWarmStartBasis::basic); 
            }
            if (model->solver()->setWarmStart(ws) == false) { 
                throw CoinError("Fail setWarmStart() after cut installation.",
                                "applyConstraints","PSCGTreeNode"); 
            }
            delete ws;
        }   
    }
    
#endif
    return status;
}

//#############################################################################

int PSCGTreeNode::
reducedCostFix(PSCGModel *model)
{ 
    int i, var;
    int status = PSCG_OK;
#if 0
    int numFixedUp = 0;
    int numFixedDown = 0;
    int numTighten = 0;

    double movement;    
    double newBound;
    double boundDistance;
    double dj;

    const double *lb = model->solver()->getColLower();
    const double *ub = model->solver()->getColUpper();
    const double *solution = model->solver()->getColSolution();
    const double *reducedCost = model->solver()->getReducedCost();
    
    double cutup = getKnowledgeBroker()->getIncumbentValue() *
        model->solver()->getObjSense();
    
    if (cutup >= ALPS_OBJ_MAX) return status;
    
    double lpObjValue = model->solver()->getObjValue() * 
        model->solver()->getObjSense();
    double epInt = 1.0e-5;
    
    int numIntegers = model->getNumIntVars();
    const int *intIndices = model->getIntVars();
    
    for (i = 0; i < numIntegers; ++i) { 
	var = intIndices[i];
        
        dj = reducedCost[var];
        
        if (fabs(dj) < epInt) continue;
        
        boundDistance = ub[var] - lb[var];
	if (boundDistance < epInt) continue;
	
        movement = floor((cutup - lpObjValue) / fabs(dj));
        
        if (solution[var] > ub[var] - epInt) {
            /* At upper bound */
            if (movement < boundDistance) {
                /* new lower bound. If movement is 0, then fix. */
                newBound = ub[var] - movement;
                newBound = std::min(newBound, ub[var]);

#ifdef PSCG_DEBUG_MORE
                printf("RED-FIX: dj %g, lb %.10g, ub %.10g, newBound %.10g, movement %g\n", dj, lb[var], ub[var], newBound, movement);
#endif
                
                if (movement <= ALPS_ZERO) {
                    ++numFixedUp;
                }
                else if (newBound < ub[var]){
                    ++numTighten;
                }
                model->solver()->setColLower(var, newBound);
            }
        }
        else if (solution[var] < lb[var] + epInt) {
            /* At lower bound */
            if (movement < boundDistance) {
                newBound = lb[var] + movement;
                newBound = std::max(newBound, lb[var]);
                
#ifdef PSCG_DEBUG_MORE
                printf("RED-FIX: dj %g, lb %g, ub %g, newBound %g, movement %g\n", dj, lb[var], ub[var], newBound, movement);
#endif
		
                if (movement <= ALPS_ZERO) {
                    ++numFixedDown;
                }
                else if(newBound > lb[var] ){
                    ++numTighten;
                }
                /* new upper bound. If movement is 0, then fix. */
                model->solver()->setColUpper(var, newBound);
            }
        }
    }
    
    //int change = numFixedUp + numFixedDown + numTighten;
    //model->reducedCostFixed_ += change;
#ifdef PSCG_DEBUG_MORE
    if (numFixedUp > 0 || numFixedDown > 0 || numTighten > 0) {
        printf("reducedCostFix: numFixedUp = %d, numFixedDown = %d, numTighten %d\n", numFixedUp, numFixedDown, numTighten);
    }
#endif
#endif
    return status; 
}

//#############################################################################

AlpsEncoded*
PSCGTreeNode::encode() const 
{
#ifdef PSCG_DEBUG
    std::cout << "PSCGTreeNode::encode()--start to encode node "
	      << index_ << std::endl;
#endif

    AlpsReturnStatus status = AlpsReturnStatusOk;
    
    // NOTE: "AlpsKnowledgeTypeNode" is used as type name.
    AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeNode);
    
    // Encode decription.
    status = desc_->encode(encoded);
    
    // Encode Alps portion.
    status = encodeAlps(encoded);
    
    // Encode Bcps portion.
    status = encodeBcps(encoded);
    
    // Nothing to encode for PSCG portion.
    
    return encoded;
}

//#############################################################################

AlpsKnowledge* 
PSCGTreeNode::decode(AlpsEncoded& encoded) const 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    PSCGTreeNode* treeNode = NULL;

    PSCGModel *model = dynamic_cast<PSCGModel*>(desc_->getModel());
    
    //------------------------------------------------------
    // Unpack decription.
    //------------------------------------------------------

    AlpsNodeDesc* nodeDesc = new PSCGNodeDesc(model);
    status = nodeDesc->decode(encoded);
    
    //------------------------------------------------------
    // Unpack node.
    //------------------------------------------------------
    
    // Unpack Alps portion.
    treeNode = new PSCGTreeNode(nodeDesc);
    treeNode->decodeAlps(encoded);
    
    // Unpack Bcps portion.
    int type = 0;
    encoded.readRep(type);	
    if (type == PSCG_BO_INT) {
	// branchObject_ is simple integer.
	PSCGBranchObjectInt *bo = new PSCGBranchObjectInt();
	status = bo->decode(encoded);

	// Set bo in treeNode.
	treeNode->setBranchObject(bo);
    }
    
    // Nothing to unpack for PSCG portion.
    
    return treeNode;
}

//#############################################################################

void 
PSCGTreeNode::convertToExplicit() 
{
#ifdef PSCG_DEBUG
    std::cout << "PSCG: convertToExplicit(); explicit_="<<explicit_ << std::endl;
#endif
#if 0
    if(!explicit_) {
	
	// Convert to explicit
	explicit_ = 1;
	
	PSCGModel* model = dynamic_cast<PSCGModel*>(desc_->getModel());
	PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc *>(desc_);
	PSCGConstraint *aCon = NULL;
	
	int numCols = model->solver()->getNumCols();

	int i, k, index;
	int tempInt;

	int numModify = 0;
	int numSoftVarLowers = 0;
	int numSoftVarUppers = 0;

	double value;
		
	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];
	
	double *fVarSoftLB = new double [numCols];
	double *fVarSoftUB = new double [numCols];
	int *fVarSoftLBInd = new int [numCols];
	int *fVarSoftUBInd = new int [numCols];

	for (k = 0; k < numCols; ++k) {
	    fVarSoftLB[k] = ALPS_DBL_MAX;
	    fVarSoftUB[k] = -ALPS_DBL_MAX;
	    fVarHardLB[k] = ALPS_DBL_MAX;
	    fVarHardUB[k] = -ALPS_DBL_MAX;
	    fVarHardLBInd[k] = k;
	    fVarHardUBInd[k] = k;
	}

	int numOldCons = 0;
	int maxOld = model->getOldConstraintsSize();
	BcpsObject ** oldConstraints = new BcpsObject* [maxOld];
	
	//--------------------------------------------------
	// Travel back to a full node, then collect diff (add/rem col/row,
	// hard/soft col/row bounds) from the node full to this node.
	//--------------------------------------------------------
    
	PSCGNodeDesc* pathDesc = NULL;
	AlpsTreeNode *parent = parent_;    
	
	model->leafToRootPath.push_back(this);
	
	while(parent) {
#ifdef PSCG_DEBUG_MORE
	    std::cout << "Parent id = " << parent->getIndex() << std::endl;
#endif     
	    model->leafToRootPath.push_back(parent);
	    if (parent->getExplicit()) {
		// Reach an explicit node, then stop.
		break;
	    }
	    else {
		parent = parent->getParent();
	    }
	}
    
#ifdef PSCG_DEBUG
	std::cout << "CONVERT TO EXP: path len = " << model->leafToRootPath.size()
		  << std::endl;
#endif
	
	//------------------------------------------------------
	// Travel back from this node to the explicit node to
	// collect full description.
	//------------------------------------------------------	


	for(i = static_cast<int> (model->leafToRootPath.size() - 1); i > -1; --i) {

	    pathDesc = dynamic_cast<PSCGNodeDesc*>((model->leafToRootPath.at(i))->
						   getDesc());

	    //--------------------------------------
	    // Full variable hard bounds.
	    //--------------------------------------
	    
	    numModify = pathDesc->getVars()->lbHard.numModify;
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->lbHard.posModify[k];
		value = pathDesc->getVars()->lbHard.entries[k];
		fVarHardLB[index] = value;
	    }
		    
	    numModify = pathDesc->getVars()->ubHard.numModify;
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->ubHard.posModify[k];
		value = pathDesc->getVars()->ubHard.entries[k];
		fVarHardUB[index] = value;
	    }
		    
	    //--------------------------------------
	    // Full variable soft bounds.
	    //--------------------------------------
	    
	    numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef PSCG_DEBUG
	    std::cout << "CONVERT: EXP: i=" << i << ", numModify soft lb="
		      << numModify << std::endl;
#endif
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->lbSoft.posModify[k];
		value = pathDesc->getVars()->lbSoft.entries[k];
		fVarSoftLB[index] = value;
	    }
		    
	    numModify = pathDesc->getVars()->ubSoft.numModify;
#ifdef PSCG_DEBUG
	    std::cout << "CONVERT: EXP: i=" << i << ", numModify soft ub="
		      << numModify << std::endl;
#endif
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->ubSoft.posModify[k];
		value = pathDesc->getVars()->ubSoft.entries[k];
		fVarSoftUB[index] = value;
	    }


            //----------------------------------------------
            // Collect all generated constraints, then remove deleted.
            //----------------------------------------------
            
            tempInt = pathDesc->getCons()->numAdd;
	    
#ifdef PSCG_DEBUG
            std::cout << "\nCONVERT: EXP: numAdd = " << tempInt << std::endl;
#endif
            
            for (k = 0; k < tempInt; ++k) {
                aCon = dynamic_cast<PSCGConstraint *>
                    (pathDesc->getCons()->objects[k]);
                
                assert(aCon);
                assert(aCon->getSize() > 0);
                assert(aCon->getSize() < 100000);
                
#ifdef PSCG_DEBUG
		std::cout << "CONVERT: EXP: k=" << k 
			  << ", len=" <<aCon->getSize() << std::endl;
#endif
                oldConstraints[numOldCons++] = aCon;
                
                if (numOldCons >= maxOld) {
                    // Need resize
#ifdef PSCG_DEBUG
                    std::cout << "CONVERT: EXP: resize, maxOld = " 
                              << maxOld << std::endl;
#endif
                    maxOld *= 2;
                    BcpsObject **tempCons = new BcpsObject* [maxOld];
		    
                    memcpy(tempCons, 
                           oldConstraints, 
                           numOldCons * sizeof(BcpsObject *));
		    
		    delete [] oldConstraints;
		    oldConstraints = tempCons;
		    tempCons = NULL;
                }
            }
	    
            //----------------------------------------------
            // Remove those deleted. 
	    // NOTE: oldConstraints stores all previously 
	    // generated active constraints at parent.
            //----------------------------------------------
	    
            tempInt = pathDesc->getCons()->numRemove;
            
            if (tempInt > 0) {
                int tempPos;
                int *tempMark = new int [numOldCons];
                CoinZeroN(tempMark, numOldCons);
                for (k = 0; k < tempInt; ++k) {
                    tempPos = pathDesc->getCons()->posRemove[k];
#ifdef PSCG_DEBUG_MORE
		    std::cout << "tempPos=" << tempPos 
			      << ", tempInt=" << tempInt 
			      << ", numOldCons=" << numOldCons << std::endl;
#endif
                    tempMark[tempPos] = 1;
                }
		
                tempInt = 0;
                for (k = 0; k < numOldCons; ++k) {
                    if (tempMark[k] != 1) {
                        // Survived.
                        oldConstraints[tempInt++] = oldConstraints[k];
                    }
                }
		if (tempInt + pathDesc->getCons()->numRemove != numOldCons) {
		    std::cout << "INSTALL: tempInt=" << tempInt
			      <<", numRemove="<<pathDesc->getCons()->numRemove
			      << ", numOldCons=" << numOldCons << std::endl;
		    
		    assert(0);
		}

                // Update number of old non-core constraints.
                numOldCons = tempInt;
                delete [] tempMark;
            }
	    
	} // EOF for (path)
        
	//------------------------------------------
	// Record hard variable bounds. FULL set.
	//------------------------------------------
	
	desc->assignVarHardBound(numCols,
				 fVarHardLBInd,
				 fVarHardLB,
				 numCols,
				 fVarHardUBInd,
				 fVarHardUB);

	//------------------------------------------
	// Recode soft variable bound. Modified.
	//------------------------------------------
	
	for (k = 0; k < numCols; ++k) {
	    if (fVarSoftLB[k] < ALPS_BND_MAX) {
		fVarSoftLBInd[numSoftVarLowers] = k;
		fVarSoftLB[numSoftVarLowers++] = fVarSoftLB[k];
	    }
	    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
		fVarSoftUBInd[numSoftVarUppers] = k;
		fVarSoftUB[numSoftVarUppers++] = fVarSoftUB[k];
	    }
	}
	// Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
	desc->assignVarSoftBound(numSoftVarLowers, 
				 fVarSoftLBInd,
				 fVarSoftLB,
				 numSoftVarUppers, 
				 fVarSoftUBInd,
				 fVarSoftUB);
	
	//------------------------------------------
	// Recode added constraints.
	//------------------------------------------

	// First make a hard copy.
	for (k = 0; k < numOldCons; ++k) {
	    aCon = dynamic_cast<PSCGConstraint *>(oldConstraints[k]);
	    assert(aCon);
	    
	    PSCGConstraint *newCon = new PSCGConstraint(*aCon);
	    oldConstraints[k] = newCon;
	}
	// Add will first delete, then add. It is safe to use in parallel.
	desc->setAddedConstraints(numOldCons, oldConstraints);

	//------------------------------------------
	// Recode deleted constraints.
	//------------------------------------------
	
	desc->delConstraints(0, NULL);
	
	//--------------------------------------------------
	// Clear path vector.
	//--------------------------------------------------
	
	model->leafToRootPath.clear();
	assert(model->leafToRootPath.size() == 0);


    } // EOF of if.
#endif
}

//#############################################################################

// Not defined yet.
void
PSCGTreeNode::convertToRelative()
{
    if(explicit_) {
	

    }
}

//#############################################################################

bool 
PSCGTreeNode::parallel(PSCGModel *model, 
		       OsiCuts *newCutSet,
		       int lastNew,
		       OsiRowCut *rowCut)
{
    bool parallel = false;
#if 0
    int k;
    double threshold = 0.999;
    
    //------------------------------------------------------
    // Compare with old cuts
    //------------------------------------------------------

    int numOldCons = model->getNumOldConstraints();
    for (k = 0; k < numOldCons; ++k) {
	PSCGConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	parallel = PSCGParallelCutCon(rowCut,
				      aCon,
				      threshold);
	if (parallel) return parallel;
    }
    

    //------------------------------------------------------
    // Compare with old cuts
    //------------------------------------------------------

    for (k = 0; k < lastNew; ++k) {
	OsiRowCut *rowCut2 = &(newCutSet->rowCut(k));
	parallel = PSCGParallelCutCut(rowCut,
				      rowCut2,
				      threshold);
	if (parallel) return parallel;
    }
#endif 
    return parallel;
}

//#############################################################################


