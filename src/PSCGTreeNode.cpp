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

//#include "CoinUtility.hpp"

#include "AlpsKnowledge.h"
#include "AlpsEnumProcessT.h"
#include "AlpsKnowledgeBroker.h"
#include "AlpsTreeNode.h"

#include "BcpsBranchStrategy.h"

#include "PSCGParams.h"
#include "PSCGTreeNode.h"

#define REMOVE_SLACK 1

//#############################################################################

AlpsTreeNode*
PSCGTreeNode::createNewTreeNode(AlpsNodeDesc *&desc) const
{
    // Create a new tree node
    PSCGTreeNode *node = new PSCGTreeNode(desc);


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

int PSCGTreeNode::process(bool isRoot, bool rampUp)
{
    int status = PSCG_OK;

    //------------------------------------------------------
    // Get model information and parameters.
    //------------------------------------------------------
    PSCGModel* model = dynamic_cast<PSCGModel*>(desc_->getModel());
    #ifndef KEEP_OMEGA
    if(model->getBestNodeQuality()==model->getBound()){model->saveOmega();}
    #endif
    bool isMPIRoot = model->getMPIRank()==0;
    int nodeDepth=getDepth();
if(isMPIRoot) {cout << "***************BEGINNING OF PROCESSING NODE*******************" << endl;}
if(isMPIRoot) {cout << "This node has depth: " << nodeDepth << endl;}
    model->printBestBounds();
    PSCGNodeDesc* desc = dynamic_cast<PSCGNodeDesc*>(desc_);
    //------------------------------------------------------
    // Check if this can be fathomed by objective cutoff.
    //------------------------------------------------------
    if(quality_ >= this->getKnowledgeBroker()->getIncumbentValue()){
	setStatus(AlpsNodeStatusFathomed);
	if(isMPIRoot){ 
	    cout << "Fathomed by bound" << endl;
	    cout << "***************END OF PROCESSING NODE*******************" << endl;
	}
        desc->freeNodeInfo();
    	return status;
    }
    
    //------------------------------------------------------
    // Extract info from this node and load subproblem into MIP solver.
    //------------------------------------------------------
    
    if(isMPIRoot) cout << "Number of nodes left: " << this->getKnowledgeBroker()->getNumNodeLeftSystem() << endl;;
    // Restore, load and solve the subproblem.
    installSubProblem(model);


    //------------------------------------------------------
    // Compute the bound.
    //------------------------------------------------------
    desc->printZBounds(model);
    bound(model); //node statuses are set here
    

    int z_status = model->getZStatus();
    if( model->getZStatus()==Z_UNKNOWN ){
      chooseBranchingObject(model);

      if(model->getNumNewNodeSPs()==0){
	setStatus(AlpsNodeStatusFathomed);
        //desc->freeNodeInfo();
	if(isMPIRoot) cout << "Fathomed by optimality" << endl;
	//model->evaluateFeasibleZ();
        model->findPrimalFeasSoln(20);
      }
      else{
        //desc->updateZ(model);
        //desc->updateBranchingZVals(model);
    	//desc->updateOmega(model); 
    	//desc->updateDesc(model);
	setStatus(AlpsNodeStatusPregnant);
	//free node information later when branching occurs
	if(isMPIRoot) cout << "Node needs to branch...with branching information:" << endl;
	if(isMPIRoot) model->printNewNodeSPInfo();
        if(isMPIRoot) cout << "Searching for feasible solution, improvement on incumbent value..." << endl;
        model->findPrimalFeasSoln(20);
      }
        //desc->freeNodeInfo();

    //------------------------------------------------------
    // Try to find a feasible solution, improve the incumbent value 
    //------------------------------------------------------
      //if(isMPIRoot) cout << "Searching for feasible solution, improvement on incumbent value..." << endl;
      //model->findPrimalFeasSoln(20);
    }
    else{
        if(model->getZStatus()==Z_INFEAS){
	    if(isMPIRoot) cout << "Fathomed due to infeasibility" << endl;
        }
	else{ //z_status==Z_BOUNDED
	    if(isMPIRoot) cout << "Fathomed by bound" << endl;
        }
        //desc->freeNodeInfo();
	setStatus(AlpsNodeStatusFathomed);
    }
    desc->freeNodeInfo();
    
    //------------------------------------------------------
    // End of process()
    //------------------------------------------------------
model->printNodeStats();
    
if(isMPIRoot) {cout << "***************END OF PROCESSING NODE*******************" << endl;}
    return status;
}

//#############################################################################

std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
PSCGTreeNode::branch()
{
//cout << "Begin branch()" << endl;
//printInstallSubproblem();
    //------------------------------------------------------
    // Change one var hard bound and record the change in nodedesc:
    // THINK: how about constraint bounds? When to update?
    // TODO: how about other SOS object, etc.?
    //------------------------------------------------------
    PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);
#if 0
    AlpsPhase phase = knowledgeBroker_->getPhase();


    //std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
	//childNodeDescs;
    
    PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);

    assert(desc->getBranchingIndex()!=-1);
    PSCGNodeDesc *childDesc = desc->createChildNodeDescUp(desc->getBranchingIndex());
    childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					    (childDesc),
					    AlpsNodeStatusCandidate,
					    quality_));
    childDesc = desc->createChildNodeDescDown(desc->getBranchingIndex());
    childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					    (childDesc),
					    AlpsNodeStatusCandidate,
					    quality_));
#endif
    status_ = AlpsNodeStatusBranched;
    
    //desc->freeNodeInfo();
//cout << "End branch()" << endl;
    return childNodeDescs;
}

//#############################################################################


//#############################################################################

int PSCGTreeNode::bound(BcpsModel *model) 
{
//cout << "Begin bound()" << endl;
    PSCGModel *m = dynamic_cast<PSCGModel *>(model);
    PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);
    int nodeDepth=getDepth();
    //updateQuality(m->computeBound(20,true));
    //if(m->getAlgorithm()==ALGDDBASIC){updateQuality(m->computeBound(100,true));}
    //updateQuality(m->computeBound(min(500,20*(nodeDepth+1)),true));
    //updateQuality(m->computeBound(min(500,2*(nodeDepth)+5),true));
    //updateQuality(m->computeBound(2*(nodeDepth)+5,true));
    //if(nodeDepth==0){
    if(false){
	updateQuality(m->computeBound(500,true));
    }
    else{
	updateQuality(m->computeBound((nodeDepth*(nodeDepth+1))/2,true));
    }
    //updateQuality(m->computeBound(5,true));
    //updateQuality(m->computeBound(5+2*nodeDepth,true));
    //if(false){updateQuality(m->computeBound(100,true));}
    //else{updateQuality(m->computeBound(5+2*nodeDepth,true));}
    //updateQuality(m->computeBound(5,true));
    //if( m->getZStatus()==Z_UNKNOWN ){
      //desc->updateOmega(m); 
    //}
//cout << "End bound()" << endl;
    return 0;
}

//#############################################################################

int PSCGTreeNode::installSubProblem(BcpsModel *m)
{
//cout << "Begin installSubProblem()" << endl;
    AlpsReturnStatus status = AlpsReturnStatusOk;

    PSCGModel *model = dynamic_cast<PSCGModel *>(m);
    assert(model);
    
    PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);
    
    desc->installSubproblemFromNodeDesc();

//cout << "End installSubProblem()" << endl;
    return status;
}

//#############################################################################

int 
PSCGTreeNode::generateConstraints(PSCGModel *model, OsiCuts & cutPool) 
{
    int status = PSCG_OPTIMAL;
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
#if 0
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
#endif 
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
#if 0
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
#endif
//#############################################################################


