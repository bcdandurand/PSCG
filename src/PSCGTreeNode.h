/*===========================================================================*
 * This file is part of the Bcps Linear Solver (BLIS).                       *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
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

#ifndef PSCGTreeNode_h_
#define PSCGTreeNode_h_

//#############################################################################

#include "AlpsNodeDesc.h"

#include "BcpsObjectPool.h"
#include "BcpsTreeNode.h"
#include "BcpsNodeDesc.h"
#include "PSCGNodeDesc.h"
#if 1
#include "PSCGModel.h"
#include "CoinUtility.hpp"

class BcpsModel;
//class PSCGModel;


//#############################################################################
/** This is the class in which we are finally able to concretely define the
    bounding procedure. Here we can assume that we have an LP solver and that
    the objects are cuts and variables, etc. */
//#############################################################################

class PSCGTreeNode : public BcpsTreeNode {
 private:

    /** No copy constructor, assignment operator. */
    PSCGTreeNode(const PSCGTreeNode&);

    PSCGTreeNode& operator=(const PSCGTreeNode&);
    
    /** Constraint pool. */
    //BcpsConstraintPool *constraintPool_;
    
    /** Variable pool. */
    //BcpsVariablePool *variablePool_;

    /** Save an explicit node description. */
    //void saveExplicit();
    std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
	childNodeDescs;
    
 public:

    /** Default constructor. */
#if 1
    PSCGTreeNode() : BcpsTreeNode() 
        { init(); }
#endif
    
    /** Useful constructor. */
    PSCGTreeNode(PSCGModel* m) : BcpsTreeNode() {
        init();
        desc_ = new PSCGNodeDesc(m);
    }

    /** Useful constructor. */
    PSCGTreeNode(AlpsNodeDesc *&desc) : BcpsTreeNode() {
        init();
        desc_ = desc;
        desc = NULL;
	//solEstimate_=dynamic_cast<PSCGNodeDesc*>(desc_)->getLB();
	quality_=dynamic_cast<PSCGNodeDesc*>(desc_)->getLB();
    }

    /** Destructor. */
    virtual ~PSCGTreeNode() {
        //delete constraintPool_;
        //delete variablePool_;
    }
    
    /** Initilize member data when constructing a node. */
    void init() {
        //constraintPool_ = new BcpsConstraintPool;
        //variablePool_ = new BcpsVariablePool;
    }
    
    /** Create a new node based on given desc. */
    AlpsTreeNode* createNewTreeNode(AlpsNodeDesc *&desc) const;

    /** Convert explicit description to difference, and vise-vesa */
    ///@{
    virtual void convertToExplicit();
    virtual void convertToRelative();
    ///@}
    
    /** intall subproblem */
    virtual int installSubProblem(BcpsModel *model);
    
    /** Performing the bounding operation. */
    virtual int process(bool isRoot = false, bool rampUp = false);
    
    /** Bounding procedure */
    virtual int bound(BcpsModel *model);
    void updateQuality(double val){
	quality_=val;
	//if(val > quality_) quality_=val;
	//else{cout << "PSCGTreeNode::updateQuality(): Warning: candidate update of quality would be diminishing: " << val << " < " << quality_ << endl;}
    	PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);
        desc->updateBd(val);
    }
    //int postProcessBound(BcpsModel *model); 
    /** Takes the explicit description of the current active node and 
        creates the children's descriptions, which contain information 
        about how the branching is to be done. The stati of the children
        are AlpsNodeStatusCandidate. */
    virtual std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
	branch();
    
    virtual int chooseBranchingObject(BcpsModel *model){
	//int branchingIndex = -1;
    	PSCGModel *m = dynamic_cast<PSCGModel*>(model);
    	PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);
        m->findBranchingIndex();
        desc->updateBranchingIndex(m);
	
	vector< vector<var_branch> > newSPInfo = m->getNewNodeSPInfo();
	for(int nn=0; nn<newSPInfo.size(); nn++){
            PSCGNodeDesc *childDesc = new PSCGNodeDesc(m);//desc->createChildNodeDescUp(desc->getBranchingIndex());
	    for(int oo=0; oo<newSPInfo[nn].size(); oo++){
		childDesc->setZLB(newSPInfo[nn][oo].index, newSPInfo[nn][oo].lb);
		childDesc->setZUB(newSPInfo[nn][oo].index, newSPInfo[nn][oo].ub);
	    }
	    //childDesc->updateOmega(m);
	    childDesc->setLB(quality_);
    	    childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
					    (childDesc),
					    AlpsNodeStatusCandidate,
					    quality_));
	}
	return m->getInfeasIndex();
    }
    /** Select a branching object based on give branching strategy. */
#if 0
    int selectBranchObject(PSCGModel *model, 
                           bool& foundSol, 
                           int numPassesLeft);
#endif

    /** To be defined. */
    //virtual int chooseBranchingObject(BcpsModel*) { return AlpsReturnStatusOk;}
    
    /** Generate constraints. */
    int generateConstraints(PSCGModel *model, OsiCuts & cutPool);

    /** Select and apply constraints. */
    int applyConstraints(PSCGModel *model,
                         OsiCuts & cutPool,
                         const double *solution); 

    /** Fix and tighten varaibles based optimality conditions. */
    int reducedCostFix(PSCGModel *model);
    
    /** Return constraint pool. */
    //BcpsConstraintPool * constraintPool() { return constraintPool_; }

    /** Return variable pool. */
    //BcpsVariablePool * variablePool() { return variablePool_; }
    
    /** Encode this node for message passing. */
    virtual AlpsEncoded* encode() const;

    /** Decode a node from an encoded object. */
    virtual AlpsKnowledge* decode(AlpsEncoded&) const;

    //For testing/debugging purposes....
    void printInstallSubproblem(){
        PSCGModel* model = dynamic_cast<PSCGModel*>(desc_->getModel());
        PSCGNodeDesc* desc = dynamic_cast<PSCGNodeDesc*>(desc_);
	cout << "(" << model->getBound() << "," << desc->getLB() << ")" << endl;
	cout << "(" << model->getInfeasIndex() << "," << desc->getBranchingIndex() << ")" << endl;
    }
};
#endif
#endif
