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
    virtual int installSubProblem(BcpsModel *mode);
    
    /** Performing the bounding operation. */
    virtual int process(bool isRoot = false, bool rampUp = false);
    
    /** Bounding procedure */
    virtual int bound(BcpsModel *model);

    /** Takes the explicit description of the current active node and 
        creates the children's descriptions, which contain information 
        about how the branching is to be done. The stati of the children
        are AlpsNodeStatusCandidate. */
    virtual std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
	branch();
    
    /** Select a branching object based on give branching strategy. */
    int selectBranchObject(PSCGModel *model, 
                           bool& foundSol, 
                           int numPassesLeft);

    /** To be defined. */
    virtual int chooseBranchingObject(BcpsModel*) { return AlpsReturnStatusOk;}
    
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
};
#endif
#endif
