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

#ifndef PSCGNodeDesc_h_
#define PSCGNodeDesc_h_

//#############################################################################

//#include "CoinWarmStartBasis.hpp"

#include "AlpsNodeDesc.h"
#include "BcpsNodeDesc.h"

//#include "PSCGHelp.h"
//class PSCGModel;
#include "PSCGModel.h"
#define CUTOFF 0.5
//#############################################################################


class PSCGNodeDesc : public BcpsNodeDesc {

 protected:
    
    double lbd_;
    int branchingIndex_;
    bool branchingIndexIsInt_;
    double branchingZVals[MAX_N_BRANCHINGS];
    int zLen_;
    //vector<double* > omega_;
    //vector<int> newBdInds_;
    //vector<int> newBdTypes_; //0 UP branch, setting LB; 1 DOWN branch, setting UB 
    //vector<double> newBds_;
    //int nNewBds_;
    double origPenalty_;
    //double(*zBounds)[2];
    double *zLBs;
    double *zUBs;

 public:
    /** Default constructor. */
    PSCGNodeDesc() : BcpsNodeDesc(),//newBdInds_(),newBdTypes_(),newBds_(),nNewBds_(0),
	lbd_(-ALPS_DBL_MAX),origPenalty_(1.0),branchingIndex_(-1),branchingIndexIsInt_(false),zLBs(NULL),zUBs(NULL){;}

    /** Useful constructor. */
    PSCGNodeDesc(PSCGModel *m) :
	BcpsNodeDesc(m),//newBdInds_(),newBdTypes_(),newBds_(),nNewBds_(0),
	lbd_(-ALPS_DBL_MAX),origPenalty_(1.0),branchingIndex_(-1),branchingIndexIsInt_(false)
    {
	    zLen_ = m->n1;
	    for(int tS=0; tS<m->nNodeSPs; tS++){
		//omega_.push_back(new double[m->n1]);
		//for(int ii=0;ii<m->n1; ii++) omega_[tS][ii]=0.0;
	    }
	origPenalty_ = m->getPenalty();
	//zBounds = new double[m->n1][2];
	zLBs = new double[m->n1];
	zUBs = new double[m->n1];
	setZBounds(m->getOrigVarLbds(),m->getOrigVarUbds());
    }
	


    /** Destructor. */
    virtual ~PSCGNodeDesc() { 
	//for(int ii=0; ii<omega_.size(); ii++){
	 //   if(omega_[ii]) delete [] omega_[ii];
	//}
	//delete [] zBounds;
	delete [] zLBs;
	delete [] zUBs;
	//if(newBdInds_) delete [] newBdInds_;
	//if(newBdTypes_) delete [] newBdTypes_;
	//if(newBds_) delete [] newBds_;
    }
    //vector<int>& getNewBdInds(){return newBdInds_;}
    //vector<int>& getNewBdTypes(){return newBdTypes_;}
    //vector<double>& getNewBds(){return newBds_;}
    //vector<double* >& getOmega(){return omega_;}
    double getLB(){return lbd_;}
    void setLB(double lb){lbd_=lb;}
    void setZBounds(const double *lbds, const double *ubds){
	for(int ii=0; ii<zLen_; ii++){
	    //zBounds[ii][0]=lbds[ii];
	    zLBs[ii]=lbds[ii];
	    //zBounds[ii][1]=ubds[ii];
	    zUBs[ii]=ubds[ii];
	}
    }
    void printZBounds(PSCGModel *m){
	if(m->getMPIRank()==0){
	  printf("\nPrinting z bounds for this node:\n[");
	  for (int ii = 0; ii < m->n1; ii++) {
		printf("(%0.2g,%0.2g) ", zLBs[ii], zUBs[ii]);
	  }
	  printf("]\n");
	}
    }
    //double* getZBounds(){return &zBounds[0][0];}
    double* getZLBs(){return zLBs;}
    double* getZUBs(){return zUBs;}
    void setZLBs(const double* lbs){memcpy(zLBs,lbs,zLen_*sizeof(double));}
    void setZUBs(const double* ubs){memcpy(zUBs,ubs,zLen_*sizeof(double));}
    int getBranchingIndex(){return branchingIndex_;}
    double getOrigPenalty(){return origPenalty_;}
    void setOrigPenalty(double p){origPenalty_=p;}
    void updateBranchingZVals(PSCGModel *m){
	memcpy(&branchingZVals[0], m->getBranchingZVals(), MAX_N_BRANCHINGS*sizeof(double));
    }
    void updateBd(PSCGModel *m){lbd_=m->getBound();}
    void updateBd(double val){lbd_=val;}
    void updateBranchingIndex(PSCGModel *m){branchingIndex_=m->getInfeasIndex();}
    void updateDesc(PSCGModel *m){
    	updateBranchingZVals(m);
    	//updateOmega(m); 
    	//updateBranchingIndex(m);
    }

    void assignNewBds(PSCGNodeDesc *parentNode){
	//vector<int> &bdInds = parentNode->getNewBdInds();
	//vector<int> &bdTypes = parentNode->getNewBdTypes();
	//vector<double> &bds = parentNode->getNewBds();
	//assignNewBds(bdInds,bdTypes,bds);
	//nNewBds_ = newBdInds_.size();
        //memcpy(zBounds,parentNode->getZBounds(),2*zLen_*sizeof(double));
        memcpy(zLBs,parentNode->getZLBs(),zLen_*sizeof(double));
        memcpy(zUBs,parentNode->getZUBs(),zLen_*sizeof(double));
    }
#if 0
    void assignNewBds(vector<int> &bdInds, vector<int> &bdTypes, vector<double> &bds){
	newBdInds_.clear();
	newBdTypes_.clear();
	newBds_.clear();

	newBdInds_=bdInds;
	newBdTypes_=bdTypes;
	newBds_=bds;
	nNewBds_ = newBdInds_.size();
    }
#endif
#if 0
    void addNewBd(int ind, int type, double bd){
	newBdInds_.push_back(ind);
	newBdTypes_.push_back(type);
	newBds_.push_back(bd);
	nNewBds_ = newBdInds_.size();
    }
#endif
    void addNewBd(int ind, int type, bool isInt, double branchingZVal){
	//newBdInds_.push_back(ind);
	//newBdTypes_.push_back(type);
	if(type==UP){ 
	    if(isInt) {
		//zBounds[ind][0]=ceil(branchingZVal);
		zLBs[ind]=ceil(branchingZVal);
	    }
	    else {
		//zBounds[ind][0]=branchingZVal;
		zLBs[ind]=branchingZVal;
	    }
  	    //newBds_.push_back(zLBs[ind]);
	}
	else if(type==DOWN){
	    if(isInt) {
		//zBounds[ind][1]=floor(branchingZVal);
		zUBs[ind]=floor(branchingZVal);
	    }
	    else {
		//zBounds[ind][1]=branchingZVal;
		zUBs[ind]=branchingZVal;
	    }
	    //newBds_.push_back(zUBs[ind]);
	}
	assert(type==UP || type==DOWN);
	//nNewBds_ = newBdInds_.size();
    }
    PSCGNodeDesc *createChildNodeDescUp(int ind){
	PSCGModel *model = dynamic_cast<PSCGModel*>(this->model_);
	PSCGNodeDesc *child = new PSCGNodeDesc(model);
	child->assignNewBds(this);
	//printZ(model);
	//child->updateBranchingZVals(model);
	//child->updateOmega(omega_,model);
	child->addNewBd(ind, UP, model->indexIsInt(ind),branchingZVals[0]);
	child->setOrigPenalty(origPenalty_);
	child->setLB(lbd_);
	//child->printZBounds(model);
	return child;
    }
    PSCGNodeDesc *createChildNodeDescDown(int ind){
	PSCGModel *model = dynamic_cast<PSCGModel*>(this->model_);
	PSCGNodeDesc *child = new PSCGNodeDesc(model);
	child->assignNewBds(this);
	//printZ(model);
	//child->updateBranchingZVals(model);
	//child->updateOmega(omega_,model);
	child->addNewBd(ind, DOWN, model->indexIsInt(ind),branchingZVals[0]);
	child->setOrigPenalty(origPenalty_);
	child->setLB(lbd_);
	//child->printZBounds(model);
	return child;
    }
    
    
 protected:

    /** Pack blis portion of node description into an encoded. */
    AlpsReturnStatus encodePSCG(AlpsEncoded *encoded) const {
	AlpsReturnStatus status = AlpsReturnStatusOk;
#if 0
	encoded->writeRep(branchedDir_);
	encoded->writeRep(branchedInd_);
	encoded->writeRep(branchedVal_);

	// Basis
	int ava = 0;
	if (basis_) {
	    ava = 1;
	    encoded->writeRep(ava);
	    PSCGEncodeWarmStart(encoded, basis_);
	}
	else {
	    encoded->writeRep(ava);
	}
#endif	
	return status;
    }

    /** Unpack blis portion of node description from an encoded. */
    AlpsReturnStatus decodePSCG(AlpsEncoded &encoded) {
	AlpsReturnStatus status = AlpsReturnStatusOk;
#if 0	
	encoded.readRep(branchedDir_);
	encoded.readRep(branchedInd_);
	encoded.readRep(branchedVal_);
	
	// Basis
	int ava;
	encoded.readRep(ava);
	if (ava == 1) {
	    basis_ = PSCGDecodeWarmStart(encoded, &status);
	}
	else {
	    basis_ = NULL;
	}
#endif	
	return status;
    }

 public:

    /** Pack node description into an encoded. */
    virtual AlpsReturnStatus encode(AlpsEncoded *encoded) const {
    	AlpsReturnStatus status = AlpsReturnStatusOk;
	
	status = encodeBcps(encoded);
	status = encodePSCG(encoded);
	
	return status;
    }

    /** Unpack a node description from an encoded. Fill member data. */
    virtual AlpsReturnStatus decode(AlpsEncoded &encoded) {
	
    	AlpsReturnStatus status = AlpsReturnStatusOk;
	
	status = decodeBcps(encoded);
	status = decodePSCG(encoded);

	return status;
    }
    void installSubproblemFromNodeDesc(){
	PSCGModel *model = dynamic_cast<PSCGModel*>(model_);
	//printZBounds(model);
	//model->installSubproblem(lbd_, newBdInds_, newBdTypes_, newBds_, nNewBds_, branchingIndex_, origPenalty_);
	model->installSubproblem(lbd_, zLBs, zUBs, branchingIndex_, origPenalty_);
    }
    
};
#endif
