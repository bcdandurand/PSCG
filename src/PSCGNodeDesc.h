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

//#############################################################################


class PSCGNodeDesc : public BcpsNodeDesc {

 protected:
    
    double lbd_;
    double *z_;
    int zLen_;
    vector<double* > omega_;
    vector<int> newBdInds_;
    vector<int> newBdTypes_; //0 UP branch, setting LB; 1 DOWN branch, setting UB 
    vector<double> newBds_;
    int nNewBds_;
 public:
    /** Default constructor. */
    PSCGNodeDesc() : BcpsNodeDesc(),z_(NULL),newBdInds_(),newBdTypes_(),newBds_(),nNewBds_(0){;}

    /** Useful constructor. */
    PSCGNodeDesc(PSCGModel *m) :
	BcpsNodeDesc(m),z_(NULL),newBdInds_(),newBdTypes_(),newBds_(),nNewBds_(0)
	{
	    zLen_ = m->n1;
	    z_ = new double[m->n1];
	    for(int tS=0; tS<m->nNodeSPs; tS++){
		omega_.push_back(new double[m->n1]);
		for(int ii=0;ii<m->n1; ii++) omega_[tS][ii]=0.0;
	    }
	}


    /** Destructor. */
    virtual ~PSCGNodeDesc() { 
	if(z_) delete [] z_;
	for(int ii=0; ii<omega_.size(); ii++){
	    if(omega_[ii]) delete [] omega_[ii];
	}
	//if(newBdInds_) delete [] newBdInds_;
	//if(newBdTypes_) delete [] newBdTypes_;
	//if(newBds_) delete [] newBds_;
    }
    vector<int>& getNewBdInds(){return newBdInds_;}
    vector<int>& getNewBdTypes(){return newBdTypes_;}
    vector<double>& getNewBds(){return newBds_;}
    void updateZ(PSCGModel *m){
	memcpy(z_,m->getZ(),(m->n1)*sizeof(double));
    }
    void updateOmega(PSCGModel *m){
	vector<double*> &omega = m->getOmega();
	for(int tS=0; tS<m->nNodeSPs; tS++){
	    memcpy(omega_[tS],omega[tS],(m->n1)*sizeof(double));
	}
    }
    void updateBd(PSCGModel *m){lbd_=m->getBound();}

    void assignNewBds(PSCGNodeDesc *parentNode){
	vector<int> &bdInds = parentNode->getNewBdInds();
	vector<int> &bdTypes = parentNode->getNewBdTypes();
	vector<double> &bds = parentNode->getNewBds();
	assignNewBds(bdInds,bdTypes,bds);
    }
    void assignNewBds(vector<int> &bdInds, vector<int> &bdTypes, vector<double> &bds){
	newBdInds_.clear();
	newBdTypes_.clear();
	newBds_.clear();

	newBdInds_=bdInds;
	newBdTypes_=bdTypes;
	newBds_=bds;
	nNewBds_ = newBdInds_.size();
    }
    void addNewBd(int ind, int type, double bd){
	newBdInds_.push_back(ind);
	newBdTypes_.push_back(type);
	newBds_.push_back(bd);
	nNewBds_++;
    }
    void addNewBd(int ind, int type){
	newBdInds_.push_back(ind);
	newBdTypes_.push_back(type);
	if(type==UP){ newBds_.push_back((int)ceil(z_[ind]));}
	else{newBds_.push_back((int)floor(z_[ind]));}
	nNewBds_++;
    }
    PSCGNodeDesc *createChildNodeDescUp(int ind){
	PSCGModel *model = dynamic_cast<PSCGModel*>(this->model_);
	PSCGNodeDesc *child = new PSCGNodeDesc(model);
	child->assignNewBds(this);
	child->updateZ(model);
	child->updateOmega(model);
	child->addNewBd(ind, UP);
	return child;
    }
    PSCGNodeDesc *createChildNodeDescDown(int ind){
	PSCGModel *model = dynamic_cast<PSCGModel*>(this->model_);
	PSCGNodeDesc *child = new PSCGNodeDesc(model);
	child->assignNewBds(this);
	child->updateZ(model);
	child->updateOmega(model);
	child->addNewBd(ind, DOWN);
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
	dynamic_cast<PSCGModel*>(model_)->installSubproblem(z_, omega_, lbd_, newBdInds_, newBdTypes_, newBds_, nNewBds_);
    }
    
};
#endif
