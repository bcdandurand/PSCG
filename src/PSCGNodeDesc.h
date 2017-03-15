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
    int *newBdInds_;
    int *newBdTypes_; //0 UP branch, setting LB; 1 DOWN branch, setting UB 
    double *newBds_;
    int nNewBds_;
 public:
#if 1
    /** Default constructor. */
    PSCGNodeDesc() : BcpsNodeDesc(),z_(NULL),newBdInds_(NULL),newBdTypes_(NULL),newBds_(NULL),nNewBds_(0){;}

    /** Useful constructor. */
    PSCGNodeDesc(PSCGModel *m) :
	BcpsNodeDesc(m),z_(NULL),newBdInds_(NULL),newBdTypes_(NULL),newBds_(NULL),nNewBds_(0)
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
	if(newBdInds_) delete [] newBdInds_;
	if(newBdTypes_) delete [] newBdTypes_;
	if(newBds_) delete [] newBds_;
    }
    
#endif
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
	dynamic_cast<PSCGModel*>(model_)->installSubproblem(z_, omega_, newBdInds_, newBdTypes_, newBds_, nNewBds_);
    }
    
};
#endif
