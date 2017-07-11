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
#include "PSCGModel.h"

//#############################################################################


class PSCGNodeDesc : public BcpsNodeDesc {

 protected:
    
    double lbd_;
    double *z_;
    double *omega_;
    int *newBdInds_;
    int *newBdTypes_; //0 UP branch, setting LB, 1 DOWN branch, setting UB 
    double *newBds_;
    int nNewBds_;
    
 public:

    /** Default constructor. */
    PSCGNodeDesc() :
	z_(NULL),omega_(NULL),newBdInds_(NULL),newBdTypes_(NULL),newBds_(NULL),nNewBds_(0),
        BcpsNodeDesc()
        {}

    /** Useful constructor. */
    PSCGNodeDesc(PSCGModel* m) :
	z_(NULL),omega_(NULL),newBdInds_(NULL),newBdTypes_(NULL),newBds_(NULL),nNewBds_(0),
	BcpsNodeDesc(m)
	{}

    /** Destructor. */
    virtual ~PSCGNodeDesc() { 
	delete [] z_;
	delete [] omega_;
	delete [] newBdInds_;
	delete [] newBdTypes_;
	delete [] newBds_;
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
	dynamic_cast<PSCGModel*>(model_)->installSubproblem(z_, omega_, newBdInds_, newBdTypes_, newBds_, nNewBds_);
    }
    
};
#endif
