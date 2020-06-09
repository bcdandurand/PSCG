// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <string>
using namespace std;

#include <cassert>
#include <iostream>
#include "time.h"

#include "CoinPragma.hpp"
#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"
#include "PSCGScen.h"

#ifdef USING_MPI
   #include <mpi.h>
#endif

//#define DATASTOCHASTICDIR "/homes/bcdandurand/COIN-OR/Data/Stochastic"
#define DATASTOCHASTICDIR "/home/bdandurand/Data/SIPLIB"

#define MIP_TOL 1e-9

//forward declarations

vector<vector<double> > getColSolutionsByStageForScn(SmiScnModel &smi, int ns);

int main(int argc, char **argv) {

    int mpiRank;
    int mpiSize;
    bool mpiHead;
    bool parallel;
    
    time_t start_t, end_t;
   //identifying whether the code is going to run in  parallel or not
    #ifdef USING_MPI
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
        parallel = (mpiSize > 1);
    #else
        mpiRank = 0;
        mpiSize = 1;
        parallel = false;
    #endif
    //Flag for the head node recognise itself as such
    mpiHead = (mpiRank == 0);
    if(mpiHead){
        time(&start_t);
        cout << "Number of processors: " << mpiSize << endl;
    }

    SmiScnModel smi;
    char p_name[256] = DATASTOCHASTICDIR;
    //strcat(p_name,"/SSLP/sslp_5_25_50");
    strcat(p_name,"/SSLP/sslp_10_50_100");
    smi.readSmps(p_name);

    if(mpiHead){
        cout << "Model generation using SMPS files for "<< p_name << " problems." << endl ;
    }

    // generate OSI solver object
    OsiCpxSolverInterface *smi_osi = new OsiCpxSolverInterface();

    // set solver object for SmiScnModel
    smi.setOsiSolverHandle(*smi_osi);

    // load solver data
    //  this step generates the deterministic equivalent
    //  and returns an OsiSolver object
    OsiSolverInterface *smi_osi_loaded = smi.loadOsiSolverData();


    // solve
    //if(mpiHead){
    if(false){
        // set some nice Hints to the OSI solver
        smi_osi_loaded->setHintParam(OsiDoPresolveInInitial,true);
        smi_osi_loaded->setHintParam(OsiDoScale,true);
        smi_osi_loaded->setHintParam(OsiDoCrash,true);
        smi_osi_loaded->messageHandler()->setLogLevel(0); //Nice hack to suppress all output
        smi_osi_loaded->initialSolve();
        smi_osi_loaded->branchAndBound();

        // print results
        cout << "Solved stochastic program " << p_name << endl;;
        cout << "Number of rows: " << smi_osi_loaded->getNumRows() << endl;
        cout << "Number of cols: " << smi_osi_loaded->getNumCols() << endl;
        cout << "Optimal value: " << smi_osi_loaded->getObjValue() << endl;
    }

    int numScenarios=smi.getNumScenarios();
    OsiCpxSolverInterface *sp_osi;
    PSCGScen *sp_wrapper;
    int cpx_status=0;
    double optval_sum=0.0;
    double optval_sum_local=0.0;
    std::vector<OsiCpxSolverInterface*> osi_cpx_scns;
    std::vector<PSCGScen*> scn_sps;
    for (int ii=mpiRank ; ii<numScenarios; ii+=mpiSize) {
        //cout << "Extracting subproblem: " << ii << endl;
        sp_osi = new OsiCpxSolverInterface( ) ;
        osi_cpx_scns.push_back( sp_osi );
        sp_wrapper = new PSCGScen(smi,*smi_osi_loaded,*sp_osi,ii);
        scn_sps.push_back( sp_wrapper );
        //cout << "Done extracting subproblem: " << ii << endl;
        sp_osi->switchToMIP();
        sp_wrapper->enforceIntegrality();
        sp_osi->messageHandler()->setLogLevel(0); //Nice hack to suppress all output
	cpx_status=CPXsetintparam( sp_osi->getEnvironmentPtr(), CPXPARAM_Threads, 1);
	cpx_status=CPXsetintparam( sp_osi->getEnvironmentPtr(), CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC); //no iteration message until solution
        sp_osi->initialSolve();
        sp_osi->branchAndBound();
#if 0
        cout << "Scenario " << ii << " optimal value: " << sp_osi->getObjValue() << endl;
        cout << "Solution is: " << endl;
        for(int nn=0; nn<sp_osi->getNumCols(); nn++){
            cout << " " << (sp_osi->getColSolution())[nn];
        } 
        cout << endl;
#endif
        optval_sum_local += sp_osi->getObjValue();
    }
    #ifdef USING_MPI
      if(mpiSize>1){
          MPI_Allreduce(&optval_sum_local, &optval_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
    #endif
    if(mpiSize==1){
        optval_sum=optval_sum_local;
    }
    if(mpiHead){
        time(&end_t);
	double time_elapsed = end_t - start_t;
        cout << "Average optimal value: " << optval_sum << endl;
    	cout << "Elapsed time is " << time_elapsed << " seconds." << endl;
    	cout << "*** Done! *** " << endl;
    }

    return 0;
}


vector<vector<double> > getColSolutionsByStageForScn(SmiScnModel &smi, int scn){
    const double * osiSoln = (smi.getOsiSolverInterface())->getColSolution();
    int numcols=0;
    vector<vector<double> > solnsByStage;
    SmiScnNode *node = smi.getLeafNode(scn);
    while (node != NULL){
        solnsByStage.push_back(vector<double>());
        node = node->getParent();
    }


    int n_stages = solnsByStage.size();
    int stg = n_stages;
    node = smi.getLeafNode(scn);
    while (node != NULL){
        stg--;
        // copy entries
        // getColStart returns the starting index of node in OSI model
        for(int j=node->getColStart(); j<node->getColStart()+node->getNumCols(); ++j){
            // getCoreColIndex returns the corresponding Core index
            // in the original (user's) ordering
            solnsByStage[stg].push_back(osiSoln[j]);
        }
        // get parent of node
        node = node->getParent();
    }
    return solnsByStage;
}

