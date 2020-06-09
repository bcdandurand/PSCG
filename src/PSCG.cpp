// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <string>
using namespace std;

#include <cassert>
#include <iostream>

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

void testingMessage(const char * const);
void SmpsIO(const char* const);
vector<vector<double> > getColSolutionsByStageForScn(SmiScnModel &smi, int ns);
vector<vector<int> > extractScnSPColsByStage(SmiScnModel &smi, int ns);
vector<vector<int> > extractScnSPRowsByStage(SmiScnModel &smi, int ns);
vector<double> extractScnProbTrace(SmiScnModel &smi, int ns);
void extractScnSubmatrix(OsiSolverInterface *smiOsi, CoinPackedMatrix &spMat, vector<vector<int> > &stg_cols, vector<vector<int> > &stg_rows);
void constructScnSubmatrix(OsiSolverInterface *smiOsi, CoinPackedMatrix &spMat, vector<vector<int> > &stg_cols, vector<vector<int> > &stg_rows);
vector<int> extractScnSP(SmiScnModel &smi, OsiSolverInterface *smiOsi, int ns, OsiSolverInterface *osi);

int main()
{

    testingMessage( "Model generation using SMPS files for Cambridge-Watson problems.\n" );
    SmpsIO(DATASTOCHASTICDIR"/SSLP/sslp_5_25_50");

    testingMessage( "*** Done! *** \n");

    return 0;
}

void SmpsIO(const char * const name )
{
        SmiScnModel smi;

        // read SMPS model from files
        //  <name>.core, <name>.time, and <name>.stoch
        smi.readSmps(name);

        // generate OSI solver object
        //  here we use OsiClp
        OsiCpxSolverInterface *clp = new OsiCpxSolverInterface();

        // set solver object for SmiScnModel
        smi.setOsiSolverHandle(*clp);

        // load solver data
        //  this step generates the deterministic equivalent
        //  and returns an OsiSolver object
        OsiSolverInterface *osiStoch = smi.loadOsiSolverData();

        // set some nice Hints to the OSI solver
        osiStoch->setHintParam(OsiDoPresolveInInitial,true);
        osiStoch->setHintParam(OsiDoScale,true);
        osiStoch->setHintParam(OsiDoCrash,true);

        // solve
#if 1
        osiStoch->initialSolve();
        osiStoch->branchAndBound();

        // print results
        printf("Solved stochastic program %s\n", name);
        printf("Number of rows: %d\n",osiStoch->getNumRows());
        printf("Number of cols: %d\n",osiStoch->getNumCols());
        printf("Optimal value: %g\n",osiStoch->getObjValue());
#endif

        // print solution to file
        int numScenarios=smi.getNumScenarios();
#if 0
        for (int i=0 ; i<numScenarios; ++i) {
            vector< vector<double> > solnsByStage = getColSolutionsByStageForScn(smi, i);
            cout << "Scenario " << i << endl;
            for(size_t stg=0; stg<solnsByStage.size(); stg++){
                cout << "\tStage " << stg + 1 << " soln: ";
                for(size_t jj=0; jj<solnsByStage[stg].size(); jj++){
                    cout << "  " << solnsByStage[stg][jj];
                }
                cout << endl;
            }
        }
#endif

        double optval_sum=0.0;
        std::vector<OsiCpxSolverInterface*> osi_cpx_scns;
        std::vector<PSCGScen*> scn_sps;
        for (int ii=0 ; ii<numScenarios; ++ii) {
            cout << "Extracting subproblem: " << ii << endl;
            osi_cpx_scns.push_back( new OsiCpxSolverInterface( ) );
            scn_sps.push_back( new PSCGScen(smi,*osiStoch,*osi_cpx_scns[ii],ii) );
            cout << "Done extracting subproblem: " << ii << endl;
            osi_cpx_scns[ii]->switchToMIP();
            //for (unsigned int nn = 0; nn < smi.getIntegerLen(); nn++) {
	    scn_sps[ii]->enforceIntegrality();
#if 0
            if( ii==0 ){
                osi_cpx_scns[ii]->writeLp("extract_model.out");
            }
#endif
            osi_cpx_scns[ii]->messageHandler()->setLogLevel(0); //Nice hack to suppress all output
            osi_cpx_scns[ii]->initialSolve();
            osi_cpx_scns[ii]->branchAndBound();
            printf("Scenario %d optimal value: %g\n",ii,osi_cpx_scns[ii]->getObjValue());
            cout << "Solution is: " << endl;
            for(int nn=0; nn<osi_cpx_scns[ii]->getNumCols(); nn++){
                cout << " " << (osi_cpx_scns[ii]->getColSolution())[nn];
            } 
            cout << endl;
            optval_sum += osi_cpx_scns[ii]->getObjValue();
        }
        cout << "Average optimal value: " << optval_sum << endl;


}


// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
//  std::cerr <<msg;
  cout <<endl <<"*****************************************"
       <<endl <<msg <<endl;
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


    size_t n_stages = solnsByStage.size();
    size_t stg = n_stages;
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

