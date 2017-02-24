/*Header for the main procedure ParallelSCG. */

#ifndef PSCGMODEL_H
#define PSCGMODEL_H

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <sys/time.h>
#include "tclap/CmdLine.h"
#include "StructureDefs.h"
#include "TssModel.h"
#include "PSCGModelScen.h"
#include "PSCGParams.h"
#include "BcpsModel.h"

#define smallNumber 0.000001
#define SSC_PARAM 0.1
#define SSC_DEN_TOL 1e-10
#define ROUNDING_TOLERANCE 1e-4
#define DEFAULT_THREADS 1
#define KIWIEL_PENALTY 0 //set 1 to use Kiwiel (2006) penalty update rule

#ifdef USING_MPI
   #include <mpi.h>
#endif

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#define MAX_INNER_LOOP 200
#define DEFAULT_PENALTY 100

// Parameters for newConvergenceCriterion

#define BG_BETA 0.05

enum COL_TYPE{
CONTINUOUS=0,
BINARY,
INTEGER
};

// If you want to add a flag, search ADDFLAG in this file and the corresponding .cpp file and follow instructions.

// ADDFLAG : Add the variable here

class PSCGModel{
public:
PSCGParams *par;
TssModel smpsModel;
vector<int> scenariosToThisModel;
vector<PSCGModelScen*> subproblemSolvers;
int nS;
int n1;
int n2;
int nNodeSPs;
vector<double> pr;
vector<double*> scaling_matrix; //glorified rho
vector<double*> omega_tilde; //dual variable
vector<double*> omega_current;

vector<double*> x_current;
vector<double*> y_current;

double* z_current;// = new double[n1];
double* z_local;// = new double[n1];
double* totalSoln_; //new double[n1+n2];

//initialising zˆ0=0
#if 0
for (int i = 0; i < n1; i++) {
    z_current[i] = 0;
}
#endif

char filepath[64];
double LagrLB_Local;
// Bounds
double LagrLB;
double currentLagrLB;
double ALVal;

// Norms
double discrepNorm;

double localReduceBuffer[3]; //0-LagrLB_Local,  1-ALVal_Local,  2-localDiscrepNorm
double reduceBuffer[3];	//0-LagrLB,  1-ALVal,  2-discrepNorm
double ALVal_Local;
double localDiscrepNorm;
double penC;
double penMult;
int maxStep;
int maxSeconds;
int fixInnerStep;
int nVerticesUsed;
int nThreads;
bool verbose;
bool debug;
bool linRelaxFirst;
bool linRelaxSecond;
bool scaling;
bool dataPathOverride;
bool LBcalc;
bool AlgorithmC;
bool disableHeuristic;
bool useVertexHistory;
int ftype;
int mpiRank;
int mpiSize;
bool mpiHead;

int numIntVars_;
int *intVar_;
char *colType_;

// Termination conditions
int totalNoGSSteps;

ProblemDataBodur pdBodur;

PSCGModel(PSCGParams *p):par(p),nNodeSPs(0),LagrLB_Local(0.0),ALVal_Local(0.0),localDiscrepNorm(1e9),discrepNorm(1e9),
	mpiRank(0),mpiSize(1),mpiHead(true),totalNoGSSteps(0){

   	//******************Read Command Line Parameters**********************
	//Params par;
	initialiseParameters();

	//******************Display Command Line Parameters**********************
	if (mpiRank==0) {
		displayParameters();
	}

	//******************Reading Problem Data**********************

	initialiseModel();
	

	if (mpiRank==0) {
		std::cout << std::endl;
		std::cout << "Problem data: " << std::endl;
		std::cout << "Number of first stage variables: " << n1 << std::endl;
		std::cout << "Number of second stage variables: " << n2 << std::endl;
		std::cout << "Number of scenarios: " << nS << std::endl;
    	}

	//******************Assign Scenarios**********************

	assignSubproblems();
	setupSolvers();
}

~PSCGModel(){
	//Clean up structures
	for (int tS = 0; tS < nNodeSPs; tS++) {
		delete [] scaling_matrix[tS];
		delete [] omega_tilde[tS];
		delete [] omega_current[tS];
		delete subproblemSolvers[tS];
	}
	delete[] z_current;
	delete[] z_local;
	delete [] totalSoln_;
}
void setMPIParams(int rank, int size){
	mpiRank=rank;
	mpiSize=size;
	mpiHead=(mpiRank==0);
}
void assignSubproblems(){
	for (int tS = 0; tS < nS; tS++) {
		//scenarioAssign[tS] = tS % mpiSize;
		if( (tS % mpiSize)==mpiRank ){
		    scenariosToThisModel.push_back(tS);
		}
	}
}


void initialiseParameters();


void initialiseModel(){
	strcpy(filepath,par->filename.c_str());
	int j=0;
	switch( ftype ){
	  case 1:
	    pdBodur.initialise(par);
	    n1 = pdBodur.get_n1();
	    n2 = pdBodur.get_n2();
	    nS = pdBodur.get_nS();
	    totalSoln_=new double[n1+n2];
	    break;
	  case 2:
	    smpsModel.readSmps(par->filename.c_str());
	    n1 = smpsModel.getNumCols(0);
	    n2 = smpsModel.getNumCols(1);
	    totalSoln_=new double[n1+n2];
	    nS = smpsModel.getNumScenarios();
		//Is there something wrong with StoModel::getNumIntegers(0) ???
	    //numIntVars_=smpsModel.getNumIntegers(0); //first-stage only
	    colType_ = new char[n1];
	    intVar_ = new int[n1];
	    numIntVars_=0;
	    memcpy(colType_,smpsModel.getCtypeCore(0),n1*sizeof(char));
	    for(int i=0; i<n1; i++){
//cout << " " << colType_[i];
		if(colType_[i]=='I' || colType_[i]=='B'){
		    intVar_[numIntVars_++]=i;
		}
	    }
	    //assert(j==numIntVars_);
	    break;
	  default:
	    throw(-1);
	    break;
	}
	z_current = new double[n1];
	z_local = new double[n1];
	//initialising zˆ0=0
	for (int i = 0; i < n1; i++) {
    	    z_current[i] = 0;
	}
}

void setupSolvers();

void upBranchAllSPsAt(int index, double bound){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->upBranchOnVar(index, bound);
    }
}

void downBranchAllSPsAt(int index, double bound){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->downBranchOnVar(index, bound);
    }
}

void initialIteration();

void regularIteration(){
	solveContinuousMPs();

	performColGenStep();
	double SSCVal = computeSSCVal();
	updateOmega(SSCVal);
	#if KIWIEL_PENALTY 
	 computeKiwielPenaltyUpdate(SSCVal);
	#endif
}

void solveContinuousMPs(){
	for(int itGS=0; itGS < fixInnerStep; itGS++) { //The inner loop has a fixed number of occurences
    	    for (int tS = 0; tS < nNodeSPs; tS++) {

		//*************************** Quadratic subproblem ***********************************
			    
		    if (useVertexHistory) { //Compute Next X, Y (With History)

					
			//Solve the quadratic master problem for x and y
			subproblemSolvers[tS]->updatePrimalVariablesHistory_OneScenario(omega_current[tS],z_current);
			
			// note: the final weight corresponds to the existing x
		    }
		    else { //might not work correctly, untested! Compute Next X, Y (Without History)

			subproblemSolvers[tS]->updatePrimalVariables_OneScenario(omega_current[tS],z_current,scaling_matrix[tS]);
		    }
	    }
	    					
	    // Update z_previous.
	    updateZ();
	    totalNoGSSteps++;
	}
}

void performColGenStep();


void updateZ(){
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += x_current[tS][i] * pr[tS];
	    }
	}
	#ifdef USING_MPI
		MPI_Allreduce(z_local, z_current, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) z_current[i] = z_local[i]; //Only one node, trivially set z_local = z_current
	#endif
}

int numIntInfeas()
{
#if 1
    int numIntegerInfs = 0;
    //int i = -1;
//printZ(); 
    for (int i = 0; i < numIntVars_; ++i) {
      if ( ! checkInteger(z_current[intVar_[i]]) ) {
//cout << "numIntVars_ " << numIntVars_ << " i = " << i << " index " << intVar_[i] << endl;
	++numIntegerInfs;
      }
    }
#endif
    return numIntegerInfs;
}

bool verifyIntegrality()
{
    return numIntInfeas()==0;
}

/** Check if a value is integer. */
bool checkInteger(double value) const {
    double integerTolerance = 1.0e-5;
    double nearest = floor(value + 0.5);
//cout << " " << fabs(value-nearest);
    if (fabs(value - nearest) <= integerTolerance) {
        return true;
    }
    else {
        return false;
    }
}

bool checkZIsFeasForScen(int tS){
    const CoinPackedMatrix *mat = dynamic_cast<PSCGModelScen_SMPS*>(subproblemSolvers[tS])->getOSI()->getMatrixByCol();
    memcpy(totalSoln_,z_current,n1*sizeof(double));
    double *ySoln = subproblemSolvers[tS]->getYVertex();
    memcpy(totalSoln_+n1,ySoln,n2*sizeof(double));
    
    return false;
}
bool checkZHasFullRecourse(){
    return false;
}


//********************** Serious Step Condition (SSC) **********************
bool shouldTerminate(){
	return (ALVal + 0.5*discrepNorm  - currentLagrLB < SSC_DEN_TOL);
}
double computeSSCVal(){
     return (LagrLB-currentLagrLB)/(ALVal + 0.5*discrepNorm - currentLagrLB);
}
		
void updateOmega(double SSCVal){
	if(SSCVal >= SSC_PARAM) {
	    for (int tS = 0; tS < nNodeSPs; tS++) {
		for (int i = 0; i < n1; i++) {
	    	    omega_current[tS][i] = omega_tilde[tS][i];
		}
	    }
	    currentLagrLB = LagrLB;
	}
	else {
	    if(mpiRank==0) cout << "Null step taken." << endl;
	}
}

		
//********************** Penalty update **********************
//This should normally be turned off, this is a rule for updating the 
//penalty from Kiwiel 2006 and Lubin et al.
void computeKiwielPenaltyUpdate(double SSCVal){	
	penC = 1.0/min( max( (2.0/penC)*(1.0-SSCVal),  max(1.0/(10.0*penC),1e-4)    ), 10.0/penC);
	for (int tS = 0; tS < nNodeSPs; tS++) {
	    subproblemSolvers[tS]->setQuadraticTerm(penC);
    	    for (int i = 0; i < n1; i++) {
		scaling_matrix[tS][i] = penC;
    	    }
	}
if(mpiRank==0) cout << "Penalty is now: " << penC << endl;
}


void displayParameters();
void printStatus(){
	printf("Lagrangian Lower Bound: %0.9g\n", currentLagrLB);
	printf("Aug. Lagrangian value: %0.9g\n", ALVal);
	printf("Norm of primal discrepancy: %0.6g\n", discrepNorm);
	std::cout << "Number of integrality infeasibilities of z: " << numIntInfeas() << std::endl;
}
void printZ(){
	printf("\nConsensus values for first-stage decision variables:\n[");
	
	for (int i = 0; i < n1; i++) {
		printf("%0.2g ", z_current[i]);
	}

	printf("]\n");
}

//*** Helper functions.
//Frobenius norm
double froNorm(double** array, int nS, int n1);
double froNorm(double* array, int n1);
//Deals with numerical imprecision from solver (MIPs)
double roundIfClose(double input) {
	double tmp = round(input);
	if (abs(tmp-input) < ROUNDING_TOLERANCE)
		{ return tmp; }
	else
		{ return input; }

}
//*** FWMM functions.
double computeNormedDiscrepancies(double** scaling_matrix, double** x, double* z, int nS, int n1);
//void weightedAverage(const std::vector<double*> &x, const std::vector<double> &p, double *localZ, double* z, int nNodeSPs, int n1, int mpiRank);
void computeFeasibleStartingPoint(int tS, double* x, double* yFeasible);
};
#endif
