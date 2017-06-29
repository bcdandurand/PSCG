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
#include "Alps.h"
#include "AlpsModel.h"
#include "AlpsKnowledge.h"
#include "BcpsSolution.h"
#include "BcpsModel.h"
#include "AlpsTreeNode.h"

//#define smallNumber 0.000001
#define SSC_PARAM 0.1
#define SSC_DEN_TOL 1e-10
//#define ROUNDING_TOLERANCE 1e-6
#define DEFAULT_THREADS 1
#define KIWIEL_PENALTY 1 //set 1 to use Kiwiel (2006) penalty update rule
//#define MAX_N_BRANCHINGS 1
//#define BRANCHING_EPS 1e-2
#define IDEPS 1e-10

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

enum BRANCH_TYPE{
UP=0,
DOWN
};

typedef struct{
  int index;
  double lb;
  double ub;
} var_branch;

#if 0
enum SolverReturnStatus {
  PSCG_OK=0,
  PSCG_OPTIMAL=0,
  PSCG_ABANDONED,
  PSCG_PRIMAL_INF_INT, //z does not satisfy integrality
  PSCG_PRIMAL_INF_REC, //z does not have recourse
  PSCG_PRIMAL_INF, //z both does not satisfy integrality and does not have recourse
  PSCG_DUAL_INF,
  PSCG_PRIMAL_LIM,
  PSCG_DUAL_LIM,
  PSCG_ITER_LIM,
  PSCG_ERR=100,
  PSCG_INF=200,
  PSCG_UNBOUND=201,
  PSCG_UNKNOWN=202
};
#endif

// If you want to add a flag, search ADDFLAG in this file and the corresponding .cpp file and follow instructions.

// ADDFLAG : Add the variable here

class PSCGModel : public BcpsModel{
public:
PSCGParams *par;
TssModel smpsModel;
int algorithm;
vector<int> scenariosToThisModel;
vector<PSCGModelScen*> subproblemSolvers;
vector< vector<var_branch> > newNodeSPInfo;
int nS;
int n1;
int n2;
int nNodeSPs;
vector<double> pr;
vector<double*> scaling_matrix; //glorified rho
vector<double*> omega_tilde; //dual variable
vector<double*> omega_current;
vector<double*> omega_saved;
bool omegaIsZero_;
bool omegaUpdated_;
bool shouldTerminate;
bool shouldFathomByOpt;

double *scaleVec_;

vector<double*> x_current;
vector<double*> y_current;

double* z_current;// = new double[n1];
double* z_old;// = new double[n1];
double* z_average;
double* z_average1;
double* z_average2;
double* z_local;// = new double[n1];
double* z_incumbent_; //this should be the last feasible z with best obj
double *z_rounded;
bool selectOnlyIntegerVars;
double *z_saved;
double* totalSoln_; //new double[n1+n2];
vector<double> constrVec_; //This is filled with the value of Ax

vector<int> spSolverStatuses_;

int modelStatus_[2];

char filepath[64];
char probname[64];
double LagrLB_Local;
// Bounds
double LagrLB;
double currentLagrLB;
double referenceLagrLB;
double ALVal;
//double incumbentVal_;
double objVal;

// Norms
double discrepNorm;

double(*recordKeeping)[4];
double localReduceBuffer[3]; //0-LagrLB_Local,  1-ALVal_Local,  2-localDiscrepNorm
double reduceBuffer[3];	//0-LagrLB,  1-ALVal,  2-discrepNorm
double ALVal_Local;
double localDiscrepNorm;
double penC;
double baselinePenalty_;
double penMult;
int maxStep;
int maxNoSteps;
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
bool AlgorithmZ;
bool disableHeuristic;
bool useVertexHistory;
int ftype;
int mpiRank;
int mpiSize;
bool mpiHead;

int numIntVars_;
int *intVar_;
char *colType_;
int infeasIndex_;
int nIntInfeas_;

double *origVarLB_;
double *origVarUB_;
double *currentVarLB_;
double *currentVarUB_;

/** Incumbent objective value. */
double incObjValue_;
/** Incumbent */

// Termination conditions
int totalNoGSSteps;
char zOptFile[128];

ProblemDataBodur pdBodur;

PSCGModel(PSCGParams *p):par(p),nNodeSPs(0),algorithm(ALGDD),referenceLagrLB(-ALPS_DBL_MAX),currentLagrLB(-ALPS_DBL_MAX),LagrLB(-ALPS_DBL_MAX),
LagrLB_Local(0.0),ALVal_Local(ALPS_DBL_MAX),ALVal(ALPS_DBL_MAX),objVal(ALPS_DBL_MAX),localDiscrepNorm(1e9),discrepNorm(1e9),
	mpiRank(0),mpiSize(1),mpiHead(true),totalNoGSSteps(0),infeasIndex_(-1),selectOnlyIntegerVars(false),maxNoSteps(20),
	nIntInfeas_(-1),omegaUpdated_(false),scaleVec_(NULL),BcpsModel(){

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
	sprintf(zOptFile,"zOpt-%s-noP%dalg%d",probname,mpiSize,algorithm);
	
}

~PSCGModel(){
	//Clean up structures
	for (int tS = 0; tS < nNodeSPs; tS++) {
		delete [] scaling_matrix[tS];
		delete [] omega_tilde[tS];
		delete [] omega_current[tS];
		delete [] omega_saved[tS];
		delete subproblemSolvers[tS];
	}
	//delete[] z_current;
	delete [] z_current;
	delete [] z_old;
	delete [] z_saved;
	delete [] z_average;
	delete [] z_average1;
	delete [] z_average2;
	delete [] z_local;
	delete [] z_incumbent_;
	delete [] z_rounded;
	delete [] origVarLB_;
	delete [] origVarUB_;
	delete [] currentVarLB_;
	delete [] currentVarUB_;
	delete [] totalSoln_;
	delete [] scaleVec_;
	delete [] recordKeeping;
	delete [] colType_;
	delete [] intVar_;
}

void setMPIParams(int rank, int size){
	mpiRank=rank;
	mpiSize=size;
	mpiHead=(mpiRank==0);
}
int getMPIRank(){return mpiRank;}

void assignSubproblems(){
	for (int tS = 0; tS < nS; tS++) {
		//scenarioAssign[tS] = tS % mpiSize;
		if( (tS % mpiSize)==mpiRank ){
		    scenariosToThisModel.push_back(tS);
		    spSolverStatuses_.push_back(SP_UNKNOWN);
		}
	}
}


void initialiseParameters();
void setAlgorithm(int algNo){
   algorithm=algNo;
}

int getAlgorithm(){return algorithm;}
void printAlgorithm(){
  if(mpiRank==0){
    switch(algorithm){
	case ALGDDBASIC:
	    cout << "Using baseline DDBB" << endl;
	    break;
	case ALGDD:
	    cout << "Using slightly enhanced DDBB" << endl;
	    break;
	case ALGDDD:
	    cout << "Using enhanced DDBB" << endl;
	    break;
	case ALGZ:
	    cout << "Using Z BB" << endl;
	    break;
    }
  }
}

void initialiseFileName(){
char str[] ="- This, a sample string.";
  char * pch, *prb;
  //printf ("Splitting string \"%s\" into tokens:\n",str);
  pch = strtok (filepath,"/");
  while (pch != NULL)
  {
    strcpy(probname,pch);
    //printf ("%s\n",pch);
    pch = strtok (NULL, "/");
  }
}

void initialiseModel(){



	int j=0;
	switch( ftype ){
	  case 1:
	    pdBodur.initialise(par);
	    n1 = pdBodur.get_n1();
	    n2 = pdBodur.get_n2();
	    nS = pdBodur.get_nS();
	    totalSoln_=new double[n1+n2];
	    colType_ = new char[n1];
	    intVar_ = new int[n1];
	    numIntVars_=0;
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
	//z_current = new double[n1];
	z_current = new double[n1];
	z_old = new double[n1];
	z_saved = new double[n1];
	z_average = new double[n1];
	z_average1 = new double[n1];
	z_average2 = new double[n1];
	z_local = new double[n1];
	z_incumbent_ = new double[n1];
	z_rounded = new double[n1];
	origVarLB_ = new double[n1];
	origVarUB_ = new double[n1];
	currentVarLB_ = new double[n1];
	currentVarUB_ = new double[n1];
	//initialising zˆ0=0
#if 0
	for (int i = 0; i < n1; i++) {
    	    z_current[i] = 0;
	}
#endif
	switch( ftype ){
	  case 1:
	    //TODO: initialise origVarBds.
	    break;
	  case 2:
	    smpsModel.copyCoreColLower(origVarLB_,0);
	    smpsModel.copyCoreColUpper(origVarUB_,0);
	    memcpy(currentVarLB_,origVarLB_,n1*sizeof(double));
	    memcpy(currentVarUB_,origVarUB_,n1*sizeof(double));
	    break;
	  default:
	    throw(-1);
	    break;
	}
}

void setupSolvers();

void fixVarAllSPsAt(int index, double fixVal){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->fixVarAt(index, fixVal);
    }
    currentVarLB_[index]=fixVal;		
    currentVarUB_[index]=fixVal;		
}

void upBranchAllSPsAt(int index, double bound){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->upBranchOnVar(index, bound);
    }
    currentVarLB_[index]=bound;		
}
void setLBsAllSPs(const double *lbs){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->setLBs(lbs, n1);
    }
    memcpy(currentVarLB_,lbs,n1*sizeof(double));		
}
void setUBsAllSPs(const double *ubs){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->setUBs(ubs, n1);
    }
    memcpy(currentVarUB_,ubs,n1*sizeof(double));		
}

void downBranchAllSPsAt(int index, double bound){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->downBranchOnVar(index, bound);
    }
    currentVarUB_[index]=bound;		
}

void restoreOriginalVarBounds(){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->unfixX(origVarLB_,origVarUB_);
    }
    memcpy(currentVarLB_,origVarLB_,n1*sizeof(double));
    memcpy(currentVarUB_,origVarUB_,n1*sizeof(double));
}
void clearSPVertexHistory(){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->clearVertexHistory();
    }
}
#if 0
void installSubproblem(double lb, vector<double*> & omega, vector<int> &indices, vector<int> &fixTypes, vector<double> &bds, int nFix, int infeasIndex, double pen){
//if(mpiRank==0){cerr << "Begin installSubproblem ";}
    restoreOriginalVarBounds();
//if(mpiRank==0){cout << "Branching at (indices,bounds,type): ";}
    for(int ii=0; ii<nFix; ii++){
//if(mpiRank==0){cout << "(" << indices[ii] << "," << bds[ii] << "," << fixTypes[ii] << ") ";}
	switch( fixTypes[ii] ){
	    case UP:
		//up-branch
		upBranchAllSPsAt(indices[ii],bds[ii]);
		break;
	    case DOWN:
		//down-branch
		downBranchAllSPsAt(indices[ii],bds[ii]);
		break;
	}
    }
//if(mpiRank==0){cout << endl;}
    //TODO: resolve ownership questions for current_z, current_omega
//if(mpiRank==0){cerr << "Reading z...";}
    //readZIntoModel(z);
//if(mpiRank==0){cerr << "Reading omega...";}
    readOmegaIntoModel(omega);
    //loadOmega();
    //zeroOmega();
    currentLagrLB=lb;
    infeasIndex_=infeasIndex;
    setPenalty(pen);
    baselinePenalty_=pen;
    printOriginalVarBds();
    printCurrentVarBds();
    //subproblemSolvers[0]->printXBounds();
//if(mpiRank==0){cerr << "End installSubproblem ";}
}
#endif


void installSubproblem(double lb, vector<double*> &omega, const double *zLBs, const double *zUBs, double pen){
//if(mpiRank==0){cerr << "Begin installSubproblem ";}
    for(int ii=0; ii<n1; ii++){
	assert(zLBs[ii] <= zUBs[ii]);
    }
    setLBsAllSPs(zLBs);
    setUBsAllSPs(zUBs);
//if(mpiRank==0){cout << "Branching at (indices,bounds,type): ";}
    readOmegaIntoModel(omega);
    //loadOmega();
    //zeroOmega();
    currentLagrLB=-ALPS_DBL_MAX;
    referenceLagrLB=lb;
    setPenalty(pen);
    printOriginalVarBds();
    printCurrentVarBds();
    //subproblemSolvers[0]->printXBounds();
//if(mpiRank==0){cerr << "End installSubproblem ";}
}
void installSubproblem(double lb, const double *zLBs, const double *zUBs, double pen){
//if(mpiRank==0){cerr << "Begin installSubproblem ";}
    for(int ii=0; ii<n1; ii++){
	assert(zLBs[ii] <= zUBs[ii]);
    }
    setLBsAllSPs(zLBs);
    setUBsAllSPs(zUBs);
//if(mpiRank==0){cout << "Branching at (indices,bounds,type): ";}
    //readOmegaIntoModel(omega);
    //loadOmega();
    zeroOmega();
    currentLagrLB=-ALPS_DBL_MAX;
    referenceLagrLB=lb;
    setPenalty(pen);
    printOriginalVarBds();
    printCurrentVarBds();
    //subproblemSolvers[0]->printXBounds();
//if(mpiRank==0){cerr << "End installSubproblem ";}
}

double *getZ(){return z_current;}
double *getZAverage(){return z_average;}
double *getZAverage1(){return z_average1;}
double *getZAverage2(){return z_average2;}
vector< vector<var_branch> >& getNewNodeSPInfo(){return newNodeSPInfo;}
int getNumNewNodeSPs(){return newNodeSPInfo.size();}


//double evaluateXDispersion(double *z=NULL, int branchingNo=0){
double findBranchingIndex(double *z=NULL, int branchingNo=0){
    //if(z==NULL) z=z_current; //Default behaviour
    newNodeSPInfo.clear();
    printTwoZDiscr();
    printAlgorithm();
#if 0
    if(AlgorithmZ){
        if(mpiRank==0){ cout << "Using the consensus vector z to determine branching..." << endl;}
    }
    else{
	if(mpiRank==0){ cout << "Using the average of scenario-specific optimal columns to determine branching..." << endl;}	
    }
#endif
    double maxDisp=ALPS_DBL_MAX;
    double zVal;
    double cVarLB,cVarUB;
    double *dispZAverLocal = new double[n1];
    double *dispZAver = new double[n1];
    double *dispZAver1Local = new double[n1];
    double *dispZAver1 = new double[n1];
    double *dispZAver2Local = new double[n1];
    double *dispZAver2 = new double[n1];
    double *dispZLocal = new double[n1];
    double *dispZ = new double[n1];
    double *dispOmegaLocal = new double[n1];
    double *dispOmega = new double[n1];
    double *dispVerticesLocal = new double[n1];
    double *dispVertices = new double[n1];
    double *coeffLocal = new double[n1];
    double *coeff = new double[n1];
    double *penSumLocal = new double[n1];
    double *penSum = new double[n1];
    double discrep;
    double branchingVal;
    int nVertices,nVerticesAdded;
    double integrDiscr;
    double integrDiscrSum=0.0;
    double integrDiscrSumLocal=0.0;
    double integrDiscrMax=0.0;
    double integrDiscrMaxLocal=0.0;
    double totalWeightLocal=0.0;
    double totalWeightLocal2=0.0;
    double totalWeight=0.0;
    double totalWeight2=0.0;
    double totalZDisp;

    for(int tS=0; tS<nNodeSPs; tS++){
	integrDiscr = subproblemSolvers[tS]->computeIntegralityDiscr();  
	integrDiscrSumLocal += integrDiscr;  
	integrDiscrMaxLocal = max(integrDiscrMaxLocal,integrDiscr);
	//if(integrDiscr > 1e-10) cout << "Node " << mpiRank << " scenario " << tS << " integrality discrepancy: " << integrDiscr << endl;  
    }
#ifdef USING_MPI
    if(mpiSize>1){
        MPI_Allreduce(&integrDiscrSumLocal, &integrDiscrSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&integrDiscrMaxLocal, &integrDiscrMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
#endif
    if(mpiSize==1){
	integrDiscrSum=integrDiscrSumLocal;
	integrDiscrMax=integrDiscrMaxLocal;
    }
    for(int ii=0; ii<n1; ii++){
	dispZLocal[ii]=0.0;
	dispZ[ii]=0.0;
	penSumLocal[ii]=0.0;
	penSum[ii]=0.0;
        for(int tS=0; tS<nNodeSPs; tS++){
	    penSumLocal[ii] += scaling_matrix[tS][ii]*pr[tS];
	    dispZLocal[ii] += scaling_matrix[tS][ii]*pr[tS]*fabs(z_current[ii] - subproblemSolvers[tS]->getX()[ii]);
	}

    }
#ifdef USING_MPI
    if(mpiSize>1){
        MPI_Allreduce(penSumLocal, penSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(dispZLocal, dispZ, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    if(mpiSize==1){
	memcpy(penSum,penSumLocal,n1*sizeof(double));
	memcpy(dispZ,dispZLocal,n1*sizeof(double));
    }
    totalZDisp = 0.0;
    for(int ii=0; ii<n1; ii++){dispZ[ii]/=penSum[ii];}
    if(mpiRank==0){
	cout << "dispZ is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << dispZ[ii];}
	cout << endl;
    }
    for(int ii=0; ii<n1; ii++){totalZDisp+=dispZ[ii];} 
    if(mpiRank==0) cout << "total integer dispersions: " << integrDiscrSum << endl;// || selectOnlyIntegerVars){ 
    if(mpiRank==0) cout << "maximum integer dispersions: " << integrDiscrMax << endl;// || selectOnlyIntegerVars){ 
    if(mpiRank==0) cout << "total z dispersions: " << totalZDisp << endl;// || selectOnlyIntegerVars){ 
    if(mpiRank==0) cout << "discrepNorm: " << discrepNorm << endl;// || selectOnlyIntegerVars){ 

  //if(integrDiscrSum > 1e-10 || totalZDisp > 1e-6){
  //if(integrDiscrSum < 1e-10 && discrepNorm < 1e-20){
  if(integrDiscrSum < 1e-10){
    //roundCurrentZ();
  }
  else{
    solveContinuousMPs(true);
    //roundCurrentZ();
    averageOfVertices();
    averageOfVertices1();
    averageOfVertices2();
    for(int ii=0; ii<n1; ii++){
	dispZAverLocal[ii]=0.0;
	dispZAver[ii]=0.0;
	dispZAver1Local[ii]=0.0;
	dispZAver1[ii]=0.0;
	dispZAver2Local[ii]=0.0;
	dispZAver2[ii]=0.0;
	dispOmegaLocal[ii]=0.0;
	dispOmega[ii]=0.0;
	dispVerticesLocal[ii]=0.0;
	dispVertices[ii]=0.0;
	coeffLocal[ii]=0.0;
	coeff[ii]=0.0;
        for(int tS=0; tS<nNodeSPs; tS++){
 	    integrDiscr = subproblemSolvers[tS]->computeIntegralityDiscr();  
	    if(!(subproblemSolvers[tS]->getXVertex()[ii] >= currentVarLB_[ii] && subproblemSolvers[tS]->getXVertex()[ii] <= currentVarUB_[ii])){
cout << "Scenario " << tS << " index " << ii << " with value " << subproblemSolvers[tS]->getXVertex()[ii] << endl;
	    assert(subproblemSolvers[tS]->getXVertex()[ii] >= currentVarLB_[ii] && subproblemSolvers[tS]->getXVertex()[ii] <= currentVarUB_[ii]);
	    }
	    coeffLocal[ii] += pr[tS]*(subproblemSolvers[tS]->getC()[ii] + omega_current[tS][ii]);
	    dispOmegaLocal[ii] += fabs(pr[tS]*omega_current[tS][ii]);
	    dispVerticesLocal[ii] += (integrDiscr)*(subproblemSolvers[tS]->getDispersions())[ii];
	    //dispZLocal[ii] += fabs(z_current[ii] - subproblemSolvers[tS]->getXVertex()[ii]);
	    dispZAverLocal[ii] += fabs(z_average[ii] - subproblemSolvers[tS]->getXVertex()[ii]);
	    dispZAver1Local[ii] += (integrDiscr)*fabs(z_average1[ii] - subproblemSolvers[tS]->getXVertex()[ii]);
	    if(ii==0){totalWeightLocal += (integrDiscr);}
	    for(int vv=0; vv<subproblemSolvers[tS]->getNVertices(); vv++){
	        dispZAver2Local[ii] += (integrDiscr)*fabs(z_average2[ii] - subproblemSolvers[tS]->getXVertexEntry(ii,vv));
	        if(ii==0){totalWeightLocal2 += (integrDiscr);}
	    }
        }//loop over tS
    }//loop over ii
#ifdef USING_MPI
    if(mpiSize>1){
	MPI_Allreduce(dispZAverLocal, dispZAver, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(dispZAver1Local, dispZAver1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(dispZAver2Local, dispZAver2, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(dispOmegaLocal, dispOmega, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(dispVerticesLocal, dispVertices, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(coeffLocal, coeff, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&totalWeightLocal, &totalWeight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&totalWeightLocal2, &totalWeight2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    if(mpiSize==1){
	memcpy(dispZAver,dispZAverLocal,n1*sizeof(double));
	memcpy(dispZAver1,dispZAver1Local,n1*sizeof(double));
	memcpy(dispZAver2,dispZAver2Local,n1*sizeof(double));
	memcpy(dispOmega,dispOmegaLocal,n1*sizeof(double));
	memcpy(dispVertices,dispVerticesLocal,n1*sizeof(double));
	memcpy(coeff,coeffLocal,n1*sizeof(double));
	totalWeight=totalWeightLocal;
	totalWeight2=totalWeightLocal2;
    }
    for(int ii=0; ii<n1; ii++){ 
	//dispVertices[ii] /= totalWeight;
	dispZAver[ii] /= ((double)nS);
	//dispZAver1[ii] /= totalWeight;
	//dispZAver2[ii] /= totalWeight2;
    }
    if(mpiRank==0){ 
	cout << "-----------------------------" << endl;
	cout << "coeff is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << coeff[ii];}
	cout << endl;
	cout << "ZAverage is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << z_average[ii];}
	cout << endl;
	cout << "ZAverage1 is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << z_average1[ii];}
	cout << endl;
	cout << "ZAverage2 is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << z_average2[ii];}
	cout << endl;
	cout << "ZCurrent is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << z_current[ii];}
	cout << endl;
	cout << "-----------------------------" << endl;
	cout << "dispZAver is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << dispZAver[ii];}
	cout << endl;
	cout << "dispZAver1 is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << dispZAver1[ii];}
	cout << endl;
	cout << "dispZAver2 is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << dispZAver2[ii];}
	cout << endl;
	cout << "dispOmega is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << dispOmega[ii];}
	cout << endl;
	cout << "dispVertices is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << dispVertices[ii];}
	cout << endl;
	cout << "Sum of dispVertices and dispZ is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << dispVertices[ii]+dispZ[ii];}
	cout << endl;
	cout << "-----------------------------" << endl;
    }



	//Now decide on branching index...
	maxDisp=0.0;
	//if(selectOnlyIntegerVars && mpiRank==0) cout << "findBranchingIndex(): only selecting integer variables" << endl;
	double minDisp=0.0;//fabs(disp[0]);
	double testDisp=0.0;
	double gapMin;
    	for(int ii=0; ii<n1; ii++){ 
#if 0
          if(selectOnlyIntegerVars){  //when there is still integer violation, address these first
	    if(colType_[ii]=='C'){continue;}
	    else{
	        testDisp = dispOmega[ii];
	    }
	    //gapMin = min( min( fabs(z[ii]-currentVarLB_[ii]), fabs(z[ii]-currentVarUB_[ii])), 1.0);
	    //gapMin = 1.0;
	  }
#endif
	    //if(selectOnlyIntegerVars && colType_[ii]=='C'){continue;}
	    switch(algorithm){
		case ALGDDBASIC:
		    testDisp = dispZAver[ii];//+dispZ[ii];
		    break;
		case ALGDD:
		    testDisp = dispZAver1[ii]+dispZ[ii];
		    break;
		case ALGDDD:
		    testDisp = dispZAver2[ii]+dispZ[ii];
		    break;
		case ALGZ:
		    testDisp = dispVertices[ii];//+dispZ[ii];
		    break;
		default:
		    testDisp = dispZAver[ii]+dispZ[ii];
		    break;
	    }
#if 0
	    if(AlgorithmZ){
		testDisp = dispVertices[ii]+dispZ[ii];
	    }
	    else{testDisp = dispZAver2[ii];}
#endif
	    //else{testDisp = gapMin*(fabs(coeff[ii])+1.0)*disp[ii];}
	    //testOptSum += testOpt[ii];
	  //}//else
	  if(testDisp > maxDisp){// && (colType_[ii]=='C' || AlgorithmZ)){
	    	maxDisp = testDisp;
	    	infeasIndex_ = ii;
	  }
    	}
        //if(maxDisp < 1e-4 && !selectOnlyIntegerVars){infeasIndex_=-1;} //z solution is regarded as feasible, should update z status to feas or opt.
	//if(infeasIndex_ != -1){
        //if(maxDisp >= 1e-4){// || selectOnlyIntegerVars){ 
        //if(fabs(totalWeight)+totalZDisp >= sqrt(SSC_DEN_TOL)){// || selectOnlyIntegerVars){ 
        //if(fabs(totalWeight)+totalZDisp >= 1e-7){// || selectOnlyIntegerVars){ 
        //if(maxDisp >= 1e-5){// || selectOnlyIntegerVars){ 
	  switch(algorithm){
		case ALGDDBASIC:
	    	    zVal=z_average[infeasIndex_];
		    break;
		case ALGDD:
	    	    zVal=z_average1[infeasIndex_];
		    break;
		case ALGDDD:
	    	    zVal=z_average2[infeasIndex_];
		    break;
		case ALGZ:
	    	    zVal=z_current[infeasIndex_];
		    break;
		default:
	    	    zVal=z_average[infeasIndex_];
		    break;
	  }
#if 0
	  if(AlgorithmZ){
	    zVal=z_current[infeasIndex_];
	  }
	  else{
	    zVal=z_average2[infeasIndex_];
	  }
#endif
	  cVarLB=currentVarLB_[infeasIndex_];
	  cVarUB=currentVarUB_[infeasIndex_];
	  
	  if(colType_[infeasIndex_]=='C'){
		//double newlb = max( (zVal+cVarLB)/2.0, zVal-BRANCHING_EPS);
		//double newub = min( (zVal+cVarUB)/2.0, zVal+BRANCHING_EPS);
		newNodeSPInfo.push_back(vector<var_branch>());
	     if(cVarLB + 1e-6 <= zVal){
		newNodeSPInfo[0].push_back({infeasIndex_, cVarLB, zVal });
	     }
	     else{
		newNodeSPInfo[0].push_back({infeasIndex_, (zVal+cVarLB)/2.0, (zVal+cVarLB)/2.0 });
	     }
		//newNodeSPInfo[0].push_back({infeasIndex_, currentVarLB_[infeasIndex_], newlb });

		//newNodeSPInfo.push_back(vector<var_branch>());
		//newNodeSPInfo[1].push_back({infeasIndex_, newlb, newub });

		newNodeSPInfo.push_back(vector<var_branch>());
	    if(zVal + 1e-6 <= cVarUB){
		newNodeSPInfo[1].push_back({infeasIndex_, zVal, cVarUB });
	    }
	    else{
		newNodeSPInfo[1].push_back({infeasIndex_, (zVal+cVarUB)/2.0, (zVal+cVarUB)/2.0 });
	    }
		//newNodeSPInfo[2].push_back({infeasIndex_, newub, currentVarUB_[infeasIndex_] });
	  }
	  else{//branching var is integer
		//if(checkInteger(zVal)){
		 if( cVarLB <= round(zVal) && round(zVal) <= cVarUB){
		   newNodeSPInfo.push_back(vector<var_branch>());
		   newNodeSPInfo[newNodeSPInfo.size()-1].push_back({infeasIndex_, round(zVal), round(zVal) });
		 }
		 if( cVarLB <= round(zVal)-1){
		    newNodeSPInfo.push_back(vector<var_branch>());
		    newNodeSPInfo[newNodeSPInfo.size()-1].push_back({infeasIndex_, cVarLB, round(zVal)-1 });
		 }
		 if( round(zVal)+1 <= cVarUB){
		    newNodeSPInfo.push_back(vector<var_branch>());
		    newNodeSPInfo[newNodeSPInfo.size()-1].push_back({infeasIndex_, round(zVal)+1, cVarUB });
		 }

	  }
	  //printNewNodeSPInfo();
	//}
    }//if(integrDiscrSum > 1e-10 || totalZDisp > 1e-6){
    delete [] dispZAverLocal;
    delete [] dispZAver;
    delete [] dispZAver1Local;
    delete [] dispZAver1;
    delete [] dispZAver2Local;
    delete [] dispZAver2;
    delete [] dispZLocal;
    delete [] dispZ;
    delete [] dispOmegaLocal;
    delete [] dispOmega;
    delete [] dispVerticesLocal;
    delete [] dispVertices;
    delete [] coeffLocal;
    delete [] coeff;
    delete [] penSumLocal;
    delete [] penSum;
    return maxDisp;
}

double* getOrigVarLbds(){return origVarLB_;}
double* getOrigVarUbds(){return origVarUB_;}
double* getCurrentVarLbds(){return currentVarLB_;}
double* getCurrentVarUbds(){return currentVarUB_;}
void printOriginalVarBds(){
  if(mpiRank==0){
    cout << "Printing original variable bounds:" << endl;
    for(int ii=0; ii<n1; ii++){
	cout << " (" << origVarLB_[ii] << "," << origVarUB_[ii] << ")";
    }
    cout << endl;
  }
}
void printCurrentVarBds(){
  if(mpiRank==0){
    cout << "Printing current variable bounds:" << endl;
    for(int ii=0; ii<n1; ii++){
	cout << " (" << currentVarLB_[ii] << "," << currentVarUB_[ii] << ")";
    }
    cout << endl;
  }
}
//void evaluateXDispersion(){
    //double disp=computeXDispersion();
    //if(mpiRank==0) cout << "X dispersion given z is: " << disp << endl;
//}

vector<double*>& getOmega(){return omega_current;}
double *getIncumbentZ(){return z_incumbent_;}
double getIncumbentVal(){return getKnowledgeBroker()->getIncumbentValue();}
double getObjVal(){return objVal;}
int getInfeasIndex(){return infeasIndex_;}
void setInfeasIndex(int ii){infeasIndex_=ii;}

void readZIntoModel(const double* z){
    memcpy(z_current,z,n1*sizeof(double));
}
void saveZ(const double* z=NULL){
    if(z==NULL){
        memcpy(z_saved,z_current,n1*sizeof(double));
    }
    else{
        memcpy(z_saved,z,n1*sizeof(double));
    }
}
void loadSavedZ(double* z=NULL){
    if(z==NULL){
        memcpy(z_current,z_saved,n1*sizeof(double));
    }
    else{
        memcpy(z,z_saved,n1*sizeof(double));
    }
}
void readOmegaIntoModel(vector<double*> &omega){
//cout << "nNodeSPs " << nNodeSPs << " size of vec " << omega.size() << endl;
    omegaIsZero_=false;
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_current[tS],omega[tS],n1*sizeof(double));
}
void loadOmega(){
    omegaIsZero_=false;
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_current[tS],omega_saved[tS],n1*sizeof(double));
}
void loadOmega(vector<double*> &omega){
    omegaIsZero_=false;
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega[tS],omega_saved[tS],n1*sizeof(double));
}
void zeroOmega(){
    omegaIsZero_=true;
    for(int tS=0; tS<nNodeSPs; tS++){ for(int ii=0; ii<n1; ii++){omega_current[tS][ii]=0.0;}}
}
void decimateOmega(double byFactor = 0.2){
    if(byFactor < 0.0) byFactor=0.0;
    else if(byFactor > 1) byFactor=1.0;
    for(int tS=0; tS<nNodeSPs; tS++){ for(int ii=0; ii<n1; ii++){omega_current[tS][ii]*=byFactor;}}
}
void saveOmega(){
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_saved[tS],omega_current[tS],n1*sizeof(double));
}

#if 1
virtual AlpsTreeNode* createRoot();
#endif

int initialIteration();

int regularIteration(bool scalePenalty=false, bool adjustPenalty=false, bool SSC=true){
//if(mpiRank==0) cout << "Begin regularIteration()" << endl;
	solveContinuousMPs();

	int SPStatus=performColGenStep();
	assert(SPStatus!=SP_INFEAS); //subproblem infeasibility should be caught in initialIteration()
	double SSCVal = computeSSCVal(); //shouldTerminate is also updated here.
    //verifyOmegaDualFeas();
    //printStatus();
	//cout << "ALVal: " << ALVal << " discrepNorm: " << discrepNorm << " currentLagrLB " << currentLagrLB << endl;
	if(SSC){updateOmega(SSCVal);}
	else{updateOmega(1.0);}
	#if KIWIEL_PENALTY 
	 if(adjustPenalty) computeKiwielPenaltyUpdate(SSCVal);
	#endif
	 if(omegaUpdated_) computePenaltyGapVec();
//printIntegralityViolations();
//if(mpiRank==0) cout << "End regularIteration()" << endl;
	return SPStatus;
}
void setMaxNoSteps(int noSteps){maxNoSteps=noSteps;}

double computeBound(int nIters, bool adjustPenalty=false){
//if(mpiRank==0) cout << "Begin computeBound()" << endl;
    modelStatus_[Z_STATUS]=Z_UNKNOWN;
    clearSPVertexHistory();
    int SPStatus=initialIteration();
    if(mpiRank==0) cout << "After initial iteration: currentLB: " << currentLagrLB << " versus refLB: " << referenceLagrLB << endl;;
    
    bool exceedingReferenceBd=true;
    bool omegaUpdatedAtLeastOnce=false;
    shouldFathomByOpt = false;
    double integralityDisc = 100.0;
    discrepNorm = 100.0;
    shouldTerminate=false;
    if(SPStatus==SP_INFEAS){
	modelStatus_[Z_STATUS]=Z_INFEAS;
if(mpiRank==0){cout << "Terminating due to subproblem infeasibility..." << endl;}
    }
    else if(currentLagrLB >= getIncumbentVal()){
	modelStatus_[Z_STATUS]=Z_BOUNDED;
if(mpiRank==0){cout << "Terminating due to exceeding cutoff..." << endl;}
      	//printStatus();
    }
    else{
      //currentLagrLB=-ALPS_DBL_MAX;
      //zeroOmega();
      bool shouldContinue=true;
      //for(int ii=0; ii<nIters || !exceedingReferenceBd || !omegaUpdatedAtLeastOnce || (discrepNorm >= 1e-20 && integralityDisc < 1e-10) ; ii++){
      for(int ii=0; shouldContinue; ii++){
if(mpiRank==0) cerr << "Regular iteration " << ii << endl;
if(mpiRank==0 && integralityDisc < 1e-10){cout << "computeBound(): Integrality satisfied: Commence with postprocessing..." << endl;}
	if(ii<10) {regularIteration(false,adjustPenalty);}
	else {regularIteration(true,false);}
	//if(omegaUpdated_) omegaUpdatedAtLeastOnce=true;
	//omegaUpdatedAtLeastOnce=true;
      	printStatusAsErr();
        if(currentLagrLB >= getIncumbentVal()){
      	    modelStatus_[Z_STATUS]=Z_BOUNDED;
if(mpiRank==0){cout << "computeBound(): Terminating due to exceeding cutoff..." << endl;}
      	    //printStatus();
	    break;
	}
        exceedingReferenceBd = currentLagrLB + SSC_DEN_TOL >= referenceLagrLB;
        //exceedingReferenceBd = true;
#if 1
	//if(shouldTerminate()){
	if(shouldTerminate){
if(mpiRank==0){cout << "computeBound(): Terminating due to reaching tolerance criterion at iteration " << ii << endl;}
      	    //updateModelStatusZ();
      	    //printStatus();
	    //return currentLagrLB;
	    break;
	}
#endif
	integralityDisc=evaluateIntegralityDiscrepancies();
	shouldContinue = (ii<nIters) || !exceedingReferenceBd || (integralityDisc < 1e-10 && discrepNorm >=1e-20); 
    //if(currentLagrLB < referenceLagrLB && mpiRank==0){
if(mpiRank==0 && ii>=(nIters-1) && exceedingReferenceBd && omegaUpdatedAtLeastOnce){cout << "computeBound(): Terminating due to reaching the maximum number of iterations." << endl;}
    }
    roundCurrentZ();
    //shouldFathomByOpt = integralityDisc < 1e-10 && (shouldTerminate || discrepNorm < 1e-20);
#if 0
integralityDisc=evaluateIntegralityDiscrepancies();
if(mpiRank==0)  cout << "computeBound(): Integrality discrepancy in total: " << integralityDisc << endl;
if(integralityDisc < 1e-10){
  if(mpiRank==0){cout << "computeBound(): Integrality satisfied: Commence with postprocessing..." << endl;}
  setPenalty(baselinePenalty_);
  for(int ii=0; (ii<200 || (discrepNorm >= 1e-20)) && integralityDisc < 1e-10; ii++){
	regularIteration(true,false);
      	printStatusAsErr();
        if(currentLagrLB >= getIncumbentVal()){
      	    modelStatus_[Z_STATUS]=Z_BOUNDED;
	    break;
	}
#if 1
	if(shouldTerminate){
	    if(mpiRank==0){cout << "computeBound(): Terminating due to reaching tolerance criterion at iteration " << ii << endl;}
	    break;
	}
#endif
	integralityDisc=evaluateIntegralityDiscrepancies();
  }
    if(mpiRank==0){cout << "New incumbent value should be: " << currentLagrLB << endl;}
    //roundCurrentZ();
}
else{
#endif
//if(mpiRank==0){cout << "computeBound(): Commencing another loop of iterations since." << endl;}
    //}while(currentLagrLB < referenceLagrLB);
    #if 0
    if(currentLagrLB < referenceLagrLB && mpiRank==0){
	cout << "computeBound(): Lagr bound does not meet reference: " << currentLagrLB << " < " << referenceLagrLB << endl;
    }
    #endif
      //assert(SPStatus!=SP_INFEAS);
      //assert(modelStatus_[Z_STATUS]==Z_UNKNOWN);
      //updateModelStatusZ();
      //loopIter++;
    //}
    //while(getZStatus()==Z_REC_INFEAS && loopIter < maxNLoopIters);
#if 0
    if(getZStatus()==Z_REC_INFEAS){
	modelStatus_[Z_STATUS]=Z_FEAS;
    }
#endif
//if(mpiRank==0) cout << "End computeBound()" << endl;
    //saveOmega(); //omega should only be saved under certain circumstances
    //saveZ();
    //setPenalty(baselinePenalty_);
//}//else integrality not satisfied
    printStatus();
    }//else
    //solveContinuousMPs();
    //for (int tS = 0; tS < nNodeSPs; tS++){subproblemSolvers[tS]->computeWeightsForCurrentSoln();}
    //loadSavedZ();
    //printOmegaProperties();
    //averageOfOptVertices();
    //double bestNodeQuality = getKnowledgeBroker()->getBestNode()->getQuality();
    return currentLagrLB;
}

double evaluateIntegralityDiscrepancies(){
    double integrDiscr,integrDiscrSumLocal=0.0,integrDiscrSum=0.0;
    for(int tS=0; tS<nNodeSPs; tS++){
	integrDiscr = subproblemSolvers[tS]->computeIntegralityDiscr();  
	integrDiscrSumLocal += integrDiscr;  
	//if(integrDiscr > 1e-10) cout << "Node " << mpiRank << " scenario " << tS << " integrality discrepancy: " << integrDiscr << endl;  
    }
#ifdef USING_MPI
    if(mpiSize>1){
        MPI_Allreduce(&integrDiscrSumLocal, &integrDiscrSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    if(mpiSize==1){
	integrDiscrSum=integrDiscrSumLocal;
    }
    return integrDiscrSum;
}
#if 0
void resetOmega(){
     currentLagrLB=-ALPS_DBL_MAX;
     decimateOmega(0.8);
     setPenalty(baselinePenalty_);
     for (int tS = 0; tS < nNodeSPs; tS++) {
	 recordKeeping[tS][0]=-ALPS_DBL_MAX;//subproblemSolvers[tS]->getLagrBd();
	 recordKeeping[tS][1]=recordKeeping[tS][0];
     }
}
#endif

//Perform more iterations, e.g., to obtain recourse when integrality met.
#if 0
double findPrimalFeasSoln2(int nIters){

  for(int ii=0; ii<10; ii++){
    setPenalty(pow(2.0,ii));
    solveQMIPsGS();
//printZ();
  }
  if(solveRecourseProblemGivenFixedZ()){
    if(mpiRank==0) {cout << "findPrimalFeasSoln(): Terminating: Found primal feas soln..." << endl;}
    evaluateFeasibleZ();
  }
  else{
    if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Primal feas soln not found." << endl;
  }
  return currentLagrLB;
}
#endif
double findPrimalFeasSoln(int nIters){
//if(mpiRank==0){cout << "Before z: " << endl;}
//printZ();
  //bool foundCandidateSoln=true;
  bool foundCandidateSoln=!selectOnlyIntegerVars;
#if 0
  double roundingDisc = 0.0;
  for(int ii=0; ii<n1; ii++){roundingDisc += fabs(z_rounded[ii]-z_current[ii]);}
  if(roundingDisc <= 1e-10){
      if(mpiRank==0) cout << "findPrimalFeasSoln(): No need for postprocessing...testing for recourse..." << endl;
      foundCandidateSoln=true;
  }
  else{//rounding discrepancy outside of tolerance

    if(mpiRank==0) cout << "findPrimalFeasSoln(): Rounding discrepancy is: " << roundingDisc << ", so there is need for postprocessing..." << endl;
    for(int ii=0; ii<numIntVars_; ii++){fixVarAllSPsAt(intVar_[ii], z_rounded[intVar_[ii]]);}
    clearSPVertexHistory();
    int SPStatus=initialIteration();

    if(SPStatus==SP_INFEAS){
	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Primal feas soln not found." << endl;
        foundCandidateSoln=false;
    }
    else if(currentLagrLB >= getIncumbentVal()){
	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Any primal feas soln to be found would not be better in value than current incumbent value." << endl;
      	printStatus();
        foundCandidateSoln=false;
    }
    else{
      currentLagrLB=-ALPS_DBL_MAX;
      zeroOmega();
      for(int ii=0; ii<nIters; ii++){
	if(mpiRank==0) cerr << "Regular iteration " << ii << endl;
	if(nIters<20) {regularIteration(false,true);}
	else {regularIteration(false,false);}
      	printStatusAsErr();

        if(currentLagrLB >= getIncumbentVal()){
	    if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Any primal feas soln to be found would not be better in value than current incumbent value." << endl;
            foundCandidateSoln=false;
	    break;
	}
#if 1
	if(shouldTerminate()){
	    if(mpiRank==0){cout << "findPrimalFeasSoln(): Terminating due to reaching tolerance criterion..." << endl;}
            foundCandidateSoln=true;
	    break;
	}
#endif
     }
     printStatus();
#if 0
     if(solveRecourseProblemGivenFixedZ()){
	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Found primal feas soln..." << endl;
	evaluateFeasibleZ();
     }
     else{
 	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Primal feas soln not found." << endl;
     }
#endif
    }//else
	foundCandidateSoln=false;
   }//else rounding discrepancy is outside of tolerance
#endif

   if(foundCandidateSoln){
     if(solveRecourseProblemGivenFixedZ()){
	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Found primal feas soln..." << endl;
	evaluateFeasibleZ();
     }
     else{
 	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Primal feas soln not found." << endl;
     }
   }
   return objVal;
}

double getBound(){return currentLagrLB;}
void setBound(double bd){currentLagrLB=bd;}

void solveContinuousMPs(bool updateDisp=false){
	bool isLastGSIt;
	double zDiff;
	//for(int itGS=0; itGS < fixInnerStep || (updateDisp && zDiff > 1e-6); itGS++) { //The inner loop has a fixed number of occurences
	for(int itGS=0; itGS < fixInnerStep; itGS++) { //The inner loop has a fixed number of occurences
	    memcpy(z_old,z_current, n1*sizeof(double));
	    zDiff=0.0;
	    isLastGSIt = (itGS==fixInnerStep-1) && updateDisp;
    	    for (int tS = 0; tS < nNodeSPs; tS++) {

		//*************************** Quadratic subproblem ***********************************
			    
		    if (useVertexHistory && subproblemSolvers[tS]->getNVertices()>2) { //Compute Next X, Y (With History)

					
			//Solve the quadratic master problem for x and y
			//subproblemSolvers[tS]->updatePrimalVariablesHistory_OneScenario(omega_current[tS],z_current,scaling_matrix[tS],isLastGSIt);
			subproblemSolvers[tS]->updatePrimalVariablesHistory_OneScenario(omega_current[tS],z_current,scaling_matrix[tS],false);
    			//if(isLastGSIt) subproblemSolvers[tS]->computeWeightsForCurrentSoln();
			
			// note: the final weight corresponds to the existing x
		    }
		    else if(!useVertexHistory || subproblemSolvers[tS]->getNVertices()==2){ //might not work correctly, untested! Compute Next X, Y (Without History)

			subproblemSolvers[tS]->updatePrimalVariables_OneScenario(omega_current[tS],z_current,scaling_matrix[tS],isLastGSIt);
		    }
	    }
	    					
	    // Update z_previous.
	    updateZ();
	    totalNoGSSteps++;
	    for(int ii=0; ii<n1; ii++){zDiff=max(zDiff,fabs(z_current[ii]-z_old[ii]));}
	}
	if (useVertexHistory && updateDisp) { 
    	    for (int tS = 0; tS < nNodeSPs; tS++) {
    		subproblemSolvers[tS]->computeWeightsForCurrentSoln(z_current);
	    }
	}

	
    	for (int tS = 0; tS < nNodeSPs; tS++) {
	    subproblemSolvers[tS]->updateALValues(omega_current[tS],z_current,scaling_matrix[tS]);
	    //if(tS==0) subproblemSolvers[tS]->printXVertices();
	}
//if(mpiRank==0){cout << "solveContinuousMPs(): max zDiff is: " << zDiff << endl;}
}
#if 0
void solveQMIPsGS(){
  for(int ii=0; ii<10; ii++){
    for (int tS = 0; tS < nNodeSPs; tS++) {
	subproblemSolvers[tS]->solveAugmentedLagrangianMIP(omega_current[tS], z_current, scaling_matrix[tS]);
	subproblemSolvers[tS]->setXToVertex();
	subproblemSolvers[tS]->setYToVertex();
    }
    updateZ(true);
  }
}
#endif

int performColGenStep();
int performColGenStepBasic();


void updateZ(bool roundZ=false){
	double *penSumLocal = new double[n1];
	double *penSum = new double[n1];
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z_current[i] = 0.0;
	    penSumLocal[i]=0.0;
            penSum[i]=0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		penSumLocal[i]+= scaling_matrix[tS][i]*pr[tS];
	    }
	}
	#ifdef USING_MPI
	    MPI_Allreduce(z_local, z_current, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(penSumLocal, penSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    for(int i=0; i<n1; i++) z_current[i] /= penSum[i];
	#else
	    for (int i=0; i<n1; i++) z_current[i] = z_local[i]/penSumLocal[i]; //Only one node, trivially set z_local = z_current
	#endif
	if(roundZ){//rounding integer restricted vars to nearest integer
	  for(int ii=0; ii< numIntVars_; ii++){
	    z_current[intVar_[ii]]=round(z_current[intVar_[ii]]);
	  }
	}
	delete [] penSumLocal;
	delete [] penSum;
}

#if 0
void averageOfVertices(double *z=NULL, bool roundZ=false){
	if(z==NULL) z=z_average;
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += subproblemSolvers[tS]->getXVertex()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
	    }
	}
	#ifdef USING_MPI
	    MPI_Allreduce(z_local, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) z[i] = z_local[i]; //Only one node, trivially set z_local = z
	#endif
	//for (int i=0; i<n1; i++) z[i] /= ((double)nS); //Only one node, trivially set z_local = z
	for (int i=0; i<n1; i++) z[i] /= ((double)nS); //Only one node, trivially set z_local = z
	if(roundZ){//rounding integer restricted vars to nearest integer
	  for(int ii=0; ii< numIntVars_; ii++){
	    z[intVar_[ii]]=round(z[intVar_[ii]]);
	  }
	}
}
#endif
void averageOfVertices(double* z=NULL,bool roundZ=false){
	if(z==NULL) z=z_average;
	double *penSumLocal = new double[n1];
	double *penSum = new double[n1];
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    penSumLocal[i]=0.0;
            penSum[i]=0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += scaling_matrix[tS][i]*pr[tS]*subproblemSolvers[tS]->getXVertex()[i];
		penSumLocal[i]+= scaling_matrix[tS][i]*pr[tS];
	    }
	}
	#ifdef USING_MPI
	    MPI_Allreduce(z_local, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(penSumLocal, penSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    for(int i=0; i<n1; i++) z[i] /= penSum[i];
	#else
	    for (int i=0; i<n1; i++) z[i] = z_local[i]/penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif
	if(roundZ){//rounding integer restricted vars to nearest integer
	  for(int ii=0; ii< numIntVars_; ii++){
	    z[intVar_[ii]]=round(z[intVar_[ii]]);
	  }
	}
	delete [] penSumLocal;
	delete [] penSum;
}

void averageOfVertices1(double *z=NULL, bool roundZ=false){
	if(z==NULL) z=z_average1;
	double localTotalWeight=0.0;
	double totalWeight=0.0;
	double weight;
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		weight = subproblemSolvers[tS]->computeIntegralityDiscr();
		z_local[i] += weight*subproblemSolvers[tS]->getXVertex()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		if(i==0){localTotalWeight += weight;}
	    }
	}
	#ifdef USING_MPI
	    MPI_Allreduce(z_local, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&localTotalWeight, &totalWeight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) z[i] = z_local[i]; //Only one node, trivially set z_local = z
	    totalWeight=localTotalWeight;
	#endif
	//for (int i=0; i<n1; i++) z[i] /= ((double)nS); //Only one node, trivially set z_local = z
	for (int i=0; i<n1; i++) z[i] /= totalWeight; //Only one node, trivially set z_local = z
	if(roundZ){//rounding integer restricted vars to nearest integer
	  for(int ii=0; ii< numIntVars_; ii++){
	    z[intVar_[ii]]=round(z[intVar_[ii]]);
	  }
	}
}
void averageOfVertices2(double *z=NULL, bool roundZ=false){
	if(z==NULL) z=z_average2;
	int numVertices;
	double localTotalWeight=0.0;
	double totalWeight=0.0;
	double weight;
	//for (int tS = 0; tS < nNodeSPs; tS++){ localTotalNVertices+=subproblemSolvers[tS]->getNVertices();}
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		numVertices = subproblemSolvers[tS]->getNVertices();
		weight = subproblemSolvers[tS]->computeIntegralityDiscr();
		for(int vv=0; vv<numVertices; vv++){
		//z_local[i] += x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		    //z_local[i] += subproblemSolvers[tS]->getXVertex()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		    z_local[i] += weight*subproblemSolvers[tS]->getXVertexEntry(i,vv);//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		    if(i==0){localTotalWeight += weight;}
		}
		//localTotalNVertices+=numVertices;
	    }
	}
	#ifdef USING_MPI
	    MPI_Allreduce(z_local, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(&localTotalWeight, &totalWeight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) z[i] = z_local[i]; //Only one node, trivially set z_local = z
	    totalWeight=localTotalWeight;
	#endif
	for (int i=0; i<n1; i++) z[i] /= totalWeight; //Only one node, trivially set z_local = z
	if(roundZ){//rounding integer restricted vars to nearest integer
	  for(int ii=0; ii< numIntVars_; ii++){
	    z[intVar_[ii]]=round(z[intVar_[ii]]);
	  }
	}
}


void averageOfOptVertices(double *z=NULL, bool roundZ=false){
	if(z==NULL) z=z_average;
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		//z_local[i] += x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		z_local[i] += subproblemSolvers[tS]->getXVertexOpt()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
	    }
	}
	#ifdef USING_MPI
//if(mpiRank==0) cerr << "updateZ(): z_local: before MPI_Allreduce()" << endl;
//cerr << "#";
	    MPI_Allreduce(z_local, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//if(mpiRank==0) cerr << endl << "updateZ(): z_local: after MPI_Allreduce()" << endl;
//if(mpiRank==0) cerr << "updateZ(): penSum: before MPI_Allreduce()" << endl;
//cerr << "#";
	#else
	    for (int i=0; i<n1; i++) z[i] = z_local[i]; //Only one node, trivially set z_local = z
	#endif
	for (int i=0; i<n1; i++) z[i] /= ((double)nS); //Only one node, trivially set z_local = z
	if(roundZ){//rounding integer restricted vars to nearest integer
	  for(int ii=0; ii< numIntVars_; ii++){
	    z[intVar_[ii]]=round(z[intVar_[ii]]);
	  }
	}
}
void printTwoZDiscr(){
if(mpiRank==0){
    cout << "Printing difference between z_current and z_average:" << endl;
    for(int ii=0; ii<n1; ii++){
	cout << "  " << z_current[ii]-z_average[ii];
    }
    cout << endl;
}
}

#if 0
int numIntInfeas(double *z=NULL, int branchingNo=0)
{
if(z==NULL){ z=z_current;} //default behaviour
#if 1
    nIntInfeas_ = 0;
    infeasIndex_=-1;
    //int i = -1;
//printZ(); 
    for (int i = 0; i < numIntVars_; ++i) {
      if ( ! checkInteger(z[intVar_[i]]) ) {
//cout << "numIntVars_ " << numIntVars_ << " i = " << i << " index " << intVar_[i] << endl;
	++nIntInfeas_;
	if(infeasIndex_==-1){
	    infeasIndex_=intVar_[i];
	    branchingZVals[branchingNo]=z[infeasIndex_];
	} //Just take the first violation index found, for now.
      }
    }
#endif
    return nIntInfeas_;
}

bool verifyIntegrality()
{
    return nIntInfeas_==0;
}
#endif

/** Check if a value is integer. */
bool checkInteger(double value) const {
    double integerTolerance = 1.0e-8;
    double nearest = floor(value + 0.5);
//cout << " " << fabs(value-nearest);
    if (fabs(value - nearest) <= integerTolerance) {
        return true;
    }
    else {
        return false;
    }
}
#if 0
bool checkZIsFeasForScen1(int tS){
    memcpy(totalSoln_,z_current,n1*sizeof(double));
    double *ySoln = subproblemSolvers[tS]->getYVertex();
    memcpy(totalSoln_+n1,ySoln,n2*sizeof(double));
    return subproblemSolvers[tS]->checkSolnForFeasibility(totalSoln_, constrVec_);
}
bool checkZIsFeasForScen2(int tS){
    subproblemSolvers[tS]->solveFeasibilityProblemWithXFixedToZ(z_current, origVarLB_, origVarUB_,colType_);
    //subproblemSolvers[tS]->solveLagrangianWithXFixedToZ(z_current, omega_current[tS], origVarLB_, origVarUB_,colType_);
    return (subproblemSolvers[tS]->getSolverStatus()==PSCG_OPTIMAL || subproblemSolvers[tS]->getSolverStatus()==PSCG_ITER_LIM);
}
bool checkZIsFeasForScen3(int tS){
    //subproblemSolvers[tS]->solveFeasibilityProblemWithXFixedToZ(z_current, origVarLB_, origVarUB_,colType_);
    subproblemSolvers[tS]->solveLagrangianWithXFixedToZ(z_current, omega_current[tS], origVarLB_, origVarUB_,colType_);
    return (subproblemSolvers[tS]->getSolverStatus()==PSCG_OPTIMAL || subproblemSolvers[tS]->getSolverStatus()==PSCG_ITER_LIM);
}
#endif

bool indexIsInt(int ii){
assert(ii>=0);
assert(ii < n1);
	if(colType_[ii]=='I' || colType_[ii]=='B'){
	  return true;
	}
	else{
	  return false;
	}

}

bool checkForCutoff(){
if(mpiRank==0){cout << getIncumbentVal() << " <= " << currentLagrLB << "???" << endl;}
    return (getIncumbentVal() <= currentLagrLB);
}

void roundCurrentZ(){
    memcpy(z_rounded,z_current,n1*sizeof(double));
    for(int ii=0; ii<numIntVars_; ii++) {z_rounded[intVar_[ii]] = round(z_current[intVar_[ii]]);}
    double roundingDisc = 0.0;
    for(int ii=0; ii<n1; ii++){roundingDisc += fabs(z_rounded[ii]-z_current[ii]);}
    selectOnlyIntegerVars = (roundingDisc > 1e-10);
    if(mpiRank==0) cout << "roundCurrentZ(): Rounding discrepancy is: " << roundingDisc << endl;
}

bool evaluateFeasibleZ(){
    if(getIncumbentVal() > objVal){ 
if(mpiRank==0){
    cout << "New incumbent value: " << objVal << " and its corresponding solution: " << endl;
    printZRounded();
}
	//getKnowledgeBroker()->setIncumbentValue(objVal);
	memcpy( z_incumbent_, z_rounded, n1*sizeof(double) );
	BcpsSolution* ksol = new BcpsSolution(n1,z_incumbent_,objVal);
        getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,ksol,objVal);
#if 0
                        new BcpsSolution(model->n1,
                                           model->getZ(),
                                           model->getIncumbentVal());
                    getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,
                                                       ksol,
                                                       model->getIncumbentVal());
#endif

	getKnowledgeBroker()->setPhase(AlpsPhaseSearch);
//cout << zOptFile << endl;
        if(mpiRank==0) getKnowledgeBroker()->printBestSolution(zOptFile);
	if(mpiRank==0) cout << "Registering solution..." << endl;
//cout << "End evaluateFeasibleZ()" << endl;
	return true;
    }
//cout << "feasSolnUpdated returns false" << endl;
//if(mpiRank==0){cout << "End evaluateFeasibleZ()" << endl;}
    return false;
}
//int solveFeasibilityProblemWithXFixedToZ(const double *z, const double *origLBs, const double *origUBs, const char *colTypes){

#if 0
bool checkZIsFeasForCurrentY(){
    for(int tS=0; tS<nNodeSPs; tS++){
	if(!checkZIsFeasForScen1(tS)){ 
	    modelStatus_[Z_STATUS]=Z_REC_INFEAS;
	    return false;
	}
    }
    return true;
}

bool checkZHasFullRecourse(){
    for(int tS=0; tS<nNodeSPs; tS++){
	if(!checkZIsFeasForScen2(tS)){ 
	    modelStatus_[Z_STATUS]=Z_REC_INFEAS;
	    return false;
	}
    }
    return true;
}
#endif

bool solveRecourseProblemGivenFixedZ();	

#if 0
void updateModelStatusSP(){
//if(mpiRank==0){cout << "Begin updateModelStatusSP()" << endl;}
enum Statuses{
    SP_STATUS=0,
    Z_STATUS
};

enum SPStatuses{
    SP_OPT=0, //PSCG termination criteria met
    SP_ITER_LIM, //otherwise feasible
    SP_INFEAS //at least one subproblem is infeasible
};

    modelStatus_[SP_STATUS] = SP_ITER_LIM;
    //modelStatus_[Z_STATUS] = Z_FEAS;
    for(int tS=0; tS<nNodeSPs; tS++){
	if(spSolverStatuses_[tS]==SP_INFEAS){
	    modelStatus_[SP_STATUS]=SP_INFEAS;
//	    modelStatus_[Z_STATUS]=Z_INFEAS;
	    return;
	}
    }
    bool PSCGTerm = shouldTerminate();
    if(modelStatus_[SP_STATUS]!=SP_INFEAS && PSCGTerm){
	modelStatus_[SP_STATUS]=SP_OPT;
    }
//if(mpiRank==0){cout << "End updateModelStatusSP()" << endl;}
}
#endif

#if 0
void updateModelStatusZ(){
//if(mpiRank==0){cout << "Begin updateModelStatusZ()" << endl;}
enum Statuses{
    SP_STATUS=0,
    Z_STATUS
};
enum Z_Statuses{
    Z_OPT=0,
    Z_FEAS, //z is feasible (has both recourse and integer feas)
    Z_REC_INFEAS, //z has recourse (but its integer feas unknown)
    Z_INT_INFEAS, //z is integer feas (but its recourse is unknown)
    Z_INFEAS //z is infeasible by both feasibility qualities
};
    //checkZIsFeasForCurrentY();
    if(!solveRecourseProblemGivenFixedZ()){
        modelStatus_[Z_STATUS]=Z_REC_INFEAS;
    }	
    else{
        modelStatus_[Z_STATUS]=Z_FEAS;
	//evaluateFeasibleZ();
    }
    //checkZHasFullRecourse();
    //if(numIntInfeas() > 0){
    if(infeasIndex_!=-1){
	if(modelStatus_[Z_STATUS]==Z_REC_INFEAS){
	    modelStatus_[Z_STATUS]=Z_INFEAS;
	}
	else{
	    modelStatus_[Z_STATUS]=Z_INT_INFEAS;
	}
    }
    if(modelStatus_[Z_STATUS]==Z_FEAS){
	//evaluateFeasibleZ();
	if(modelStatus_[SP_STATUS]==SP_OPT){
	    modelStatus_[Z_STATUS]=Z_OPT;
	}
    }
//if(mpiRank==0){cout << "End updateModelStatusZ()" << endl;}
}
#endif

int getSPStatus(){return modelStatus_[SP_STATUS];}
void setSPStatus(int stat){modelStatus_[SP_STATUS]=stat;}
int getZStatus(){return modelStatus_[Z_STATUS];}
void setZStatus(int stat){modelStatus_[Z_STATUS]=stat;}


//********************** Serious Step Condition (SSC) **********************
#if 0
bool shouldTerminate(){
//if(mpiRank==0){cout << "shouldTerminate(): " << (ALVal + 0.5*discrepNorm  - currentLagrLB) << endl;}
    return (ALVal + 0.5*discrepNorm  - currentLagrLB < SSC_DEN_TOL);
}
#endif
double computeSSCVal(){
    if(ALVal - currentLagrLB < -1e-6){
        if(mpiRank==0){cout << " (ALVal,currentLagrLB) (" << ALVal << "," << currentLagrLB << ")" << endl;}
	shouldTerminate = true;
	return 0.0;
    }
    if(ALVal + 0.5*discrepNorm  - LagrLB < -1e-6){
        if(mpiRank==0){cout << mpiRank << " (ALVal,LagrLB) (" << ALVal << "," << LagrLB << ")" << endl;	
	cout << "regularIteration(): Something probably went wrong with the last computation of LagrLB, returning..." << endl;}
	shouldTerminate = true;
	return 0.0;
    }
    shouldTerminate = (ALVal + 0.5*discrepNorm  - currentLagrLB < SSC_DEN_TOL);
    if(shouldTerminate){return 0.0;}
    else{return (LagrLB-currentLagrLB)/(ALVal + 0.5*discrepNorm - currentLagrLB);}
}
		
bool updateOmega(double SSCVal){
	if(SSCVal >= SSC_PARAM) {
	    for (int tS = 0; tS < nNodeSPs; tS++) {
		memcpy(omega_current[tS],omega_tilde[tS],n1*sizeof(double));
#if 0
		for (int i = 0; i < n1; i++) {
	    	    omega_current[tS][i] = omega_tilde[tS][i];
		}
#endif
	        recordKeeping[tS][1]=recordKeeping[tS][0];
	        subproblemSolvers[tS]->updateOptSoln();
	    }
	    currentLagrLB = LagrLB;
	    //recordKeeping[0]=LagrLB;
    	    omegaUpdated_ = true;
    	    omegaIsZero_=false;
//if(mpiRank==0){cout << "Updating omega..." << endl;}
	}
	else {
    	    omegaUpdated_ = false;
	    //if(mpiRank==0) cout << "Null step taken..." << endl;
	}
	return omegaUpdated_;
}
#if 1
bool printOmegaProperties(){
    double *omegaLocal = new double[n1];
    double *omegaSum = new double[n1]; 
    double omegaFroNormLocal = 0.0;
    double omegaFroNorm = 0.0; 
    double cFroNormLocal = 0.0;
    double cFroNorm = 0.0; 
    for (int ii = 0; ii < n1; ii++) {
	omegaLocal[ii]=0.0;
    	for (int tS = 0; tS < nNodeSPs; tS++) {
    	    omegaLocal[ii] +=pr[tS]*omega_current[tS][ii];
	    omegaFroNormLocal += omega_current[tS][ii]*omega_current[tS][ii];
	    cFroNormLocal += (subproblemSolvers[tS]->getC()[ii])*(subproblemSolvers[tS]->getC()[ii]);
	}
    }
#ifdef USING_MPI
    if(mpiSize>1){
	MPI_Allreduce(omegaLocal, omegaSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&omegaFroNormLocal, &omegaFroNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&cFroNormLocal, &cFroNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    if(mpiSize==1){
	memcpy(omegaSum,omegaLocal,n1*sizeof(double));
	omegaFroNorm = omegaFroNormLocal;
	cFroNorm = cFroNormLocal;
    }
if(mpiRank==0){
cout << "Printing omega weighted sum..." << endl;
    for(int ii=0; ii<n1; ii++){
	cout << " " << omegaSum[ii];	
    }
cout << endl;
cout << "...and the Frobenius norm of omega: " << sqrt(omegaFroNorm) << endl;
cout << "...versus the Frobenius norm of c: " << sqrt(cFroNorm) << endl;
}
    delete [] omegaLocal;
    delete [] omegaSum;
}
#endif

void printNewNodeSPInfo(){
  if(mpiRank==0){
    cout << "Printing new node SP info: (index, LB, UB)" << endl;
    for(int nn=0; nn<newNodeSPInfo.size(); nn++){
	for(int oo=0; oo<newNodeSPInfo[nn].size(); oo++){
	    cout << " (" << newNodeSPInfo[nn][oo].index << "," << newNodeSPInfo[nn][oo].lb << "," << newNodeSPInfo[nn][oo].ub << ")";
	}
	cout << endl;
    }
  }
}

		
//********************** Penalty update **********************
double getPenalty(){return penC;}
double getBaselinePenalty(){return baselinePenalty_;}
void setPenalty(double p){
    penC=p;
    for (int tS = 0; tS < nNodeSPs; tS++) {
	subproblemSolvers[tS]->setQuadraticTerm(penC);
        for (int i = 0; i < n1; i++) {
	    scaling_matrix[tS][i] = penC;
    	}
    }
//if(mpiRank==0) cout << "Penalty is now: " << penC << endl;
}
//This should normally be turned off, this is a rule for updating the 
//penalty from Kiwiel 2006 and Lubin et al.
void computeKiwielPenaltyUpdate(double SSCVal){	
	penC = 1.0/min( max( (2.0/penC)*(1.0-SSCVal),  max(1.0/(10.0*penC),1e-4)    ), 10.0/penC);
	setPenalty(penC);
}
void computeScalingPenaltyUpdate(double scaling){	
	//penC = 1.0/min( max( (2.0/penC)*(1.0-SSCVal),  max(1.0/(10.0*penC),1e-4)    ), 10.0/penC);
	setPenalty(scaling*penC);
}
void computePenaltyGapVec(){
    for (int tS = 0; tS < nNodeSPs; tS++) {
        for (int i = 0; i < n1; i++) {
	    scaling_matrix[tS][i] *= scaleVec_[tS];
    	}
	subproblemSolvers[tS]->setQuadraticTerm(scaling_matrix[tS]);
    }
}


void displayParameters();
double getBestNodeQuality(){
    double bestNodeQuality = -ALPS_DBL_MAX;
    if(getKnowledgeBroker()->getBestNode()){
        bestNodeQuality = getKnowledgeBroker()->getBestNode()->getQuality();
    }
    return bestNodeQuality;
}
void printStatus(){
  if(mpiRank==0){
cout << "________________________________________________________________________" << endl;
	printf("LagrLB: %0.9g, ALVal: %0.9g, sqrDiscNorm: %0.9g\n", currentLagrLB, ALVal, discrepNorm);
	printf("Best node: %0.9g, Incumbent value: %0.9g\n", getBestNodeQuality(), getIncumbentVal());
cout << "________________________________________________________________________" << endl;
	//printf("Aug. Lagrangian value: %0.9g\n", ALVal);
	//printf("Norm of primal discrepancy: %0.6g\n", discrepNorm);
	//printf("Current penalty: %0.2g\n",penC);
	//std::cout << "Number of integrality infeasibilities of z: " << nIntInfeas_ << std::endl;
  }
}
void printBestBounds(){
  if(mpiRank==0){
cout << "________________________________________________________________________" << endl;
	printf("Best node: %0.9g, Incumbent value: %0.9g\n", getBestNodeQuality(), getIncumbentVal());
cout << "________________________________________________________________________" << endl;
  }
}
void printStatusAsErr(){
  if(mpiRank==0){
	fprintf(stderr, "LagrLB: %0.9g, ALVal: %0.9g, sqrDiscNorm: %0.9g\n", currentLagrLB, ALVal, discrepNorm);
	double bestNodeQuality = -1e20;
	if(getKnowledgeBroker()->getBestNode()){
	  bestNodeQuality = getKnowledgeBroker()->getBestNode()->getQuality();
	}
	fprintf(stderr,"Best node: %0.9g, Incumbent value: %0.9g\n", bestNodeQuality, getIncumbentVal());
	//printf("Aug. Lagrangian value: %0.9g\n", ALVal);
	//printf("Norm of primal discrepancy: %0.6g\n", discrepNorm);
	//printf("Current penalty: %0.2g\n",penC);
	//std::cout << "Number of integrality infeasibilities of z: " << nIntInfeas_ << std::endl;
  }
}
void printZ(double *z=NULL){
  if(mpiRank==0){
    if(z==NULL) z=z_current;
	printf("\nConsensus values for first-stage decision variables:\n[");
	
	for (int i = 0; i < n1; i++) {
		printf("%0.8g ", z[i]);
	}

	printf("]\n");
  }
}
void printIncumbentValue(){
  if(mpiRank==0){
	printf("\nIncumbent value: %0.6g\n", getIncumbentVal());
  }
}
void printZRounded(){
  if(mpiRank==0){
	printf("\nRounded consensus values for first-stage decision variables:\n[");
	
	for (int i = 0; i < n1; i++) {
		printf("%0.6g ", z_rounded[i]);
	}

	printf("]\n");
  }
}
void printIntegralityViolations(){
  if(mpiRank==0){
    for(int ii=0; ii<numIntVars_; ii++){
	if(!checkInteger(z_current[intVar_[ii]])){
	    cout << " (" << intVar_[ii] << "," << z_current[intVar_[ii]] << ")";
	}
    }
    cout << endl;
  }
}

void printNodeStats(){
  if(mpiRank==0){
    cout << "Number of nodes processed: " << getKnowledgeBroker()->getNumNodesProcessed() << endl;;
    cout << "Number of nodes branched: " << getKnowledgeBroker()->getNumNodesBranched() << endl;;
    cout << "Depth of tree: " << getKnowledgeBroker()->getTreeDepth() << endl;
    cout << "Number of nodes discarded: " << getKnowledgeBroker()->getNumNodesDiscarded() << endl;;
    cout << "Number of pregnant nodes: " << getKnowledgeBroker()->getNumNodesPartial() << endl;;
    cout << "Number of nodes left: " << getKnowledgeBroker()->getNumNodeLeftSystem() << endl;;
  }
}

void finalPrintout(){
  if(mpiRank==0){
    cout << "************************************************************************" << endl;
	printf("Final incumbent value: %0.12g\n", getIncumbentVal());
	printZ(z_incumbent_);
    cout << "************************************************************************" << endl;
  }
}


//*** Helper functions.
//Frobenius norm
double froNorm(double** array, int nS, int n1);
double froNorm(double* array, int n1);
//Deals with numerical imprecision from solver (MIPs)
#if 0
double roundIfClose(double input) {
	double tmp = round(input);
	if (abs(tmp-input) < ROUNDING_TOLERANCE)
		{ return tmp; }
	else
		{ return input; }

}
#endif
//*** FWMM functions.
double computeNormedDiscrepancies(double** scaling_matrix, double** x, double* z, int nS, int n1);
//void weightedAverage(const std::vector<double*> &x, const std::vector<double> &p, double *localZ, double* z, int nNodeSPs, int n1, int mpiRank);
void computeFeasibleStartingPoint(int tS, double* x, double* yFeasible);
};
#endif
