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

#define smallNumber 0.000001
#define SSC_PARAM 0.1
#define SSC_DEN_TOL 1e-4
#define ROUNDING_TOLERANCE 1e-6
#define DEFAULT_THREADS 1
#define KIWIEL_PENALTY 1 //set 1 to use Kiwiel (2006) penalty update rule
#define MAX_N_BRANCHINGS 1

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
vector<double*> omega_saved;

double *scaleVec_;

vector<double*> x_current;
vector<double*> y_current;

double* z_current;// = new double[n1];
double* z_average;
double* z_local;// = new double[n1];
double* z_incumbent_; //this should be the last feasible z with best obj
double *z_rounded;
double* totalSoln_; //new double[n1+n2];
vector<double> constrVec_; //This is filled with the value of Ax

vector<int> spSolverStatuses_;

int modelStatus_[2];
bool omegaUpdated_;

char filepath[64];
char probname[64];
double LagrLB_Local;
// Bounds
double LagrLB;
double currentLagrLB;
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
double branchingZVals[MAX_N_BRANCHINGS];

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

PSCGModel(PSCGParams *p):par(p),nNodeSPs(0),currentLagrLB(-ALPS_DBL_MAX),LagrLB(-ALPS_DBL_MAX),
LagrLB_Local(0.0),ALVal_Local(ALPS_DBL_MAX),ALVal(ALPS_DBL_MAX),objVal(ALPS_DBL_MAX),localDiscrepNorm(1e9),discrepNorm(1e9),
	mpiRank(0),mpiSize(1),mpiHead(true),totalNoGSSteps(0),infeasIndex_(-1),nIntInfeas_(-1),omegaUpdated_(false),scaleVec_(NULL),BcpsModel(){

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
	sprintf(zOptFile,"zOpt-%s-noP%dalgZ%d",probname,mpiSize,AlgorithmZ);
	
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
	delete [] z_average;
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
	z_average = new double[n1];
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
	    memcpy(origVarLB_,origVarLB_,n1*sizeof(double));
	    memcpy(origVarUB_,origVarUB_,n1*sizeof(double));
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
void installSubproblem(double lb, vector<int> &indices, vector<int> &fixTypes, vector<double> &bds, int nFix, int infeasIndex, double pen){
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
    //readOmegaIntoModel(omega);
    loadOmega();
    currentLagrLB=lb;
    infeasIndex_=infeasIndex;
    setPenalty(pen);
    printOriginalVarBds();
    printCurrentVarBds();
    //subproblemSolvers[0]->printXBounds();
//if(mpiRank==0){cerr << "End installSubproblem ";}
}


void installSubproblem(double lb, const double *zLBs, const double *zUBs, int infeasIndex, double pen){
//if(mpiRank==0){cerr << "Begin installSubproblem ";}
    setLBsAllSPs(zLBs);
    setUBsAllSPs(zUBs);
    //restoreOriginalVarBounds();
//if(mpiRank==0){cout << "Branching at (indices,bounds,type): ";}
#if 0
    for(int ii=0; ii<n1; ii++){
//if(mpiRank==0){cout << "(" << indices[ii] << "," << bds[ii] << "," << fixTypes[ii] << ") ";}
	setLBsAllSPs(const double *lbs){
	upBranchAllSPsAt(ii,bds[ii]);
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
#endif
//if(mpiRank==0){cout << endl;}
    //TODO: resolve ownership questions for current_z, current_omega
//if(mpiRank==0){cerr << "Reading z...";}
    //readZIntoModel(z);
//if(mpiRank==0){cerr << "Reading omega...";}
    //readOmegaIntoModel(omega);
    loadOmega();
    currentLagrLB=lb;
    //infeasIndex_=infeasIndex;
    setPenalty(pen);
    printOriginalVarBds();
    printCurrentVarBds();
    //subproblemSolvers[0]->printXBounds();
//if(mpiRank==0){cerr << "End installSubproblem ";}
}

double *getZ(){return z_current;}
double *getZAverage(){return z_average;}
double findBranchingIndex(double *z=NULL, int branchingNo=0){
//if(z==NULL) z=z_current;
printTwoZDiscr();
    if(AlgorithmZ){
        if(mpiRank==0){ cout << "Using the consensus vector z to determine branching..." << endl;}
	if(z==NULL) z=z_current;
        //solveContinuousMPs(); //already done
    }
    else{
	if(mpiRank==0){ cout << "Using the average of scenario-specific optimal columns to determine branching..." << endl;}	
	if(z==NULL) z=z_average;
	//averageOfOptVertices();
    }
    infeasIndex_=-1;
    numIntInfeas(z,branchingNo);
    double maxDisp=ALPS_DBL_MAX;
    if(infeasIndex_==-1){
	maxDisp = evaluateXDispersion(z,branchingNo);
    }	
    if(maxDisp < 1e-5){infeasIndex_=-1;} //z solution is regarded as feasible, should update z status to feas or opt.
    return maxDisp;
}


double evaluateXDispersion(double *z=NULL, int branchingNo=0){
    if(z==NULL) z=z_current; //Default behaviour
    double maxDisp=ALPS_DBL_MAX;
    double *dispLocal = new double[n1];
    double *disp = new double[n1];
    double wSum;
    double cii;
    double discrep;
    int nVertices,nVerticesAdded;
    for(int ii=0; ii<n1; ii++){
	dispLocal[ii]=0.0;
	disp[ii]=0.0;
	//wSum=0.0;
        for(int tS=0; tS<nNodeSPs; tS++){
	    //dispLocal[ii] += pr[tS]*(subproblemSolvers[tS]->getXVertex()[ii] - z[ii])*scaling_matrix[tS][ii]*(subproblemSolvers[tS]->getXVertex()[ii] - z[ii]);
	    //dispLocal[ii] += pr[tS]*(subproblemSolvers[tS]->getXVertexOpt()[ii] - z[ii])*scaling_matrix[tS][ii]*(subproblemSolvers[tS]->getXVertexOpt()[ii] - z[ii]);
	    if(!(subproblemSolvers[tS]->getXVertexOpt()[ii] >= currentVarLB_[ii] && subproblemSolvers[tS]->getXVertexOpt()[ii] <= currentVarUB_[ii])){
cout << "Scenario " << tS << " index " << ii << " with value " << subproblemSolvers[tS]->getXVertexOpt()[ii] << endl;
	    assert(subproblemSolvers[tS]->getXVertexOpt()[ii] >= currentVarLB_[ii] && subproblemSolvers[tS]->getXVertexOpt()[ii] <= currentVarUB_[ii]);
	    }
	//cout << " (" << currentVarLB_[ii] << "," << currentVarUB_[ii] << ")";
	    cii = subproblemSolvers[tS]->getC()[ii];
#if 0
	    discrep = 0.0;
	    nVertices = subproblemSolvers[tS]->getNVertices();
	    nVerticesAdded=0;

	    for(int vv=nVertices-1; vv>=0 && vv>=nVertices-20; vv--){ 
		discrep += fabs(subproblemSolvers[tS]->getXVertexEntry(ii,vv) -z[ii]); 
		nVerticesAdded++;
	    }
	    discrep /= ((double)nVerticesAdded);
	    //discrep = subproblemSolvers[tS]->getXVertexOpt()[ii] - z[ii];
//double getXVertexEntry(int entry, int vertex){
//int getNVertices(){return nVertices;}
#endif
	    //dispLocal[ii] += pr[tS]*fabs(cii*discrep);
	    dispLocal[ii] += pr[tS]*(cii+omega_current[tS][ii]);
	    //wSum += pr[tS]*scaling_matrix[tS][ii];
        }
	//cout << endl;
//cerr << "Rank: " << mpiRank << " wSum: " << wSum << endl;
	//dispLocal[ii] /= wSum;
    }
#ifdef USING_MPI
    if(mpiSize>1){
//if(mpiRank==0) cerr << "evaluateXDispersion(): before MPI_Allreduce()" << endl;
//cerr << "#";
	MPI_Allreduce(dispLocal, disp, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//if(mpiRank==0) cerr << "evaluateXDispersion(): after MPI_Allreduce()" << endl;
    }
#endif
    if(mpiSize==1){
	memcpy(disp,dispLocal,n1*sizeof(double));;
    }
    for(int ii=0; ii<n1; ii++){disp[ii] = fabs(disp[ii]);}
    if(mpiRank==0){ 
	cout << "Dispersion given z is: " << endl;
	for(int ii=0; ii<n1; ii++){ cout << " " << disp[ii];}
	cout << endl;
    }
	maxDisp=0.0;
    	for(int ii=0; ii<n1; ii++){ 
	    if(disp[ii] > maxDisp && colType_[ii]=='C'){
	    	maxDisp = disp[ii];
	    	infeasIndex_ = ii;
		branchingZVals[branchingNo]=z[infeasIndex_];
	    }
    	}
    delete [] dispLocal;
    delete [] disp;
    return maxDisp;
}

double* getBranchingZVals(){return &branchingZVals[0];}
double* getOrigVarLbds(){return origVarLB_;}
double* getOrigVarUbds(){return origVarUB_;}
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
void readOmegaIntoModel(vector<double*> &omega){
//cout << "nNodeSPs " << nNodeSPs << " size of vec " << omega.size() << endl;
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_current[tS],omega[tS],n1*sizeof(double));
}
void loadOmega(){
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_current[tS],omega_saved[tS],n1*sizeof(double));
}
void saveOmega(){
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_saved[tS],omega_current[tS],n1*sizeof(double));
}

#if 1
virtual AlpsTreeNode* createRoot();
#endif

int initialIteration();

int regularIteration(bool scalePenalty=false, bool adjustPenalty=false){
//if(mpiRank==0) cout << "Begin regularIteration()" << endl;
	solveContinuousMPs();

	int SPStatus=performColGenStep();
	assert(SPStatus!=SP_INFEAS); //subproblem infeasibility should be caught in initialIteration()
  if(mpiRank==0){
    if(ALVal - currentLagrLB < -1e-6){
        cout << " (ALVal,currentLagrLB) (" << ALVal << "," << currentLagrLB << ")" << endl;	
    }
    if(ALVal + 0.5*discrepNorm  - LagrLB < -1e-6){
        cout << mpiRank << " (ALVal,LagrLB) (" << ALVal << "," << LagrLB << ")" << endl;	
    }
  }
    //verifyOmegaDualFeas();
    //printStatus();
	//cout << "ALVal: " << ALVal << " discrepNorm: " << discrepNorm << " currentLagrLB " << currentLagrLB << endl;
	double SSCVal = computeSSCVal();
	updateOmega(SSCVal);
	#if KIWIEL_PENALTY 
	 if(adjustPenalty) computeKiwielPenaltyUpdate(SSCVal);
	#endif
	 if(omegaUpdated_) computePenaltyGapVec();
//printIntegralityViolations();
//if(mpiRank==0) cout << "End regularIteration()" << endl;
	return SPStatus;
}

double computeBound(int nIters, bool adjustPenalty=false){
//if(mpiRank==0) cout << "Begin computeBound()" << endl;
    modelStatus_[Z_STATUS]=Z_UNKNOWN;
    clearSPVertexHistory();
    int SPStatus=initialIteration();
    if(SPStatus==SP_INFEAS){
	modelStatus_[Z_STATUS]=Z_INFEAS;
if(mpiRank==0){cout << "Terminating due to subproblem infeasibility..." << endl;}
	//return currentLagrLB;
    }
    else if(currentLagrLB >= getIncumbentVal()){
        //updateModelStatusZ();
	modelStatus_[Z_STATUS]=Z_BOUNDED;
if(mpiRank==0){cout << "Terminating due to exceeding cutoff..." << endl;}
      	printStatus();
	//return currentLagrLB;
    }
    else{
    //int maxNLoopIters = 1;
    //int loopIter = 0;
    //do
    //{
      //if(loopIter>0) {computeScalingPenaltyUpdate(10.0);}
      for(int ii=0; ii<nIters; ii++){
if(mpiRank==0) cerr << "Regular iteration " << ii << endl;
	if(ii<10) {regularIteration(false,adjustPenalty);}
	else {regularIteration(true,false);}
      	printStatusAsErr();
	//verifyOmegaDualFeas();
	//printZ();
        if(currentLagrLB >= getIncumbentVal()){
      	    modelStatus_[Z_STATUS]=Z_BOUNDED;
if(mpiRank==0){cout << "Terminating due to exceeding cutoff..." << endl;}
      	    printStatus();
	    //return currentLagrLB;
	    break;
	}
#if 1
	if(shouldTerminate()){
if(mpiRank==0){cout << "Terminating due to reaching tolerance criterion..." << endl;}
      	    //updateModelStatusZ();
      	    printStatus();
	    //return currentLagrLB;
	    break;
	}
#endif
      }
      //assert(SPStatus!=SP_INFEAS);
      //assert(modelStatus_[Z_STATUS]==Z_UNKNOWN);
      printStatus();
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
    }//else
    saveOmega();
    solveContinuousMPs();
    averageOfOptVertices();
    roundCurrentZ();
    return currentLagrLB;
}

//Perform more iterations, e.g., to obtain recourse when integrality met.
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
double findPrimalFeasSoln(int nIters){
//if(mpiRank==0){cout << "Before z: " << endl;}
//printZ();
    for(int ii=0; ii<numIntVars_; ii++){fixVarAllSPsAt(intVar_[ii], z_rounded[intVar_[ii]]);}
    clearSPVertexHistory();
    int SPStatus=initialIteration();

    if(SPStatus==SP_INFEAS){
	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Primal feas soln not found." << endl;
	//return currentLagrLB;
    }
    else if(currentLagrLB >= getIncumbentVal()){
	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Any primal feas soln to be found would not be better in value than current incumbent value." << endl;
      	printStatus();
	//return currentLagrLB;
    }
    else{
     for(int ii=0; ii<nIters; ii++){
if(mpiRank==0) cerr << "Regular iteration " << ii << endl;
	if(nIters<20) {regularIteration(false,true);}
	else {regularIteration(false,false);}
      	printStatusAsErr();

        if(currentLagrLB >= getIncumbentVal()){
	    if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Any primal feas soln to be found would not be better in value than current incumbent value." << endl;
      	    //printStatus();
	    //return currentLagrLB;
	    break;
	}
#if 1
	if(shouldTerminate()){
	    if(mpiRank==0){cout << "findPrimalFeasSoln(): Terminating due to reaching tolerance criterion..." << endl;}
            //printStatus();
      	    //updateModelStatusZ();
      	    #if 0
	    if(solveRecourseProblemGivenFixedZ()){
	 	if(mpiRank==0) cout << "\tfindPrimalFeasSoln(): Terminating: Found primal feas soln..." << endl;
		evaluateFeasibleZ();
	    }
	    else{
	 	if(mpiRank==0) cout << "\tfindPrimalFeasSoln(): Terminating: Primal feas soln not found." << endl;
	    }
	    #endif
	    //return currentLagrLB;
	    break;
	}
#endif
     }
     //assert(SPStatus!=SP_INFEAS);
     printStatus();
     if(solveRecourseProblemGivenFixedZ()){
	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Found primal feas soln..." << endl;
	evaluateFeasibleZ();
     }
     else{
 	if(mpiRank==0) cout << "findPrimalFeasSoln(): Terminating: Primal feas soln not found." << endl;
     }
    }//else
    return currentLagrLB;
}

double getBound(){return currentLagrLB;}
void setBound(double bd){currentLagrLB=bd;}

void solveContinuousMPs(){
	for(int itGS=0; itGS < fixInnerStep; itGS++) { //The inner loop has a fixed number of occurences
    	    for (int tS = 0; tS < nNodeSPs; tS++) {

		//*************************** Quadratic subproblem ***********************************
			    
		    if (useVertexHistory) { //Compute Next X, Y (With History)

					
			//Solve the quadratic master problem for x and y
			subproblemSolvers[tS]->updatePrimalVariablesHistory_OneScenario(omega_current[tS],z_current,scaling_matrix[tS]);
			
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
    	for (int tS = 0; tS < nNodeSPs; tS++) {
	    subproblemSolvers[tS]->updateALValues(omega_current[tS],z_current,scaling_matrix[tS]);
	    //if(tS==0) subproblemSolvers[tS]->printXVertices();
	}
}
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

int performColGenStep();


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
//if(mpiRank==0) cerr << "updateZ(): z_local: before MPI_Allreduce()" << endl;
//cerr << "#";
	    MPI_Allreduce(z_local, z_current, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//if(mpiRank==0) cerr << endl << "updateZ(): z_local: after MPI_Allreduce()" << endl;
//if(mpiRank==0) cerr << "updateZ(): penSum: before MPI_Allreduce()" << endl;
//cerr << "#";
	    MPI_Allreduce(penSumLocal, penSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//if(mpiRank==0) cerr << endl <<  "updateZ(): penSum: after MPI_Allreduce()" << endl;
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

void averageOfOptVertices(bool roundZ=false){
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z_average[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		//z_local[i] += x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		z_local[i] += subproblemSolvers[tS]->getXVertexOpt()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
	    }
	}
	#ifdef USING_MPI
//if(mpiRank==0) cerr << "updateZ(): z_local: before MPI_Allreduce()" << endl;
//cerr << "#";
	    MPI_Allreduce(z_local, z_average, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//if(mpiRank==0) cerr << endl << "updateZ(): z_local: after MPI_Allreduce()" << endl;
//if(mpiRank==0) cerr << "updateZ(): penSum: before MPI_Allreduce()" << endl;
//cerr << "#";
	#else
	    for (int i=0; i<n1; i++) z_average[i] = z_local[i]; //Only one node, trivially set z_local = z_average
	#endif
	for (int i=0; i<n1; i++) z_average[i] /= ((double)nS); //Only one node, trivially set z_local = z_average
	if(roundZ){//rounding integer restricted vars to nearest integer
	  for(int ii=0; ii< numIntVars_; ii++){
	    z_average[intVar_[ii]]=round(z_average[intVar_[ii]]);
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

/** Check if a value is integer. */
bool checkInteger(double value) const {
    double integerTolerance = 1.0e-6;
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
}

bool evaluateFeasibleZ(){
//if(mpiRank==0){cout << "Begin evaluateFeasibleZ()" << endl;}
    //assert(modelStatus_[Z_STATUS]==Z_FEAS || modelStatus_[Z_STATUS]==Z_OPT || modelStatus_[Z_STATUS]==Z_INFEAS);
    //assert(modelStatus_[Z_STATUS]!=Z_REC_INFEAS);
//cout << "incumbent >??? currentLagrLB " << incumbentVal_ << " ??? " << currentLagrLB << endl;
    //if(incumbentVal_ > currentLagrLB){ 
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
bool shouldTerminate(){
//if(mpiRank==0){cout << "shouldTerminate(): " << (ALVal + 0.5*discrepNorm  - currentLagrLB) << endl;}
    return (ALVal + 0.5*discrepNorm  - currentLagrLB < SSC_DEN_TOL);
}
double computeSSCVal(){
    return (LagrLB-currentLagrLB)/(ALVal + 0.5*discrepNorm - currentLagrLB);
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
//if(mpiRank==0){cout << "Updating omega..." << endl;}
	}
	else {
    	    omegaUpdated_ = false;
	    //if(mpiRank==0) cout << "Null step taken..." << endl;
	}
	return omegaUpdated_;
}
#if 0
bool verifyOmegaDualFeas(){
    double *omegaLocal = new double[n1];
    double *omegaSum = new double[n1]; 
    for (int ii = 0; ii < n1; ii++) {
	omegaLocal[ii]=0.0;
    	for (int tS = 0; tS < nNodeSPs; tS++) {
    	    omegaLocal[ii] +=pr[tS]*omega_current[tS][ii];
	}
    }
#ifdef USING_MPI
    if(mpiSize>1){
	MPI_Allreduce(omegaLocal, omegaSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    if(mpiSize==1){
	memcpy(omegaSum,omegaLocal,n1*sizeof(double));;
    }
if(mpiRank==0){
cout << "Printing omega weighted sum..." << endl;
    for(int ii=0; ii<n1; ii++){
	cout << " " << roundIfClose(omegaSum[ii]);	
    }
cout << endl;
}
    delete [] omegaLocal;
    delete [] omegaSum;
}
#endif

		
//********************** Penalty update **********************
double getPenalty(){return penC;}
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
void printStatus(){
  if(mpiRank==0){
	printf("LagrLB: %0.9g, ALVal: %0.9g, sqrDiscNorm: %0.9g\n", currentLagrLB, ALVal, discrepNorm);
	double bestNodeQuality = -1e20;
	if(getKnowledgeBroker()->getBestNode()){
	  bestNodeQuality = getKnowledgeBroker()->getBestNode()->getQuality();
	}
	printf("Best node: %0.9g, Incumbent value: %0.9g\n", bestNodeQuality, getIncumbentVal());
	//printf("Aug. Lagrangian value: %0.9g\n", ALVal);
	//printf("Norm of primal discrepancy: %0.6g\n", discrepNorm);
	//printf("Current penalty: %0.2g\n",penC);
	//std::cout << "Number of integrality infeasibilities of z: " << nIntInfeas_ << std::endl;
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
void printZ(){
  if(mpiRank==0){
	printf("\nConsensus values for first-stage decision variables:\n[");
	
	for (int i = 0; i < n1; i++) {
		printf("%0.5g ", roundIfClose(z_current[i]));
	}

	printf("]\n");
  }
}
void printZRounded(){
  if(mpiRank==0){
	printf("\nRounded consensus values for first-stage decision variables:\n[");
	
	for (int i = 0; i < n1; i++) {
		printf("%0.5g ", z_rounded[i]);
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
    cout << "Number of nodes discarded: " << getKnowledgeBroker()->getNumNodesDiscarded() << endl;;
    cout << "Number of pregnant nodes: " << getKnowledgeBroker()->getNumNodesPartial() << endl;;
    cout << "Number of nodes left: " << getKnowledgeBroker()->getNumNodeLeftSystem() << endl;;
  }
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
