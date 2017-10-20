
/*Header for the main procedure ParallelSCG. */

#ifndef PSCGMODEL_H
#define PSCGMODEL_H

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <string>
#include <memory>
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
#include "PSCGScen.h"
#include "PSCGParams.h"

#define SSC_DEN_TOL 1e-10
#define DEFAULT_THREADS 1
#define KIWIEL_PENALTY 1 //set 1 to use Kiwiel (2006) penalty update rule
#define MIN_PEN 0.0 

#ifdef USING_MPI
   #include <mpi.h>
#endif

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#define MAX_INNER_LOOP 200
#define DEFAULT_PENALTY 100
#define KEEP_LOG

// Parameters for newConvergenceCriterion

#define BG_BETA 0.05

enum COL_TYPE{
CONTINUOUS=0,
BINARY,
INTEGER
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

class PSCG {
public:
PSCGParams *par;
TssModel smpsModel;
int algorithm;
vector<int> scenariosToThisModel;
vector<PSCGScen*> subproblemSolvers;
#ifdef KEEP_LOG
  vector< shared_ptr<ofstream> > logFiles;
#endif
//vector< vector<var_branch> > newNodeSPInfo;
int nS;
int n1;
int n2;
int nNodeSPs;
int currentIter_;
int noConseqNullSteps;
int noSeriousSteps;
int phase;
bool shouldContinue;
vector<double> pr;
vector<double*> scaling_matrix; //This should not change after the initial iteration
vector<double*> omega_tilde; //dual variable
vector<double*> omega_current;
vector<double*> omega_centre;
vector<double*> omega_sp;
vector<double*> omega_saved;
bool omegaIsZero_;
bool omegaUpdated_;
bool shouldTerminate;

double SSCVal;
double SSCParam;
//double *scaleVec_;

vector<double*> x_current;
vector<double*> y_current;

double* z_current;// = new double[n1];
double* z_old;// = new double[n1];
double* z_local;// = new double[n1];
double* z_incumbent_; //this should be the last feasible z with best obj
double *z_rounded;
double *z_saved;
double* totalSoln_; //new double[n1+n2];
vector<double> constrVec_; //This is filled with the value of Ax
double *penSumLocal;// = new double[n1];
double *penSum;// = new double[n1];
double *weights_;// = new double[nS];
double *integrDiscr_;

vector<int> spSolverStatuses_;

int modelStatus_[2];

char filepath[64];
char probname[64];
double LagrLB_Local;
double innerLagrLB_Local;
// Bounds
double trialLagrLB;
double currentLagrLB;
double innerLagrLB;
double centreLagrLB;
double referenceLagrLB;
double cutoffLagrLB;
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
double rho;
double baselineRho;
//double penMult;
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
//bool LBcalc;
//bool AlgorithmZ;
//bool disableHeuristic;
int ftype;
int mpiRank;
int mpiSize;

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
//double incObjValue_;
/** Incumbent */

// Termination conditions
int totalNoGSSteps;
char zOptFile[128];

ProblemDataBodur pdBodur;

PSCG(PSCGParams *p):par(p),nNodeSPs(0),algorithm(ALGDD),referenceLagrLB(-COIN_DBL_MAX),cutoffLagrLB(COIN_DBL_MAX),currentLagrLB(-COIN_DBL_MAX),centreLagrLB(-COIN_DBL_MAX),trialLagrLB(-COIN_DBL_MAX),
LagrLB_Local(0.0),ALVal_Local(COIN_DBL_MAX),ALVal(COIN_DBL_MAX),objVal(COIN_DBL_MAX),localDiscrepNorm(1e9),discrepNorm(1e9),
	mpiRank(0),mpiSize(1),totalNoGSSteps(0),infeasIndex_(-1),maxNoSteps(1e6),
	nIntInfeas_(-1),omegaUpdated_(false),SSCParam(0.0),phase(0){

   	//******************Read Command Line Parameters**********************
	//Params par;
	initialiseParameters();

	//******************Display Command Line Parameters**********************
	if (mpiRank==0) {
		displayParameters();
	}

	//******************Reading Problem Data**********************

	initialiseModel();
	


	//******************Assign Scenarios**********************
	assignSubproblems();
	setupSolvers();
	sprintf(zOptFile,"zOpt-%s-noP%dalg%d",probname,mpiSize,algorithm);
	
}

~PSCG(){
	//Clean up structures
	for (int tS = 0; tS < nNodeSPs; tS++) {
		delete [] scaling_matrix[tS];
		delete [] omega_tilde[tS];
		delete [] omega_current[tS];
		delete [] omega_centre[tS];
		delete [] omega_sp[tS];
		delete [] omega_saved[tS];
		delete subproblemSolvers[tS];
	    #ifdef KEEP_LOG
	        if(logFiles.size() > tS) logFiles[tS]->close();
	    #endif
	}
	//delete[] z_current;
	delete [] z_current;
#if 0
	delete [] z_intdisp;
	delete [] z_old;
	delete [] z_saved;
	delete [] z_average;
	delete [] z_vertex_average;
	delete [] z_average1;
	delete [] z_average2;
#endif
	delete [] z_local;
	delete [] z_incumbent_;
	delete [] z_rounded;
	delete [] penSumLocal;
	delete [] penSum;
	delete [] origVarLB_;
	delete [] origVarUB_;
	delete [] currentVarLB_;
	delete [] currentVarUB_;
	delete [] totalSoln_;
	delete [] recordKeeping;
	delete [] colType_;
	delete [] intVar_;
	delete [] weights_;
	delete [] integrDiscr_;
}

void setMPIParams(int rank, int size){
	mpiRank=rank;
	mpiSize=size;
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

	z_current = new double[n1];
	z_old = new double[n1];
	z_saved = new double[n1];
	z_local = new double[n1];
	z_incumbent_ = new double[n1];
	z_rounded = new double[n1];
	penSumLocal = new double[n1];
	penSum = new double[n1];
	origVarLB_ = new double[n1];
	origVarUB_ = new double[n1];
	currentVarLB_ = new double[n1];
	currentVarUB_ = new double[n1];
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
	if (mpiRank==0) {
		std::cout << std::endl;
		std::cout << "Problem data: " << std::endl;
		std::cout << "Number of first stage variables: " << n1 << std::endl;
		std::cout << "Number of second stage variables: " << n2 << std::endl;
		std::cout << "Number of scenarios: " << nS << std::endl;
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
    currentLagrLB=-COIN_DBL_MAX;
    centreLagrLB=-COIN_DBL_MAX;
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
    //zeroOmega();
    currentLagrLB=-COIN_DBL_MAX;
    centreLagrLB=-COIN_DBL_MAX;
    referenceLagrLB=lb;
    setPenalty(pen);
    printOriginalVarBds();
    printCurrentVarBds();
    //subproblemSolvers[0]->printXBounds();
//if(mpiRank==0){cerr << "End installSubproblem ";}
}

double *getZ(){return z_current;}

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

vector<double*>& getOmegaCurrent(){return omega_current;}
vector<double*>& getOmegaCentre(){return omega_centre;}
double *getIncumbentZ(){return z_incumbent_;}
double getObjVal(){return objVal;}

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
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_centre[tS],omega[tS],n1*sizeof(double));
}
void saveOmega(){
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_saved[tS],omega_current[tS],n1*sizeof(double));
}
void loadOmega(){
    omegaIsZero_=false;
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_current[tS],omega_saved[tS],n1*sizeof(double));
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_centre[tS],omega_saved[tS],n1*sizeof(double));
}
void loadOmega(vector<double*> &omega){
    omegaIsZero_=false;
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega[tS],omega_saved[tS],n1*sizeof(double));
}
void zeroOmega(){
    omegaIsZero_=true;
    for(int tS=0; tS<nNodeSPs; tS++){ for(int ii=0; ii<n1; ii++){omega_current[tS][ii]=0.0;}}
    for(int tS=0; tS<nNodeSPs; tS++){ for(int ii=0; ii<n1; ii++){omega_centre[tS][ii]=0.0;}}
}

int initialIteration();

void updateSSC(int k){
    SSCParam = 0.5;
    //SSCParam = k/(100.0+k);
    //if(SSCParam > 0.5) SSCParam=0.5;
}

void updatePenalty(int k){
    //updateSSC(k);
    //if(k % 50 == 0 && omegaUpdated_) computeScalingPenaltyUpdate(2.0);
    //if(currentIter_<=10){
    computeKiwielPenaltyUpdate();
    //}
#if 0
    switch(phase){
	case 0:
	  computeKiwielPenaltyUpdate();
	  break;
	case 1:
#if 0
	  SSCParam = 0.10;
	  for (int tS = 0; tS < nNodeSPs; tS++) {
	    computeScalingPenaltyUpdate(tS,scaleVec_[tS]);
	  }
#endif
	  break;
	case 2:
	  penaltyScale = (double)k;
	  penaltyScale /= (penaltyScale-1.0);
	  //penaltyScale = 1.05;
	  penalty = (double)k;
	  //penalty = pow(1.05,k);
    	  //setPenalty(penalty);
	  computeScalingPenaltyUpdate(penaltyScale);
	  //penalty = pow(1.618, (currentIter_-100.0));
	  //SSCParam = penalty / (10.0*baselineRho + penalty);
    	  if(mpiRank==0) cout << "Penalty: " << penalty << " and SSCParam: " << SSCParam << endl;
	case 3:
	  break; //do nothing
	default:
	  cerr << "Invalid phase: " << phase << endl;
	  break;

    }
#endif
#if 0
    if(currentIter_ < 10) {computeKiwielPenaltyUpdate();}
    else if(omegaUpdated_){ 
	for (int tS = 0; tS < nNodeSPs; tS++) {
	    computeScalingPenaltyUpdate(tS,scaleVec_[tS]);
	}
    }
#endif
}

//Update penalty and SSC parameter as motivated by convergence analysis
#if 0
void updateParams(int k){
    setPenalty(baselineRho*k);
    SSCParam = baselineRho*k / ((9+k)*baselineRho  );
    if(mpiRank==0) cout << "Penalty: " << baselineRho*k << " and SSCParam: " << SSCParam << endl;
}
#endif
void updateVertexHistory(int tS){
	if(subproblemSolvers[tS]->getNumVertices() < nVerticesUsed || nVerticesUsed==0){subproblemSolvers[tS]->addVertex();}
	else{// if(nVerticesUsed > 0){// && subproblemSolvers[tS]->getNumVertices() >= nVerticesUsed) {
		subproblemSolvers[tS]->replaceOldestVertex();
	}
}

int regularIteration(bool adjustPenalty=false, bool SSC=true){
//if(mpiRank==0) cout << "Begin regularIteration()" << endl;
	//if(currentIter_ > 10 && currentIter_<200) preSolveMP();
        double integralityDisc = 100.0;
	if(phase==2) preSolveMP();
	else{
	//else{solveContinuousMPs();}
        int numInnerSolves=1;	
	while(solveContinuousMPs()){
	    //if(mpiRank==0) cout << "***************************************Need to continue with solveMP..." << endl;
	    numInnerSolves++;
	}
	if(mpiRank==0) cout << "Number of inner solve MP calls: " << numInnerSolves << endl;
	}
#if 0
        if(numInnerSolves > 100){ 
	    computeScalingPenaltyUpdate(0.5);
	    if(mpiRank==0) cout << "Reducing penalty by half." << endl;
	}
#endif

	int SPStatus=performColGenStep();
	assert(SPStatus!=SP_INFEAS); //subproblem infeasibility should be caught in initialIteration()
	SSCVal = computeSSCVal(); //shouldTerminate is also updated here.
	if(mpiRank==0) cout << "SSC value is: " << SSCVal << endl;
    //verifyOmegaDualFeas();
    //printStatus();
	//cout << "ALVal: " << ALVal << " discrepNorm: " << discrepNorm << " currentLagrLB " << currentLagrLB << endl;
	updateOmega(SSC);
	if(adjustPenalty) updatePenalty(noSeriousSteps);
//if(mpiRank==0) cout << "End regularIteration()" << endl;
	printStatus();

        if(currentLagrLB >= cutoffLagrLB){
      	    modelStatus_[Z_STATUS]=Z_BOUNDED;
	    if(mpiRank==0){cout << "computeBound(): Terminating due to exceeding cutoff..." << endl;}
            shouldContinue=false;
	}

	integralityDisc=evaluateIntegralityDiscrepancies();
	if((shouldTerminate) || (currentIter_ >= maxNoSteps) ){
	    if(mpiRank==0 && (currentIter_ >= maxNoSteps) ){
		cout << "computeBound(): Terminating due to reaching tolerance criterion at iteration " << currentIter_ << endl;
	    }
	    else{
   		if(mpiRank==0){cout << "computeBound(): Terminating due to reaching maximum number of iterations: " << currentIter_ << endl;}
	    }
            shouldContinue=false;
	}
	if(mpiRank==0) cerr << "Omega updated " << noSeriousSteps << " times, " 
		<< " integrality disc: " << integralityDisc << endl;
	return SPStatus;
}


void setMaxNoSteps(int noSteps){maxNoSteps=noSteps;}

double computeBound(int maxNoConseqNullSteps, bool adjustPenalty=false){
//if(mpiRank==0) cout << "Begin computeBound()" << endl;

    modelStatus_[Z_STATUS]=Z_UNKNOWN;
    clearSPVertexHistory();
    noSeriousSteps=0;
    noConseqNullSteps=0;
    //SSCParam = 0.10;
    #ifdef KEEP_LOG
	for (int tS = 0; tS < nNodeSPs; tS++) {*(logFiles[tS]) << "Initial iteration: " << mpiRank << " scen " << tS << endl;}
    #endif
    int SPStatus=initialIteration();
    if(mpiRank==0) cout << "After initial iteration: currentLB: " << trialLagrLB << " versus refLB: " << referenceLagrLB << endl;;
    
    bool omegaUpdatedAtLeastOnce=false;
    //shouldFathomByOpt = false;
    discrepNorm = 100.0;
    //shouldTerminate=false;
    if(SPStatus==SP_INFEAS){
	modelStatus_[Z_STATUS]=Z_INFEAS;
        shouldContinue=false;
	if(mpiRank==0){cout << "Terminating due to subproblem infeasibility..." << endl;}
    }
    else if(currentLagrLB >= cutoffLagrLB){
	modelStatus_[Z_STATUS]=Z_BOUNDED;
        shouldContinue=false;
	if(mpiRank==0){cout << "Terminating due to exceeding cutoff..." << endl;}
      	//printStatus();
    }
    else{
        shouldContinue=true;
    }//else
    for(currentIter_=0; shouldContinue; currentIter_++){
	if(mpiRank==0) cout << "Regular iteration " << currentIter_ << endl;
	#ifdef KEEP_LOG
	    for(int tS=0; tS<nNodeSPs; tS++){*(logFiles[tS]) << endl << "Regular iteration: " << currentIter_ << endl;}
	#endif

        if(currentIter_ < 10) phase=0;
        else if(currentIter_ < 20) phase=1;
        else phase=2;

	if(phase==0) regularIteration(true,true);
	else if(phase==1) regularIteration(false,true);
	else regularIteration(false,false);


    }// main loop over currentIter_

    solveContinuousMPs();
    printStatus();
    //integralityDisc=evaluateIntegralityDiscrepancies();
    objVal=findPrimalFeasSolnWith(z_current);
    return currentLagrLB;
}

double evaluateIntegralityDiscrepancies(){
    double integrDiscrSum=0.0;
    for(int tS=0; tS<nNodeSPs; tS++){
	integrDiscr_[tS] = subproblemSolvers[tS]->computeIntegralityDiscr();  
	integrDiscrSum += integrDiscr_[tS];  
	//if(integrDiscr > 1e-10) cout << "Node " << mpiRank << " scenario " << tS << " integrality discrepancy: " << integrDiscr << endl;  
    }
    return integrDiscrSum;
}
#if 0
void resetOmega(){
     currentLagrLB=-ALPS_DBL_MAX;
     decimateOmega(0.8);
     setPenalty(baselineRho);
     for (int tS = 0; tS < nNodeSPs; tS++) {
	 recordKeeping[tS][0]=-ALPS_DBL_MAX;//subproblemSolvers[tS]->getLagrBd();
	 recordKeeping[tS][1]=recordKeeping[tS][0];
     }
}
#endif

//Perform more iterations, e.g., to obtain recourse when integrality met.
double findPrimalFeasSolnWith(double *z){
     roundCurrentZ(z);
     if(solveRecourseProblemGivenFixedZ()){
	evaluateFeasibleZ();
     }
     else{
 	if(mpiRank==0) cerr << "findPrimalFeasSoln(): Terminating: Primal feas soln not found." << endl;
     }
     return objVal;
}


double getBound(){return currentLagrLB;}
void setBound(double bd){currentLagrLB=bd;}

void solveForWeights(){
   for (int tS = 0; tS < nNodeSPs; tS++) {
        if (subproblemSolvers[tS]->getNVertices()>2) { //Compute Next X, Y (With History)
	     //subproblemSolvers[tS]->solveMPHistory(omega_current[tS],z_current,scaling_matrix[tS],true);
	     subproblemSolvers[tS]->computeWeightsForCurrentSoln(NULL);
	}
	//if(mpiRank==0){subproblemSolvers[tS]->printWeights();}
   }

}

void preSolveMP(){
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_tilde[tS],omega_centre[tS],n1*sizeof(double));
    //int maxNoGSIts = 100+currentIter_;//max(20,currentIter_/2);
    //int maxNoGSIts = 1;//max(20,currentIter_/2);
    //int maxNoGSIts = min(100*currentIter_,10000000);//max(20,currentIter_/2);
    int maxNoGSIts = 100;//max(20,currentIter_/2);
    for(int itGS=0; itGS < maxNoGSIts; itGS++) { //The inner loop has a fixed number of occurences
	for (int tS = 0; tS < nNodeSPs; tS++) {
		//if(subproblemSolvers[tS]->getNVertices()>0){subproblemSolvers[tS]->solveMPVertices(omega_tilde[tS],z_current,scaling_matrix[tS]);}
		if(subproblemSolvers[tS]->getNVertices()>0){subproblemSolvers[tS]->solveMPHistory(omega_tilde[tS],z_current,NULL,NULL,scaling_matrix[tS],false);}
	}
	updateZ();
    	for (int tS = 0; tS < nNodeSPs; tS++) {
	    for (int i = 0; i < n1; i++) {
		omega_tilde[tS][i] += scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
      	    subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]); //prepares next call of solveMPLineSearch(omega,z,scaling_vector)
	}
    }

}
bool solveContinuousMPs(){
	bool isLastGSIt;
	bool needToContinue = false;
	double zDiff;
        //for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_tilde[tS],omega_centre[tS],n1*sizeof(double));
	//double integrDiscr;
	assert(currentIter_ >= 0);
	//int maxNoGSIts = 20+currentIter_;
	//int maxNoGSIts = 100 + 10*currentIter_;//max(20,currentIter_/2);
	//int maxNoGSIts = min(100*currentIter_,10000000);//max(20,currentIter_/2);
	int maxNoGSIts = 100;//max(20,currentIter_/2);
	for(int itGS=0; itGS < maxNoGSIts; itGS++) { //The inner loop has a fixed number of occurences
	    memcpy(z_old,z_current, n1*sizeof(double));
	    zDiff=0.0;
    	    for (int tS = 0; tS < nNodeSPs; tS++) {
		//*************************** Quadratic subproblem ***********************************
			    
		//if(subproblemSolvers[tS]->getNVertices()>0){subproblemSolvers[tS]->solveMPVertices(omega_centre[tS],z_current,scaling_matrix[tS]);}
		if(subproblemSolvers[tS]->getNVertices()>0){subproblemSolvers[tS]->solveMPHistory(omega_centre[tS],z_current,NULL,NULL,scaling_matrix[tS],false);}
		#ifdef KEEP_LOG
		    //if(itGS==0){subproblemSolvers[tS]->printWeights(logFiles[tS]);}
		#endif
	    }
	    					
	    // Update z_previous.
	    updateZ();
	    totalNoGSSteps++;
	    for(int ii=0; ii<n1; ii++){zDiff=max(zDiff,fabs(z_current[ii]-z_old[ii]));}
	}
#if 1
	innerLagrLB_Local = 0.0;
	ALVal_Local = 0.0;
	localDiscrepNorm = 0.0;
	reduceBuffer[0]=0.0;
	reduceBuffer[1]=0.0;
	reduceBuffer[2]=0.0;
	double refVal, LagrLB_tS, ALVal_tS,sqrDiscrNorm_tS;
    	for (int tS = 0; tS < nNodeSPs; tS++) {
		#ifdef KEEP_LOG
		    //subproblemSolvers[tS]->printWeights(logFiles[tS]);
		#endif
		for (int i = 0; i < n1; i++) {
		    omega_tilde[tS][i] = omega_centre[tS][i] + scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
		}
		refVal = subproblemSolvers[tS]->evaluateSolution(omega_tilde[tS]);
		LagrLB_tS = subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]);
		innerLagrLB_Local += pr[tS]*LagrLB_tS;
	        ALVal_tS = subproblemSolvers[tS]->updateALValues(omega_centre[tS],z_current,scaling_matrix[tS]);
		ALVal_Local += pr[tS]*ALVal_tS;
		sqrDiscrNorm_tS = subproblemSolvers[tS]->getSqrNormDiscr();
		localDiscrepNorm += pr[tS]*sqrDiscrNorm_tS;
		#ifdef KEEP_LOG
		//	*(logFiles[tS]) << "  gap val (master problem): " << LagrLB_tS - refVal;
		//*(logFiles[tS]) << endl;
		#endif
	}
	#ifdef USING_MPI
	if (mpiSize > 1) {
		localReduceBuffer[0]=innerLagrLB_Local;
		localReduceBuffer[1]=ALVal_Local;
		localReduceBuffer[2]=localDiscrepNorm;
		MPI_Allreduce(localReduceBuffer, reduceBuffer, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		innerLagrLB = reduceBuffer[0];
		ALVal = reduceBuffer[1];
		discrepNorm = reduceBuffer[2];
	}
	#endif
	if (mpiSize == 1) {
		innerLagrLB = LagrLB_Local;
		ALVal = ALVal_Local;
		discrepNorm = localDiscrepNorm;
	}
	double innerSSCVal;
    	//if(mpiRank==0) cout << "**************  " << innerLagrLB - centreLagrLB << endl;
    	//if(mpiRank==0) cout << "**************  " << ALVal + 0.5*discrepNorm  - centreLagrLB << endl;
    	if( (ALVal + 0.5*discrepNorm  - centreLagrLB > SSC_DEN_TOL) ){
	    innerSSCVal = (innerLagrLB-centreLagrLB)/(ALVal + 0.5*discrepNorm - centreLagrLB);
	    //if(mpiRank==0) cout << "innerSSCVal: " << innerSSCVal << "  vs.  " << SSCParam << endl;
	    needToContinue = (innerLagrLB-centreLagrLB)/(ALVal + 0.5*discrepNorm - centreLagrLB) < SSCParam;
	}
#endif

//if(mpiRank==0){cout << "solveContinuousMPs(): max zDiff is: " << zDiff << endl;}

    return needToContinue;
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


void updateZ(double *z=NULL,double *dispZ=NULL){
	if(z==NULL) z=z_current;
	for(int tS=0; tS<nNodeSPs; tS++){
	    weights_[tS] = pr[tS];
	}
        averageOfX(z, n1, weights_);
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
void computeAverageAcrossProcs(double *partialAveVec, double* retVec, double weight, int vecSize){
    double weightG=0.0;
    #ifdef USING_MPI
    if(mpiSize>1){
	MPI_Allreduce(partialAveVec, retVec, vecSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&weight, &weightG, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for(int i=0; i<vecSize; i++) retVec[i] /= weightG;
    }
    #endif
    if(mpiSize==1){
        assert(weight > 0.0);
	for (int i=0; i<vecSize; i++) retVec[i] = partialAveVec[i]/weight; //Only one node, trivially set z_local = z
    }
}

//Weighted average of the current_x
void averageOfX(double *retVec, int vecSize, double *weights=NULL){
    double weightSum=0.0;
    for (int i = 0; i < vecSize; i++){
	z_local[i] = 0.0;
    }
    if(weights!=NULL){
      for (int tS = 0; tS < nNodeSPs; tS++){
         for (int i = 0; i < vecSize; i++){
	    z_local[i] += weights[tS]*subproblemSolvers[tS]->getX()[i];
	 }
	 weightSum += weights[tS];
      }
    }
    else{
      for (int tS = 0; tS < nNodeSPs; tS++){
         for (int i = 0; i < n1; i++){
	    z_local[i] += subproblemSolvers[tS]->getX()[i];
	 }
	 weightSum += 1.0;
      }
    }
    computeAverageAcrossProcs(z_local, retVec, weightSum, vecSize);
}
void dispOfX(double *aveVec, double *retVec, int vecSize, double *weights=NULL){
    double weightSum=0.0;
    for(int ii=0; ii<n1; ii++){
	z_local[ii]=0.0;
    }
    if(weights!=NULL){
      for(int tS=0; tS<nNodeSPs; tS++){
        for(int ii=0; ii<n1; ii++){
    	     z_local[ii] += weights[tS]*fabs(aveVec[ii] - x_current[tS][ii]);
	}
	weightSum += weights[tS];
      }
    }
    else{
      for(int tS=0; tS<nNodeSPs; tS++){
        for(int ii=0; ii<n1; ii++){
    	     z_local[ii] += fabs(aveVec[ii] - x_current[tS][ii]);
	}
	weightSum += 1.0;
      }
    }
    computeAverageAcrossProcs(z_local, retVec, weightSum, vecSize);
}
#if 0
void averageOfVertices(double* retVal, double *zDisp=NULL){
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    penSumLocal[i]=0.0;
            penSum[i]=0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += subproblemSolvers[tS]->getXVertex()[i];
		penSumLocal[i]+= 1.0;
	    }
	}
	#ifdef USING_MPI
	    MPI_Allreduce(z_local, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(penSumLocal, penSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    for(int i=0; i<n1; i++) z[i] /= penSum[i];
	#else
	    for (int i=0; i<n1; i++) z[i] = z_local[i]/penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif

	if(zDisp!=NULL){
	   double *zDispLocal = new double[n1]; 
	   for (int i = 0; i < n1; i++){
	     zDispLocal[i] = 0.0;
	     zDisp[i] = 0.0;
	     for (int tS = 0; tS < nNodeSPs; tS++)
	     {
		zDispLocal[i] += fabs(subproblemSolvers[tS]->getXVertex()[i] - z[i]);
	     }
	   }
	#ifdef USING_MPI
	    MPI_Allreduce(zDispLocal, zDisp, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    //for(int i=0; i<n1; i++) zDisp[i] /= penSum[i];
	#else
	    for (int i=0; i<n1; i++) zDisp[i] = zDispLocal[i];///penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif
	   
	   delete [] zDispLocal; 
	}
}
#endif
#if 0
void averageOfVerticesWeightedByRho(double* z=NULL, double *zDisp=NULL){
	//if(z==NULL) z=z_average;
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

	if(zDisp!=NULL){
	   double *zDispLocal = new double[n1]; 
	   for (int i = 0; i < n1; i++){
	     zDispLocal[i] = 0.0;
	     zDisp[i] = 0.0;
	     for (int tS = 0; tS < nNodeSPs; tS++)
	     {
		zDispLocal[i] += scaling_matrix[tS][i]*pr[tS]*fabs(subproblemSolvers[tS]->getXVertex()[i] - z[i]);
	     }
	   }
	#ifdef USING_MPI
	    MPI_Allreduce(zDispLocal, zDisp, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) zDisp[i] = zDispLocal[i];///penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif
	    for(int i=0; i<n1; i++) zDisp[i] /= penSum[i];
	   
	   delete [] zDispLocal; 
	}
}
#endif

#if 0
void averageOfVertices1(double *z=NULL, double *zDisp=NULL){
	//if(z==NULL) z=z_average1;
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
		for(int vv=0; vv<numVertices; vv++){
		//z_local[i] += x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		    //z_local[i] += subproblemSolvers[tS]->getXVertex()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		    z_local[i] += subproblemSolvers[tS]->getXVertexEntry(i,vv);//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		    if(i==0){localTotalWeight += 1.0;}
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


	if(zDisp!=NULL){
	   double *zDispLocal = new double[n1]; 
	   for (int i = 0; i < n1; i++){
	     zDispLocal[i] = 0.0;
	     zDisp[i] = 0.0;
	     for (int tS = 0; tS < nNodeSPs; tS++)
	     {
		numVertices = subproblemSolvers[tS]->getNVertices();
		for(int vv=0; vv<numVertices; vv++){
		    zDispLocal[i] += fabs(subproblemSolvers[tS]->getXVertexEntry(i,vv)-z[i]);//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		}
	     }
	   }
	#ifdef USING_MPI
	    MPI_Allreduce(zDispLocal, zDisp, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) zDisp[i] = zDispLocal[i];///penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif
	   for(int i=0; i<n1; i++) zDisp[i] /= totalWeight;
	   
	   delete [] zDispLocal; 
	}
}
#endif
#if 0
void averageOfVertices2(double *z=NULL, double *zDisp=NULL){
	//if(z==NULL) z=z_average2;
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

	if(zDisp!=NULL){
	   double *zDispLocal = new double[n1]; 
	   for (int i = 0; i < n1; i++){
	     zDispLocal[i] = 0.0;
	     zDisp[i] = 0.0;
	     for (int tS = 0; tS < nNodeSPs; tS++)
	     {
		numVertices = subproblemSolvers[tS]->getNVertices();
		weight = subproblemSolvers[tS]->computeIntegralityDiscr();
		for(int vv=0; vv<numVertices; vv++){
		    zDispLocal[i] += weight*fabs(subproblemSolvers[tS]->getXVertexEntry(i,vv)-z[i]);//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		}
	     }
	   }
	#ifdef USING_MPI
	    MPI_Allreduce(zDispLocal, zDisp, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    //for(int i=0; i<n1; i++) zDisp[i] /= totalWeight;
	#else
	    for (int i=0; i<n1; i++) zDisp[i] = zDispLocal[i];///penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif
	   
	   delete [] zDispLocal; 
	}

}
#endif
#if 0
void averageOfVerticesZ(double *z=NULL, double *zDisp=NULL){
	//if(z==NULL) z=z_intdisp;
        solveForWeights();
	int numVertices;
	//double localTotalWeight=0.0;
	//double totalWeight=0.0;
	double intDiscr;
	double weight;
	//for (int tS = 0; tS < nNodeSPs; tS++){ localTotalNVertices+=subproblemSolvers[tS]->getNVertices();}
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    penSumLocal[i]=0.0;
            penSum[i]=0.0;
	    
	    
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		numVertices = subproblemSolvers[tS]->getNVertices();
		intDiscr = subproblemSolvers[tS]->computeIntegralityDiscr();
		for(int vv=0; vv<numVertices; vv++){
		//zDispLocal[i] += scaling_matrix[tS][i]*pr[tS]*fabs(subproblemSolvers[tS]->getXVertex()[i] - z[i]);
		    //weight = intDiscr*(subproblemSolvers[tS]->getWeight(vv))*pr[tS]*scaling_matrix[tS][i];
		    //weight = intDiscr*(subproblemSolvers[tS]->getWeight(vv));
		    weight = (subproblemSolvers[tS]->getWeight(vv))*pr[tS]*scaling_matrix[tS][i];
		    z_local[i] += weight*subproblemSolvers[tS]->getXVertexEntry(i,vv);//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		    penSumLocal[i] += weight;
		}
	    }
	}
	#ifdef USING_MPI
	    MPI_Allreduce(z_local, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(penSumLocal, penSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) z[i] = z_local[i]; //Only one node, trivially set z_local = z
	    for (int i=0; i<n1; i++) penSum[i] = penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif
	for (int i=0; i<n1; i++) z[i] /= penSum[i]; 

	if(zDisp!=NULL){
	   double *zDispLocal = new double[n1]; 
	   for (int i = 0; i < n1; i++){
	     zDispLocal[i] = 0.0;
	     zDisp[i] = 0.0;
	     for (int tS = 0; tS < nNodeSPs; tS++)
	     {
		numVertices = subproblemSolvers[tS]->getNVertices();
		intDiscr = subproblemSolvers[tS]->computeIntegralityDiscr();
		for(int vv=0; vv<numVertices; vv++){
		    //weight = intDiscr*(subproblemSolvers[tS]->getWeight(vv))*pr[tS]*scaling_matrix[tS][i];
		    //weight = intDiscr*(subproblemSolvers[tS]->getWeight(vv));
		    weight = (subproblemSolvers[tS]->getWeight(vv))*pr[tS]*scaling_matrix[tS][i];
		    zDispLocal[i] += weight*fabs(subproblemSolvers[tS]->getXVertexEntry(i,vv)-z[i]);//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		}
	     }
	   }
	#ifdef USING_MPI
	    MPI_Allreduce(zDispLocal, zDisp, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) zDisp[i] = zDispLocal[i];///penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif
	    //for(int i=0; i<n1; i++) zDisp[i] /= totalWeight;
	   
	   delete [] zDispLocal; 
	}

}
#endif

#if 0
void computeVertexDispAroundX(double *zDisp){
        solveForWeights();
	int numVertices;
	double intDiscr;
	double weight;
	if(zDisp!=NULL){
	   double *zDispLocal = new double[n1]; 
	   for (int i = 0; i < n1; i++){
	     zDispLocal[i] = 0.0;
	     zDisp[i] = 0.0;
	     for (int tS = 0; tS < nNodeSPs; tS++)
	     {
		numVertices = subproblemSolvers[tS]->getNVertices();
		intDiscr = subproblemSolvers[tS]->computeIntegralityDiscr();
		for(int vv=0; vv<numVertices; vv++){
		    //weight = pr[tS]*scaling_matrix[tS][i]*(subproblemSolvers[tS]->getWeight(vv));
		    weight = (subproblemSolvers[tS]->getWeight(vv));
		    zDispLocal[i] += intDiscr*weight*fabs(subproblemSolvers[tS]->getXVertexEntry(i,vv)-x_current[tS][i]);//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		}
	     }
	   }
	#ifdef USING_MPI
	    MPI_Allreduce(zDispLocal, zDisp, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) zDisp[i] = zDispLocal[i];///penSumLocal[i]; //Only one node, trivially set z_local = z
	#endif
	    //for(int i=0; i<n1; i++) zDisp[i] /= totalWeight;
	   
	   delete [] zDispLocal; 
	}

}
#endif

#if 0
void sumWeightedVertexDisp(double *zDisp){
    double intDiscr;
    double *zDispLocal = new double[n1]; 
    for(int ii=0; ii<n1; ii++){zDispLocal[ii]=0.0;}
    for (int tS = 0; tS < nNodeSPs; tS++){
	intDiscr = subproblemSolvers[tS]->computeIntegralityDiscr();
	for(int ii=0; ii<n1; ii++){
	    zDispLocal[ii] += intDiscr*subproblemSolvers[tS]->getDispersions2()[ii];
	}
    }
    #ifdef USING_MPI
        MPI_Allreduce(zDispLocal, zDisp, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #else
        for (int i=0; i<n1; i++) zDisp[i] = zDispLocal[i];///penSumLocal[i]; //Only one node, trivially set z_local = z
    #endif
    delete [] zDispLocal; 

}
#endif

void computeOmegaDisp(double *dispOmega){
    double *dispOmegaLocal = new double[n1];
    for(int ii=0; ii<n1; ii++){
	dispOmegaLocal[ii]=0.0;
	dispOmega[ii]=0.0;
        for(int tS=0; tS<nNodeSPs; tS++){
	    dispOmegaLocal[ii] += fabs(pr[tS]*omega_centre[tS][ii]);

	}
    }
#ifdef USING_MPI
    if(mpiSize>1){
	MPI_Allreduce(dispOmegaLocal, dispOmega, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    if(mpiSize==1){
	memcpy(dispOmega,dispOmegaLocal,n1*sizeof(double));
    }
    delete [] dispOmegaLocal;
}



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
int checkZIsInfeasForScen(int tS){
    memcpy(totalSoln_,z_rounded,n1*sizeof(double));
    double *ySoln = subproblemSolvers[tS]->getY();
    memcpy(totalSoln_+n1,ySoln,n2*sizeof(double));
    bool feas = subproblemSolvers[tS]->checkSolnForFeasibility(totalSoln_, constrVec_);
    if(feas) {return 0;}
    else{return 1;}
}

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

#if 0
bool checkForCutoff(){
if(mpiRank==0){cout << getIncumbentVal() << " <= " << currentLagrLB << "???" << endl;}
    return (getIncumbentVal() <= currentLagrLB);
}
#endif

void roundCurrentZ(double *z=NULL){
    if(z==NULL){z=z_current;}
//double *origVarLB_;
    for(int ii=0; ii<n1; ii++){
	if(z[ii] < origVarLB_[ii]){z[ii]=origVarLB_[ii];}
	if(z[ii] > origVarUB_[ii]){z[ii]=origVarUB_[ii];}
    }
    memcpy(z_rounded,z,n1*sizeof(double));
    for(int ii=0; ii<numIntVars_; ii++) {z_rounded[intVar_[ii]] = round(z[intVar_[ii]]);}
    double roundingDisc = 0.0;
    for(int ii=0; ii<n1; ii++){roundingDisc += fabs(z_rounded[ii]-z[ii]);}
}

bool evaluateFeasibleZ(){
    //if(getIncumbentVal() > objVal){ 
    if(true){ 
if(mpiRank==0){
    cout << "New incumbent value: " << objVal << " and its corresponding solution: " << endl;
    printZRounded();
}
	memcpy( z_incumbent_, z_rounded, n1*sizeof(double) );
	//BcpsSolution* ksol = new BcpsSolution(n1,z_incumbent_,objVal);
        //getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,ksol,objVal);
#if 0
                        new BcpsSolution(model->n1,
                                           model->getZ(),
                                           model->getIncumbentVal());
                    getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution,
                                                       ksol,
                                                       model->getIncumbentVal());
#endif

	//getKnowledgeBroker()->setPhase(AlpsPhaseSearch);
//cout << zOptFile << endl;
        //if(mpiRank==0) getKnowledgeBroker()->printBestSolution(zOptFile);
	//if(mpiRank==0) cout << "Registering solution..." << endl;
//cout << "End evaluateFeasibleZ()" << endl;
	return true;
    }
//cout << "feasSolnUpdated returns false" << endl;
//if(mpiRank==0){cout << "End evaluateFeasibleZ()" << endl;}
    return false;
}
//int solveFeasibilityProblemWithXFixedToZ(const double *z, const double *origLBs, const double *origUBs, const char *colTypes){


bool solveRecourseProblemGivenFixedZ();	
int getSPStatus(){return modelStatus_[SP_STATUS];}
void setSPStatus(int stat){modelStatus_[SP_STATUS]=stat;}
int getZStatus(){return modelStatus_[Z_STATUS];}
void setZStatus(int stat){modelStatus_[Z_STATUS]=stat;}


//********************** Serious Step Condition (SSC) **********************
double computeSSCVal(){
    if(ALVal - centreLagrLB < -SSC_DEN_TOL){
        if(mpiRank==0){cout << " (ALVal,centreLagrLB) (" << setprecision(10) << ALVal << "," << setprecision(10) << centreLagrLB << ")" << endl;}
    }
    if(ALVal + 0.5*discrepNorm  - trialLagrLB < -SSC_DEN_TOL){
        if(mpiRank==0){cout << mpiRank << " (ALVal,trialLagrLB) (" << setprecision(10) << ALVal << "," << setprecision(10) << trialLagrLB << ")" << endl;	
	cout << "regularIteration(): Something probably went wrong with the last computation of trialLagrLB, returning..." << endl;}
    }
    shouldTerminate = (ALVal + 0.5*discrepNorm  - centreLagrLB < SSC_DEN_TOL) && (discrepNorm < 1e-20);
    if(shouldTerminate){return 0.0;}
    else{
	return (trialLagrLB-centreLagrLB)/(ALVal + 0.5*discrepNorm - centreLagrLB);
    }
}
		
bool updateOmega(bool useSSC){
	if(SSCVal >= SSCParam || !useSSC) {
	    for (int tS = 0; tS < nNodeSPs; tS++) {
		memcpy(omega_centre[tS],omega_tilde[tS],n1*sizeof(double));
		memcpy(omega_current[tS],omega_tilde[tS],n1*sizeof(double));
	        recordKeeping[tS][1]=recordKeeping[tS][0];
	        //subproblemSolvers[tS]->updateOptSoln();
	    }
	    centreLagrLB = trialLagrLB;
	    //currentLagrLB = centreLagrLB;
	    //recordKeeping[0]=LagrLB;
    	    omegaUpdated_ = true;
    	    omegaIsZero_=false;
//if(mpiRank==0){cout << "Updating omega..." << endl;}
	}
	else {
    	    omegaUpdated_ = false;
	    //if(mpiRank==0) cout << "Null step taken..." << endl;
	}
	if(omegaUpdated_){ 
            noConseqNullSteps=0;
 	    noSeriousSteps++;
	}
	else{noConseqNullSteps++;}

	if(trialLagrLB > currentLagrLB){
	    currentLagrLB = trialLagrLB;
	    for (int tS = 0; tS < nNodeSPs; tS++) {
		memcpy(omega_current[tS],omega_tilde[tS],n1*sizeof(double));
	    }
	}
	return omegaUpdated_;
}
#if 1

void repairOmega(){
    double *omegaLocal = new double[n1];
    double *omegaSum = new double[n1]; 
    for (int ii = 0; ii < n1; ii++) {
	omegaLocal[ii]=0.0;
    	for (int tS = 0; tS < nNodeSPs; tS++) {
    	    omegaLocal[ii] +=pr[tS]*omega_centre[tS][ii];
	}
    }
#ifdef USING_MPI
    if(mpiSize>1){
	MPI_Allreduce(omegaLocal, omegaSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
    if(mpiSize==1){
	memcpy(omegaSum,omegaLocal,n1*sizeof(double));
    }
#if 0
cout << "Printing omega weighted sum..." << endl;
    for(int ii=0; ii<n1; ii++){
	cout << " " << omegaSum[ii];	
    }
cout << endl;
#endif
//if(mpiRank==0){cout << "Repairing dual feasibility of omega..." << endl;}
    for (int tS = 0; tS < nNodeSPs; tS++) {
      for (int ii = 0; ii < n1; ii++) {
        omega_centre[tS][ii] -= (1.0/((double)nS*pr[tS]))*omegaSum[ii];
      }
    }


}


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
    	    omegaLocal[ii] +=pr[tS]*omega_centre[tS][ii];
	    omegaFroNormLocal += omega_centre[tS][ii]*omega_centre[tS][ii];
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

#if 0
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
#endif

		
//********************** Penalty update **********************
double getPenalty(){return rho;}
double getBaselinePenalty(){return baselineRho;}
void setPenalty(double p){
    rho=max(p,baselineRho);
    for (int tS = 0; tS < nNodeSPs; tS++) {
	setPenalty(tS,rho);
    }
//if(mpiRank==0) cout << "Penalty is now: " << rho << endl;
}
void setPenalty(int tS, double p){
    double penalty = max(p,baselineRho);
    subproblemSolvers[tS]->setQuadraticTerm(penalty);
    for (int i = 0; i < n1; i++) {
        scaling_matrix[tS][i] = penalty;
    }
//if(mpiRank==0) cout << "Penalty is now: " << rho << endl;
}
//This should normally be turned off, this is a rule for updating the 
//penalty from Kiwiel 2006 and Lubin et al.
void computeKiwielPenaltyUpdate(){	
	//rho = 1.0/min( max( (2.0/rho)*(1.0-SSCVal),  max(1.0/(10.0*rho),1e-4)    ), 10.0/rho);
	//setPenalty(rho);
    	//if( SSCVal >= 0.5 ){ computeScalingPenaltyUpdate(min( 0.5/(1-SSCVal),2.0));}
	double limit = 100.0;
    	if(fabs(1.0-SSCVal) > 1e-10){
	     computeScalingPenaltyUpdate( max( min( 0.5/(1.0-SSCVal),limit), 1.0/limit));
	}
	else{
    	    computeScalingPenaltyUpdate( limit );
	}
}
void computeScalingPenaltyUpdate(double scaling){	
    for (int tS = 0; tS < nNodeSPs; tS++) {
	computeScalingPenaltyUpdate(tS,scaling);
#if 0
        for (int i = 0; i < n1; i++) {
	    scaling_matrix[tS][i] *= scaling;
	    scaling_matrix[tS][i] = max( scaling_matrix[tS][i], MIN_PEN);
    	}
	subproblemSolvers[tS]->setQuadraticTerm(scaling_matrix[tS]);
#endif
    }
}
void computeScalingPenaltyUpdate(int tS, double scaling){	
    for (int i = 0; i < n1; i++) {
	scaling_matrix[tS][i] *= scaling;
	scaling_matrix[tS][i] = max( scaling_matrix[tS][i], baselineRho);
    }
    subproblemSolvers[tS]->setQuadraticTerm(scaling_matrix[tS]);
}
#if 0
void computePenaltyGapVec(){
    for (int tS = 0; tS < nNodeSPs; tS++) {
    computeScalingPenaltyUpdate(tS,scaleVec_[tS]);
#if 0
        for (int i = 0; i < n1; i++) {
	    scaling_matrix[tS][i] *= scaleVec_[tS];
	    scaling_matrix[tS][i] = max( scaling_matrix[tS][i], MIN_PEN);
    	}
	subproblemSolvers[tS]->setQuadraticTerm(scaling_matrix[tS]);
#endif
    }
}
#endif


void displayParameters();
#if 0
double getBestNodeQuality(){
    double bestNodeQuality = -ALPS_DBL_MAX;
    if(getKnowledgeBroker()->getBestNode()){
        bestNodeQuality = getKnowledgeBroker()->getBestNode()->getQuality();
    }
    return bestNodeQuality;
}
#endif
void printStatus(){
  if(mpiRank==0){
cout << "________________________________________________________________________" << endl;
	printf("currentBestLagrLB: %0.14g, ALVal: %0.14g, sqrDiscNorm: %0.14g\n", currentLagrLB, ALVal, discrepNorm);
//	printf("Best node: %0.14g, Incumbent value: %0.14g\n", getBestNodeQuality(), getIncumbentVal());
cout << "________________________________________________________________________" << endl;
	//printf("Aug. Lagrangian value: %0.9g\n", ALVal);
	//printf("Norm of primal discrepancy: %0.6g\n", discrepNorm);
	//printf("Current penalty: %0.2g\n",rho);
	//std::cout << "Number of integrality infeasibilities of z: " << nIntInfeas_ << std::endl;
  }
}
void printZ(double *z=NULL){
  if(mpiRank==0){
    if(z==NULL) z=z_current;
	printf("\nConsensus values for first-stage decision variables:\n[");
	
	for (int i = 0; i < n1; i++) {
		printf("%0.14g ", z[i]);
	}

	printf("]\n");
  }
}
void printZRounded(){
  if(mpiRank==0){
	printf("\nRounded consensus values for first-stage decision variables:\n[");
	
	for (int i = 0; i < n1; i++) {
		printf("%0.14g ", z_rounded[i]);
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


//*** Helper functions.
//Frobenius norm
double froNorm(double** array, int nS, int n1);
double froNorm(double* array, int n1);
//Deals with numerical imprecision from solver (MIPs)
double computeNormedDiscrepancies(double** scaling_matrix, double** x, double* z, int nS, int n1);
void computeFeasibleStartingPoint(int tS, double* x, double* yFeasible);
};
#endif
/*Header for the main procedure ParallelSCG. */


