
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
#include "PSCGScen.h"
#include "PSCGParams.h"

#define SSC_PARAM 0.50
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

class PSCGModel {
public:
PSCGParams *par;
TssModel smpsModel;
int algorithm;
vector<int> scenariosToThisModel;
vector<PSCGModelScen*> subproblemSolvers;
//vector< vector<var_branch> > newNodeSPInfo;
int nS;
int n1;
int n2;
int nNodeSPs;
int currentIter_;
vector<double> pr;
vector<double*> scaling_matrix; //glorified rho
vector<double*> omega_tilde; //dual variable
vector<double*> omega_current;
vector<double*> omega_sp;
vector<double*> omega_saved;
bool omegaIsZero_;
bool omegaUpdated_;
bool shouldTerminate;

double SSCVal;
double *scaleVec_;

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

PSCGModel(PSCGParams *p):par(p),nNodeSPs(0),algorithm(ALGDD),referenceLagrLB(-COIN_DBL_MAX),currentLagrLB(-COIN_DBL_MAX),LagrLB(-COIN_DBL_MAX),
LagrLB_Local(0.0),ALVal_Local(COIN_DBL_MAX),ALVal(COIN_DBL_MAX),objVal(COIN_DBL_MAX),localDiscrepNorm(1e9),discrepNorm(1e9),
	mpiRank(0),mpiSize(1),mpiHead(true),totalNoGSSteps(0),infeasIndex_(-1),maxNoSteps(20),
	nIntInfeas_(-1),omegaUpdated_(false),scaleVec_(NULL){

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

~PSCGModel(){
	//Clean up structures
	for (int tS = 0; tS < nNodeSPs; tS++) {
		delete [] scaling_matrix[tS];
		delete [] omega_tilde[tS];
		delete [] omega_current[tS];
		delete [] omega_sp[tS];
		delete [] omega_saved[tS];
		delete subproblemSolvers[tS];
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
	z_old = new double[n1];
	z_saved = new double[n1];
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

vector<double*>& getOmega(){return omega_current;}
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
}
void saveOmega(){
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_saved[tS],omega_current[tS],n1*sizeof(double));
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

int initialIteration();

void updatePenalty(){
    //if((omegaUpdated_ || currentIter_ < 10)) computeKiwielPenaltyUpdate(SSCVal);
    //if(omegaUpdated_) computeKiwielPenaltyUpdate();
    //if(omegaUpdated_) computeScalingPenaltyUpdate( 2.0 );
    if(omegaUpdated_) computePenaltyGapVec();
}

int regularIteration(bool adjustPenalty=false, bool SSC=true){
//if(mpiRank==0) cout << "Begin regularIteration()" << endl;
	solveContinuousMPs();

	int SPStatus=performColGenStep();
	assert(SPStatus!=SP_INFEAS); //subproblem infeasibility should be caught in initialIteration()
	SSCVal = computeSSCVal(); //shouldTerminate is also updated here.
	if(mpiRank==0) cout << "SSC value is: " << SSCVal << endl;
    //verifyOmegaDualFeas();
    //printStatus();
	//cout << "ALVal: " << ALVal << " discrepNorm: " << discrepNorm << " currentLagrLB " << currentLagrLB << endl;
	updateOmega(SSC);
	//#if KIWIEL_PENALTY 
	 //if(adjustPenalty) computeKiwielPenaltyUpdate();
	 if(adjustPenalty) updatePenalty();
	//#endif
	 //if(omegaUpdated_ && scalePenalty) computePenaltyGapVec();
	 //if(scalePenalty) computePenaltyGapVec();
//printIntegralityViolations();
//if(mpiRank==0) cout << "End regularIteration()" << endl;
	return SPStatus;
}
void setMaxNoSteps(int noSteps){maxNoSteps=noSteps;}

double computeBound(int maxNoConseqNullSteps, bool adjustPenalty=false){
//if(mpiRank==0) cout << "Begin computeBound()" << endl;
    //fixInnerStep = 20;
    modelStatus_[Z_STATUS]=Z_UNKNOWN;
    clearSPVertexHistory();
    int maxNoIts = 100;
    int SPStatus=initialIteration();
    if(mpiRank==0) cout << "After initial iteration: currentLB: " << LagrLB << " versus refLB: " << referenceLagrLB << endl;;
    
    //bool exceedingReferenceBd = currentLagrLB + SSC_DEN_TOL >= referenceLagrLB;
    bool exceedingReferenceBd = false;
    bool omegaUpdatedAtLeastOnce=false;
    //shouldFathomByOpt = false;
    double integralityDisc = 100.0;
    discrepNorm = 100.0;
    shouldTerminate=false;
    if(SPStatus==SP_INFEAS){
	modelStatus_[Z_STATUS]=Z_INFEAS;
if(mpiRank==0){cout << "Terminating due to subproblem infeasibility..." << endl;}
    }
    //else if(currentLagrLB >= getIncumbentVal()){
    else if(false){
	modelStatus_[Z_STATUS]=Z_BOUNDED;
if(mpiRank==0){cout << "Terminating due to exceeding cutoff..." << endl;}
      	//printStatus();
    }
    else{
      //objVal=findPrimalFeasSolnWith(z_current);
      //currentLagrLB=-ALPS_DBL_MAX;
      //zeroOmega();
      int noTimesOmegaUpdated=0;
      int noConseqTimesOmegaNotUpdated=0;
      bool shouldContinue=true;
      bool postprocessing=false;
      bool reachedMaxNoIts = false;
      //for(int ii=0; ii<nIters || !exceedingReferenceBd || !omegaUpdatedAtLeastOnce || (discrepNorm >= 1e-20 && integralityDisc < 1e-10) ; ii++){
      for(int ii=0; shouldContinue; ii++){
	currentIter_=ii;
if(mpiRank==0) cerr << "Regular iteration " << ii << endl;
	//if(postprocessing || ii > 100){
#if 0
	if(postprocessing){
    	     //if(mpiRank==0)cout << "computeBound(): Integrality satisfied: Commence with postprocessing..." << endl;
    	     computeScalingPenaltyUpdate(1.05);	
    	     regularIteration(false,false);
	}
#endif
	{regularIteration(true,true);}
	//{regularIteration(false,adjustPenalty);}
	//if(ii<20) {regularIteration(false,adjustPenalty);}
	//else{regularIteration(true,false);}

	//regularIteration(true,false);
	//if(omegaUpdated_) omegaUpdatedAtLeastOnce=true;
	if(omegaUpdated_){ 
            noConseqTimesOmegaNotUpdated=0;
	    if(exceedingReferenceBd){
    	        //objVal=findPrimalFeasSolnWith(z_current);
		noTimesOmegaUpdated++;
	    }
	}
	else{noConseqTimesOmegaNotUpdated++;}
	//omegaUpdatedAtLeastOnce=true;
      	//printStatusAsErr();
	printStatus();
#if 0
        if(currentLagrLB >= getIncumbentVal()){
      	    modelStatus_[Z_STATUS]=Z_BOUNDED;
if(mpiRank==0){cout << "computeBound(): Terminating due to exceeding cutoff..." << endl;}
      	    //printStatus();
	    break;
	}
#endif
        exceedingReferenceBd = currentLagrLB + SSC_DEN_TOL >= referenceLagrLB;
	integralityDisc=evaluateIntegralityDiscrepancies();
	//postprocessing = ((integralityDisc < 1e-10) && (ii>20)) || postprocessing;
	//postprocessing = ((integralityDisc < 1e-10) && (ii>20)) && (noConseqTimesOmegaNotUpdated < 2*maxNoConseqNullSteps);
        //exceedingReferenceBd = true;
#if 1
	//if(shouldTerminate()){
        reachedMaxNoIts = (ii>=maxNoIts);
	if((shouldTerminate) || reachedMaxNoIts){// || noTimesOmegaUpdated > 50){
if(mpiRank==0 && reachedMaxNoIts){cout << "computeBound(): Terminating due to reaching tolerance criterion at iteration " << ii << endl;}
else{
   if(mpiRank==0){cout << "computeBound(): Terminating due to reaching maximum number of iterations: " << ii << endl;}
}
      	    //updateModelStatusZ();
      	    //printStatus();
	    //return currentLagrLB;
	    break;
	}
#endif
	//shouldContinue = (ii<nIters) || !exceedingReferenceBd || (integralityDisc < 1e-10 && discrepNorm >=1e-20); 
	//shouldContinue = (noTimesOmegaUpdated < 8) || !exceedingReferenceBd || (integralityDisc < 1e-10 && discrepNorm >=1e-20); 
	//shouldContinue = (noConseqTimesOmegaNotUpdated < maxNoConseqNullSteps || noTimesOmegaUpdated < 10) || !exceedingReferenceBd || (integralityDisc < 1e-10 && discrepNorm >=1e-20); 
	shouldContinue = (noConseqTimesOmegaNotUpdated < maxNoConseqNullSteps || noTimesOmegaUpdated < 10);// || !exceedingReferenceBd || postprocessing; 
if(mpiRank==0) cerr << "Omega updated " << noTimesOmegaUpdated << " times, max no. consec null steps: " << maxNoConseqNullSteps << " integrality disc: " << integralityDisc << endl;

    //if(currentLagrLB < referenceLagrLB && mpiRank==0){
//if(mpiRank==0 && ii>=(nIters-1) && exceedingReferenceBd && omegaUpdatedAtLeastOnce){cout << "computeBound(): Terminating due to reaching the maximum number of iterations." << endl;}
    }
    solveContinuousMPs();
    printStatus();
    //integralityDisc=evaluateIntegralityDiscrepancies();
    objVal=findPrimalFeasSolnWith(z_current);
#if 1
    if(objVal - currentLagrLB >= 1e-6 && integralityDisc < 1e-10 && modelStatus_[Z_STATUS]!=Z_BOUNDED && (shouldTerminate) ){
    //if(false){
    //if(true){
	//setPenalty(1.0);
	int remoteNNodeSPs;
	double *zBuffer = new double[n1];
	for(int proc=0; proc<mpiSize; proc++){
	     if(proc==mpiRank){ remoteNNodeSPs = nNodeSPs; }
             MPI_Bcast(&remoteNNodeSPs, 1, MPI_INT, proc, MPI_COMM_WORLD);
             for(int tS=0; tS<remoteNNodeSPs; tS++){
		if(proc==mpiRank){memcpy(zBuffer,subproblemSolvers[tS]->getX(),n1*sizeof(double));}
                MPI_Bcast(zBuffer, n1, MPI_DOUBLE, proc, MPI_COMM_WORLD);
		objVal=findPrimalFeasSolnWith(zBuffer);
	     }
	//cerr << "My rank " << mpiRank << " remoteNNodeSPs " << remoteNNodeSPs << endl;     
	}
	delete [] zBuffer;
    }
#endif
    //printStatus();
    }//else
    //for (int tS = 0; tS < nNodeSPs; tS++){subproblemSolvers[tS]->computeWeightsForCurrentSoln();}
    //loadSavedZ();
    //printOmegaProperties();
    //averageOfOptVertices();
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
        if (useVertexHistory && subproblemSolvers[tS]->getNVertices()>2) { //Compute Next X, Y (With History)
	     //subproblemSolvers[tS]->updatePrimalVariablesHistory_OneScenario(omega_current[tS],z_current,scaling_matrix[tS],true);
	     subproblemSolvers[tS]->computeWeightsForCurrentSoln(NULL);
	}
	//if(mpiRank==0){subproblemSolvers[tS]->printWeights();}
   }

}
void solveContinuousMPs(){
	bool isLastGSIt;
	double zDiff;
        for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_tilde[tS],omega_current[tS],n1*sizeof(double));
	//double integrDiscr;
	//for(int itGS=0; itGS < fixInnerStep || (updateDisp && zDiff > 1e-6); itGS++) { //The inner loop has a fixed number of occurences
	assert(currentIter_ >=0 && currentIter_ < 1000);
	//int maxNoGSIts = 20+currentIter_;
	int maxNoGSIts = 20;
	//averageOfVertices(z_vertex_average);
	for(int itGS=0; itGS < maxNoGSIts; itGS++) { //The inner loop has a fixed number of occurences
	    memcpy(z_old,z_current, n1*sizeof(double));
	    zDiff=0.0;
    	    for (int tS = 0; tS < nNodeSPs; tS++) {
		//*************************** Quadratic subproblem ***********************************
			    
		    if (useVertexHistory && subproblemSolvers[tS]->getNVertices()>2) { //Compute Next X, Y (With History)

					
			//Solve the quadratic master problem for x and y
			subproblemSolvers[tS]->updatePrimalVariablesHistory_OneScenario(omega_tilde[tS],z_current,currentVarLB_, currentVarUB_, scaling_matrix[tS],false);
			
		    }
		    else if(!useVertexHistory || subproblemSolvers[tS]->getNVertices()==2){ //might not work correctly, untested! Compute Next X, Y (Without History)
			//subproblemSolvers[tS]->updatePrimalVariables_OneScenario(omega_tilde[tS],z_current,scaling_matrix[tS],z_vertex_average);
			subproblemSolvers[tS]->updatePrimalVariables_OneScenario(omega_tilde[tS],z_current,scaling_matrix[tS],NULL);
		    }
	    }
	    					
	    // Update z_previous.
	    updateZ();
#if 0
    	    for (int tS = 0; tS < nNodeSPs; tS++) {
		for (int i = 0; i < n1; i++) {
		    omega_tilde[tS][i] += scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
		}
	    }
#endif
	    totalNoGSSteps++;
	    for(int ii=0; ii<n1; ii++){zDiff=max(zDiff,fabs(z_current[ii]-z_old[ii]));}
	    //if(zDiff < SSC_DEN_TOL){ break; }
	}
#if 1
    	    for (int tS = 0; tS < nNodeSPs; tS++) {
		for (int i = 0; i < n1; i++) {
		    omega_tilde[tS][i] += scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
		}
	    }
#endif

if(mpiRank==0){cout << "solveContinuousMPs(): max zDiff is: " << zDiff << endl;}
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
		z_local[i] += x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
		penSumLocal[i]+= scaling_matrix[tS][i]*pr[tS];
	    }
	}
	#ifdef USING_MPI
	    MPI_Allreduce(z_local, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    MPI_Allreduce(penSumLocal, penSum, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    for(int i=0; i<n1; i++) z[i] /= penSum[i];
	#else
	    for (int i=0; i<n1; i++) z[i] = z_local[i]/penSumLocal[i]; //Only one node, trivially set z_local = z_current
	#endif

	if(dispZ!=NULL){
    	    double *dispZLocal = new double[n1];
    	    for(int ii=0; ii<n1; ii++){
		dispZLocal[ii]=0.0;
		dispZ[ii]=0.0;
        	for(int tS=0; tS<nNodeSPs; tS++){
	    	     dispZLocal[ii] += scaling_matrix[tS][ii]*pr[tS]*fabs(z[ii] - x_current[tS][ii]);
		}
    	    }
	#ifdef USING_MPI
    	    if(mpiSize>1){
		MPI_Allreduce(dispZLocal, dispZ, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    	    }
	#endif
    	    if(mpiSize==1){
		memcpy(dispZ,dispZLocal,n1*sizeof(double));
    	    }
    	    for(int ii=0; ii<n1; ii++){dispZ[ii]/=penSum[ii];}
	
	    delete [] dispZLocal;
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
void averageOfXUniformWeights(double* z=NULL, double *dispZ=NULL){
	//if(z==NULL) z=z_average;
	double *penSumLocal = new double[n1];
	double *penSum = new double[n1];
	double integrDiscr;
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    penSumLocal[i]=0.0;
            penSum[i]=0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += subproblemSolvers[tS]->getX()[i];
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
	if(dispZ!=NULL){
    	    double *dispZLocal = new double[n1];
    	    for(int ii=0; ii<n1; ii++){
		dispZLocal[ii]=0.0;
		dispZ[ii]=0.0;
        	for(int tS=0; tS<nNodeSPs; tS++){
	    	     dispZLocal[ii] += fabs(z[ii] - x_current[tS][ii]);
		}
    	    }
	#ifdef USING_MPI
    	    if(mpiSize>1){
		MPI_Allreduce(dispZLocal, dispZ, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    	    }
	#endif
    	    if(mpiSize==1){
		memcpy(dispZ,dispZLocal,n1*sizeof(double));
    	    }
    	    for(int ii=0; ii<n1; ii++){dispZ[ii]/=penSum[ii];}
	
	    delete [] dispZLocal;
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
	delete [] penSumLocal;
	delete [] penSum;
}
#endif
void averageOfXweightedByIntDisp(double* z=NULL,bool roundZ=false){
	//if(z==NULL) z=z_intdisp;
	double *penSumLocal = new double[n1];
	double *penSum = new double[n1];
	double integrDiscr;
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    penSumLocal[i]=0.0;
            penSum[i]=0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		integrDiscr = subproblemSolvers[tS]->computeIntegralityDiscr();  
		z_local[i] += integrDiscr*subproblemSolvers[tS]->getX()[i];
		penSumLocal[i]+= integrDiscr;
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
void averageOfVertices(double* z=NULL, double *zDisp=NULL){
	//if(z==NULL) z=z_average;
#if 0
	for(int tS=0; tS<nNodeSPs; tS++){
	  if(!(subproblemSolvers[tS]->getXVertex()[ii] >= currentVarLB_[ii] && subproblemSolvers[tS]->getXVertex()[ii] <= currentVarUB_[ii])){
cout << "Scenario " << tS << " index " << ii << " with value " << subproblemSolvers[tS]->getXVertex()[ii] << endl;
	    assert(subproblemSolvers[tS]->getXVertex()[ii] >= currentVarLB_[ii] && subproblemSolvers[tS]->getXVertex()[ii] <= currentVarUB_[ii]);
	  }
	}
#endif
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
	delete [] penSumLocal;
	delete [] penSum;
}
void averageOfVerticesWeightedByRho(double* z=NULL, double *zDisp=NULL){
	//if(z==NULL) z=z_average;
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
	delete [] penSumLocal;
	delete [] penSum;
}
#if 0
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
#endif
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
void averageOfVerticesZ(double *z=NULL, double *zDisp=NULL){
	//if(z==NULL) z=z_intdisp;
        solveForWeights();
	int numVertices;
	//double localTotalWeight=0.0;
	//double totalWeight=0.0;
	double *penSumLocal = new double[n1];
	double *penSum = new double[n1];
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
	delete [] penSumLocal;
	delete [] penSum;

}

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

void computeOmegaDisp(double *dispOmega){
    double *dispOmegaLocal = new double[n1];
    for(int ii=0; ii<n1; ii++){
	dispOmegaLocal[ii]=0.0;
	dispOmega[ii]=0.0;
        for(int tS=0; tS<nNodeSPs; tS++){
	    dispOmegaLocal[ii] += fabs(pr[tS]*omega_current[tS][ii]);

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

#if 0
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
    if(ALVal - currentLagrLB < -SSC_DEN_TOL){
        if(mpiRank==0){cout << " (ALVal,currentLagrLB) (" << setprecision(10) << ALVal << "," << setprecision(10) << currentLagrLB << ")" << endl;}
	//shouldTerminate = true;
	shouldTerminate = true;
	return 0.0;
    }
    if(ALVal + 0.5*discrepNorm  - LagrLB < -SSC_DEN_TOL){
        if(mpiRank==0){cout << mpiRank << " (ALVal,LagrLB) (" << setprecision(10) << ALVal << "," << setprecision(10) << LagrLB << ")" << endl;	
	cout << "regularIteration(): Something probably went wrong with the last computation of LagrLB, returning..." << endl;}
	//shouldTerminate = true;
	shouldTerminate = true;
	return 0.0;
    }
    shouldTerminate = (ALVal + 0.5*discrepNorm  - currentLagrLB < SSC_DEN_TOL) && (discrepNorm < 1e-20);
    if(shouldTerminate){return 0.0;}
    else{
	return (LagrLB-currentLagrLB)/(ALVal + 0.5*discrepNorm - currentLagrLB);
    }
}
		
bool updateOmega(bool useSSC){
	if(SSCVal >= SSC_PARAM || !useSSC) {
	    for (int tS = 0; tS < nNodeSPs; tS++) {
		memcpy(omega_current[tS],omega_tilde[tS],n1*sizeof(double));
#if 0
		for (int i = 0; i < n1; i++) {
	    	    omega_current[tS][i] = omega_tilde[tS][i];
		}
#endif
	        recordKeeping[tS][1]=recordKeeping[tS][0];
	        //subproblemSolvers[tS]->updateOptSoln();
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

void repairOmega(){
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
	memcpy(omegaSum,omegaLocal,n1*sizeof(double));
    }
#if 0
cout << "Printing omega weighted sum..." << endl;
    for(int ii=0; ii<n1; ii++){
	cout << " " << omegaSum[ii];	
    }
cout << endl;
#endif
if(mpiRank==0){cout << "Repairing dual feasibility of omega..." << endl;}
    for (int tS = 0; tS < nNodeSPs; tS++) {
      for (int ii = 0; ii < n1; ii++) {
        omega_current[tS][ii] -= (1.0/((double)nS*pr[tS]))*omegaSum[ii];
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
double getPenalty(){return penC;}
double getBaselinePenalty(){return baselinePenalty_;}
void setPenalty(double p){
    penC=max(p,MIN_PEN);
    for (int tS = 0; tS < nNodeSPs; tS++) {
	setPenalty(tS,penC);
    }
//if(mpiRank==0) cout << "Penalty is now: " << penC << endl;
}
void setPenalty(int tS, double p){
    double penalty = max(p,MIN_PEN);
    subproblemSolvers[tS]->setQuadraticTerm(penalty);
    for (int i = 0; i < n1; i++) {
        scaling_matrix[tS][i] = penalty;
    }
//if(mpiRank==0) cout << "Penalty is now: " << penC << endl;
}
//This should normally be turned off, this is a rule for updating the 
//penalty from Kiwiel 2006 and Lubin et al.
void computeKiwielPenaltyUpdate(){	
	//penC = 1.0/min( max( (2.0/penC)*(1.0-SSCVal),  max(1.0/(10.0*penC),1e-4)    ), 10.0/penC);
	//setPenalty(penC);
    	//if( SSCVal >= 0.5 ){ computeScalingPenaltyUpdate(min( 0.5/(1-SSCVal),2.0));}
    	computeScalingPenaltyUpdate( max( min( 0.5/(1.0-SSCVal),10.0), 0.1));
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
	scaling_matrix[tS][i] = max( scaling_matrix[tS][i], MIN_PEN);
    }
    subproblemSolvers[tS]->setQuadraticTerm(scaling_matrix[tS]);
}
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
	printf("LagrLB: %0.14g, ALVal: %0.14g, sqrDiscNorm: %0.14g\n", currentLagrLB, ALVal, discrepNorm);
//	printf("Best node: %0.14g, Incumbent value: %0.14g\n", getBestNodeQuality(), getIncumbentVal());
cout << "________________________________________________________________________" << endl;
	//printf("Aug. Lagrangian value: %0.9g\n", ALVal);
	//printf("Norm of primal discrepancy: %0.6g\n", discrepNorm);
	//printf("Current penalty: %0.2g\n",penC);
	//std::cout << "Number of integrality infeasibilities of z: " << nIntInfeas_ << std::endl;
  }
}
#if 0
void printBestBounds(){
  if(mpiRank==0){
cout << "________________________________________________________________________" << endl;
	printf("Best node: %0.14g, Incumbent value: %0.14g\n", getBestNodeQuality(), getIncumbentVal());
//	printf("Best node: %0.9g, Incumbent value: %0.9g\n", getBestNodeQuality(), getIncumbentVal());
cout << "________________________________________________________________________" << endl;
  }
}
void printStatusAsErr(){
  if(mpiRank==0){
	fprintf(stderr, "LagrLB: %0.14g, ALVal: %0.14g, sqrDiscNorm: %0.14g\n", currentLagrLB, ALVal, discrepNorm);
	double bestNodeQuality = -1e20;
	if(getKnowledgeBroker()->getBestNode()){
	  bestNodeQuality = getKnowledgeBroker()->getBestNode()->getQuality();
	}
	fprintf(stderr,"Best node: %0.14g, Incumbent value: %0.14g, Candidate value: %0.14g \n", bestNodeQuality, getIncumbentVal(), objVal);
	//printf("Aug. Lagrangian value: %0.9g\n", ALVal);
	//printf("Norm of primal discrepancy: %0.6g\n", discrepNorm);
	//printf("Current penalty: %0.2g\n",penC);
	//std::cout << "Number of integrality infeasibilities of z: " << nIntInfeas_ << std::endl;
  }
}
#endif
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
#if 0
void printIncumbentValue(){
  if(mpiRank==0){
	printf("\nIncumbent value: %0.14g\n", getIncumbentVal());
  }
}
#endif
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

#if 0
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
#endif


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
/*Header for the main procedure ParallelSCG. */


