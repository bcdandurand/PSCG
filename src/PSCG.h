
/*Header for the main procedure ParallelSCG. */

#ifndef PSCG_H
#define PSCG_H

#include <stdio.h>
//#include <math.h>
#include <cmath>
#include <string>
//#include <memory>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
//#include <list>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <sys/time.h>
#include "PSCGParams.h"
#include <utility>
#include "DecTssModel.h"
#include "PSCGScen.h"

using namespace std;

#define SSC_DEN_TOL 1e-10
#define DEFAULT_NO_THREADS 1
#define KIWIEL_PENALTY 1 //set 1 to use Kiwiel (2006) penalty update rule
#define MIN_PEN 0.0 
#define MAX_NO_INNERSTEPS 1000

#ifdef USING_MPI
   #include <mpi.h>
#endif

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#define MAX_INNER_LOOP 200
#define DEFAULT_PENALTY 100
//#define KEEP_LOG

// Parameters for newConvergenceCriterion

#define BG_BETA 0.05

enum COL_TYPE{
CONTINUOUS=0,
BINARY,
INTEGER
};

typedef struct{
    //double a,b,c,d;
    double rank,scen,index,disp,brVal,brLBUp,brUBUp,brLBDn,brUBDn;
} BranchingVarInfo;

typedef struct{
    vector<int> ranks;
    vector<int> sps;
    vector<int> inds;
    vector<double> lbs;
    vector<double> ubs;
} BranchingBDs;
#ifdef USING_MPI
#if 0
static void compareBranchInfoOld(BranchingVarInfo *in, BranchingVarInfo *inout, int *len, MPI_Datatype *dptr){
  for(int ii=0; ii< *len; ii++){
    if(in[ii].intDiscr > inout[ii].intDiscr){
	inout[ii].rank = in[ii].rank;
	inout[ii].scen = in[ii].scen;
	inout[ii].index = in[ii].index;
	inout[ii].disp = in[ii].disp;
	inout[ii].brVal = in[ii].brVal;
	inout[ii].brLBUp = in[ii].brLBUp;
	inout[ii].brUBUp = in[ii].brUBUp;
	inout[ii].brLBDn = in[ii].brLBDn;
	inout[ii].brUBDn = in[ii].brUBDn;
    }
  }
}
#endif
static void compareBranchInfo(BranchingVarInfo *in, BranchingVarInfo *inout, int *len, MPI_Datatype *dptr){
  for(int ii=0; ii< *len; ii++){
    if(in[ii].disp > inout[ii].disp){
	inout[ii].rank = in[ii].rank;
	inout[ii].scen = in[ii].scen;
	inout[ii].index = in[ii].index;
	inout[ii].disp = in[ii].disp;
	inout[ii].brVal = in[ii].brVal;
	inout[ii].brLBUp= in[ii].brLBUp;
	inout[ii].brUBUp = in[ii].brUBUp;
	inout[ii].brLBDn = in[ii].brLBDn;
	inout[ii].brUBDn = in[ii].brUBDn;
    }
  }
}
#endif


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
//PSCGParams *par;
DecTssModel &smpsModel;
IloEnv env;
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
int noInnerSolves;
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


double innerSSCVal;
double SSCVal;
double SSCParam;
double innerSSCParam;
double tCritVal;
double tCritParam;
//double *scaleVec_;

vector<double*> x_current;
vector<double*> y_current;

double* z_current;// = new double[n1];
double* z_old;// = new double[n1];
double* z_local;// = new double[n1];
double* z_incumbent_; //this should be the last feasible z with best obj
double *z_rounded;
double *z_saved;
double *z_average;
double *z_dispersions;
vector<double> constrVec_; //This is filled with the value of Ax
double *penSumLocal;// = new double[n1];
double *penSum;// = new double[n1];
double *weights_;// = new double[nS];
double *integrDiscr_;

BranchingBDs branchingBDs_;
BranchingVarInfo brVarInfo;

vector<int> spSolverStatuses_;

int modelStatus_[2];

char filepath[256];
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

double objVal;
double incumbentVal;

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
int maxNoSteps;
int maxNoInnerSteps;
int maxNoGSSteps;
int maxNoConseqNullSteps;
int noGSIts;
int maxSeconds;
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
#ifdef USING_MPI
MPI_Comm comm_;
#endif

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
//char zOptFile[128];

ProblemDataBodur pdBodur;

PSCG(PSCGParams *p);
PSCG(DecTssModel &model);
#ifdef USING_MPI
PSCG(DecTssModel &model, MPI_Comm comm);
#endif

~PSCG();

void setMPIParams(int rank, int size){
	mpiRank=rank;
	mpiSize=size;
}
int getMPIRank(){return mpiRank;}
int getN1(){return n1;}
int getN2(){return n2;}
int getNS(){return nS;}

void assignSubproblems();
IloEnv& getEnv(){return env;}

void initialiseParameters();


void initialiseFileName(){
  char filepathcpy[256];
  strcpy(filepathcpy,filepath);
  char * pch, *prb;
  //printf ("Splitting string \"%s\" into tokens:\n",str);
  pch = strtok (filepathcpy,"/");
  while (pch != NULL)
  {
    strcpy(probname,pch);
    //printf ("%s\n",pch);
    pch = strtok (NULL, "/");
  }
}

void initialiseModel();


void setupSolvers();

#if 0
void fixVarAllSPsAt(int index, double fixVal){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->fixVarAt(index, fixVal);
    }
    currentVarLB_[index]=fixVal;		
    currentVarUB_[index]=fixVal;		
}
#endif
#if 0
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
#endif
void setBdForSP(const int sp, const int ind, const double lb, const double ub){
    subproblemSolvers[sp]->setBound(ind, lb, ub);
    currentVarLB_[ind]=lb;
    currentVarUB_[ind]=ub;
}
void setBdAllSPs(const int ind, const double lb, const double ub){
    for(int tS=0; tS<nNodeSPs; tS++){
	setBdForSP(tS,ind,lb,ub);	
    }
}

void setCutoffLagrBD(double val){cutoffLagrLB = val;}
#if 0
void setBdsAllSPs(const vector<int> &inds, const vector<double> &lbs, const vector<double> &ubs){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->setBounds(inds, lbs, ubs);
    }
    for(int ii=0; ii<inds.size(); ii++){
	currentVarLB_[inds[ii]]=lbs[ii];
	currentVarUB_[inds[ii]]=ubs[ii];
    }
}
#endif

void clearBranchingBounds(){
    branchingBDs_.ranks.clear();
    branchingBDs_.sps.clear();
    branchingBDs_.inds.clear();
    branchingBDs_.lbs.clear();
    branchingBDs_.ubs.clear();
}
void restoreOriginalVarBounds(){
    clearBranchingBounds();
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->restoreBounds();//unfixX(origVarLB_,origVarUB_);
    }
    memcpy(currentVarLB_,origVarLB_,n1*sizeof(double));
    memcpy(currentVarUB_,origVarUB_,n1*sizeof(double));
}
void clearSPVertexHistory(){
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->clearVertexHistory();
    }
}

void addBranchVarBd(int br_rank, int br_SP, int br_index, double br_lb, double br_ub){
    branchingBDs_.ranks.push_back(br_rank);
    branchingBDs_.sps.push_back(br_SP);
    branchingBDs_.inds.push_back(br_index);
    branchingBDs_.lbs.push_back(br_lb);
    branchingBDs_.ubs.push_back(br_ub);
    if(br_index >=0){
        if(br_rank < 0 && br_SP < 0){
	  setBdAllSPs(br_index, br_lb, br_ub);
        }
        else if(br_rank==mpiRank && br_SP>=0 && br_SP < nNodeSPs){
	    setBdForSP(br_SP, br_index, br_lb, br_ub);
        }
	//else this processor and/or subproblem does nothing
    }
    else{
	cout << "No branch found, this node should be solved to optimality." << endl;
    }

}

void installSubproblem(double lb, vector<double*> &omega, const double *zLBs, const double *zUBs);
void installSubproblem(double lb, vector<double*> &omega, const vector<int> &indices, const vector<double> &zLBs, const vector<double> &zUBs);


void installSubproblem(double lb, const double *zLBs, const double *zUBs);
void installSubproblem(double lb, const vector<int> &indices, const vector<double> &zLBs, const vector<double> &zUBs);


double *getZ(){return z_current;}
double *getZRounded(){return z_rounded;}
double *getZIncumbent(){return z_incumbent_;}
double *getZAverage(){return z_average;}
double getSqrDiscrNorm(){return discrepNorm;}
double getTCritVal(){return tCritVal;}
void printTCritVal(){if(mpiRank==0) cout << "tCritVal: " << tCritVal << endl;}

double* getOrigVarLbds(){return origVarLB_;}
double* getOrigVarUbds(){return origVarUB_;}
double* getCurrentVarLbds(){return currentVarLB_;}
double* getCurrentVarUbds(){return currentVarUB_;}

void cloneCurrentVarBds(double*& lbs, double*& ubs){
printCurrentVarBds();
    lbs = new double[n1];
    memcpy(lbs,currentVarLB_,n1*sizeof(double));
    ubs = new double[n1];
    memcpy(ubs,currentVarUB_,n1*sizeof(double));
}
void readCurrentVarBds(const double *lbs, const double *ubs){
    memcpy(currentVarLB_,lbs,n1*sizeof(double));
    memcpy(currentVarUB_,ubs,n1*sizeof(double));
}
void writeCurrentVarBds(double *lbs, double *ubs){
    memcpy(lbs,currentVarLB_,n1*sizeof(double));
    memcpy(ubs,currentVarUB_,n1*sizeof(double));
}

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
double getPrimalObjVal(){return objVal;}
double getIncumbentVal(){return incumbentVal;}

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
void loadOmegaCurrent(){
    omegaIsZero_=false;
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_current[tS],omega_saved[tS],n1*sizeof(double));
}
void loadOmegaCentre(){
    omegaIsZero_=false;
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
void writeOmegaCurrent(vector<double*> &writeOmega){
    //assert(writeOmega.size() >= nNodeSPs); 
    if(writeOmega.size()==0){
        for(int tS=0; tS<nNodeSPs; tS++){
	    writeOmega.push_back(new double[n1]);
	}
    }
    for(int tS=0; tS<nNodeSPs; tS++){
	assert(writeOmega[tS]!=NULL);
        memcpy(writeOmega[tS],omega_current[tS],n1*sizeof(double));
    }
}

int initialIteration();

#if 0
void updateSSC(int k){
    SSCParam = 0.5;
    //SSCParam = k/(100.0+k);
    //if(SSCParam > 0.5) SSCParam=0.5;
}
#endif

void updatePenalty(){
    //updateSSC(k);
    //if(k % 50 == 0 && omegaUpdated_) computeScalingPenaltyUpdate(2.0);
    //if(currentIter_<=10){
    //computeKiwielPenaltyUpdate();
    //double penaltyScale = 1.05;
    //computeScalingPenaltyUpdate(penaltyScale);
    setPenalty(rho+baselineRho);
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
}

//Update penalty and SSC parameter as motivated by convergence analysis
void updateVertexHistory(int tS){
	if(subproblemSolvers[tS]->getNumVertices() < nVerticesUsed || nVerticesUsed==0){subproblemSolvers[tS]->addVertex();}
	else{// if(nVerticesUsed > 0){// && subproblemSolvers[tS]->getNumVertices() >= nVerticesUsed) {
		subproblemSolvers[tS]->replaceOldestVertex();
	}
}



void setMaxNoSteps(int noSteps){maxNoSteps=noSteps;}
void setMaxNoGSSteps(int noSteps){maxNoGSSteps=noSteps;}
void setMaxNoConseqNullSteps(int noNullSteps){maxNoConseqNullSteps=noNullSteps;}
void setSSCParam(double ssc){SSCParam=ssc;}
void setTCritParam(double tcrit){tCritParam=tcrit;}
void setPhase(int ph){phase=ph;}

void preSolveMP(){
    for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_tilde[tS],omega_centre[tS],n1*sizeof(double));
    //int maxNoGSIts = 100+currentIter_;//max(20,currentIter_/2);
    //int maxNoGSIts = 1;//max(20,currentIter_/2);
    //int maxNoGSIts = min(100*currentIter_,10000000);//max(20,currentIter_/2);
    int maxNoGSIts = 100;//max(20,currentIter_/2);
    for(int itGS=0; itGS < maxNoGSIts; itGS++) { //The inner loop has a fixed number of occurences
	for (int tS = 0; tS < nNodeSPs; tS++) {
		//if(subproblemSolvers[tS]->getNVertices()>0){subproblemSolvers[tS]->solveMPVertices(omega_tilde[tS],z_current,rho,scaling_matrix[tS]);}
		if(subproblemSolvers[tS]->getNVertices()>0){subproblemSolvers[tS]->solveMPHistory(omega_tilde[tS],z_current,NULL,NULL,rho,scaling_matrix[tS],false);}
	}
	updateZ();
    	for (int tS = 0; tS < nNodeSPs; tS++) {
	    for (int i = 0; i < n1; i++) {
		omega_tilde[tS][i] += rho*scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
      	    subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]); //prepares next call of solveMPLineSearch(omega,z,scaling_vector)
	}
    }

}

bool solveContinuousMPs(bool adjustPenalty){
	bool isLastGSIt;
	bool needToContinue = false;
	double zDiff;
        //for(int tS=0; tS<nNodeSPs; tS++) memcpy(omega_tilde[tS],omega_centre[tS],n1*sizeof(double));
	//double integrDiscr;
	assert(currentIter_ >= 0);
	//int maxNoGSIts = 20+currentIter_;
	//int maxNoGSIts = 100 + 10*currentIter_;//max(20,currentIter_/2);
	//int maxNoGSIts = min(100*currentIter_,10000000);//max(20,currentIter_/2);
#if 0
    	for (int tS = 0; tS < nNodeSPs; tS++) {
	    for (int i = 0; i < n1; i++) {
	        omega_tilde[tS][i] = omega_centre[tS][i] + rho*scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
	    subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]);
	}
#endif
	for(int itGS=0; itGS < noGSIts; itGS++) { //The inner loop has a fixed number of occurences
	    memcpy(z_old,z_current, n1*sizeof(double));
	    zDiff=0.0;
    	    for (int tS = 0; tS < nNodeSPs; tS++) {
		//*************************** Quadratic subproblem ***********************************
			    
		if(subproblemSolvers[tS]->getNVertices()>0){
		    for(int iii=0; iii<10; iii++){
	    		for (int i = 0; i < n1; i++) {
	        	    omega_tilde[tS][i] = omega_centre[tS][i] + rho*scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    		}
	    	    	subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]);
		    	subproblemSolvers[tS]->solveMPVertices(omega_centre[tS],z_current,rho,scaling_matrix[tS]);
		    }
		}
		//if(subproblemSolvers[tS]->getNVertices()>0){subproblemSolvers[tS]->solveMPHistory(omega_centre[tS],z_current,NULL,NULL,rho,scaling_matrix[tS],false);}
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
		    omega_tilde[tS][i] = omega_centre[tS][i] + rho*scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
		}
		refVal = subproblemSolvers[tS]->evaluateSolution(omega_tilde[tS]);
		LagrLB_tS = subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]);
		innerLagrLB_Local += pr[tS]*LagrLB_tS;
	        ALVal_tS = subproblemSolvers[tS]->updateALValues(omega_centre[tS],z_current,rho,scaling_matrix[tS]);
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
		MPI_Allreduce(localReduceBuffer, reduceBuffer, 3, MPI_DOUBLE, MPI_SUM, comm_);
		innerLagrLB = reduceBuffer[0];
		ALVal = reduceBuffer[1];
		discrepNorm = reduceBuffer[2];
	}
	#endif
	if (mpiSize == 1) {
		innerLagrLB = innerLagrLB_Local;
		ALVal = ALVal_Local;
		discrepNorm = localDiscrepNorm;
	}
    	if( (ALVal + 0.5*discrepNorm  - centreLagrLB > SSC_DEN_TOL) ){
	    innerSSCVal =    (innerLagrLB-centreLagrLB)/(ALVal + 0.5*discrepNorm - centreLagrLB);
	    needToContinue = innerSSCVal < innerSSCParam;
	}
#endif
//if(mpiRank==0){cout << "solveContinuousMPs(): max zDiff is: " << zDiff << endl;}

    return needToContinue;
}

int regularIteration(bool adjustPenalty=false, bool SSC=true){
//if(mpiRank==0) cout << "Begin regularIteration()" << endl;

	int SPStatus=performColGenStep();
	assert(SPStatus!=SP_INFEAS); //subproblem infeasibility should be caught in initialIteration()
	if(currentIter_>0)SSCVal = computeSSCVal(); //shouldTerminate is also updated here.
	else{
	     SSCVal=1.0;
	     shouldTerminate=false;
	}
    //verifyOmegaDualFeas();
    //printStatus();
	//cout << "ALVal: " << ALVal << " discrepNorm: " << discrepNorm << " currentLagrLB " << currentLagrLB << endl;
	updateOmega(SSC);
	if(adjustPenalty) updatePenalty();

	if(phase==2) preSolveMP();
	else{
          noInnerSolves=0;	
	  do{
	    solveContinuousMPs(adjustPenalty);
	    //if(mpiRank==0) cout << "***************************************Need to continue with solveMP..." << endl;
	    noInnerSolves++;
	  }while(innerSSCVal < innerSSCParam && noInnerSolves<maxNoInnerSteps);
	}
//if(mpiRank==0) cout << "End regularIteration()" << endl;
	if(currentIter_%10==0) printStatus();


        if(currentLagrLB >= cutoffLagrLB){
      	    modelStatus_[Z_STATUS]=Z_BOUNDED;
	    if(mpiRank==0){cout << "computeBound(): Terminating due to exceeding cutoff..." << endl;}
            shouldContinue=false;
	}

	double integralityDisc=evaluateIntegralityDiscrepancies();
	if((shouldTerminate) || (currentIter_ >= maxNoSteps) ){
	    if(mpiRank==0 && (currentIter_ >= maxNoSteps) ){
		cout << "computeBound(): Terminating due to reaching tolerance criterion at iteration " << currentIter_ << endl;
	    }
	    else{
   		if(mpiRank==0){cout << "computeBound(): Terminating due to reaching maximum number of iterations: " << currentIter_ << endl;}
	    }
            shouldContinue=false;
	}
#if 0
	if(mpiRank==0) cerr << "Omega updated " << noSeriousSteps << " times, " 
		<< " integrality disc: " << integralityDisc << endl;
if(mpiRank==0){cout << "Rho is: " << rho << endl;}
#endif
	return SPStatus;
}

double computeBound(){
//if(mpiRank==0) cout << "Begin computeBound()" << endl;

    //SSCParam = 0.10;
    if(mpiRank==0) printCurrentVarBds();
#if 1
    for (int tS = 0; tS < nNodeSPs; tS++) {
	if(subproblemSolvers[tS]->printModifiedYBounds()){
	    cout << "[" << mpiRank << "," << tS << "]" << endl;
	}
    }
#endif
    
    int SPStatus=initialIteration();
    //if(mpiRank==0) cout << "After initial iteration: currentLB: " << trialLagrLB << " versus refLB: " << referenceLagrLB << endl;;
    if(modelStatus_[SP_STATUS]==SP_INFEAS){
	currentLagrLB=COIN_DBL_MAX;
        brVarInfo={-1.0,-1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0};
	return currentLagrLB;
    }
    else{return processBound();}
}

double processBound(){
    bool updatePenalty=true;
    for(currentIter_=0; shouldContinue; currentIter_++){
	#ifdef KEEP_LOG
	    for(int tS=0; tS<nNodeSPs; tS++){*(logFiles[tS]) << endl << "Regular iteration: " << currentIter_ << endl;}
	#endif
        //if(currentIter_ < 1) phase=0;
        //else if(currentIter_ < 20) phase=1;
        //else phase=1;

	phase=0;

	if(currentIter_>200) updatePenalty=false;
	//if(currentIter_ < 100) regularIteration(true,true);
	//else regularIteration(false,true);
    //printTCritVal();
        //if(currentIter_ < 10 || (currentIter_%10==0)) objVal=findPrimalFeasSolnWith(z_current);
#if 1
	if(phase==0) regularIteration(updatePenalty,true);
	else if(phase==1) regularIteration(updatePenalty,true);
	else regularIteration(updatePenalty,false);
#endif
    }// main loop over currentIter_
    printStatus();

    //integralityDisc=evaluateIntegralityDiscrepancies();
    objVal=findPrimalFeasSolnWith(z_current);
    brVarInfo = findBranchingIndex();
    //if(mpiRank==0) cout << "Deciding whether to proceed to solve current node to higher precision..." << endl;
    while(brVarInfo.index<0 && incumbentVal - currentLagrLB > 1e-6 && currentIter_<5000){
      if(mpiRank==0) cout << "No branching, proceeding to solve current node to higher precision..." << endl;
      for(int kk=0; kk<100; kk++){
	regularIteration(updatePenalty,true);
	currentIter_++;
      }
      objVal=findPrimalFeasSolnWith(z_current);
      brVarInfo = findBranchingIndex();
    }
    return currentLagrLB;
}

double evaluateIntegralityDiscrepancies(){
    double integrDiscrSumLocal=0.0,integrDiscrSum=0.0;
    for(int tS=0; tS<nNodeSPs; tS++){
	integrDiscr_[tS] = subproblemSolvers[tS]->computeIntegralityDiscr();  
	integrDiscrSumLocal += integrDiscr_[tS];  
	//if(integrDiscr > 1e-10) cout << "Node " << mpiRank << " scenario " << tS << " integrality discrepancy: " << integrDiscr << endl;  
    }
    #ifdef USING_MPI
    if(mpiSize>1){
        MPI_Allreduce(&integrDiscrSumLocal, &integrDiscrSum, 1, MPI_DOUBLE, MPI_SUM, comm_);
    }
    #endif
    if(mpiSize==1){
        integrDiscrSum=integrDiscrSumLocal;
    }

    return integrDiscrSum;
}

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


double getLagrBound(){return currentLagrLB;}
void setBound(double bd){currentLagrLB=bd;}
void setReferenceLB(double bd){referenceLagrLB=bd;}

#if 0
void solveForWeights(){
   for (int tS = 0; tS < nNodeSPs; tS++) {
        if (subproblemSolvers[tS]->getNVertices()>2) { //Compute Next X, Y (With History)
	     //subproblemSolvers[tS]->solveMPHistory(omega_current[tS],z_current,rho,scaling_matrix[tS],true);
	     subproblemSolvers[tS]->computeWeightsForCurrentSoln(NULL);
	}
	//if(mpiRank==0){subproblemSolvers[tS]->printWeights();}
   }
}
#endif




#if 0
void solveQMIPsGS(){
  for(int ii=0; ii<10; ii++){
    for (int tS = 0; tS < nNodeSPs; tS++) {
	subproblemSolvers[tS]->solveAugmentedLagrangianMIP(omega_current[tS], z_current, rho,scaling_matrix[tS]);
	subproblemSolvers[tS]->setXToVertex();
	subproblemSolvers[tS]->setYToVertex();
    }
    updateZ(true);
  }
}
#endif

int performColGenStep();
int performColGenStepBasic();


void updateZ(){
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
        delete [] penSumLocal;
        delete [] penSum;
}

#if 1
void averageOfBestVertices(double *z=NULL, double *weightVec=NULL){
	if(z==NULL) z=z_average;
	double weight=0.0;
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += subproblemSolvers[tS]->getXVertex()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
	    }
	}
	if(weightVec!=NULL){
	  for (int tS = 0; tS < nNodeSPs; tS++) weight += weightVec[tS];
	}
	else{
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += subproblemSolvers[tS]->getXVertex()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
	    }
	}
	  for (int tS = 0; tS < nNodeSPs; tS++) weight += 1.0;
	}
	computeAverageAcrossProcs(z_local, z, weight, n1);
}

void averageOfAll(double *z=NULL, double *weightVec=NULL){
	if(z==NULL) z=z_average;
	double weight=0.0;
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += subproblemSolvers[tS]->getXVertex()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
	    }
	}
	if(weightVec!=NULL){
	  for (int tS = 0; tS < nNodeSPs; tS++) weight += weightVec[tS];
	}
	else{
	for (int i = 0; i < n1; i++)
	{
	    z_local[i] = 0.0;
	    z[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		z_local[i] += subproblemSolvers[tS]->getXVertex()[i];//x_current[tS][i] * scaling_matrix[tS][i]*pr[tS];
	    }
	}
	  for (int tS = 0; tS < nNodeSPs; tS++) weight += 1.0;
	}
	computeAverageAcrossProcs(z_local, z, weight, n1);
}

void dispOfBestXVertices(double *retVec=NULL, double *aveVec=NULL, double *weightVec=NULL){
    double weightSum=0.0;
    if(retVec==NULL){retVec=z_dispersions;}
    if(aveVec==NULL){aveVec=z_average;}
    for(int ii=0; ii<n1; ii++){
	z_local[ii]=0.0;
    }
    if(weightVec!=NULL){
      for(int tS=0; tS<nNodeSPs; tS++){
        for(int ii=0; ii<n1; ii++){
    	     z_local[ii] += weightVec[tS]*fabs(aveVec[ii] - subproblemSolvers[tS]->getXVertex()[ii]);
	}
	weightSum += weightVec[tS];
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
    computeAverageAcrossProcs(z_local, retVec, weightSum, n1);
}

void dispOfAllXVertices(double *retVec=NULL){
    double *dispVec;
    if(retVec==NULL){retVec=z_dispersions;}
    for(int ii=0; ii<n1; ii++){
	z_local[ii]=0.0;
    }
    for(int tS=0; tS<nNodeSPs; tS++){
	dispVec = subproblemSolvers[tS]->computeDispersions();
        for(int ii=0; ii<n1; ii++){
    	     z_local[ii] += dispVec[ii];
	}
    }
    #ifdef USING_MPI
    if(mpiSize>1){
	MPI_Allreduce(z_local, retVec, n1, MPI_DOUBLE, MPI_SUM, comm_);
    }
    #endif
    if(mpiSize==1){
	for (int i=0; i<n1; i++) retVec[i] = z_local[i];
    }
    //computeAverageAcrossProcs(z_local, retVec, weightSum, n1);
}

#if 0
BranchingVarInfo findBranchingIndexOld(){
    double *dispVec;
    BranchingVarInfo localBranchInfo={-1.0,-1.0,-1.0,0.0,0.0,0.0};
    BranchingVarInfo branchInfo={-1.0,-1.0,-1.0,0.0,0.0,0.0};
    
    for(int tS=0; tS<nNodeSPs; tS++){
	integrDiscr_[tS] = subproblemSolvers[tS]->computeIntegralityDiscr();  
	if(integrDiscr_[tS] < 1e-10) integrDiscr_[tS]=0.0;
	dispVec = subproblemSolvers[tS]->computeDispersions();
	if(integrDiscr_[tS] > localBranchInfo.intDiscr){
	    localBranchInfo.intDiscr=integrDiscr_[tS];
	    localBranchInfo.index=-1;
	    localBranchInfo.disp=0.0;
	    for(int ii=0; ii<n1; ii++){
	      if(dispVec[ii] > localBranchInfo.disp){ 
		localBranchInfo.disp = dispVec[ii];
	        localBranchInfo.index = ii;
	        localBranchInfo.brVal = subproblemSolvers[tS]->getX()[ii];
	      }
	    }
	    
	}
    }
#ifdef USING_MPI
    MPI_Op myOp;
    MPI_Datatype mpitype;
    MPI_Type_contiguous(8, MPI_DOUBLE, &mpitype);
    MPI_Type_commit( &mpitype);
    MPI_Op_create( (MPI_User_function *) compareBranchInfoOld, true, &myOp ); 
    if(mpiSize>1){
	MPI_Allreduce(&localBranchInfo, &branchInfo, 1, mpitype, myOp, comm_);
    }
#endif
    if(mpiSize==1){
	branchInfo = localBranchInfo;
    }
    if(mpiRank==0){
        cout << "Branching on index: " << branchInfo.index << " disps: " << branchInfo.disp << " intDiscr " << branchInfo.intDiscr << " with value " << branchInfo.brVal << endl;
	cout << "Branching determined on process: " << branchInfo.rank << " subproblem " << branchInfo.scen << endl;
    }
    return branchInfo;
}
#endif
BranchingVarInfo getBranchingVarInfo(){return brVarInfo;}

BranchingVarInfo findBranchingIndex(){

    brVarInfo={-1.0,-1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double totalIntDiscr = evaluateIntegralityDiscrepancies();
    if(mpiRank==0) cout << "Total integrality discrepancy is: " << totalIntDiscr << endl;
    if(totalIntDiscr < 1e-10){
	if(mpiRank==0) cout << "Total integrality discrepancy is negligible, no branching..." << endl;
	return brVarInfo;
    }

    double *dispVec;
    BranchingVarInfo localBranchInfo={-1.0,-1.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0};
    dispOfAllXVertices(); //sets z_dispersions
    localBranchInfo.index=-1;
    localBranchInfo.disp=0.0;
    for(int ii=0; ii<n1; ii++){
	if(z_dispersions[ii] > localBranchInfo.disp){ 
	    localBranchInfo.index = ii;
	    localBranchInfo.disp = z_dispersions[ii];
	    localBranchInfo.brVal = z_current[ii];
	    if(subproblemSolvers[0]->getColTypes()[ii]==0){
	        localBranchInfo.brLBUp = localBranchInfo.brVal;
	        localBranchInfo.brUBUp = currentVarUB_[ii];
	        localBranchInfo.brLBDn = currentVarLB_[ii];
	        localBranchInfo.brUBDn = localBranchInfo.brVal;
	    }
	    else{
	        localBranchInfo.brLBUp = ceil(localBranchInfo.brVal);
	        localBranchInfo.brUBUp = currentVarUB_[ii];
	        localBranchInfo.brLBDn = currentVarLB_[ii];
	        localBranchInfo.brUBDn = floor(localBranchInfo.brVal);
	    }
        }
    }
    
    for(int tS=0; tS<nNodeSPs; tS++){
	integrDiscr_[tS] = subproblemSolvers[tS]->computeIntegralityDiscr();  
	if(integrDiscr_[tS] < 1e-10) integrDiscr_[tS]=0.0;
	    dispVec = subproblemSolvers[tS]->computeDispersions();
	    //localBranchInfo.intDiscr=integrDiscr_[tS];
	    for(int ii=n1; ii<n1+n2; ii++){
	      //if(!indexIsInt(ii)){continue;}
	      if(subproblemSolvers[tS]->getColTypes()[ii]==0){continue;}
	      if(dispVec[ii] > localBranchInfo.disp){ 
		localBranchInfo.rank = mpiRank;
		localBranchInfo.scen = tS;
		localBranchInfo.disp = dispVec[ii];
	        localBranchInfo.index = ii;
	        localBranchInfo.brVal = subproblemSolvers[tS]->getY()[ii-n1];
	    if(subproblemSolvers[tS]->getColTypes()[ii]==0){
	        localBranchInfo.brLBUp = localBranchInfo.brVal;
	        localBranchInfo.brUBUp = subproblemSolvers[tS]->getUB(ii);
	        localBranchInfo.brLBDn = subproblemSolvers[tS]->getLB(ii);
	        localBranchInfo.brUBDn = localBranchInfo.brVal;
	    }
	    else{
	        localBranchInfo.brLBUp = ceil(localBranchInfo.brVal);
	        localBranchInfo.brUBUp = subproblemSolvers[tS]->getUB(ii);
	        localBranchInfo.brLBDn = subproblemSolvers[tS]->getLB(ii);
	        localBranchInfo.brUBDn = floor(localBranchInfo.brVal);
	    }
	        //localBranchInfo.brLB = subproblemSolvers[tS]->getLB(ii);
	        //localBranchInfo.brUB = subproblemSolvers[tS]->getUB(ii);
	      }
	    }
	    
    }
#ifdef USING_MPI
    MPI_Op myOp;
    MPI_Datatype mpitype;
    MPI_Type_contiguous(9, MPI_DOUBLE, &mpitype);
    MPI_Type_commit( &mpitype);
    MPI_Op_create( (MPI_User_function *) compareBranchInfo, true, &myOp ); 
    if(mpiSize>1){
	MPI_Allreduce(&localBranchInfo, &brVarInfo, 1, mpitype, myOp, comm_);
    }
#endif
    if(mpiSize==1){
	brVarInfo = localBranchInfo;
    }
    if(mpiRank==0){
        cout << "Branching on (rank, sp, index): (" << brVarInfo.rank <<","<<brVarInfo.scen<<","<<brVarInfo.index<<") " << " disps: " << brVarInfo.disp << " with value " << brVarInfo.brVal << " and branches " 
	<< "(" << brVarInfo.brLBUp << "," << brVarInfo.brUBUp << ")" << " and "
	<< "(" << brVarInfo.brLBDn << "," << brVarInfo.brUBDn << ")" << endl;
    }
    return brVarInfo;
}

#endif
void computeAverageAcrossProcs(double *partialAveVec, double* retVec, double weight, int vecSize){
    double weightG=0.0;
    #ifdef USING_MPI
    if(mpiSize>1){
	MPI_Allreduce(partialAveVec, retVec, vecSize, MPI_DOUBLE, MPI_SUM, comm_);
	MPI_Allreduce(&weight, &weightG, 1, MPI_DOUBLE, MPI_SUM, comm_);
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
void dispOfX(double *retVec, double *aveVec=NULL, double *weights=NULL){
    double weightSum=0.0;
    if(aveVec==NULL) aveVec=z_current;
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
    computeAverageAcrossProcs(z_local, retVec, weightSum, n1);
}




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
	MPI_Allreduce(dispOmegaLocal, dispOmega, n1, MPI_DOUBLE, MPI_SUM, comm_);
    }
#endif
    if(mpiSize==1){
	memcpy(dispOmega,dispOmegaLocal,n1*sizeof(double));
    }
    delete [] dispOmegaLocal;
}



/** Check if a value is integer. */
bool checkInteger(double value) const {
    double integerTolerance = 1.0e-10;
    double nearest = floor(value + 0.5);
//cout << " " << fabs(value-nearest);
    if (fabs(value - nearest) <= integerTolerance) {
        return true;
    }
    else {
        return false;
    }
}

#if 1
bool firstStageIndexIsInt(int ii){
    const char *cTypes = subproblemSolvers[0]->getColTypes();  //Using only first-stage, so it doesn't matter which scenario is chosen.
    assert(ii>=0);
    assert(ii < n1);
    if(cTypes[ii]!=0){
	return true;
    }
    else{
	return false;
    }
}
#endif

#if 0
bool checkForCutoff(){
if(mpiRank==0){cout << getIncumbentVal() << " <= " << currentLagrLB << "???" << endl;}
    return (getIncumbentVal() <= currentLagrLB);
}
#endif

void roundCurrentZ(double *z=NULL){
    const char *cTypes = subproblemSolvers[0]->getColTypes(); //Using only first-stage, so it doesn't matter which scenario is chosen.
    if(z==NULL){z=z_current;}
//double *origVarLB_;
    for(int ii=0; ii<n1; ii++){
	if(z[ii] < origVarLB_[ii]){z[ii]=origVarLB_[ii];}
	if(z[ii] > origVarUB_[ii]){z[ii]=origVarUB_[ii];}
    }
    memcpy(z_rounded,z,n1*sizeof(double));
    double roundingDisc = 0.0;
    for(int ii=0; ii<n1; ii++) {
	if(cTypes[ii]!=0) {
	    z_rounded[ii] = round(z[ii]);
    	    roundingDisc += fabs(z_rounded[ii]-z[ii]);
	}
    }
}

bool evaluateFeasibleZ(){
    if(incumbentVal > objVal){ 
	incumbentVal = objVal;
if(mpiRank==0){
    cout << "New incumbent value: " << incumbentVal << " and its corresponding solution: " << endl;
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
bool statusIsFeasible(){return modelStatus_[SP_STATUS]!=SP_INFEAS;}
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
    tCritVal = ALVal + 0.5*discrepNorm  - centreLagrLB;
    shouldTerminate = (tCritVal<tCritParam);// && (discrepNorm < 1e-20);
    if(shouldTerminate){return 0.0;}
    else{
	return (trialLagrLB-centreLagrLB)/(ALVal + 0.5*discrepNorm - centreLagrLB);
    }
}
		
bool updateOmega(bool useSSC){

#if 0
	double stepLength = exp(SSCVal-1.0);
	//double centreLagrLBLocal = 0.0;
	for (int tS = 0; tS < nNodeSPs; tS++) {
	  for (int i = 0; i < n1; i++) {
	    omega_centre[tS][i] += stepLength*rho*scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	  }
	  //centreLagrLBLocal += pr[tS]*subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_centre[tS]);
	  recordKeeping[tS][1]=recordKeeping[tS][0];
	}
	centreLagrLB = (1.0-stepLength)*centreLagrLB + stepLength*trialLagrLB;
    	omegaUpdated_ = true;
#endif
#if 0
	#ifdef USING_MPI
	if (mpiSize > 1) {
		MPI_Allreduce(&centreLagrLBLocal, &centreLagrLB, 1, MPI_DOUBLE, MPI_SUM, comm_);
	}
	#endif
	if (mpiSize == 1) {
		centreLagrLB = centreLagrLBLocal;
	}
#endif

#if 1
	if(SSCVal >= SSCParam || !useSSC || currentIter_==0) {
	    for (int tS = 0; tS < nNodeSPs; tS++) {
		//memcpy(omega_centre[tS],omega_tilde[tS],n1*sizeof(double));
		//memcpy(omega_current[tS],omega_tilde[tS],n1*sizeof(double));
		memcpy(omega_centre[tS],omega_tilde[tS],n1*sizeof(double));
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
#endif
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
	MPI_Allreduce(omegaLocal, omegaSum, n1, MPI_DOUBLE, MPI_SUM, comm_);
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
	MPI_Allreduce(omegaLocal, omegaSum, n1, MPI_DOUBLE, MPI_SUM, comm_);
	MPI_Allreduce(&omegaFroNormLocal, &omegaFroNorm, 1, MPI_DOUBLE, MPI_SUM, comm_);
	MPI_Allreduce(&cFroNormLocal, &cFroNorm, 1, MPI_DOUBLE, MPI_SUM, comm_);
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
    //rho=max(p,baselineRho);
    rho=max(p,1e-6);
#if 0
    for (int tS = 0; tS < nNodeSPs; tS++) {
	//setPenalty(tS,rho);
    	//subproblemSolvers[tS]->setQuadraticTerm(rho);
        subproblemSolvers[tS]->setQuadraticTerm(rho,scaling_matrix[tS]);
    }
#endif
//if(mpiRank==0) cout << "Penalty is now: " << rho << endl;
}
#if 0
void setPenalty(int tS, double p){
    //double penalty = max(p,baselineRho);
    rho = max(p,1e-6);
    subproblemSolvers[tS]->setQuadraticTerm(rho);
#if 0
    for (int i = 0; i < n1; i++) {
        scaling_matrix[tS][i] = 1.0;
    }
#endif
//if(mpiRank==0) cout << "Penalty is now: " << rho << endl;
}
#endif
//This should normally be turned off, this is a rule for updating the 
//penalty from Kiwiel 2006 and Lubin et al.
double computeKiwielScalingFactor(){
    double limit = 10.0;
    if(fabs(1.0-SSCVal) > 1e-10){
        return (max( min( 0.5/(1.0-SSCVal),limit), 1.0/limit));
    }
    else{
	return limit;
    }
}
double computeKiwielPenaltyUpdate(){	
	double scalingFactor = computeKiwielScalingFactor();
	computeScalingPenaltyUpdate( scalingFactor );
	return scalingFactor;
}
void computeScalingPenaltyUpdate(double scaling){	
    rho *= scaling;
#if 0
    for (int tS = 0; tS < nNodeSPs; tS++) {
	//computeScalingPenaltyUpdate(tS,scaling);
        subproblemSolvers[tS]->setQuadraticTerm(rho,scaling_matrix[tS]);
    }
#endif
}
#if 0
void computeScalingPenaltyUpdate(int tS, double scaling){	
    rho *= scaling;
   // for (int i = 0; i < n1; i++) {
	//scaling_matrix[tS][i] *= scaling;
	//scaling_matrix[tS][i] = max( scaling_matrix[tS][i], baselineRho);
    //}
    subproblemSolvers[tS]->setQuadraticTerm(rho,scaling_matrix[tS]);
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
  int wid = 20;
  if(mpiRank==0){
	//cout << setfill( ' ' );
	//cout << "Iter " << currentIter_ << setw(wid) << setprecision(10) << "LagrLB " << currentLagrLB << setw(wid) << "ALVal " << ALVal << setw(wid) << setprecision(5) << "primDiscr " << discrepNorm
	//	<< setw(wid) << "SSCVal " << SSCVal << setw(wid) << "innSSCVal " << innerSSCVal << setw(wid) << "MP iters: " << noInnerSolves << endl;
	printf("Iter %6d: LagrLB: %-10.8gALVal: %-10.8gbestPr: %-10.8gpDisc: %-8.3gSSCVal: %-8.3ginnSSCVal: %-8.3gtcrit: %-8.3gnoInnerMP: %-5d\n", currentIter_, currentLagrLB, ALVal, incumbentVal,discrepNorm, SSCVal, innerSSCVal, tCritVal, noInnerSolves);
//	printf("Best node: %0.14g, Incumbent value: %0.14g\n", getBestNodeQuality(), getIncumbentVal());
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
#if 1
void printIntegralityViolations(){
  const char *cTypes = subproblemSolvers[0]->getColTypes(); //Using only first-stage, so it doesn't matter which scenario is chosen.
  if(mpiRank==0){
    for(int ii=0; ii<n1; ii++){
      if(cTypes[ii]!=0){
	if(!checkInteger(z_current[ii])){
	    cout << " (" << ii << "," << z_current[ii] << ")";
	}
      }
    }
    cout << endl;
  }
}
#endif


//*** Helper functions.
//Frobenius norm
double froNorm(double** array, int nS, int n1);
double froNorm(double* array, int n1);
//Deals with numerical imprecision from solver (MIPs)
//double computeNormedDiscrepancies(double rho, double** scaling_matrix, double** x, double* z, int nS, int n1);
//void computeFeasibleStartingPoint(int tS, double* x, double* yFeasible);
};
#endif
/*Header for the main procedure ParallelSCG. */


