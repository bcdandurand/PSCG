#ifndef PSCGMODELSCEN_H
#define PSCGMODELSCEN_H

#include "StructureDefs.h"
#include <list>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "TssModel.h"
#include "OsiCpxSolverInterface.hpp"
#include "PSCGParams.h"
#include "ProblemDataBodur.h"

#define ptrModel IloModel*
#define ptrRange IloRange*
#define ptrptrRange IloRange**
#define ptrObjective IloObjective*
#define ptrNumVarArray IloNumVarArray*

#define ptrDouble double*

ILOSTLBEGIN

#if 1
enum SolverReturnStatus {
  PSCG_OK=0,
  PSCG_OPTIMAL=0,
  PSCG_OPT_INT_FEAS,
  PSCG_INT_FEAS,
  PSCG_ABANDONED,
  PSCG_PRIMAL_INF_INT, //z does not satisfy integrality
  PSCG_PRIMAL_INF_REC, //z does not have recourse
  PSCG_PRIMAL_INF_INT_REC, //z both does not satisfy integrality and does not have recourse
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

class PSCGModelScen{
public:
PSCGModelScen():n1(0),n2(0),nS(0),tS(-1),initialised(false),env(),disableHeuristic(false),nThreads(0),solverStatus_(0),
x(NULL),y(NULL),c(NULL),d(NULL),cplexMP(env),x_vertex(NULL),y_vertex(NULL),weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),nVertices(0),LagrBd(0.0),
mpModel(env),mpObjective(env),mpWeightConstraints(env),mpVertexConstraints(env),mpWeightVariables(env),mpWeight0(env,0.0,1.0),mpAuxVariables(env),pr(0.0){;}

//copy constructor
PSCGModelScen(const PSCGModelScen &other):n1(0),n2(0),nS(0),tS(-1),initialised(false),env(),disableHeuristic(false),nThreads(0),solverStatus_(0),
x(NULL),y(NULL),c(NULL),d(NULL),cplexMP(env),x_vertex(NULL),y_vertex(NULL),weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),nVertices(0),LagrBd(0.0),
mpModel(env),mpObjective(env),mpWeightConstraints(env),mpVertexConstraints(env),mpWeightVariables(env),mpWeight0(env,0.0,1.0),mpAuxVariables(env),pr(0.0){;}

//void initialiseBodur(ProblemDataBodur &pdBodur, SMIP_fileRequest* probSpecs, int scenario);
//int initialiseSMPS(SMIP_fileRequest* probSpecs, int scenario);
//int initialiseSMPS(SMIP_fileRequest* probSpecs, TssModel &smpsModel, int scenario);
void finishInitialisation();

~PSCGModelScen(){
  delete [] x;
  delete [] y;
  delete [] c;
  delete [] d;
  delete [] x_vertex;
  delete [] y_vertex;
  mpModel.end();
  mpObjective.end();
  mpWeightConstraints.end();
  mpVertexConstraints.end();
  mpWeight0.end();
  mpWeightVariables.end();
  mpAuxVariables.end();
  cplexMP.end();
  weightSoln.end();
}

virtual int solveLagrangianProblem(const double* omega)=0;
virtual int solveFeasibilityProblem()=0;
int getSolverStatus(){return solverStatus_;}
virtual int solveLagrangianWithXFixedToZ(const double *z, const double *omega, const double *origLBs, const double *origUBs, const char *colTypes)=0;
virtual int solveFeasibilityProblemWithXFixedToZ(const double *z, const double *origLBs, const double *origUBs, const char *colTypes)=0;
virtual void setSolverStatus()=0;
virtual void upBranchOnVar(int varIndex, double bound)=0;
virtual void downBranchOnVar(int varIndex, double bound)=0;
virtual void printColTypesFirstStage(){;}
virtual void printColBds(){;}

void updateVertexHistory(){
	for(int i=0; i<n1; i++) {
	   xVertices[i].push_back(x_vertex[i]);
	}

	for(int i=0; i<n2; i++) { 
	   yVertices[i].push_back(y_vertex[i]);
	}

	//upvhWeightVariables.add(IloNumVar(env, 0.0, 1.0));
	//upvhWeightConstraints[0].setLinearCoef(upvhWeightVariables[nVertices], 1.0);	
	mpWeightVariables.add(IloNumVar(mpWeightConstraints[0](1.0)));
	mpWeightVariables[nVertices].setBounds(0.0,1.0);
	
	for(int i=0; i<n1; i++) {
	    mpVertexConstraints[i].setLinearCoef(mpWeightVariables[nVertices], x_vertex[i]);
	}

	baseWeightObj.push_back(0.0);
	weightSoln.add(0.0);
	weightObjective.add(0.0);

	for (int i = 0; i < n1; i++) {
		baseWeightObj[nVertices] += x_vertex[i] * c[i];
	}

	for (int j = 0; j < n2; j++) {
		baseWeightObj[nVertices] += y_vertex[j] * d[j];
	}

	nVertices++;
}
void clearVertexHistory(){
  while(nVertices >0){
    removeBackVertex();
  }
#if 0
	for(int i=0; i<n1; i++) {
	   xVertices[i].clear();
	}

	for(int i=0; i<n2; i++) { 
	   yVertices[i].clear();
	}
	baseWeightObj.clear();
	mpWeightVariables//.clear();
	mpVertexConstraints[i]
	weightSoln//.clear();
	weightObjective//.clear();
	nVertices=0;

		mpWeightVariables[0].end();
		mpWeightVariables.remove(0);
		weightObjective.remove(0);
#endif
}
virtual bool updateSolnInfo()=0;

virtual void fixXToZ(const double *z)=0;
virtual void unfixX(const double* origLBs, const double* origUBs)=0;
virtual void unfixX(const double* origLBs, const double* origUBs, const char* colTypes)=0;

void fixWeightToZero(int index){
    if(index >=0 && index < nVertices){
	mpWeightVariables[index].setBounds(0.0,0.0);
    }
}
void unfixWeight(int index){
    if(index >=0 && index < nVertices){
	mpWeightVariables[index].setBounds(0.0,1.0);
    }
}


int getNVertices(){return nVertices;}

//Not tested, purpose is to modify quadratic master problem formulation to reflect removing of vertices.
void removeBackVertex() {
    if(nVertices>0) {
		for(int i=0;i<n1;i++) xVertices[i].erase(xVertices[i].begin());

		for(int j=0;j<n2;j++) yVertices[j].erase(yVertices[j].begin());

		//mpModel.remove(mpWeightVariables[0]);
		mpWeightVariables[0].end();
		mpWeightVariables.remove(0);
		weightObjective.remove(0);
		baseWeightObj.erase(baseWeightObj.begin());
		nVertices--;
    }
}

void setQuadraticTerm(const double scaling_const) {
	for (int i = 0; i < n1; i++) {
		quadraticTerm.setQuadCoef(mpAuxVariables[i], mpAuxVariables[i], 0.5 * scaling_const);
	}
}

void setQuadraticTerm(const double *scaling_vector) {
	for (int i = 0; i < n1; i++) {
		quadraticTerm.setQuadCoef(mpAuxVariables[i], mpAuxVariables[i], 0.5 * scaling_vector[i]);
	}
}

void updatePrimalVariables_OneScenario(const double *omega, const double *z, const double *scaling_vector); 
void updatePrimalVariablesHistory_OneScenario(const double *omega, const double *z);

//void getLagrangianGradient(SMIP_qu_getLagrangianGradient* question, SMIP_ans_getLagrangianGradient* answer);

double getDefaultPenaltyParameter();

int getNumVertices(){return nVertices;}
double * getX(){return x;}
double * getY(){return y;}
double * getXVertex(){return x_vertex;}
double * getYVertex(){return y_vertex;}

void setXToVertex(){for(int i=0;i<n1;i++) x[i]=x_vertex[i];}
void setYToVertex(){for(int i=0;i<n2;i++) y[i]=y_vertex[i];}

double getLagrBd(){return LagrBd;}
double getProbabilities(){return pr;}

void updateALValues(const double *omega, const double *z, const double *scaling_vector){
    ALVal = 0.0;
    sqrNormDiscr = 0.0;
    for(int ii=0; ii<n1; ii++) {
		sqrNormDiscr += scaling_vector[ii]*(x[ii]-z[ii])*(x[ii]-z[ii]);
		ALVal += (c[ii]+omega[ii])*x[ii];
    }
    ALVal += 0.5*sqrNormDiscr;
    
    for(int jj=0; jj<n2; jj++) {
		ALVal += d[jj]*y[jj];
    }
}
double getALVal(){return ALVal;}
double getSqrNormDiscr(){return sqrNormDiscr;}
virtual bool checkSolnForFeasibility(const double *soln, vector<double> &constrVec){
printColTypesFirstStage();
    return false;
}

protected:

int n1;
int n2;
int nS;
int tS; //scenario
bool initialised;

IloEnv env;
bool disableHeuristic;
int nThreads;
double pr;
double *c;
double *d;

double LagrBd;

double *x;
double *y;

int nVertices;
double* x_vertex;
double* y_vertex;
vector< vector<double> > xVertices; //row ordered
vector< vector<double> > yVertices;

IloCplex cplexMP;
IloModel mpModel;
IloNumArray weightSoln;//(env, nVertices);
IloNumArray weightObjective;
vector<double> baseWeightObj;
IloExpr quadraticTerm;
IloObjective mpObjective;
IloRangeArray mpWeightConstraints;
IloRangeArray mpVertexConstraints;
IloNumVar mpWeight0;
IloNumVarArray mpWeightVariables;
IloNumVarArray mpAuxVariables;

double ALVal;
double sqrNormDiscr;
int solverStatus_;

};



class PSCGModelScen_SMPS : public PSCGModelScen{
public:
PSCGModelScen_SMPS():PSCGModelScen(),LagrMIPInterface_(NULL){;}

//copy constructor
PSCGModelScen_SMPS(const PSCGModelScen_SMPS &other):PSCGModelScen(other),LagrMIPInterface_(NULL){;}

int initialiseSMPS(PSCGParams *par, TssModel &smpsModel, int scenario);

virtual int solveLagrangianProblem(const double* omega);
virtual int solveFeasibilityProblem();
virtual void setSolverStatus(){
	OsiCpxSolverInterface *osi = LagrMIPInterface_;

    	if (osi->isAbandoned()) {
#if 0
        std::cout << "BOUND: is abandoned" << std::endl;
#endif
        	solverStatus_ = PSCG_ABANDONED;
    	}
    	else if (osi->isProvenOptimal()) {
#if 0
        std::cout << "BOUND: is lp optimal" << std::endl;
#endif
        	solverStatus_ = PSCG_OPTIMAL;
        	//PSCGNodeDesc *desc = dynamic_cast<PSCGNodeDesc*>(desc_);


    	}
	else if (osi->isProvenPrimalInfeasible()) {
#if 0
        std::cout << "BOUND: is primal inf" << std::endl;
#endif
        	solverStatus_ = PSCG_PRIMAL_INF;
    	}
    	else if (osi->isProvenDualInfeasible()) {
#if 0
        std::cout << "BOUND: is dual inf" << std::endl;
#endif
        	solverStatus_ = PSCG_DUAL_INF;
    	}
    	else if (osi->isPrimalObjectiveLimitReached()) {
#if 0
        	std::cout << "BOUND: is primal limit" << std::endl;
#endif
        	solverStatus_ = PSCG_PRIMAL_LIM;
    	}
    	else if (osi->isDualObjectiveLimitReached()) {
#if 0
        std::cout << "BOUND: is dual limit" << std::endl;
#endif
        	solverStatus_ = PSCG_DUAL_LIM;
    	}
    	else if (osi->isIterationLimitReached()) {
#if 0
        std::cout << "BOUND: is iter limit" << std::endl;
#endif
        	solverStatus_ = PSCG_ITER_LIM;
    	}
    	else {//handle other cases
		int errCode = getCPLEXErrorStatus();
		switch(errCode){
		  case CPXMIP_INFEASIBLE:
		      solverStatus_ = PSCG_PRIMAL_INF; 
		      break;
		  default:
        	    std::cout << "UNKNOWN SOLVER STATUS" << std::endl;
		    std::cout << "CPLEX Error Code is: " << getCPLEXErrorStatus() << endl;
        	    assert(0);
		}
    	}
}

virtual void upBranchOnVar(int varIndex, double bound){
    if(varIndex >= 0 && varIndex < n1+n2){
	LagrMIPInterface_->setColLower(varIndex,bound);
	//cout << "Col bds at: " << varIndex << " are " << LagrMIPInterface_->getColLower() << " and " << LagrMIPInterface_->getColUpper() << endl;
    }
}
virtual void downBranchOnVar(int varIndex, double bound){
    if(varIndex >= 0 && varIndex < n1+n2){
	LagrMIPInterface_->setColUpper(varIndex,bound);
    }
}

const char* getColTypes(){
    return LagrMIPInterface_->getColType();
} 

int getCPLEXErrorStatus(){
	OsiCpxSolverInterface *osi = LagrMIPInterface_;
	return CPXgetstat(osi->getEnvironmentPtr(), osi->getLpPtr());
}

virtual void printColTypesFirstStage(){
    for(int i=0; i<n1; i++){
	printColTypeFirstStage(i);
    }
    cout << endl;
}
virtual void printColTypeFirstStage(int i){
	if(getColTypes()[i]==0)
	    cout << " C";
	else if(getColTypes()[i]==1)
	    cout << " B";
	else if(getColTypes()[i]==2)
	    cout << " I";
}
virtual void printColumnBound(int i){
	cout << " (" << LagrMIPInterface_->getColLower()[i] << "," << LagrMIPInterface_->getColUpper()[i] << ")";
}
virtual void printColBds(){
	for(int i=0; i<n1; i++) printColumnBound(i);
	cout << endl;
}

OsiCpxSolverInterface* getOSI(){return LagrMIPInterface_;}

virtual bool checkSolnForFeasibility(const double *soln, vector<double> &constrVec){
    const CoinPackedMatrix *mat = LagrMIPInterface_->getMatrixByCol();
    if(mat->getNumRows() > constrVec.size()) constrVec.resize(mat->getNumRows());
    double *constrVals = constrVec.data();
    const double* rowLHS = LagrMIPInterface_->getRowLower();
    const double* rowRHS = LagrMIPInterface_->getRowUpper();
    mat->times(soln,constrVals);
#if 0
for(int ii=0; ii<mat->getNumRows(); ii++) cout << rowLHS[ii] << " <= " << constrVals[ii] << " <= " << rowRHS[ii] << endl;
#endif
    for(int ii=0; ii<mat->getNumRows(); ii++){
	if( !( (rowLHS[ii] <=  constrVals[ii]+1.0e-6) && (constrVals[ii]-1.0e-6 <= rowRHS[ii]) ) ) 
	{
#if 1
for(int ii=0; ii<mat->getNumRows(); ii++) cout << rowLHS[ii] << " <= " << constrVals[ii] << " <= " << rowRHS[ii] << endl;
#endif
	    return false;
	}
    }
    return true;
}

virtual void fixXToZ(const double *z){
   for(int ii=0; ii<n1; ii++){
      LagrMIPInterface_->setContinuous(ii);
      LagrMIPInterface_->setColBounds(ii,z[ii],z[ii]);
   }
   return; 
}

//undo the fixing of fixXToZ(), for when there is no need to reset variable types
virtual void unfixX(const double* origLBs, const double* origUBs){
   for(int ii=0; ii<n1; ii++){
      LagrMIPInterface_->setColBounds(ii,origLBs[ii],origUBs[ii]);
   }
   //cout << endl;
   return;

}
//undo the fixing of fixXToZ(), includes resetting variable types
virtual void unfixX(const double* origLBs, const double* origUBs, const char* colTypes){
   for(int ii=0; ii<n1; ii++){
      //cout << "  " << colTypes[ii];
      if(colTypes[ii]=='I' || colTypes[ii]=='B'){
	 LagrMIPInterface_->setInteger(ii);
      }
      else{
	 LagrMIPInterface_->setContinuous(ii);
      }
      LagrMIPInterface_->setColBounds(ii,origLBs[ii],origUBs[ii]);
   }
   //cout << endl;
   return;
}

virtual int solveLagrangianWithXFixedToZ(const double *z, const double *omega, const double *origLBs, const double *origUBs, const char *colTypes){
//cout << "*********" << endl;
//printColTypesFirstStage();
//printColBds();
   fixXToZ(z);
//printColTypesFirstStage();
//printColBds();
   solverStatus_ = solveLagrangianProblem(omega);
   unfixX(origLBs, origUBs, colTypes); 
//printColTypesFirstStage();
//printColBds();
   return solverStatus_;
}
virtual int solveFeasibilityProblemWithXFixedToZ(const double *z, const double *origLBs, const double *origUBs, const char *colTypes){
   fixXToZ(z);
   solverStatus_ = solveFeasibilityProblem();
   unfixX(origLBs, origUBs, colTypes); 
   return solverStatus_;
}

virtual bool updateSolnInfo(){
	OsiCpxSolverInterface *osi = LagrMIPInterface_;
	if(solverStatus_==PSCG_OPTIMAL || solverStatus_==PSCG_ITER_LIM){	
		if(solverStatus_==PSCG_ITER_LIM) cerr << "Flagging: SMPS MIP solver iteration limit reached." << endl;
		const double* solution = osi->getColSolution();
		LagrBd = osi->getObjValue()*osi->getObjSense();
		memcpy(x_vertex,solution,n1*sizeof(double));
		memcpy(y_vertex,solution+n1,n2*sizeof(double));
	}
	else{
		cerr << "Flagging: SMPS MIP solver indicated isProvenOptimal() == false." << endl;
		cerr << "CPLEX error code: " << getCPLEXErrorStatus() << endl;
	        assert(solverStatus_==PSCG_OPTIMAL || solverStatus_==PSCG_ITER_LIM);	
	}
}

private:
OsiCpxSolverInterface *LagrMIPInterface_;

};

class PSCGModelScen_Bodur : public PSCGModelScen{
public:
PSCGModelScen_Bodur():PSCGModelScen(),
cplexMIP(env),xVariables(env),yVariables(env),slpModel(env),c_vec(env),d_vec(env),slpObjective(env){;}
PSCGModelScen_Bodur(const PSCGModelScen_Bodur &other):PSCGModelScen(other), 
cplexMIP(env),xVariables(env),yVariables(env),slpModel(env),c_vec(env),d_vec(env),slpObjective(env){;}

void initialiseBodur(PSCGParams *par, ProblemDataBodur &pdBodur, int scenario);
virtual int solveLagrangianProblem(const double* omega);
virtual int solveFeasibilityProblem();
virtual void setSolverStatus(){
    solverStatus_ = PSCG_OPTIMAL;
}
virtual int solveLagrangianWithXFixedToZ(const double *z, const double *omega, const double *origLBs, const double *origUBs, const char *colTypes){
   fixXToZ(z);
   int solveStatus = solveLagrangianProblem(omega);
   unfixX(origLBs, origUBs, colTypes); 
}
virtual int solveFeasibilityProblemWithXFixedToZ(const double *z, const double *origLBs, const double *origUBs, const char *colTypes){
   fixXToZ(z);
   int solveStatus = solveFeasibilityProblem();
   unfixX(origLBs, origUBs, colTypes); 
}
virtual void fixXToZ(const double *z){
//TODO: set the xVariables to be continuous
    for(int ii=0; ii<n1; ii++){
	xVariables[ii].setBounds(z[ii],z[ii]);
    }
}
virtual void unfixX(const double* origLBs, const double* origUBs){
    for(int ii=0; ii<n1; ii++){
	xVariables[ii].setBounds(origLBs[ii],origUBs[ii]);
    }
}
virtual void unfixX(const double* origLBs, const double* origUBs, const char* colTypes){
//TODO: set the xVariables back to their original types
    for(int ii=0; ii<n1; ii++){
	xVariables[ii].setBounds(origLBs[ii],origUBs[ii]);
    }
}



virtual void upBranchOnVar(int varIndex, double bound){
    if(varIndex >=0 && varIndex < n1){
	xVariables[varIndex].setLB(bound);
    }
    else if(varIndex >= n1 && varIndex < n1+n2){
	yVariables[varIndex-n1].setLB(bound);
    }
}

virtual void downBranchOnVar(int varIndex, double bound){
    if(varIndex >=0 && varIndex < n1){
	xVariables[varIndex].setUB(bound);
    }
    else if(varIndex >= n1 && varIndex < n1+n2){
	yVariables[varIndex-n1].setUB(bound);
    }
}

virtual bool updateSolnInfo(){
    assert(solverStatus_==PSCG_OPTIMAL || solverStatus_==PSCG_ITER_LIM);	
    if(solverStatus_==PSCG_OPTIMAL || solverStatus_==PSCG_ITER_LIM){	
	if(solverStatus_==PSCG_ITER_LIM) cerr << "Flagging: SMPS MIP solver indicated iteration limit reached." << endl;
	for(int ii=0; ii<n1; ii++) x_vertex[ii] = cplexMIP.getValue(xVariables[ii]);
	//for(int ii=0; ii<n1; ii++) cout << " (" << cplexMIP.getValue(xVariables[ii]) << ","<<x_vertex[ii] << ")";
	//cout << endl;
	for(int jj=0; jj<n2; jj++) y_vertex[jj] = cplexMIP.getValue(yVariables[jj]);
	LagrBd = cplexMIP.getObjValue();
    }
    else{
	cerr << "In updateSolnInfo(): Flagging: SMPS MIP solver indicated isProvenOptimal() == false." << endl;
	assert(0);
    }
}


private:
IloCplex cplexMIP;
IloNumVarArray xVariables;
IloNumVarArray yVariables;
IloModel slpModel;
IloNumArray c_vec;
IloNumArray d_vec;
IloObjective slpObjective;
};
#endif
