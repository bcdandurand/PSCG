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
//#include "CoinUtils.h"

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
double objVal;
double gapVal;

double *x;
double *dispersions;
double *y;

int nVertices;
int maxNVertices;
int oldestVertexIndex;
double* x_vertex;
double* y_vertex;
double* x_vertex_opt;
double* y_vertex_opt;
vector< vector<double> > xVertices; //row ordered
vector< vector<double> > yVertices;

IloCplex cplexMP;
IloModel mpModel;
IloNum weight0;
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


public:
PSCGModelScen():
n1(0),n2(0),nS(0),tS(-1),initialised(false),env(),disableHeuristic(false),nThreads(0),solverStatus_(0),
x(NULL),y(NULL),c(NULL),d(NULL),cplexMP(env),x_vertex(NULL),y_vertex(NULL),x_vertex_opt(NULL),y_vertex_opt(NULL),oldestVertexIndex(-1),
weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),nVertices(0),maxNVertices(0),LagrBd(-COIN_DBL_MAX),objVal(-COIN_DBL_MAX),
mpModel(env),mpObjective(env),mpWeightConstraints(env),mpVertexConstraints(env),
mpWeightVariables(env),mpWeight0(env,0.0,1.0),mpAuxVariables(env),pr(0.0){;}

//copy constructor
PSCGModelScen(const PSCGModelScen &other):
n1(0),n2(0),nS(0),tS(-1),initialised(false),env(),disableHeuristic(false),nThreads(0),solverStatus_(0),
x(NULL),y(NULL),c(NULL),d(NULL),cplexMP(env),x_vertex(NULL),y_vertex(NULL),oldestVertexIndex(-1),
weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),nVertices(0),maxNVertices(0),LagrBd(-COIN_DBL_MAX),objVal(-COIN_DBL_MAX),
mpModel(env),mpObjective(env),mpWeightConstraints(env),mpVertexConstraints(env),
mpWeightVariables(env),mpWeight0(env,0.0,1.0),mpAuxVariables(env),pr(0.0){;}

//void initialiseBodur(ProblemDataBodur &pdBodur, SMIP_fileRequest* probSpecs, int scenario);
//int initialiseSMPS(SMIP_fileRequest* probSpecs, int scenario);
//int initialiseSMPS(SMIP_fileRequest* probSpecs, TssModel &smpsModel, int scenario);
void finishInitialisation();

~PSCGModelScen(){
  delete [] x;
  delete [] dispersions;
  delete [] y;
  delete [] c;
  delete [] d;
  delete [] x_vertex;
  delete [] y_vertex;
  delete [] x_vertex_opt;
  delete [] y_vertex_opt;

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

virtual int solveLagrangianProblem(const double* omega=NULL)=0;
virtual int solveAugmentedLagrangianMIP(const double* omega, const double* z, const double* scal)=0;
virtual int solveFeasibilityProblem()=0;
virtual int getCPLEXErrorStatus(){
    cerr << "getCPLEXErrorStatus(): default implementation, does nothing, returning 0" << endl;
    return 0;
}
int getSolverStatus(){return solverStatus_;}
double* getC(){return c;}
double* getD(){return d;}
double* getDispersions(){return dispersions;}
virtual int solveLagrangianWithXFixedToZ(const double *z, const double *omega, const double *origLBs, const double *origUBs, const char *colTypes)=0;
//virtual int solveFeasibilityProblemWithXFixedToZ(const double *z, const double *origLBs, const double *origUBs, const char *colTypes)=0;
virtual void setSolverStatus()=0;
virtual void upBranchOnVar(int varIndex, double bound)=0;
virtual void downBranchOnVar(int varIndex, double bound)=0;
virtual void setLBs(const double *lbs, int nLBs){
cerr << "setLBs(): Default implementation does nothing." << endl;
}
virtual void setUBs(const double *ubs, int nUBs){
cerr << "setLBs(): Default implementation does nothing." << endl;
}
void updateOptSoln(){
   if(x_vertex_opt==NULL) x_vertex_opt = new double[n1];
   if(y_vertex_opt==NULL) y_vertex_opt = new double[n2];
   memcpy(x_vertex_opt,x_vertex,n1*sizeof(double));
   memcpy(y_vertex_opt,y_vertex,n2*sizeof(double));
}
double evaluateVertexOptSolution(const double* omega){
    double retVal=0.0;
    if(omega==NULL) {for(int ii=0; ii<n1; ii++) retVal+= (c[ii])*x_vertex_opt[ii];}
    else {for(int ii=0; ii<n1; ii++) retVal+= (c[ii]+omega[ii])*x_vertex_opt[ii];}
    for(int jj=0; jj<n2; jj++) retVal+= d[jj]*y_vertex_opt[jj];
    return retVal;
}
double evaluateVertexSolution(const double* omega){
    double retVal=0.0;
    if(omega==NULL){for(int ii=0; ii<n1; ii++) retVal+= (c[ii])*x_vertex[ii];}
    else{for(int ii=0; ii<n1; ii++) retVal+= (c[ii]+omega[ii])*x_vertex[ii];}
    for(int jj=0; jj<n2; jj++) retVal+= d[jj]*y_vertex[jj];
    return retVal;
}
double evaluateSolution(const double* omega){
    double retVal=0.0;
    if(omega==NULL){for(int ii=0; ii<n1; ii++) retVal+= (c[ii])*x[ii];}
    else{for(int ii=0; ii<n1; ii++) retVal+= (c[ii]+omega[ii])*x[ii];}
    for(int jj=0; jj<n2; jj++) retVal+= d[jj]*y[jj];
    return retVal;
}
double evaluateSolution(const double* omega, const double *z, const double *scaling_vector){
//TODO
    double retVal=0.0;
    if(omega==NULL){for(int ii=0; ii<n1; ii++) retVal+= (c[ii])*x[ii];}
    else{for(int ii=0; ii<n1; ii++) retVal+= (c[ii]+omega[ii])*x[ii];}
    for(int jj=0; jj<n2; jj++) retVal+= d[jj]*y[jj];
    return retVal;
}
void evaluateVertexHistory(const double *omega){
    double retVal;
    for(int vv=0; vv<nVertices; vv++){
    	retVal=0.0;
    	if(omega==NULL){for(int ii=0; ii<n1; ii++) retVal+= (c[ii])*xVertices[ii][vv];}
    	else{for(int ii=0; ii<n1; ii++) retVal+= (c[ii]+omega[ii])*xVertices[ii][vv];}
    	for(int jj=0; jj<n2; jj++) retVal+= d[jj]*yVertices[jj][vv];
	cout << "\tVertex " << vv << " results in a value of " << retVal << endl;
    }
//virtual void setColSolution(const double *colsol) = 0;
}

double getXVertexEntry(int entry, int vertex){
  return xVertices[entry][vertex];
}

void compareLagrBds(const double* omega){
   cout << LagrBd << " ***versus*** " << evaluateVertexSolution(omega) << endl;
}
virtual void printColTypesFirstStage(){;}
virtual void printColBds(){;}
void printC(){
cout << "Printing c: " << endl;
 for(int ii=0; ii<n1; ii++){
  cout << " " << c[ii];
 }
 cout << endl;
}
virtual void printLagrSoln(){
return;
}
#if 1
bool roundIfClose(double &toBeRounded){
    if(fabs(toBeRounded-floor(toBeRounded)) < 1e-6){ 
	toBeRounded = floor(toBeRounded);
	return true;
    }
    if(fabs(ceil(toBeRounded)-toBeRounded) < 1e-6){ 
	toBeRounded = ceil(toBeRounded);
	return true;
    }
    return false;
}
#endif

void printWeights(){
  cout << "Printing weight solutions in continuous MP: " << endl;
try{
  cout << weight0 << " ";
  for(int ii=0; ii<nVertices; ii++){
    cout << "  " << weightSoln[ii];
  }
  cout << endl;
}
catch(IloException& e){
  cout << "printWeights() error: " << e.getMessage() << endl;
  //cout << "Exception caught...Refreshing solution..." << endl;
  //refresh();
  e.end();
}

}

//x,y,x_vertex, and y_vertex should all be set to something meaningful
double updateGapVal(const double *omega){
  gapVal=0.0;
  if(omega==NULL){for(int ii=0; ii<n1; ii++) {gapVal-= (c[ii])*(x_vertex[ii]-x[ii]);}}
  else{for(int ii=0; ii<n1; ii++) {gapVal-= (c[ii]+omega[ii])*(x_vertex[ii]-x[ii]);}}
  for(int jj=0; jj<n2; jj++) {gapVal-= d[jj]*(y_vertex[jj]-y[jj]);}
  return gapVal;
}

double getGapVal(){return gapVal;}

#if 0
double computeXDispersion(const double *z){
    double disp = 0.0;
    for(int ii=0; ii<n1; ii++){
	disp+= (x_vertex[ii]-z[ii])*(x_vertex[ii]-z[ii]);
    }
    return disp;
}
#endif

void addVertex(){
    if(nVertices==maxNVertices){
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
	maxNVertices++;
    }
    else{ //nVertices<maxNVertices
	assert(nVertices < maxNVertices);
	replaceVertexAtIndex(nVertices);
	nVertices++;
    }
    if(oldestVertexIndex==-1) oldestVertexIndex=0;
}

void replaceVertexAtIndex(int iii){

	//mpWeightVariables.add(IloNumVar(mpWeightConstraints[0](1.0)));
	//mpWeightVariables[nVertices].setBounds(0.0,1.0);
	assert(iii>=0);
	assert(iii<maxNVertices);
	for(int i=0; i<n1; i++) {
	   xVertices[i][iii] = x_vertex[i];
	}
	for(int i=0; i<n2; i++) { 
	   yVertices[i][iii] = y_vertex[i];
	}

	for(int i=0; i<n1; i++) {
	    mpVertexConstraints[i].setLinearCoef(mpWeightVariables[iii], x_vertex[i]);
	}

	baseWeightObj[iii]=0.0; //.push_back(0.0);
	//weightSoln.add(0.0);
	weightObjective[iii]=0.0;//.add(0.0);

	for (int i = 0; i < n1; i++) {
		baseWeightObj[iii] += x_vertex[i] * c[i];
	}

	for (int j = 0; j < n2; j++) {
		baseWeightObj[iii] += y_vertex[j] * d[j];
	}
	mpWeightVariables[iii].setBounds(0.0,1.0);

}

void replaceOldestVertex(){
//cout << "Oldest vertex index was: " << oldestVertexIndex;
    replaceVertexAtIndex(oldestVertexIndex);
    oldestVertexIndex = (++oldestVertexIndex)%nVertices;
//cout << " ...but is now " << oldestVertexIndex << endl;
}

void zeroOutVertexAtIndex(int iii){
#if 0
	for(int i=0; i<n1; i++) {
	    mpVertexConstraints[i].setLinearCoef(mpWeightVariables[iii], 0.0);
	}
#endif
	//baseWeightObj[iii]=0.0 //.push_back(0.0);
	//weightSoln.add(0.0);
	//weightObjective[iii]=0.0//.add(0.0);
	mpWeightVariables[iii].setBounds(0.0,0.0);
        weightSoln[iii]=0.0;//(env, nVertices);
}

void clearVertexHistory(){
  for(int v=0; v<maxNVertices; v++) zeroOutVertexAtIndex(v);
  oldestVertexIndex=-1;
  nVertices=0;
  for(int ii=0; ii<n1; ii++){ dispersions[ii]=0.0;}
}

virtual void printLinCoeffs(){
return;
}

void printVertices(){
  assert(nVertices>0);
cout << "Printing Vertices: " << endl;
  for(int ii=0; ii<n1; ii++){
      for(int vv=0; vv<nVertices; vv++){
	cout << " " << xVertices[ii][vv];
      }
      cout << endl;
  }
}

void printX(){
 assert(n1>0);
 for(int ii=0; ii<n1; ii++){
  cout << " " << x[ii];
 }
 cout << endl;
}

void printY(){
 assert(n2>0);
 for(int jj=0; jj<n2; jj++){
  cout << " " << y[jj];
 }
 cout << endl;
}

virtual void printXBounds(){
    cout << "printXBounds() default: doing nothing..." << endl;
}
virtual void printYBounds(){
    cout << "printYBounds() default: doing nothing..." << endl;
}

#if 0
void clearVertexHistory(){
  while(nVertices >0){
    removeBackVertex();
  }
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
}
#endif
virtual void updateSolnInfo()=0;

virtual void fixVarAt(int index, double fixVal){
    cout << "fixVarAt() is doing nothing..." << endl;
}
//virtual void fixXToZ(const double *z)=0;
virtual void fixXToZ(const double *z, const char* colTypes)=0;
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

void updatePrimalVariables_OneScenario(const double *omega, const double *z, const double *scaling_vector, bool updateDisp=false); 
void updatePrimalVariablesHistory_OneScenario(const double *omega, const double *z, const double *scaling_vector, bool updateDisp=false);

//void getLagrangianGradient(SMIP_qu_getLagrangianGradient* question, SMIP_ans_getLagrangianGradient* answer);

double getDefaultPenaltyParameter();

int getNumVertices(){return nVertices;}
double * getX(){return x;}
double * getY(){return y;}
double * getXVertex(){return x_vertex;}
double * getXVertexOpt(){return x_vertex_opt;}
double * getYVertex(){return y_vertex;}
double * getYVertexOpt(){return y_vertex_opt;}

void setXToVertex(){for(int i=0;i<n1;i++) x[i]=x_vertex[i];}
void setYToVertex(){for(int i=0;i<n2;i++) y[i]=y_vertex[i];}
void setXToOptVertex(){for(int i=0;i<n1;i++) x[i]=x_vertex_opt[i];}
void setYToOptVertex(){for(int i=0;i<n2;i++) y[i]=y_vertex_opt[i];}

void refresh(){
    clearVertexHistory();
    memcpy(x_vertex,x_vertex_opt,n1*sizeof(double));
    memcpy(y_vertex,y_vertex_opt,n2*sizeof(double));
    addVertex();
    
    setXToOptVertex();
    setYToOptVertex();
    weightSoln[0]=1.0;//(env, nVertices);
    for(int ii=1; ii<maxNVertices; ii++){weightSoln[ii]=0.0;}
}
void refresh(const double *omega, const double *z, const double *scaling_vector){
    
    //setXToOptVertex();
    //setYToOptVertex();
    updatePrimalVariables_OneScenario(omega, z, scaling_vector);

    clearVertexHistory();
    memcpy(x_vertex,x_vertex_opt,n1*sizeof(double));
    memcpy(y_vertex,y_vertex_opt,n2*sizeof(double));
    addVertex();

    //weightSoln[0]=1.0;//(env, nVertices);
    //for(int ii=1; ii<maxNVertices; ii++){weightSoln[ii]=0.0;}
}

virtual void polishSolution(){
cerr << "polishSolution(): Doing nothing." << endl;
}

void polishWeights(){
//assert(weight0 >=0 && weight0 <=1);
	double weightSum=0.0;
	if(weight0 >=0 && weight0 <=1){
    	    weightSum += weight0;
	}
	else if(weight0 > 1.0){
//cout << "Flagging: weight0 = " << weight0 << endl;
    	    weight0=1.0;
    	    weightSum += weight0;
	}
	else{
//cout << "Flagging: weight0 = " << weight0 << endl;
	    weight0=0.0;
	}

	for(int wI=0; wI<nVertices; wI++) {
	    if(weightSoln[wI] >=0 && weightSoln[wI] <=1){
    	        weightSum += weightSoln[wI];
	    }
	    else if(weightSoln[wI] > 1.0){
//cout << "Flagging: weight[" << wI << "] = " << weightSoln[wI] << endl;
    	        weightSoln[wI]=1.0;
    	        weightSum += weightSoln[wI];
	    }
	    else{
//cout << "Flagging: weight[" << wI << "] = " << weightSoln[wI] << endl;
		weightSoln[wI]=0.0;
	    }
	}
	assert(weightSum > 0);
	weight0 /= weightSum;
	for(int wI=0; wI<nVertices; wI++) {
	    weightSoln[wI] /=weightSum;
	}
}

double getLagrBd(){return LagrBd;}
double getProbabilities(){return pr;}

void updateALValues(const double *omega, const double *z, const double *scaling_vector){
    ALVal = 0.0;
    sqrNormDiscr = 0.0;
    for(int ii=0; ii<n1; ii++) {
		sqrNormDiscr += scaling_vector[ii]*(x[ii]-z[ii])*(x[ii]-z[ii]);
		if(omega==NULL){ALVal += (c[ii])*x[ii];}
		else{ALVal += (c[ii]+omega[ii])*x[ii];}
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


};



class PSCGModelScen_SMPS : public PSCGModelScen{
public:
PSCGModelScen_SMPS():PSCGModelScen(),LagrMIPInterface_(NULL){;}

//copy constructor
PSCGModelScen_SMPS(const PSCGModelScen_SMPS &other):PSCGModelScen(other),LagrMIPInterface_(NULL){;}

int initialiseSMPS(PSCGParams *par, TssModel &smpsModel, int scenario);

virtual void fixVarAt(int index, double fixVal){
      LagrMIPInterface_->setColBounds(index,fixVal,fixVal);
}

virtual int solveLagrangianProblem(const double* omega=NULL);
virtual int solveAugmentedLagrangianMIP(const double* omega, const double* z, const double* scal);
virtual int solveFeasibilityProblem();
virtual void polishSolution(){
  for(int ii=0; ii<n1; ii++){
      if(LagrMIPInterface_->getColLower()[ii] > x_vertex[ii]){
	x_vertex[ii] = LagrMIPInterface_->getColLower()[ii];
      }
      if(LagrMIPInterface_->getColUpper()[ii] < x_vertex[ii]){
	x_vertex[ii] = LagrMIPInterface_->getColUpper()[ii];
      }
      if(x_vertex[ii]==-0) x_vertex[ii]=0.0;
  }
  for(int ii=0; ii<n2; ii++){
      if(LagrMIPInterface_->getColLower()[n1+ii] > y_vertex[ii]){
	y_vertex[ii] = LagrMIPInterface_->getColLower()[n1+ii];
      }
      if(LagrMIPInterface_->getColUpper()[n1+ii] < y_vertex[ii]){
	y_vertex[ii] = LagrMIPInterface_->getColUpper()[n1+ii];
      }
      if(y_vertex[ii]==-0) y_vertex[ii]=0.0;
  }
}
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
virtual void setLBs(const double *lbs, int nLBs){
    for(int ii=0; ii<nLBs; ii++){
	LagrMIPInterface_->setColLower(ii,lbs[ii]);
    }
}
virtual void downBranchOnVar(int varIndex, double bound){
    if(varIndex >= 0 && varIndex < n1+n2){
	LagrMIPInterface_->setColUpper(varIndex,bound);
    }
}
virtual void setUBs(const double *ubs, int nUBs){
    for(int ii=0; ii<nUBs; ii++){
	LagrMIPInterface_->setColUpper(ii,ubs[ii]);
    }
}

const char* getColTypes(){
    return LagrMIPInterface_->getColType();
} 

virtual int getCPLEXErrorStatus(){
	OsiCpxSolverInterface *osi = LagrMIPInterface_;
	return CPXgetstat(osi->getEnvironmentPtr(), osi->getLpPtr());
}

int setGapTolerances(double absGap, double relGap){
	OsiCpxSolverInterface *osi = LagrMIPInterface_;
	CPXsetdblparam( osi->getEnvironmentPtr(), CPX_PARAM_EPAGAP, absGap );
	return CPXsetdblparam( osi->getEnvironmentPtr(), CPX_PARAM_EPGAP, relGap );
}

virtual void printXBounds(){
cout << "Printing subproblem bounds: " << endl;
  for(int ii=0; ii<n1; ii++){
      cout << " (" << LagrMIPInterface_->getColLower()[ii] << "," << LagrMIPInterface_->getColUpper()[ii] << ")";
  }
  cout << endl;
}
virtual void printYBounds(){
cout << "Printing subproblem bounds: " << endl;
  for(int ii=0; ii<n2; ii++){
      cout << " (" << LagrMIPInterface_->getColLower()[n1+ii] << "," << LagrMIPInterface_->getColUpper()[n1+ii] << ")";
  }
  cout << endl;
}
virtual void printLinCoeffs(){
cout << "Printing linear coefficients: " << endl;
  for(int ii=0; ii<n1+n2; ii++){
    cout << " " << LagrMIPInterface_->getObjCoefficients()[ii];
  }
  cout << endl;
}
virtual void printLagrSoln(){
cout << "Printing col solns: " << endl;
  const double* solution = LagrMIPInterface_->getColSolution();
  for(int ii=0; ii<n1+n2; ii++){
    cout << " " << solution[ii];
  }
  cout << endl;
}
double computeMIPVal(const double *omega){
  const double* solution = LagrMIPInterface_->getColSolution();
  const double* coeff = LagrMIPInterface_->getObjCoefficients();
  double retVal=0.0;
 if(omega==NULL){
  for(int ii=0; ii<n1+n2; ii++){
    retVal +=solution[ii]*coeff[ii];
  }
 }
 else{
  for(int ii=0; ii<n1; ii++){
    retVal +=solution[ii]*(coeff[ii]+omega[ii]);
  }
  for(int ii=n1; ii<n2; ii++){
    retVal +=solution[ii]*coeff[ii];
  }
 }
  return retVal;
}

int changeFromMIQPToMILP(){
	OsiCpxSolverInterface *osi = LagrMIPInterface_;
        if(CPXgetprobtype(osi->getEnvironmentPtr(),osi->getLpPtr()) == CPXPROB_MIQP){
	    return CPXchgprobtype(osi->getEnvironmentPtr(), osi->getLpPtr(), CPXPROB_MILP);
	}
	else{return 0;}

}
int changeFromMILPToMIQP(){
	OsiCpxSolverInterface *osi = LagrMIPInterface_;
        if(CPXgetprobtype(osi->getEnvironmentPtr(),osi->getLpPtr()) == CPXPROB_MILP){
	    return CPXchgprobtype(osi->getEnvironmentPtr(), osi->getLpPtr(), CPXPROB_MIQP);
	}
	else{
	    return 0;
	}
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
#if 0
for(int ii=0; ii<mat->getNumRows(); ii++) cout << rowLHS[ii] << " <= " << constrVals[ii] << " <= " << rowRHS[ii] << endl;
#endif
	    return false;
	}
    }
    return true;
}

#if 0
virtual void fixXToZ(const double *z){
   for(int ii=0; ii<n1; ii++){
      LagrMIPInterface_->setContinuous(ii);
      LagrMIPInterface_->setColBounds(ii,z[ii],z[ii]);
   }
   return; 
}
#endif
virtual void fixXToZ(const double *z, const char* colTypes){
   for(int ii=0; ii<n1; ii++){
      //LagrMIPInterface_->setContinuous(ii);
      if(colTypes[ii]=='I' || colTypes[ii]=='B'){
          LagrMIPInterface_->setColBounds(ii,round(z[ii]),round(z[ii]));
      }
      else{
          LagrMIPInterface_->setColBounds(ii,z[ii],z[ii]);
      }
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
   fixXToZ(z,colTypes);
//printColTypesFirstStage();
//printColBds();
   solverStatus_ = solveLagrangianProblem(omega);
   //unfixX(origLBs, origUBs, colTypes); 
   unfixX(origLBs, origUBs); 
//printColTypesFirstStage();
//printColBds();
   return solverStatus_;
}
#if 0
virtual int solveFeasibilityProblemWithXFixedToZ(const double *z, const double *origLBs, const double *origUBs, const char *colTypes){
   fixXToZ(z,colTypes);
   solverStatus_ = solveFeasibilityProblem();
   unfixX(origLBs, origUBs); 
   return solverStatus_;
}
#endif

virtual void updateSolnInfo(){
	OsiCpxSolverInterface *osi = LagrMIPInterface_;
	if(solverStatus_==PSCG_OPTIMAL || solverStatus_==PSCG_ITER_LIM){	
		if(solverStatus_==PSCG_ITER_LIM) cerr << "Flagging: SMPS MIP solver iteration limit reached." << endl;
		const double* solution = osi->getColSolution();
		//LagrBd = osi->getObjValue()*osi->getObjSense();
		memcpy(x_vertex,solution,n1*sizeof(double));
		memcpy(y_vertex,solution+n1,n2*sizeof(double));
		//for(int ii=0; ii<n1; ii++){ roundIfClose(x_vertex[ii]);}
		//for(int ii=0; ii<n2; ii++){ roundIfClose(y_vertex[ii]);}
		//for(int ii=0; ii<n1; ii++){ (x_vertex[ii]);}
		//for(int ii=0; ii<n2; ii++){ (y_vertex[ii]);}
		polishSolution();
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
virtual int solveLagrangianProblem(const double* omega=NULL);
virtual int solveAugmentedLagrangianMIP(const double* omega, const double* z, const double* scal);
virtual int solveFeasibilityProblem();
virtual void setSolverStatus(){
    solverStatus_ = PSCG_OPTIMAL;
}
virtual int solveLagrangianWithXFixedToZ(const double *z, const double *omega, const double *origLBs, const double *origUBs, const char *colTypes){
   fixXToZ(z,colTypes);
   int solveStatus = solveLagrangianProblem(omega);
   unfixX(origLBs, origUBs); 
   return solveStatus;
}
#if 0
virtual int solveFeasibilityProblemWithXFixedToZ(const double *z, const double *origLBs, const double *origUBs, const char *colTypes){
   fixXToZ(z,colTypes);
   int solveStatus = solveFeasibilityProblem();
   unfixX(origLBs, origUBs); 
   return solveStatus;
}
virtual void fixXToZ(const double *z){
//TODO: set the xVariables to be continuous
    for(int ii=0; ii<n1; ii++){
	xVariables[ii].setBounds(z[ii],z[ii]);
    }
}
#endif
virtual void fixXToZ(const double *z, const char* colTypes){;
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

virtual void updateSolnInfo(){
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
