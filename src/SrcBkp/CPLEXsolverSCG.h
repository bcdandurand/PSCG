#ifndef CPLEXSOLVER_SCG_H
#define CPLEXSOLVER_SCG_H

#include "StructureDefs.h"
#include "ProblemDataBodur.h"
//#include "ProblemDataSMPS.h"
#include "TssModel.h"
#include "OsiCpxSolverInterface.hpp"
#include <list>
#include <iostream>
#include <ilcplex/ilocplex.h>


#define ptrModel IloModel*
#define ptrRange IloRange*
#define ptrptrRange IloRange**
#define ptrObjective IloObjective*
#define ptrNumVarArray IloNumVarArray*

#define ptrDouble double*

ILOSTLBEGIN

//extern "C" {
class CPLEXsolverSCG{
public:
CPLEXsolverSCG():n1(0),n2(0),nS(0),c(NULL),d(NULL),tS(-1),initialised(false),env(),
LagrMIPInterface(NULL),filetype(0),disableHeuristic(false),nThreads(0),
x(NULL),y(NULL),cplexMIP(env),cplexQP(env),weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),xVariables(env),yVariables(env),x_vertex(env),y_vertex(env),
xVertices(),yVertices(),nVertices(0),LagrBd(0.0),
slpModel(env),p_times_c(env),p_times_d(env),slpObjective(env),
sfssModels(NULL),sfssObjectives(NULL),
upvhModel(env),upvhObjective(env),upvhWeightConstraints(env),upvhVertexConstraints(env),upvhWeightVariables(env),upvhWeight0(env,0.0,1.0),upvhDummyVariables(env),pr(0.0){;}

//copy constructor
CPLEXsolverSCG(const CPLEXsolverSCG &other):n1(0),n2(0),nS(0),c(NULL),d(NULL),tS(-1),initialised(false),env(),
LagrMIPInterface(NULL),filetype(0),disableHeuristic(false),nThreads(0),
x(NULL),y(NULL),cplexMIP(env),cplexQP(env),weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),xVariables(env),yVariables(env),x_vertex(env),y_vertex(env),
xVertices(),yVertices(),nVertices(0),LagrBd(0.0),
slpModel(env),p_times_c(env),p_times_d(env),slpObjective(env),
sfssModels(NULL),sfssObjectives(NULL),
upvhModel(env),upvhObjective(env),upvhWeightConstraints(env),upvhVertexConstraints(env),upvhWeightVariables(env),upvhWeight0(env,0.0,1.0),upvhDummyVariables(env),pr(0.0){;}

void initialiseBodur(ProblemDataBodur &pdBodur, SMIP_fileRequest* probSpecs, int scenario);
//int initialiseSMPS(SMIP_fileRequest* probSpecs, int scenario);
int initialiseSMPS(SMIP_fileRequest* probSpecs, TssModel &smpsModel, int scenario);
void finishInitialisation();

~CPLEXsolverSCG(){
    x_vertex.end();
    y_vertex.end();
    slpModel.end();
    slpObjective.end();
    p_times_c.end();
    p_times_d.end();
  
  if(initialised) {
    delete [] x;
    delete [] y;
    delete [] c;
    delete [] d;
  }

  for(int i=0; i<n1; i++) xVertices[i].end();

  upvhModel.end();
  upvhObjective.end();
  upvhWeightConstraints.end();
  upvhVertexConstraints.end();
  upvhWeight0.end();
  upvhWeightVariables.end();
  upvhDummyVariables.end();
  cplexMIP.end();
  cplexQP.end();
  weightSoln.end();
}

int solveLagrangianProblem(const double* dual_var);
int solveLagrangianProblem_Bodur(const double* dual_var);
int solveLagrangianProblem_SMPS(const double* dual_var);

void updateVertexHistory(){
	for(int i=0; i<n1; i++) {
	   xVertices[i].add(x_vertex[i]);
	}

	for(int i=0; i<n2; i++) { 
	   yVertices[i].push_back(y_vertex[i]);
	}

	//upvhWeightVariables.add(IloNumVar(env, 0.0, 1.0));
	//upvhWeightConstraints[0].setLinearCoef(upvhWeightVariables[nVertices], 1.0);	
	upvhWeightVariables.add(IloNumVar(upvhWeightConstraints[0](1.0)));
	upvhWeightVariables[nVertices].setBounds(0.0,1.0);
	
	for(int i=0; i<n1; i++) {
	    upvhVertexConstraints[i].setLinearCoef(upvhWeightVariables[nVertices], x_vertex[i]);
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

//Not tested, purpose is to modify quadratic master problem formulation to reflect removing of vertices.
void removeBackVertex() {
    if(nVertices>0) {
		for(int i=0;i<n1;i++) xVertices[i].remove(0);

		for(int j=0;j<n2;j++) yVertices[j].erase(yVertices[j].begin());

		upvhWeightVariables.remove(0);
		weightObjective.remove(0);
		baseWeightObj.erase(baseWeightObj.begin());
		nVertices--;
    }
}

void setQuadraticTerm(const double scaling_const) {
	for (int i = 0; i < n1; i++) {
		quadraticTerm.setQuadCoef(upvhDummyVariables[i], upvhDummyVariables[i], 0.5 * scaling_const);
	}
}

void setQuadraticTerm(const double *scaling_vector) {
	for (int i = 0; i < n1; i++) {
		quadraticTerm.setQuadCoef(upvhDummyVariables[i], upvhDummyVariables[i], 0.5 * scaling_vector[i]);
	}
}

void updatePrimalVariables_OneScenario(const double *dvar, const double *z, const double *scaling_vector); 
void updatePrimalVariablesHistory_OneScenario(const double *omega, const double *z);

void getLagrangianGradient(SMIP_qu_getLagrangianGradient* question, SMIP_ans_getLagrangianGradient* answer);

double getDefaultPenaltyParameter();

int getNumVertices(){return nVertices;}
double * getX(){return x;}
double * getY(){return y;}

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

private:

int n1;
int n2;
int nS;
int tS; //scenario
bool initialised;

IloEnv env;
int filetype;
bool disableHeuristic;
int nThreads;
double pr;
double *c;
double *d;

IloCplex cplexMIP;
IloNumVarArray xVariables;
IloNumVarArray yVariables;
IloNumArray x_vertex;
IloNumArray y_vertex;
IloModel slpModel;
IloNumArray p_times_c;
IloNumArray p_times_d;
IloObjective slpObjective;
OsiCpxSolverInterface *LagrMIPInterface;
double LagrBd;

// Related to solving 2nd stage problem given fixed 1st stage solution.
IloModel* sfssModels;
IloObjective* sfssObjectives;

IloCplex cplexQP;
IloNumArray weightSoln;//(env, nVertices);
IloNumArray weightObjective;
vector<double> baseWeightObj;
IloExpr quadraticTerm;

double *x;
double *y;
vector<IloNumArray> xVertices; //row ordered
vector< vector<double> > yVertices;

int nVertices;
IloModel upvhModel;
IloObjective upvhObjective;
IloRangeArray upvhWeightConstraints;
IloRangeArray upvhVertexConstraints;
IloNumVar upvhWeight0;
IloNumVarArray upvhWeightVariables;
IloNumVarArray upvhDummyVariables;
double ALVal;
double sqrNormDiscr;

};
#endif
