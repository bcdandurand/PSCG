#ifndef CPLEXSCGSOLVER_H
#define CPLEXSCGSOLVER_H

#include "StructureDefs.h"
#include <list>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "TssModel.h"
#include "OsiCpxSolverInterface.hpp"
#include "ProblemDataBodur.h"

#define ptrModel IloModel*
#define ptrRange IloRange*
#define ptrptrRange IloRange**
#define ptrObjective IloObjective*
#define ptrNumVarArray IloNumVarArray*

#define ptrDouble double*

ILOSTLBEGIN

class CPLEXSCGSolver{
public:
CPLEXSCGSolver():n1(0),n2(0),nS(0),tS(-1),initialised(false),env(),disableHeuristic(false),nThreads(0),
x(NULL),y(NULL),c(NULL),d(NULL),cplexMP(env),x_vertex(env),y_vertex(env),weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),nVertices(0),LagrBd(0.0),
mpModel(env),mpObjective(env),mpWeightConstraints(env),mpVertexConstraints(env),mpWeightVariables(env),mpWeight0(env,0.0,1.0),mpAuxVariables(env),pr(0.0){;}

//copy constructor
CPLEXSCGSolver(const CPLEXSCGSolver &other):n1(0),n2(0),nS(0),tS(-1),initialised(false),env(),disableHeuristic(false),nThreads(0),
x(NULL),y(NULL),c(NULL),d(NULL),cplexMP(env),x_vertex(env),y_vertex(env),weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),nVertices(0),LagrBd(0.0),
mpModel(env),mpObjective(env),mpWeightConstraints(env),mpVertexConstraints(env),mpWeightVariables(env),mpWeight0(env,0.0,1.0),mpAuxVariables(env),pr(0.0){;}

//void initialiseBodur(ProblemDataBodur &pdBodur, SMIP_fileRequest* probSpecs, int scenario);
//int initialiseSMPS(SMIP_fileRequest* probSpecs, int scenario);
//int initialiseSMPS(SMIP_fileRequest* probSpecs, TssModel &smpsModel, int scenario);
void finishInitialisation();

~CPLEXSCGSolver(){
  delete [] x;
  delete [] y;
  delete [] c;
  delete [] d;
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

//Not tested, purpose is to modify quadratic master problem formulation to reflect removing of vertices.
void removeBackVertex() {
    if(nVertices>0) {
		for(int i=0;i<n1;i++) xVertices[i].erase(xVertices[i].begin());

		for(int j=0;j<n2;j++) yVertices[j].erase(yVertices[j].begin());

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
IloNumArray x_vertex;
IloNumArray y_vertex;
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

};



class CPLEXSCGSolver_SMPS : public CPLEXSCGSolver{
public:
CPLEXSCGSolver_SMPS():CPLEXSCGSolver(),LagrMIPInterface_(NULL){;}

//copy constructor
CPLEXSCGSolver_SMPS(const CPLEXSCGSolver_SMPS &other):CPLEXSCGSolver(other),LagrMIPInterface_(NULL){;}

int initialiseSMPS(SMIP_fileRequest *request, TssModel &smpsModel, int scenario);

virtual int solveLagrangianProblem(const double* omega);

OsiCpxSolverInterface* getOSI(){return LagrMIPInterface_;}

private:
OsiCpxSolverInterface *LagrMIPInterface_;

};

class CPLEXSCGSolver_Bodur : public CPLEXSCGSolver{
public:
CPLEXSCGSolver_Bodur():CPLEXSCGSolver(),
cplexMIP(env),xVariables(env),yVariables(env),slpModel(env),c_vec(env),d_vec(env),slpObjective(env){;}
CPLEXSCGSolver_Bodur(const CPLEXSCGSolver_Bodur &other):CPLEXSCGSolver(other), 
cplexMIP(env),xVariables(env),yVariables(env),slpModel(env),c_vec(env),d_vec(env),slpObjective(env){;}

void initialiseBodur(ProblemDataBodur &pdBodur, SMIP_fileRequest* probSpecs, int scenario);
virtual int solveLagrangianProblem(const double* omega);

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
