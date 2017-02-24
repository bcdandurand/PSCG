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

class PSCGModelScen{
public:
PSCGModelScen():n1(0),n2(0),nS(0),tS(-1),initialised(false),env(),disableHeuristic(false),nThreads(0),
x(NULL),y(NULL),c(NULL),d(NULL),cplexMP(env),x_vertex(NULL),y_vertex(NULL),weightSoln(env),weightObjective(env),quadraticTerm(env,0.0),nVertices(0),LagrBd(0.0),
mpModel(env),mpObjective(env),mpWeightConstraints(env),mpVertexConstraints(env),mpWeightVariables(env),mpWeight0(env,0.0,1.0),mpAuxVariables(env),pr(0.0){;}

//copy constructor
PSCGModelScen(const PSCGModelScen &other):n1(0),n2(0),nS(0),tS(-1),initialised(false),env(),disableHeuristic(false),nThreads(0),
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
virtual void upBranchOnVar(int varIndex, double bound)=0;
virtual void downBranchOnVar(int varIndex, double bound)=0;
virtual void printColTypesFirstStage(){;}

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

};



class PSCGModelScen_SMPS : public PSCGModelScen{
public:
PSCGModelScen_SMPS():PSCGModelScen(),LagrMIPInterface_(NULL){;}

//copy constructor
PSCGModelScen_SMPS(const PSCGModelScen_SMPS &other):PSCGModelScen(other),LagrMIPInterface_(NULL){;}

int initialiseSMPS(PSCGParams *par, TssModel &smpsModel, int scenario);

virtual int solveLagrangianProblem(const double* omega);

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
virtual void printColTypesFirstStage(){
    for(int i=0; i<n1; i++){
	if(LagrMIPInterface_->getColType()[i]==0)
	    cout << " C";
	else if(LagrMIPInterface_->getColType()[i]==1)
	    cout << " B";
	else if(LagrMIPInterface_->getColType()[i]==2)
	    cout << " I";
    }
    cout << endl;
}

OsiCpxSolverInterface* getOSI(){return LagrMIPInterface_;}

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
