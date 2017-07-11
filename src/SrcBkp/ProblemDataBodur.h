/**
*
*
*
*
*/

/* A SMIP program has the form:
min  c*x + d_s*y_s
s.t. A*x <= b
	 Aeq*x == beq
	 T_s*x + W_s*y_s <= h_s
	 Teq_s*x + Weq_s*y_s == heq_s
	 x in X
	 y_s in Y_s
	 
	 The _s indicates dependence on scenario.
*/

#ifndef ProblemDataBodur_H
#define ProblemDataBodur_H

#include <ilcplex/ilocplex.h>
#include "StructureDefs.h"
#include <iostream>
#include <fstream>
#include <string>

ILOSTLBEGIN

typedef IloArray<IloNumArray> TwoDMatrix;

typedef IloArray<IloNumVarArray> NumVarMatrix;

typedef TwoDMatrix* ptrTwoDMatrix;
typedef IloNumArray* ptrIloNumArray;
typedef IloNumVarArray* ptrIloNumVarArray;

class ProblemDataBodur {
public:
	void initialise( SMIP_fileRequest* request, SMIP_fileReply *reply);
	void cleanup(){
	    if(initialised){
		initialised=false;
		env.end();
	    }
	}
	
	//IloNumVarArray& get_x();
	IloNumArray& get_c(){return c;}
	TwoDMatrix& get_A(){return A;}
	IloNumArray& get_b(){return b;}
	
	//IloNumVarArray& get_y(int scenIndex);
	IloNumArray& get_d(int scenIndex){return d;}
	TwoDMatrix& get_T(int scenIndex){return T;}
	TwoDMatrix& get_W(int scenIndex){return W;}
	IloNumArray& get_h(int scenIndex){return hArray[scenIndex];}
	IloNum& get_p(int scenIndex);
	
	int readBodur(SMIP_fileRequest* request, SMIP_fileReply* reply);

	int get_n1(){return n1;}
	int get_n2(){return n2;}
	int get_nS(){return nS;}
	
//protected:
	
	ProblemDataBodur(int noScenarios):initialised(false),n1(0),n2(0),nS(noScenarios),env(),c(env),A(env),b(env),d(env),T(env),W(env),hArray(env),p(0){;}
	~ProblemDataBodur(){
	    if(initialised){
		initialised=false;
		env.end();
	    }
	}
	
	
	
private:
	bool initialised;

	int n1;
	int n2;
	int nS;
	
	IloEnv env;
	IloNumArray c;
	TwoDMatrix A;
	IloNumArray  b;
	IloNumArray  d;
	TwoDMatrix  T;
	TwoDMatrix  W;
	TwoDMatrix hArray;
	IloNum p;
	
};

#endif
