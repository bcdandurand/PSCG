#ifndef ProblemDataSMPS_H
#define ProblemDataSMPS_H

#define ILOUSESTL

#include "StructureDefs.h"

#include <string>
#include <cassert>
#include <iostream>

#include "CoinPragma.hpp"
#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"


class ProblemDataSMPS {
public:
	static ProblemDataSMPS* getSingleton();
	
	static void initialise(SMIP_fileRequest* request, SMIP_fileReply *reply);
	static void cleanup();
	
	int get_n1();
	int get_n2();
	int get_nS();
	double* get_p();
	
	OsiCpxSolverInterface* instantiateScenarioSolverInterface(unsigned int tS, double &prVal);
	//OsiCpxSolverInterface* getScenarioSolverInterface(unsigned int tS);
	//OsiCpxSolverInterface* getSubproblemSolverInterface(unsigned int tS);
	//double* getBasicObjective(int tS);
	
protected:
	
	ProblemDataSMPS();
	~ProblemDataSMPS();
	
private:

	static ProblemDataSMPS* singleton;
	static bool initialised;
	
	int n1;
	int n2;
	int nS;
	//double* pr;
	//double** basicObjective;
	double total_prob;
	
	SmiScnModel* smi;
	SmiCoreData * core;
	//std::vector<OsiCpxSolverInterface*> *lagrangianInterfaces;
	//std::vector<OsiCpxSolverInterface*> *subproblemInterfaces;
	
	
	int readSMPS(SMIP_fileRequest* request, SMIP_fileReply* reply);
};


#endif
