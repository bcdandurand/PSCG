/**
*
*
*
*
*/

#include "ProblemDataSMPS.h"
//#include "TssModel.h"

#include <ilcplex/ilocplex.h>

typedef std::vector<OsiCpxSolverInterface*> InterfacePtrVector;

//InterfacePtrVector *lagrangianInterfaces = 0;
ProblemDataSMPS* ProblemDataSMPS::singleton = 0;
bool ProblemDataSMPS::initialised = false;

// Returns the singleton ProblemData. In theory this could become a collection of singletons, all to be solved in one execution.
ProblemDataSMPS* ProblemDataSMPS::getSingleton() {
	return ProblemDataSMPS::singleton;
}

ProblemDataSMPS::ProblemDataSMPS() {
	smi = new SmiScnModel();
}

ProblemDataSMPS::~ProblemDataSMPS() {
	delete smi;
	delete core;
}

// Creates the singleton and reads an SMPS instance into it.
void ProblemDataSMPS::initialise(SMIP_fileRequest* request, SMIP_fileReply* reply) {
	ProblemDataSMPS::initialised = true;

	ProblemDataSMPS::singleton = new ProblemDataSMPS();
	
	switch( request->ftype ) 
	{
	case 2:
		ProblemDataSMPS::singleton->readSMPS(request, reply);
		break;
	default:
		throw(-1);
		break;
	}
}

int ProblemDataSMPS::readSMPS(SMIP_fileRequest* request, SMIP_fileReply* reply) {

	std::string filename(request->filename);
	std::string filepath(request->filepath);
	
	std::string full_filename;
	full_filename += filepath;
	full_filename += filename;

	cout << full_filename.c_str() << endl;
	
	smi->readSmps(full_filename.c_str());
	int FileNoScenarios = smi->getNumScenarios();
	
	if (FileNoScenarios < request->totalScenarios)
	{
		cout << "More scenarios requested than available in problem; aborting" << endl;
		throw(-1);
	}
	
	if (request->totalScenarios > 0) 
	{
		nS = request->totalScenarios;
	}
	else
	{
		nS = FileNoScenarios;
	}
	//pr = new double[nS];
	total_prob = 0;
	
	
	//Unclear what of the following is necessary
	///////////////////////////////////
	int *intInds = smi->getIntegerInd();
	int nInts = smi->getIntegerLen();
	int *binInds = smi->getBinaryInd();
	int nBins = smi->getBinaryLen();
	core = smi->getCore();
	const int nStages = core->getNumStages();
	vector<int> stageBegins(nStages);
	vector<int> stageEnds(nStages);
	for(int s=0; s<nStages-1; s++){
		stageBegins[s]=core->getColStart(s);
		stageEnds[s]=core->getColStart(s+1);
	}
	stageBegins[nStages-1]=core->getColStart(nStages-1);
	stageEnds[nStages-1]=core->getNumCols();


	n1 = stageEnds[0] - stageBegins[0];
	n2 = stageEnds[1] - stageBegins[1];
	///////////////////////////////////
	
	reply->n1 = n1;
	reply->n2 = n2;
	reply->nS = nS;
	
	reply->fileLoaded = true;
	
	return 1;
}

OsiCpxSolverInterface* ProblemDataSMPS::instantiateScenarioSolverInterface(unsigned int tS, double &prVal){
	OsiCpxSolverInterface *cpxOsi = new OsiCpxSolverInterface();
	cpxOsi = dynamic_cast<OsiCpxSolverInterface*>(smi->addSingleScenarioSubproblem(cpxOsi, tS, 1.0, prVal)); 
	total_prob += prVal;
	return cpxOsi;
}


int ProblemDataSMPS::get_n1() {return n1;}
int ProblemDataSMPS::get_n2() {return n2;}
int ProblemDataSMPS::get_nS() {return nS;}
//double* ProblemDataSMPS::get_p() {return pr;}

void ProblemDataSMPS::cleanup() {
	// TODO: a lot of stuff is not deleted!
	if (ProblemDataSMPS::initialised) {
		delete ProblemDataSMPS::singleton;
		ProblemDataSMPS::initialised = false;
	}
		
}
