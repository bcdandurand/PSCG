/**
*
*
*
*
*/

#include "ProblemDataBodur.h"

ILOSTLBEGIN

//ProblemData* ProblemData::singleton = 0;
//bool ProblemData::initialised = false;

// Returns the singleton ProblemData. In theory this could become a collection of singletons, all to be solved in one execution.
#if 0
ProblemData* ProblemData::getSingleton() {
	return ProblemData::singleton;
}

// Constructor.
ProblemData::ProblemData(int noScenarios) {
	nS = noScenarios;
	
	y = new StochasticParameter<IloNumVarArray*>(noScenarios);
	d = new StochasticParameter<IloNumArray*>(noScenarios);
	T = new StochasticParameter<TwoDMatrix*>(noScenarios);
	W = new StochasticParameter<TwoDMatrix*>(noScenarios);
	h = new StochasticParameter<IloNumArray*>(noScenarios);
	p = new StochasticParameter<double*>(noScenarios);
}

// Destructor.
ProblemData::~ProblemData() {
	delete y;
	delete d;
	delete T;
	delete W;
	delete h;
}
#endif

// Creates the singleton and reads a Bodur instance into it. (Could be other types at some point.)
void ProblemDataBodur::initialise(SMIP_fileRequest* request, SMIP_fileReply* reply) {
	//env = iloenv;
	
	if (request->totalScenarios > 0)
	{
		reply->nS = request->totalScenarios;
		reply->p = new double[reply->nS];
		readBodur(request, reply);
	}
	else
	{
		std::cout << "Cannot read Bodur instance without no. of scenarios specified" << endl;
		throw(-1);
	}
	initialised = true;
}

// Deletes the singleton.
#if 0
void ProblemData::cleanup() {
	//TODO: delete problem stuff! currently leaking like a sieve... or is it all dealt with by delete env?
	if (ProblemData::initialised) {
		delete ProblemData::singleton;
		ProblemData::initialised = false;
	}
		
}
#endif

// Reads a CAP problem from Bodur et al. (2014)
int ProblemDataBodur::readBodur(SMIP_fileRequest* request, SMIP_fileReply* reply) {
	
	std::string filename(request->filename);
	std::string filepath(request->filepath);
	
	std::string full_filename = filepath + filename;
	
	int i, j, k;
	
	IloNum btmp;
	IloInt I, J;
	IloInt K = nS;
	IloInt KMax;
	
	// Read in the constraints and objective function.
	
	p = (1.0 / nS);
	
	ifstream file(full_filename.c_str());
	
	if (!file) {
		cerr << "ERROR: could not open file '" << filename
		<< "' for reading" << endl;
		throw(-1);
	}
	
	file >> c >> d >> T >> hArray >> A;
	
	file.close();
	
	I = c.getSize();
	J = T.getSize()-I;
	
	KMax = hArray.getSize();
	
	IloNumArray summands(env,K);
	for(k=0; k<K; k++) {
		summands[k] = IloSum(hArray[k]);
	}
	b.add(IloMax(summands));
	
	for(j=0; j<J; j++){
		IloNumArray rowJ(env,I*J);
		for(i=0; i<I; i++){
			rowJ[i*J+j]=1;
		}
		W.add(rowJ);
	}
	for(i=0; i<I; i++){
		IloNumArray rowI(env,I*J);
		for(j=0; j<J; j++){
			rowI[i*J+j]=-1;
		}
		W.add(rowI);
	}
	
	// Create the variables.
	
	// Fill in the ProblemData object.
	n1 = c.getSize();
	n2 = d.getSize();
	reply->n1 = n1;
	reply->n2 = n2;
	for (int i = 0; i < nS; i++) {
		reply->p[i] = 1.0 / nS;
		
	}
	reply->fileLoaded = true;
	
	return 1;
}

// Standard get functions.
#if 0
IloNumVarArray* ProblemData::get_x() {return this->x;}
IloNumArray* ProblemData::get_c() {return this->c;}
TwoDMatrix* ProblemData::get_A() {return this->A;}
IloNumArray* ProblemData::get_b() {return this->b;}
	
IloNumVarArray* ProblemData::get_y(int scenIndex) {return this->y->getParameter(scenIndex);}
IloNumArray* ProblemData::get_d(int scenIndex) {return this->d->getParameter(scenIndex);}
TwoDMatrix* ProblemData::get_T(int scenIndex) {return this->T->getParameter(scenIndex);}
TwoDMatrix* ProblemData::get_W(int scenIndex) {return this->W->getParameter(scenIndex);}
IloNumArray* ProblemData::get_h(int scenIndex) {return this->h->getParameter(scenIndex);}

IloNum* ProblemData::get_p(int scenIndex) {return this->p->getParameter(scenIndex);}

int ProblemData::get_n1() {return n1;}
int ProblemData::get_n2() {return n2;}
int ProblemData::get_nS() {return nS;};
#endif
