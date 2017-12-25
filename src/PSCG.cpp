/*
Finds the optimal solution to the Lagrangian dual of a two-stage stochastic
integer program using a Frank-Wolfe-based Method of Multipliers approach.
*/

#include "PSCG.h"
//#include "tclap/CmdLine.h"
//#include "StructureDefs.h"



using namespace std;

#if 0
PSCG::PSCG(PSCGParams *p):par(p),env(),nNodeSPs(0),referenceLagrLB(-COIN_DBL_MAX),cutoffLagrLB(COIN_DBL_MAX),currentLagrLB(-COIN_DBL_MAX),centreLagrLB(-COIN_DBL_MAX),trialLagrLB(-COIN_DBL_MAX),
LagrLB_Local(0.0),ALVal_Local(COIN_DBL_MAX),ALVal(COIN_DBL_MAX),objVal(COIN_DBL_MAX),
	incumbentVal(COIN_DBL_MAX),localDiscrepNorm(1e9),discrepNorm(1e9),mpiRank(0),mpiSize(1),
	totalNoGSSteps(0),infeasIndex_(-1),maxNoSteps(1e6),maxNoGSSteps(1),maxNoConseqNullSteps(1e6),noGSIts(1),
	omegaUpdated_(false),SSCParam(0.1),innerSSCParam(0.5),phase(0),tCritParam(1e-10){

   	//******************Read Command Line Parameters**********************
	//Params par;
	initialiseParameters();

	//******************Display Command Line Parameters**********************
	if (mpiRank==0) {
		displayParameters();
	}

	//******************Reading Problem Data**********************

	initialiseModel();
	


	//******************Assign Scenarios**********************
	assignSubproblems();
	setupSolvers();
	//sprintf(zOptFile,"zOpt-%s-noP%dalg%d",probname,mpiSize,algorithm);
	
}
#endif

PSCG::PSCG(DecTssModel &model):smpsModel(model),env(),nNodeSPs(0),referenceLagrLB(-COIN_DBL_MAX),cutoffLagrLB(COIN_DBL_MAX),currentLagrLB(-COIN_DBL_MAX),centreLagrLB(-COIN_DBL_MAX),trialLagrLB(-COIN_DBL_MAX),
LagrLB_Local(0.0),ALVal_Local(COIN_DBL_MAX),ALVal(COIN_DBL_MAX),objVal(COIN_DBL_MAX),
	incumbentVal(COIN_DBL_MAX),localDiscrepNorm(1e9),discrepNorm(1e9),mpiRank(0),mpiSize(1),
	totalNoGSSteps(0),infeasIndex_(-1),maxNoSteps(1000000),maxNoGSSteps(1),maxNoInnerSteps(MAX_NO_INNERSTEPS),maxNoConseqNullSteps(1e6),noGSIts(1),
	baselineRho(1.0),rho(1.0),nThreads(1),nVerticesUsed(100),
	nS(-1),ftype(2),omegaUpdated_(false),SSCParam(0.1),innerSSCParam(0.5),phase(0),tCritVal(1e10),tCritParam(1e-10){

	//smpsModel.readSmps(par->filename.c_str());
	
	n1 = smpsModel.getNumCols(0);
	n2 = smpsModel.getNumCols(1);
	if(nS<=0) nS = smpsModel.getNumScenarios();
	if (mpiRank==0) {
		std::cout << std::endl;
		std::cout << "Problem data: " << std::endl;
		std::cout << "Number of first stage variables: " << n1 << std::endl;
		std::cout << "Number of second stage variables: " << n2 << std::endl;
		std::cout << "Number of scenarios: " << nS << std::endl;
    	}

	//******************Reading Problem Data**********************
	initialiseModel();

	//******************Assign Scenarios**********************
	assignSubproblems();
	setupSolvers();
}
#ifdef USING_MPI
PSCG::PSCG(DecTssModel &model, MPI_Comm comm):smpsModel(model),env(),comm_(comm),nNodeSPs(0),referenceLagrLB(-COIN_DBL_MAX),cutoffLagrLB(COIN_DBL_MAX),currentLagrLB(-COIN_DBL_MAX),centreLagrLB(-COIN_DBL_MAX),trialLagrLB(-COIN_DBL_MAX),
LagrLB_Local(0.0),ALVal_Local(COIN_DBL_MAX),ALVal(COIN_DBL_MAX),objVal(COIN_DBL_MAX),
	incumbentVal(COIN_DBL_MAX),localDiscrepNorm(1e9),discrepNorm(1e9),
	totalNoGSSteps(0),infeasIndex_(-1),maxNoSteps(1000000),maxNoGSSteps(1),maxNoInnerSteps(MAX_NO_INNERSTEPS),maxNoConseqNullSteps(1e6),noGSIts(1),
	baselineRho(1.0),rho(1.0),nThreads(1),nVerticesUsed(100),
	nS(-1),ftype(2),omegaUpdated_(false),SSCParam(0.1),innerSSCParam(0.5),phase(0),tCritVal(1e10),tCritParam(1e-10){

	//smpsModel.readSmps(par->filename.c_str());
	
	MPI_Comm_rank(comm_,&mpiRank);
	MPI_Comm_size(comm_,&mpiSize);
	n1 = smpsModel.getNumCols(0);
	n2 = smpsModel.getNumCols(1);
	if(nS<=0) nS = smpsModel.getNumScenarios();
	if (mpiRank==0) {
		std::cout << std::endl;
		std::cout << "Problem data: " << std::endl;
		std::cout << "Number of first stage variables: " << n1 << std::endl;
		std::cout << "Number of second stage variables: " << n2 << std::endl;
		std::cout << "Number of scenarios: " << nS << std::endl;
    	}

	//******************Reading Problem Data**********************
	initialiseModel();

	//******************Assign Scenarios**********************
	assignSubproblems();
	setupSolvers();
}
#endif

PSCG::~PSCG(){
	//Clean up structures
	for (int tS = 0; tS < nNodeSPs; tS++) {
		delete [] scaling_matrix[tS];
		delete [] omega_tilde[tS];
		delete [] omega_current[tS];
		delete [] omega_centre[tS];
		delete [] omega_sp[tS];
		delete [] omega_saved[tS];
		delete subproblemSolvers[tS];
	    #ifdef KEEP_LOG
	        if(logFiles.size() > tS) logFiles[tS]->close();
	    #endif
	}

	delete [] z_current;
	delete [] z_average;
	delete [] z_dispersions;
	delete [] z_old;
#if 0
	delete [] z_intdisp;
	delete [] z_vertex_average;
	delete [] z_average1;
	delete [] z_average2;
#endif
	delete [] z_saved;
	delete [] z_local;
	delete [] z_incumbent_;
	delete [] z_rounded;
	delete [] penSumLocal;
	delete [] penSum;
	delete [] origVarLB_;
	delete [] origVarUB_;
	delete [] currentVarLB_;
	delete [] currentVarUB_;
	delete [] recordKeeping;
	delete [] weights_;
	delete [] integrDiscr_;
	env.end();
}

#if 0
void PSCG::initialiseParameters(){
	//par->readParameters(argc, argv);
	setMPIParams(par->mpiRank,par->mpiSize);

	std::string str_filename(par->filename);
	std::string str_outputFileID(par->outputFilename);
	strcpy(filepath,par->filename.c_str());
	initialiseFileName();
	nS = par->noScenarios;
	baselineRho = 1.0;//par->penalty;
	rho = baselineRho;
	//penMult = par->penaltyMult;
	if (rho <= 0) {
	    rho = DEFAULT_PENALTY;

	    if (mpiRank==0) {
		std::cout << "Penalty parameter chosen by heuristic: " << rho << std::endl;
	    }
#if 0
	    if (penMult > 0) {
		rho = rho * penMult;
		
		if (mpiRank==0) {
			std::cout << "Penalty parameter altered by multiplier: " << rho << std::endl;
		}
	    }
#endif
    	}

	maxNoSteps = par->maxStep;
	maxSeconds = par->maxSeconds;
	if(par->fixInnerStep > 1) noGSIts = par->fixInnerStep;
	nVerticesUsed = par->UseVertexHistory;
	if (par->threads < 0) { nThreads = DEFAULT_NO_THREADS; }
	else{nThreads = par->threads;}
	verbose = par->verbose;
	debug = par->debug;
	linRelaxFirst = par->linRelaxFirst;
	linRelaxSecond = par->linRelaxSecond;
	scaling = par->scaling;
	dataPathOverride = par->dataPathOverride;
	//LBcalc = par->LBcalc;
	//AlgorithmZ = par->AlgorithmZ; 
	//disableHeuristic = par->disableHeuristic;
	ftype = par->filetype;

}
#endif

void PSCG::assignSubproblems(){
	if(mpiRank==0) { cout << "Assigning " << nS << " subproblems with " << mpiSize << " processors. " << endl;}
	for (int tS = 0; tS < nS; tS++) {
		//scenarioAssign[tS] = tS % mpiSize;
		if( (tS % mpiSize)==mpiRank ){
		    scenariosToThisModel.push_back(tS);
		    spSolverStatuses_.push_back(SP_UNKNOWN);
		}
	}
	if(mpiRank==0) { cout << "Done assigning " << nS << " subproblems. " << endl;}
}

void PSCG::setupSolvers(){
	nNodeSPs = scenariosToThisModel.size();
cout << "Begin setting up " << nNodeSPs << " solvers at process " << mpiRank << endl;
	recordKeeping = new double[nNodeSPs][4];
	integrDiscr_ = new double[nNodeSPs];
	weights_ = new double[nNodeSPs];
	//scaleVec_ = new double[nNodeSPs];
	char logFileName[256];
	for (int tS = 0; tS < nNodeSPs; tS++) {
	  #ifdef KEEP_LOG
	    shared_ptr<ofstream> logfile(new ofstream);
	    sprintf(logFileName,"Logs/logProc%dScen%d.txt",mpiRank,tS);
	    logfile->open(logFileName);
	    logFiles.push_back(logfile);
	    //logFiles.emplace_back(ofstream(logFileName));
	  #endif
	    //CPLEXsolverSCG &currentSPSolver = subproblemSolvers[tS];
#if 1
	    switch( ftype ){
		    case 1: //CAP Problem
#if 0
			subproblemSolvers.push_back( new PSCGScen_Bodur() );
		     	dynamic_cast<PSCGScen_Bodur*>(subproblemSolvers[tS])->initialiseBodur(par,pdBodur,scenariosToThisModel[tS]);
#endif
		      	break;
		    case 2: //SIPLIB problems
			subproblemSolvers.push_back( new PSCGScen_SMPS(env) );
		      	dynamic_cast<PSCGScen_SMPS*>(subproblemSolvers[tS])->initialiseSMPS(smpsModel,scenariosToThisModel[tS]); 
			//subproblemSolvers[tS]->setNThreads(par->threads);
			subproblemSolvers[tS]->setNThreads(nThreads);
		      	break;
		    default:
			throw(-1);
			break;
	    }
#endif
	    subproblemSolvers[tS]->finishInitialisation(); 
	    
	    x_current.push_back(subproblemSolvers[tS]->getX());
	    y_current.push_back(subproblemSolvers[tS]->getY());
	    pr.push_back(subproblemSolvers[tS]->getProbabilities());
	    scaling_matrix.push_back(new double[n1]);
	    omega_tilde.push_back(new double[n1]);
	    omega_current.push_back(new double[n1]);
	    omega_centre.push_back(new double[n1]);
	    omega_sp.push_back(new double[n1]);
	    omega_saved.push_back(new double[n1]);
	    
	    for (int i = 0; i < n1; i++) {
		//scaling_matrix[tS][i] = rho ;
		scaling_matrix[tS][i] = max( fabs(subproblemSolvers[tS]->getC()[i]),1.0);
		//scaling_matrix[tS][i] = 1.0;
		//omega_saved[tS][i] = 0.0; //omega will be initialised from the node.
		omega_centre[tS][i] = 0.0; 
	    }
	    //subproblemSolvers[tS]->setQuadraticTerm(rho,scaling_matrix[tS]);
	}
	zeroOmega();
cout << "End setting up solvers at process " << mpiRank << endl;
}

void PSCG::initialiseModel(){

	z_current = new double[n1];
	z_old = new double[n1];
	z_saved = new double[n1];
	z_local = new double[n1];
	z_incumbent_ = new double[n1];
	z_rounded = new double[n1];
	z_average = new double[n1];
	z_dispersions = new double[n1];
	penSumLocal = new double[n1];
	penSum = new double[n1];
	origVarLB_ = new double[n1];
	origVarUB_ = new double[n1];
	currentVarLB_ = new double[n1];
	currentVarUB_ = new double[n1];

	if(ftype==2){
#if 0
	    //memcpy(colType_,smpsModel.getCtypeCore(0),n1*sizeof(char));
	    //memcpy(colType_+n1,smpsModel.getCtypeCore(1),n2*sizeof(char));
	    for(int i=0; i<n1; i++){
		colType_[i] = smpsModel.getCtypeCore(0)[i];
		if(colType_[i]=='I' || colType_[i]=='B'){
		    intVar_[numIntVars_++]=i;
		}
	    }
	    for(int i=n1; i<n1+n2; i++){
	    //for(int i=0; i<n1; i++){
		colType_[i] = smpsModel.getCtypeCore(1)[i-n1];
		if(colType_[i]=='I' || colType_[i]=='B'){
		    intVar_[numIntVars_++]=i;
		}
	    }
#endif
	    smpsModel.copyCoreColLower(origVarLB_,0);
	    smpsModel.copyCoreColUpper(origVarUB_,0);
	    memcpy(currentVarLB_,origVarLB_,n1*sizeof(double));
	    memcpy(currentVarUB_,origVarUB_,n1*sizeof(double));
	    cout << "******************" << endl;
	    printCurrentVarBds();
	    cout << "******************" << endl;
	}
}

#if 0
void PSCG::installSubproblem(double lb, const double *zLBs, const double *zUBs){
    for(int ii=0; ii<n1; ii++){
	if(zLBs[ii] > zUBs[ii]) cerr << "Bounds inconsistent: " << zLBs[ii] << " > " << zUBs[ii] << endl;
	assert(zLBs[ii] <= zUBs[ii]);
    }
    setLBsAllSPs(zLBs);
    setUBsAllSPs(zUBs);
    currentLagrLB=-COIN_DBL_MAX;
    centreLagrLB=-COIN_DBL_MAX;
    referenceLagrLB=lb;
    printOriginalVarBds();
    printCurrentVarBds();
}
void PSCG::installSubproblem(double lb, vector<double*> &omega, const double *zLBs, const double *zUBs){
    installSubproblem(lb, zLBs, zUBs);
    readOmegaIntoModel(omega);
}

void PSCG::installSubproblem(double lb, const vector<int> &indices, const vector<double> &zLBs, const vector<double> &zUBs){
    for(int ii=0; ii<indices.size(); ii++){
	if(zLBs[ii] > zUBs[ii]) cerr << "Bounds inconsistent: " << zLBs[ii] << " > " << zUBs[ii] << endl;
	assert(zLBs[ii] <= zUBs[ii]);
    }
    restoreOriginalVarBounds();
    setBdsAllSPs(indices, zLBs, zUBs);
    currentLagrLB=-COIN_DBL_MAX;
    centreLagrLB=-COIN_DBL_MAX;
    referenceLagrLB=lb;
    printOriginalVarBds();
    printCurrentVarBds();
}
void PSCG::installSubproblem(double lb, vector<double*> &omega, const vector<int> &indices, const vector<double> &zLBs, const vector<double> &zUBs){
    installSubproblem(lb, indices, zLBs, zUBs);
    readOmegaIntoModel(omega);
}
#endif
void PSCG::restoreOriginalVarBounds(){
    clearBranchingBounds();
    for(int tS=0; tS<nNodeSPs; tS++){
	subproblemSolvers[tS]->restoreBounds();//unfixX(origVarLB_,origVarUB_);
    }
    memcpy(currentVarLB_,origVarLB_,n1*sizeof(double));
    memcpy(currentVarUB_,origVarUB_,n1*sizeof(double));
}
void PSCG::addBranchVarBd(int br_rank, int br_SP, int br_index, double br_lb, double br_ub){
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
//cout <<br_SP <<"*"<<br_index<<"*"<<br_lb<<"*"<<br_ub<<endl;
	    setBdForSP(br_SP, br_index, br_lb, br_ub);
        }
	//else this processor and/or subproblem does nothing
    }
    else{
	cout << "No branch found, this node should be solved to optimality." << endl;
    }

}
void PSCG::setBdForSP(const int sp, const int ind, const double lb, const double ub){
    subproblemSolvers[sp]->setBound(ind, lb, ub);
    //currentVarLB_[ind]=lb;
    //currentVarUB_[ind]=ub;
}
void PSCG::setBdAllSPs(const int ind, const double lb, const double ub){
    for(int tS=0; tS<nNodeSPs; tS++){
	setBdForSP(tS,ind,lb,ub);	
    }
}

int PSCG::initialIteration(){
if(mpiRank==0){cerr << "Begin initialIteration()" << endl;}
	// This is part of the initialisation - initial Lagrangian subproblem computations
        modelStatus_[Z_STATUS]=Z_UNKNOWN;
	modelStatus_[SP_STATUS]=SP_ITER_LIM;
	cumNoInnerSolves=0;
try{
    	clearSPVertexHistory();
}
catch(std::exception &e){
    std::cout << "exception caught while updating vertex history in initialIteration: " << e.what() << endl;
    exit(1);
}
    	noSeriousSteps=0;
    	noConseqNullSteps=0;
	innerSSCVal = 0.5;
	double LagrLB_Local = 0.0;
        //setPenalty(1.0);
	rho=baselineRho;
	if(mpiRank==0){cout << "Starting rho is: " << rho << endl;}
    	#ifdef KEEP_LOG
	    for (int tS = 0; tS < nNodeSPs; tS++) {*(logFiles[tS]) << "Initial iteration: " << mpiRank << " scen " << tS << endl;}
    	#endif
	for (int tS = 0; tS < nNodeSPs; tS++) {
 	    //updateTimer.start();
    	    //subproblemSolvers[tS]->setInitialSolution(NULL);

	   subproblemSolvers[tS]->resetDispersionsToZero();
    	   try{
#if 0
	    spSolverStatuses_[tS] = subproblemSolvers[tS]->initialLPSolve(omega_centre[tS]);
            if(spSolverStatuses_[tS]==PSCG_PRIMAL_INF || spSolverStatuses_[tS]==PSCG_DUAL_INF){
	    	    LagrLB_Local = COIN_DBL_MAX;
		    cerr << "initialIteration(): Subproblem " << tS << " infeasible on proc with initialLPSolve" << mpiRank; 
		    cerr << " with CPLEX status: " << subproblemSolvers[tS]->getCPLEXErrorStatus() << endl;
	    	    continue;

	    }
#endif

	    spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_centre[tS],true);
            if(spSolverStatuses_[tS]==PSCG_PRIMAL_INF || spSolverStatuses_[tS]==PSCG_DUAL_INF){
	    	    LagrLB_Local = COIN_DBL_MAX;
		    cerr << "initialIteration(): Subproblem " << tS << " infeasible on proc " << mpiRank; 
		    cerr << " with CPLEX status: " << subproblemSolvers[tS]->getCPLEXErrorStatus() << endl;
	    	    continue;
	    }

	   }
	   catch(CoinError &e){
		cerr << "Exception thrown during MIP solve phase." << endl;
	   }

try{
	    subproblemSolvers[tS]->updateSolnInfo();
	    subproblemSolvers[tS]->addVertex();
	    subproblemSolvers[tS]->setXToVertex();
	    subproblemSolvers[tS]->setYToVertex();
	    //subproblemSolvers[tS]->updateOptSoln();
	    //recordKeeping[tS][0]=-ALPS_DBL_MAX;
	    recordKeeping[tS][0]=subproblemSolvers[tS]->getLagrBd();
	    recordKeeping[tS][1]=recordKeeping[tS][0];

	    //updateTimer.stop();
	    //updateTimer.addTime(updateTimeThisStep);
	    LagrLB_Local += pr[tS]*( subproblemSolvers[tS]->getLagrBd() );//obj;

	    updateVertexHistory(tS);
}
catch(std::exception &e){
    std::cout << "exception caught while updating vertex history in initialIteration: " << e.what() << endl;
    exit(1);
}
#if 0
	    for (int i = 0; i < n1; i++) {
	    	omega_current[tS][i] += rho*scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
#endif
	}//for tS

	
//if(mpiRank==0){cout << "After z: " << endl;}
//printZ();
	//weightedAverage(x_current, pr, z_local, z_current, nNodeSPs, n1, mpiRank);

	#ifdef USING_MPI
	if (mpiSize > 1) {
		MPI_Allreduce(&LagrLB_Local, &trialLagrLB, 1, MPI_DOUBLE, MPI_SUM, comm_);
	}
	#endif
	if(mpiSize == 1){
	    trialLagrLB = LagrLB_Local;
	}
	

	//currentLagrLB =-ALPS_DBL_MAX;
	currentLagrLB = trialLagrLB;
	centreLagrLB = trialLagrLB;

	//recordKeeping[0]=trialLagrLB;
	if(trialLagrLB > COIN_DBL_MAX/10.0){ 
	    modelStatus_[SP_STATUS]=SP_INFEAS;
	    modelStatus_[Z_STATUS]=Z_INFEAS;
            shouldContinue=false;
	    if(mpiRank==0){cout << "Terminating due to subproblem infeasibility..." << endl;}
	} //for parallel
	else{
	    //Update of z
	    //currentLagrLB=-ALPS_DBL_MAX;
	    modelStatus_[SP_STATUS]=SP_OPT;
	    updateZ();

	    ALVal_Local = 0.0;
	    localDiscrepNorm = 0.0;
	    reduceBuffer[0]=0.0;
	    reduceBuffer[1]=0.0;
	    double ALVal_tS,sqrDiscrNorm_tS;
	    for (int tS = 0; tS < nNodeSPs; tS++) {
	        ALVal_tS = subproblemSolvers[tS]->updateALValues(omega_centre[tS],z_current,rho,scaling_matrix[tS]);
		ALVal_Local += pr[tS]*ALVal_tS;
		sqrDiscrNorm_tS = subproblemSolvers[tS]->getSqrNormDiscr();
		localDiscrepNorm += pr[tS]*sqrDiscrNorm_tS;
		localDiscrepNorm += 0.5*pr[tS]*subproblemSolvers[tS]->getSqrNormDiscr();
		for (int i = 0; i < n1; i++) {
		    omega_tilde[tS][i] = omega_centre[tS][i] + rho*scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
		}
	    }
	    #ifdef USING_MPI
	    if (mpiSize > 1) {
		localReduceBuffer[0]=ALVal_Local;
		localReduceBuffer[1]=localDiscrepNorm;
		MPI_Allreduce(localReduceBuffer, reduceBuffer, 2, MPI_DOUBLE, MPI_SUM, comm_);
		ALVal = reduceBuffer[0];
		discrepNorm = reduceBuffer[1];
	    }
	    #endif
	    if(mpiSize == 1){
		ALVal = ALVal_Local;
		discrepNorm = localDiscrepNorm;
	    }
            shouldContinue=true;

#if 0
	    if(discrepNorm >= 1e-20){
		if(mpiRank==0) cout << "Initial penalty value: " << min(1e10,fabs(centreLagrLB)/discrepNorm) << endl;
		baselineRho = min(1e10,fabs(centreLagrLB)/discrepNorm);
		setPenalty( baselineRho );
	    }
	    else{
        	shouldContinue=false;
		if(mpiRank==0){cout << "Terminating due to optimality at initial iteration..." << endl;}
	    }
    	    if(currentLagrLB >= cutoffLagrLB){
		modelStatus_[Z_STATUS]=Z_BOUNDED;
        	shouldContinue=false;
		if(mpiRank==0){cout << "Terminating due to exceeding cutoff..." << endl;}
      		//printStatus();
    	    }
#else
	    //rho=baselineRho;
#endif
	}
	if(mpiRank==0){cerr << "End initialIteration()" << endl;}
	return modelStatus_[SP_STATUS];
}

bool PSCG::solveRecourseProblemGivenFixedZ(){	
if(mpiRank==0){cout << "Begin solveRecourseProblemGivenFixedZ()" << endl;}
	double localObjVal=0.0;
	bool isFeas=true;
	for (int tS = 0; tS < nNodeSPs; tS++) {

		//****************** Compute Next Vertex **********************
		// Find next vertex and Lagrangian lower bound
		
		//Solve Lagrangian MIP
		int solverStatus=subproblemSolvers[tS]->solveLagrangianWithXFixedToZ(z_rounded, NULL, currentVarLB_, currentVarUB_);

                if(solverStatus==PSCG_PRIMAL_INF || solverStatus==PSCG_DUAL_INF)
		{
		    objVal = COIN_DBL_MAX;
		    isFeas=false;
		    localObjVal = COIN_DBL_MAX;
		    break;
		}

		//Acumulate the values for the Lagrangian subproblems
		localObjVal += pr[tS]*subproblemSolvers[tS]->getLagrBd();//FWSPoptvalk;}
	}
	#ifdef USING_MPI
	if (mpiSize > 1) {

		localReduceBuffer[0]=localObjVal;
		MPI_Allreduce(&localObjVal, &objVal, 1, MPI_DOUBLE, MPI_SUM, comm_);
	}
	#endif

	if (mpiSize == 1) {
		objVal=localObjVal;
	}
if(mpiRank==0){cout << "End solveRecourseProblemGivenFixedZ() with value: " << objVal << endl;}
	return isFeas;
}

int PSCG::performColGenStep(){	
//{cout << "Proc: " << mpiRank << " Begin performColGenStep(): "  << endl;}
//if(true){cerr << "Proc: " << mpiRank << " Begin performColGenStep(): "  << endl;}
	LagrLB_Local = 0.0;
	ALVal_Local = 0.0;
	localDiscrepNorm = 0.0;
	reduceBuffer[0]=0.0;
	reduceBuffer[1]=0.0;
	reduceBuffer[2]=0.0;

	double gapVal = 0.0, sqrDiscrNorm_tS = 0.0;

	double ALVal_tS;
	double LagrLB_tS;
	double lhsCritVal;
	
	for (int tS = 0; tS < nNodeSPs; tS++) {
		//ALVal_tS = subproblemSolvers[tS]->getALVal();
		//scaleVec_[tS] = 1.0;
		lhsCritVal=0.0;
	        ALVal_tS = subproblemSolvers[tS]->updateALValues(omega_centre[tS],z_current,rho,scaling_matrix[tS]);
		ALVal_Local += pr[tS]*ALVal_tS;
		sqrDiscrNorm_tS = subproblemSolvers[tS]->getSqrNormDiscr();
		localDiscrepNorm += pr[tS]*sqrDiscrNorm_tS;

		//****************** Compute Next Vertex **********************
	// Find next vertex and Lagrangian lower bound
//if(tS==0) subproblemSolvers[tS]->printC();
		recordKeeping[tS][2]=ALVal_tS;
		recordKeeping[tS][3]=sqrDiscrNorm_tS;

		//Solve Lagrangian MIP


		
		dynamic_cast<PSCGScen_SMPS*>(subproblemSolvers[tS])->getOSI()->resolve();

		//LagrLB_tS = subproblemSolvers[tS]->getLagrBd();	
		spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_tilde[tS]);
		if(spSolverStatuses_[tS]==PSCG_PRIMAL_INF || spSolverStatuses_[tS]==PSCG_DUAL_INF){
	    	    LagrLB_Local = COIN_DBL_MAX;
	    	    ALVal_Local = COIN_DBL_MAX;
cerr << "performColGenStep(): Subproblem " << tS << " infeasible on proc " << mpiRank << endl;
	    	    break;
		}

		LagrLB_tS = subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]);
#ifdef KEEP_LOG
    	    *(logFiles[tS]) << "old LagrLB_ts: " << LagrLB_tS << "\tnew LagrLB_ts" << subproblemSolvers[tS]->getLagrBd() << "\tALVal_ts: " << ALVal_tS << "\tsqrDiscrNorm: " << sqrDiscrNorm_tS << endl;
    	    *(logFiles[tS]) << "Testing whether vertex is redundant: " << -(LagrLB_tS - subproblemSolvers[tS]->getLagrBd()) << endl;;
#endif
#if 0
	if(LagrLB_tS <= subproblemSolvers[tS]->getLagrBd()+1e-10){
#ifdef KEEP_LOG
    	    *(logFiles[tS]) << "  FLAGGING: Redundant vertex found" << endl;;
	    *(logFiles[tS]) << " Is vertex really redundant? " << subproblemSolvers[tS]->checkWhetherVertexIsRedundant() << endl;
#endif
	}
	else{
	    updateVertexHistory(tS);
	}
#endif
	updateVertexHistory(tS);
#if 1
		LagrLB_tS = subproblemSolvers[tS]->getLagrBd();	
		lhsCritVal = ALVal_tS+0.5*sqrDiscrNorm_tS;
#ifdef KEEP_LOG
		*(logFiles[tS]) << "Printing penalities: ";
#endif
		for(int ii=0; ii<n1; ii++){
		    lhsCritVal += rho*scaling_matrix[tS][ii]*z_current[ii]*(x_current[tS][ii]-z_current[ii]);
#ifdef KEEP_LOG
		    *(logFiles[tS]) << " " << rho*scaling_matrix[tS][ii];
#endif
		}
#ifdef KEEP_LOG
		*(logFiles[tS]) << endl;
		*(logFiles[tS]) << setprecision(10) << lhsCritVal << " should be >= " << setprecision(10) << LagrLB_tS << 
			", a difference of: " << lhsCritVal-LagrLB_tS << ", compared with sqrDiscrNorm: " << sqrDiscrNorm_tS << endl;
#endif
		if(lhsCritVal  < LagrLB_tS){
		  #ifdef KEEP_LOG
		    if(lhsCritVal + 1e-6 < LagrLB_tS){*(logFiles[tS]) << "performColGenStep():  iteration: " << currentIter_  << ": lhsCritVal condition not met: " 
			<< setprecision(10) << lhsCritVal << " should be >= " << setprecision(10) << LagrLB_tS << endl;
		        //subproblemSolvers[tS]->setMIPPrintLevel(1, 5, false);
			//spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_tilde[tS]);
			//dynamic_cast<PSCGScen_SMPS*>(subproblemSolvers[tS])->getOSI()->resolve();
		        //subproblemSolvers[tS]->setMIPPrintLevel(0, 0, false);
		    }
		  #endif
		    //subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]);
		    subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]);
		    //subproblemSolvers[tS]->refresh(omega_current[tS], z_current, scaling_matrix[tS]);
		    //updateZ();
		    //updateOmegaTilde();
    		    if(lhsCritVal + 1e-6 < LagrLB_tS){
		        //subproblemSolvers[tS]->printLinCoeffs();
			//cout << "Repairing Lagrangian subproblem solution: opt value was: " << setprecision(10) << LagrLB_tS << " and is now " << subproblemSolvers[tS]->getLagrBd() << endl;
			//subproblemSolvers[tS]->evaluateVertexHistory(omega_tilde[tS]);
		    }
			//subproblemSolvers[tS]->evaluateVertexHistory(omega_tilde[tS]);
		    //if(mpiRank==10) subproblemSolvers[tS]->printLagrSoln();

#if 0
    		    subproblemSolvers[tS]->setInitialSolution(NULL);
		    cout << "performColGenStep(): trying again with scenario " << tS << " of node " << mpiRank << endl; 
		    spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_tilde[tS]);
#endif
		}

		//else{
		//}

#endif
		//if(mpiRank==10) cout << "************************************************************" << endl;



		//subproblemSolvers[tS]->compareLagrBds(omega_tilde[tS]);
		gapVal=subproblemSolvers[tS]->updateGapVal(omega_tilde[tS]);
		#ifdef KEEP_LOG
		*(logFiles[tS]) << "Condition for scaling penalty: " << gapVal << "  versus  " << sqrDiscrNorm_tS << endl;
		subproblemSolvers[tS]->evaluateVertexHistory(omega_tilde[tS], logFiles[tS]);
		#endif
		//updateTimer.addTime(updateTimeThisStep);
#if 0
		if(sqrDiscrNorm_tS >= SSC_DEN_TOL*SSC_DEN_TOL){ 
		   //scaleVec_[tS] = 0.618+ 1.0/(exp(lhsCritVal - LagrLB_tS));
		   //scaleVec_[tS] = 0.5+ 1.5/((1.0/sqrt(sqrDiscrNorm_tS))*exp(lhsCritVal - LagrLB_tS));
		   //scaleVec_[tS] = pow(scaleVec_[tS], 1.0/(0.2*currentIter_+1.0));
		   //scaleVec_[tS] = (gapVal <= 0.5*sqrDiscrNorm_tS) ? 2.0:1.0;//2.0:0.5;
		   //if(gapVal <= 0.1*sqrDiscrNorm_tS){scaleVec_[tS] *=2.0;}//2.0:0.5;
		   //if(gapVal >= 10.0*sqrDiscrNorm_tS){scaleVec_[tS] *=0.5;}//2.0:0.5;
		   //scaleVec_[tS] = (gapVal <= 0.5*sqrDiscrNorm_tS) ? 1.618:0.618;//2.0:0.5;
		   //scaleVec_[tS] = min(max(0.5,sqrt(pow( gapVal/sqrDiscrNorm_tS , -1.0))),1.5) + 0.5; 
		}
#endif
#if 0
		if(gapVal <= 0.25*sqrDiscrNorm_tS){
		    scaleVec_[tS] = 2.0;
		}
		else if(gapVal <= 0.5*sqrDiscrNorm_tS){
		    scaleVec_[tS] = 1.414;
		}
		else if(gapVal > 4*sqrDiscrNorm_tS){
		    scaleVec_[tS] = 0.5;
		}
		else if(gapVal > 2*sqrDiscrNorm_tS){
		    scaleVec_[tS] = 0.7071;
		}
		else{
		    scaleVec_[tS] = 1.0;
		}
#endif
		//scaleVec_[tS] = (gapVal <= 0.5*sqrDiscrNorm_tS) ? 2.0:1.0;//2.0:0.5;
		//scaleVec_[tS] = pow(1.001,0.5*sqrDiscrNorm_tS - gapVal);//2.0:0.5;
//cerr << "ScaleVec: " << scaleVec_[tS] << endl;

		//Acumulate the values for the Lagrangian subproblems
	        recordKeeping[tS][0]=subproblemSolvers[tS]->getLagrBd();
	    	//cout << recordKeeping[tS][0] << " versus " << subproblemSolvers[tS]->evaluateVertexSolution(omega_tilde[tS]) << endl;
		LagrLB_Local += pr[tS]*subproblemSolvers[tS]->getLagrBd();//FWSPoptvalk;}
//cerr << "Proc: " << mpiRank << ": Subproblem " << tS << " solve finished." << endl;
	}
	#ifdef USING_MPI
//if(mpiRank==0) cerr << "PerformColGenStep(): MPI before MPI_Barrier" << endl;;
//cerr << "*";
//MPI_Barrier(comm_);
	if (mpiSize > 1) {
		localReduceBuffer[0]=LagrLB_Local;
		localReduceBuffer[1]=ALVal_Local;
		localReduceBuffer[2]=localDiscrepNorm;
		MPI_Allreduce(localReduceBuffer, reduceBuffer, 3, MPI_DOUBLE, MPI_SUM, comm_);

		trialLagrLB = reduceBuffer[0];
		ALVal = reduceBuffer[1];
		discrepNorm = reduceBuffer[2];
		}
	#endif

	if (mpiSize == 1) {
		trialLagrLB = LagrLB_Local;
		ALVal = ALVal_Local;
		discrepNorm = localDiscrepNorm;
	}
	for(int tS = 0; tS<nNodeSPs; tS++){
	    //cout << recordKeeping[tS][1] << " <=??? " << recordKeeping[tS][2] << "  (discr: " << recordKeeping[tS][3] << ")" << endl;
	    //cout << recordKeeping[tS][1] << " versus " << subproblemSolvers[tS]->evaluateVertexOptSolution(omega_current[tS]) << endl;;
	    //cout << " and " << recordKeeping[tS][2] << " versus " << subproblemSolvers[tS]->evaluateSolution(omega_current[tS])+0.5*recordKeeping[tS][3] << endl;
	    //if(recordKeeping[tS][1] > recordKeeping[tS][2]){
	    if(recordKeeping[tS][1] >  1e-2 + recordKeeping[tS][2] + 0.5*recordKeeping[tS][3]){
#if 0
	        cout << recordKeeping[tS][1] << " <=??? " << recordKeeping[tS][2] << "  (discr: " << recordKeeping[tS][3] << ")" << endl;
	        cout << recordKeeping[tS][1] << " <=??? " << recordKeeping[tS][2]+0.5*recordKeeping[tS][3] << "  (discr: " << recordKeeping[tS][3] << ")" << endl;
		cout << "Current solution value: " << subproblemSolvers[tS]->evaluateSolution(omega_centre[tS]) << endl;;
#endif

		//subproblemSolvers[tS]->evaluateVertexHistory(omega_current[tS]);
	        //subproblemSolvers[tS]->printWeights();	

	    }
	    //assert(recordKeeping[tS][1] <=  1e-2 + recordKeeping[tS][2] + 0.5*recordKeeping[tS][3]);
	}
	//updateModelStatusSP();
//if(mpiRank==0) cerr << endl;
	if(trialLagrLB > COIN_DBL_MAX/10.0) modelStatus_[SP_STATUS]=SP_INFEAS; //for parallel
	else{modelStatus_[SP_STATUS]=SP_OPT;}
//if(mpiRank==0){cerr << "End performColGenStep(): "  << endl;}
	return modelStatus_[SP_STATUS];
}


//******************Command line paramater handler functions**********************



void PSCG::displayParameters(){
  if(mpiRank==0){
	std::cout << "Parameters for PSCG:" << endl;
	#ifdef USING_MPI
	std::cout << "USINGMPI:TRUE" << endl;
	#else
	std::cout << "USINGMPI:FALSE" << endl;
	#endif

	std::cout << "ALGORITHM: PSCG" << std::endl;
	std::cout << "Number of processors: " << mpiSize << std::endl;
	std::cout << "Data file path: " << filepath << std::endl;
	//std::cout << "Problem filename: " << str_filename << std::endl;

  	switch (ftype) {
		case 0:
			std::cout << "Problem file type: NONE" << std::endl;
			break;
		case 1:
			std::cout << "Problem file type: CAP" << std::endl;
			break;
		case 2:
			std::cout << "Problem file type: SMPS" << std::endl;
			break;
		default:
			std::cout << "Problem file type: SOMETHING HAS GONE HORRIBLY WRONG" << std::endl;
			break;
	}

	if (nS > 0) {
		std::cout << "Manually chosen number of scenarios: " << nS << std::endl;
	}
	else {
		std::cout << "Manually chosen number of scenarios: " << "Maximum available in file" << std::endl;
	}

	std::cout << "Maximum outer step: " << maxNoSteps << std::endl;
	std::cout << "Maximum seconds spent on main updates: " << maxSeconds << std::endl;

	if (noGSIts > 0) {
		std::cout << "Number of inner loop iterations: " << noGSIts << std::endl;
	}
	else {
		std::cout << "Number of inner loop iterations: " << "1" << std::endl;
	}

	if (rho > 0) {
		std::cout << "Manually chosen starting penalty parameter: " << rho << std::endl;
	}
	else {
		std::cout << "Manually chosen starting penalty parameter: " << "Chosen by heuristic" << std::endl;
#if 0
		if (penMult > 0) {
			std::cout << "Multiplier applied to heuristic parameter: " << penMult << std::endl;
		}
#endif
	}
	std::cout << "Verbose output: " << verbose << std::endl;
	//std::cout << "Disabled heuristic: " << disableHeuristic << std::endl;
	std::cout << "Number of threads for CPLEX: " << nThreads << std::endl;
	std::cout << "First stage linear relaxation: " << linRelaxFirst << std::endl;
	std::cout << "Second stage linear relaxation: " << linRelaxSecond << std::endl;
	std::cout << "Scaling: " << scaling << std::endl;
	std::cout << "Debug output: " << debug << std::endl;
	std::cout << "Finding feasible point for first vertex: Yes" << std::endl;
	//std::cout << "Using Algorithm variant to branch using z: " << AlgorithmZ << std::endl;
	//printAlgorithm();

	if (nVerticesUsed < 0) {
		std::cout << "Number of vertexes used: 1"<< std::endl;
	} 
	else if (nVerticesUsed == 0) {
		std::cout << "Number of vertexes used: All"<< std::endl;
	} 
	else {
		std::cout << "Number of vertexes used: " << nVerticesUsed << std::endl;
	}

	std::cout << std::endl << std::endl;
    	std::cout << std::endl << std::endl;
  }//if mpiRank==0
}

