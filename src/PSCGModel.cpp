/*
Finds the optimal solution to the Lagrangian dual of a two-stage stochastic
integer program using a Frank-Wolfe-based Method of Multipliers approach.
*/

#include "PSCGModel.h"
//#include "CPLEXsolverSCG.h"
#include "ProblemDataBodur.h"
#include "PSCGModelScen.h"
#include "TssModel.h"
#include "PSCGTreeNode.h"



using namespace std;

void PSCGModel::initialiseParameters(){
	//par->readParameters(argc, argv);
	setMPIParams(par->mpiRank,par->mpiSize);

	std::string str_filename(par->filename);
	std::string str_outputFileID(par->outputFilename);
	nS = par->noScenarios;
	penC = par->penalty;
	penMult = par->penaltyMult;
	if (penC <= 0) {
	    penC = DEFAULT_PENALTY;

	    if (mpiHead) {
		std::cout << "Penalty parameter chosen by heuristic: " << penC << std::endl;
	    }

	    if (penMult > 0) {
		penC = penC * penMult;
		
		if (mpiHead) {
			std::cout << "Penalty parameter altered by multiplier: " << penC << std::endl;
		}
	    }
    	}

	maxStep = par->maxStep;
	maxSeconds = par->maxSeconds;
	fixInnerStep = par->fixInnerStep;
	nVerticesUsed = par->UseVertexHistory;
	nThreads = par->threads;
	if (nThreads < 0) { nThreads = DEFAULT_THREADS; }
	verbose = par->verbose;
	debug = par->debug;
	linRelaxFirst = par->linRelaxFirst;
	linRelaxSecond = par->linRelaxSecond;
	scaling = par->scaling;
	dataPathOverride = par->dataPathOverride;
	LBcalc = par->LBcalc;
	AlgorithmC = par->AlgorithmC; //obsolete, get rid of.
	disableHeuristic = par->disableHeuristic;
	ftype = par->filetype;

	if (nVerticesUsed < 0) { useVertexHistory = false; }
	if (nVerticesUsed == 0) { useVertexHistory = true; }
	if (nVerticesUsed > 0) { useVertexHistory = true; }
}

void PSCGModel::setupSolvers(){
	nNodeSPs = scenariosToThisModel.size();
	for (int tS = 0; tS < nNodeSPs; tS++) {
	    //CPLEXsolverSCG &currentSPSolver = subproblemSolvers[tS];
#if 1
	    switch( ftype ){
		    case 1: //CAP Problem
#if 1
			subproblemSolvers.push_back( new PSCGModelScen_Bodur() );
		     	dynamic_cast<PSCGModelScen_Bodur*>(subproblemSolvers[tS])->initialiseBodur(par,pdBodur,scenariosToThisModel[tS]);
#endif
		      	break;
		    case 2: //SIPLIB problems
			subproblemSolvers.push_back( new PSCGModelScen_SMPS() );
		      	dynamic_cast<PSCGModelScen_SMPS*>(subproblemSolvers[tS])->initialiseSMPS(par,smpsModel,scenariosToThisModel[tS]); 
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
	    omega_current.push_back(NULL);
	    
	    for (int i = 0; i < n1; i++) {
		scaling_matrix[tS][i] = penC ;
		//omega_current[tS][i] = 0.0;
	    }
	    subproblemSolvers[tS]->setQuadraticTerm(scaling_matrix[tS]);
	}
}

int PSCGModel::initialIteration(){
	// This is part of the initialisation - initial Lagrangian subproblem computations
	double LagrLB_Local = 0.0;
	modelStatus_[SP_STATUS]=SP_ITER_LIM;
	for (int tS = 0; tS < nNodeSPs; tS++) {
 	    //updateTimer.start();
	    spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_current[tS]);
    	    if(spSolverStatuses_[tS]==PSCG_OPTIMAL || spSolverStatuses_[tS]==PSCG_ITER_LIM){	
		subproblemSolvers[tS]->updateSolnInfo();
	    	subproblemSolvers[tS]->setXToVertex();
	    	subproblemSolvers[tS]->setYToVertex();
	    }
	    else{
	    	    LagrLB = ALPS_DBL_MAX;	    
		    currentLagrLB = LagrLB;
	    	    modelStatus_[SP_STATUS]=SP_INFEAS;
	    	    return SP_INFEAS;
	    }

	    //updateTimer.stop();
	    //updateTimer.addTime(updateTimeThisStep);

	    LagrLB_Local += pr[tS]*( subproblemSolvers[tS]->getLagrBd() );//obj;

	    if (useVertexHistory) {
		subproblemSolvers[tS]->updateVertexHistory();
		if (nVerticesUsed > 0 && subproblemSolvers[tS]->getNumVertices() > nVerticesUsed) {
		    subproblemSolvers[tS]->removeBackVertex();
		    //subproblemSolvers[tS]->fixWeightToZero( subproblemSolvers[tS]->getNumVertices() - nVerticesUsed-1);
		}
	    }
#if 0
	    for (int i = 0; i < n1; i++) {
	    	omega_current[tS][i] += scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
#endif
	}
	
	//Update of z
	updateZ();
	//weightedAverage(x_current, pr, z_local, z_current, nNodeSPs, n1, mpiRank);

	#ifdef USING_MPI
	if (mpiSize > 1) {
			// Send decision variables to other processes.
		MPI_Allreduce(&LagrLB_Local, &LagrLB, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			//commsTimer.stop();
			//commsTimer.addTime(commsTimeThisStep);
	}
	#endif
	if(mpiSize == 1){
	    LagrLB = LagrLB_Local;
	}
	
	currentLagrLB = LagrLB;
#if 0
	for (int tS = 0; tS < nNodeSPs; tS++) {
	    for (int i = 0; i < n1; i++) {
	    	omega_current[tS][i] += scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
	}
#endif
	return modelStatus_[SP_STATUS];
}
bool PSCGModel::solveRecourseProblemGivenFixedZ(){	
	LagrLB_Local = 0.0;
	bool isFeas=true;
	for (int tS = 0; tS < nNodeSPs; tS++) {

		//****************** Compute Next Vertex **********************
		// Find next vertex and Lagrangian lower bound
		
		//Solve Lagrangian MIP
		subproblemSolvers[tS]->solveLagrangianWithXFixedToZ(z_current, NULL, origVarLB_, origVarUB_, colType_);
		//subproblemSolvers[tS]->updateSolnInfo();
                if(!(subproblemSolvers[tS]->getSolverStatus()==PSCG_OPTIMAL || subproblemSolvers[tS]->getSolverStatus()==PSCG_ITER_LIM))
		{
	    	    LagrLB = ALPS_DBL_MAX;	    
		    isFeas=false;
		    return isFeas;
		}
		//updateTimer.addTime(updateTimeThisStep);

		//Acumulate the values for the Lagrangian subproblems
		LagrLB_Local += pr[tS]*subproblemSolvers[tS]->getLagrBd();//FWSPoptvalk;}
//cout << " (" << pr[tS] << "," << subproblemSolvers[tS]->getLagrBd() << ")";
#if 0
	 	if(useVertexHistory){
		    subproblemSolvers[tS]->updateVertexHistory();
		    if (nVerticesUsed > 0 && subproblemSolvers[tS]->getNumVertices() > nVerticesUsed) {
			subproblemSolvers[tS]->removeBackVertex();
		        //subproblemSolvers[tS]->fixWeightToZero( subproblemSolvers[tS]->getNumVertices() - nVerticesUsed-1);
		    }
		}
#endif
	}
cout << endl;
cout << "LagrLB_Local******************************: " << LagrLB_Local << endl;
	#ifdef USING_MPI
	if (mpiSize > 1) {

		localReduceBuffer[0]=LagrLB_Local;

		MPI_Allreduce(localReduceBuffer, reduceBuffer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		LagrLB = reduceBuffer[0];
	}
	#endif

	if (mpiSize == 1) {
		LagrLB = LagrLB_Local;
		//ALVal = ALVal_Local;
		//discrepNorm = localDiscrepNorm;
		//cout << "Testing: " << ALVal+ discrepNorm  << " >= " << currentLagrLB << endl;
	}
        if(isFeas){
	if(incumbentVal_ > LagrLB){ 
	    incumbentVal_=LagrLB;
	    memcpy( z_incumbent_, z_current, n1*sizeof(double) );
        }
	}
	return isFeas;
}

int PSCGModel::performColGenStep(){	
	LagrLB_Local = 0.0;
	ALVal_Local = 0.0;
	localDiscrepNorm = 0.0;
	
	for (int tS = 0; tS < nNodeSPs; tS++) {

		//****************** Compute Next Vertex **********************
		// Find next vertex and Lagrangian lower bound
		
		subproblemSolvers[tS]->updateALValues(omega_current[tS],z_current,scaling_matrix[tS]);
		ALVal_Local += pr[tS]*subproblemSolvers[tS]->getALVal();
		localDiscrepNorm += pr[tS]*subproblemSolvers[tS]->getSqrNormDiscr();
	
		for (int i = 0; i < n1; i++) {
		    omega_tilde[tS][i] = omega_current[tS][i] + scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
		}

		//Solve Lagrangian MIP
		spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_tilde[tS]);
    	        if(spSolverStatuses_[tS]==PSCG_OPTIMAL || spSolverStatuses_[tS]==PSCG_ITER_LIM){	
		    subproblemSolvers[tS]->updateSolnInfo();
		}
		else{
	    	    LagrLB = ALPS_DBL_MAX;	    
	    	    modelStatus_[SP_STATUS]=SP_INFEAS;
	    	    return SP_INFEAS;
		}
		//updateTimer.addTime(updateTimeThisStep);

		//Acumulate the values for the Lagrangian subproblems
		LagrLB_Local += pr[tS]*subproblemSolvers[tS]->getLagrBd();//FWSPoptvalk;}

	 	if(useVertexHistory){
		    subproblemSolvers[tS]->updateVertexHistory();
		    if (nVerticesUsed > 0 && subproblemSolvers[tS]->getNumVertices() > nVerticesUsed) {
			subproblemSolvers[tS]->removeBackVertex();
		        //subproblemSolvers[tS]->fixWeightToZero( subproblemSolvers[tS]->getNumVertices() - nVerticesUsed-1);
		    }
		}
	}
	
	#ifdef USING_MPI
	if (mpiSize > 1) {
		//commsTimer.start();
		// Send decision variables to other processes.

		localReduceBuffer[0]=LagrLB_Local;
		localReduceBuffer[1]=ALVal_Local;
		localReduceBuffer[2]=localDiscrepNorm;

		MPI_Allreduce(localReduceBuffer, reduceBuffer, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		LagrLB = reduceBuffer[0];
		ALVal = reduceBuffer[1];
		discrepNorm = reduceBuffer[2];
		//if(mpiRank==0){ cout << "Testing: " << ALVal+ discrepNorm  << " >= " << currentLagrLB << endl;}

		//commsTimer.stop();
		//commsTimer.addTime(commsTimeThisStep);
		}
	#endif

	if (mpiSize == 1) {
		LagrLB = LagrLB_Local;
		ALVal = ALVal_Local;
		discrepNorm = localDiscrepNorm;
		//cout << "Testing: " << ALVal+ discrepNorm  << " >= " << currentLagrLB << endl;
	}
	updateModelStatusSP();
	return modelStatus_[SP_STATUS];
}
//******************FWPH functions**********************


//NOT CORRECTLY IMPLEMENTED
void computeFeasibleStartingPoint(int tS, double* x, double* yFeasible) {
	SMIP_qu_secondStage question(tS, x);
	SMIP_ans_secondStage answer(yFeasible);
	//CPLEXsolverSCG subproblemSolvers;
	//subproblemSolvers.solveForSecondStage(&question, &answer);
}

AlpsTreeNode* PSCGModel::createRoot(){
    PSCGTreeNode *root = new PSCGTreeNode;
    PSCGNodeDesc *desc = new PSCGNodeDesc(this);
    root->setDesc(desc);
    root->setExplicit(1);
    return root;
}
//******************Command line paramater handler functions**********************



void PSCGModel::displayParameters(){
	#ifdef USING_MPI
	std::cout << "USINGMPI:TRUE" << endl;
	#else
	std::cout << "USINGMPI:FALSE" << endl;
	#endif

	std::cout << "ALGORITHM:FWPH" << std::endl;
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

	std::cout << "Maximum outer step: " << maxStep << std::endl;
	std::cout << "Maximum seconds spent on main updates: " << maxSeconds << std::endl;

	if (fixInnerStep > 0) {
		std::cout << "Number of inner loop iterations: " << fixInnerStep << std::endl;
	}
	else {
		std::cout << "Number of inner loop iterations: " << "1" << std::endl;
	}

	if (penC > 0) {
		std::cout << "Manually chosen starting penalty parameter: " << penC << std::endl;
	}
	else {
		std::cout << "Manually chosen starting penalty parameter: " << "Chosen by heuristic" << std::endl;
		if (penMult > 0) {
			std::cout << "Multiplier applied to heuristic parameter: " << penMult << std::endl;
		}
	}
	std::cout << "Verbose output: " << verbose << std::endl;
	std::cout << "Disabled heuristic: " << disableHeuristic << std::endl;
	std::cout << "Number of threads for CPLEX: " << nThreads << std::endl;
	std::cout << "First stage linear relaxation: " << linRelaxFirst << std::endl;
	std::cout << "Second stage linear relaxation: " << linRelaxSecond << std::endl;
	std::cout << "Scaling: " << scaling << std::endl;
	std::cout << "Debug output: " << debug << std::endl;
	std::cout << "Finding feasible point for first vertex: Yes" << std::endl;
	std::cout << "Using Algorithm C variant: " << AlgorithmC << std::endl;

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
}

