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
	strcpy(filepath,par->filename.c_str());
	initialiseFileName();
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
	AlgorithmZ = par->AlgorithmZ; 
	disableHeuristic = par->disableHeuristic;
	ftype = par->filetype;

	if (nVerticesUsed < 0) { useVertexHistory = false; }
	if (nVerticesUsed == 0) { useVertexHistory = true; }
	if (nVerticesUsed > 0) { useVertexHistory = true; }
}

void PSCGModel::setupSolvers(){
	nNodeSPs = scenariosToThisModel.size();
	recordKeeping = new double[nNodeSPs][4];
	scaleVec_ = new double[nNodeSPs];
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
	    omega_current.push_back(new double[n1]);
	    omega_saved.push_back(new double[n1]);
	    
	    for (int i = 0; i < n1; i++) {
		scaling_matrix[tS][i] = penC ;
		omega_current[tS][i] = 0.0; //omega will be initialised from the node.
		omega_saved[tS][i] = 0.0; //omega will be initialised from the node.
	    }
	    subproblemSolvers[tS]->setQuadraticTerm(scaling_matrix[tS]);
	}
}

int PSCGModel::initialIteration(){
if(mpiRank==0){cerr << "Begin initialIteration()" << endl;}
	// This is part of the initialisation - initial Lagrangian subproblem computations
	double LagrLB_Local = 0.0;
	modelStatus_[SP_STATUS]=SP_ITER_LIM;
	for (int tS = 0; tS < nNodeSPs; tS++) {
 	    //updateTimer.start();
	    spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_current[tS]);
            if(spSolverStatuses_[tS]==PSCG_PRIMAL_INF || spSolverStatuses_[tS]==PSCG_DUAL_INF){
	    	    LagrLB_Local = ALPS_DBL_MAX;
cerr << "Subproblem " << tS << " infeasible on proc " << mpiRank << endl;
	    	    continue;

	    }
	    subproblemSolvers[tS]->updateSolnInfo();
	    subproblemSolvers[tS]->setXToVertex();
	    subproblemSolvers[tS]->setYToVertex();
	    subproblemSolvers[tS]->updateOptSoln();
	    recordKeeping[tS][0]=subproblemSolvers[tS]->getLagrBd();
	    recordKeeping[tS][1]=recordKeeping[tS][0];

	    //updateTimer.stop();
	    //updateTimer.addTime(updateTimeThisStep);
	    LagrLB_Local += pr[tS]*( subproblemSolvers[tS]->getLagrBd() );//obj;

	    if(useVertexHistory){
	      if(subproblemSolvers[tS]->getNumVertices() < nVerticesUsed || nVerticesUsed==0){subproblemSolvers[tS]->addVertex();}
	      else if (nVerticesUsed > 0){// && subproblemSolvers[tS]->getNumVertices() >= nVerticesUsed) {
		//subproblemSolvers[tS]->removeBackVertex();
		subproblemSolvers[tS]->replaceOldestVertex();
	        //subproblemSolvers[tS]->fixWeightToZero( subproblemSolvers[tS]->getNumVertices() - nVerticesUsed-1);
	      }
	    }
#if 0
	    if (useVertexHistory) {
		subproblemSolvers[tS]->updateVertexHistory();
		if (nVerticesUsed > 0 && subproblemSolvers[tS]->getNumVertices() > nVerticesUsed) {
		    subproblemSolvers[tS]->removeBackVertex();
		    //subproblemSolvers[tS]->fixWeightToZero( subproblemSolvers[tS]->getNumVertices() - nVerticesUsed-1);
		}
	    }
#endif
#if 0
	    for (int i = 0; i < n1; i++) {
	    	omega_current[tS][i] += scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
#endif
	}
	
//if(mpiRank==0){cout << "After z: " << endl;}
//printZ();
	//weightedAverage(x_current, pr, z_local, z_current, nNodeSPs, n1, mpiRank);

	#ifdef USING_MPI
	if (mpiSize > 1) {
			// Send decision variables to other processes.
//if(mpiRank==0)cerr << "InitialIteration: before MPI_Allreduce "  << endl;
//cerr << "#";
		MPI_Allreduce(&LagrLB_Local, &LagrLB, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//if(mpiRank==0)cerr << "InitialIteration: after MPI_Allreduce "  << endl;

			//commsTimer.stop();
			//commsTimer.addTime(commsTimeThisStep);
	}
	#endif
	if(mpiSize == 1){
	    LagrLB = LagrLB_Local;
	}
	
	currentLagrLB = LagrLB;
	//recordKeeping[0]=LagrLB;
	assert(ALPS_DBL_MAX > ALPS_INFINITY);
	if(LagrLB > ALPS_INFINITY){ modelStatus_[SP_STATUS]=SP_INFEAS;} //for parallel
	else{
	    //Update of z
	    modelStatus_[SP_STATUS]=SP_OPT;
	    updateZ();
	}
#if 0
	for (int tS = 0; tS < nNodeSPs; tS++) {
	    for (int i = 0; i < n1; i++) {
	    	omega_current[tS][i] += scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
	}
#endif
if(mpiRank==0){cerr << "End initialIteration()" << endl;}
	return modelStatus_[SP_STATUS];
}
bool PSCGModel::solveRecourseProblemGivenFixedZ(){	
if(mpiRank==0){cout << "Begin solveRecourseProblemGivenFixedZ()" << endl;}
//printOriginalVarBds();
//printCurrentVarBds();
	double localObjVal=0.0;
	//objVal = 0.0;
	bool isFeas=true;
	memcpy(z_rounded,z_current,n1*sizeof(double));
	for(int ii=0; ii<numIntVars_; ii++) {z_rounded[intVar_[ii]] = round(z_current[intVar_[ii]]);}
//printZRounded();
	for (int tS = 0; tS < nNodeSPs; tS++) {

		//****************** Compute Next Vertex **********************
		// Find next vertex and Lagrangian lower bound
		
		//Solve Lagrangian MIP
		subproblemSolvers[tS]->solveLagrangianWithXFixedToZ(z_rounded, NULL, currentVarLB_, currentVarUB_, colType_);
		//if(tS==0) subproblemSolvers[tS]->printLinCoeffs();
		//if(tS==0) subproblemSolvers[tS]->printC();
		//if(tS==0) subproblemSolvers[tS]->printLagrSoln();
		//subproblemSolvers[tS]->updateSolnInfo();
                if(subproblemSolvers[tS]->getSolverStatus()==PSCG_PRIMAL_INF || subproblemSolvers[tS]->getSolverStatus()==PSCG_DUAL_INF)
		{
	    	    //LagrLB = ALPS_DBL_MAX;	    
		    objVal = ALPS_DBL_MAX;
		    if(mpiRank==0){cout << "Recourse problem given fixed z is infeasible due to subproblem " << tS << endl;}
		    if(mpiRank==0){cout << "End solveRecourseProblemGivenFixedZ() with value: " << objVal << endl;}
		    isFeas=false;
		    localObjVal = ALPS_DBL_MAX;
		    //return isFeas;
		    break;
		}
		//updateTimer.addTime(updateTimeThisStep);

		//Acumulate the values for the Lagrangian subproblems
		//cout << "  " << subproblemSolvers[tS]->getLagrBd();
		localObjVal += pr[tS]*subproblemSolvers[tS]->getLagrBd();//FWSPoptvalk;}
	}
//cout << "LagrLB_Local******************************: " << objVal << endl;
	#ifdef USING_MPI
	if (mpiSize > 1) {

		localReduceBuffer[0]=localObjVal;
//if(mpiRank==0) cerr << "solveRecourseProblemGivenFixedZ(): before MPI_Allreduce()" << endl;
//cerr << "#";
		MPI_Allreduce(localReduceBuffer, reduceBuffer, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//if(mpiRank==0) cerr << "solveRecourseProblemGivenFixedZ(): end MPI_Allreduce()" << endl;
		if(mpiRank==0){cout << "Recourse problem given fixed z has value: " << reduceBuffer[0] << " versus current LB " << currentLagrLB << endl;}
		objVal=reduceBuffer[0];
	}
	#endif

	if (mpiSize == 1) {
		objVal=localObjVal;
		if(mpiRank==0){cout << "Recourse problem given fixed z has value: " << objVal << " versus current LB " << currentLagrLB << endl;}
	}
if(mpiRank==0){cout << "End solveRecourseProblemGivenFixedZ() with value: " << objVal << endl;}
	return isFeas;
}

int PSCGModel::performColGenStep(){	
//if(true){cerr << "Proc: " << mpiRank << " Begin performColGenStep(): "  << endl;}
	LagrLB_Local = 0.0;
	ALVal_Local = 0.0;
	localDiscrepNorm = 0.0;
	reduceBuffer[0]=0.0;
	reduceBuffer[1]=0.0;
	reduceBuffer[2]=0.0;
	double gapVal = 0.0, sqrDiscrNorm = 0.0;
	
	for (int tS = 0; tS < nNodeSPs; tS++) {
		//****************** Compute Next Vertex **********************
	// Find next vertex and Lagrangian lower bound
//if(tS==0) subproblemSolvers[tS]->printC();
		recordKeeping[tS][2]=subproblemSolvers[tS]->getALVal();
		recordKeeping[tS][3]=subproblemSolvers[tS]->getSqrNormDiscr();
		ALVal_Local += pr[tS]*subproblemSolvers[tS]->getALVal();
		sqrDiscrNorm = subproblemSolvers[tS]->getSqrNormDiscr();
		localDiscrepNorm += pr[tS]*subproblemSolvers[tS]->getSqrNormDiscr();
	
		for (int i = 0; i < n1; i++) {
		    omega_tilde[tS][i] = omega_current[tS][i] + scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
		}

		//Solve Lagrangian MIP
		spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_tilde[tS]);
		if(spSolverStatuses_[tS]==PSCG_PRIMAL_INF || spSolverStatuses_[tS]==PSCG_DUAL_INF){
	    	    //LagrLB = ALPS_DBL_MAX;	    
	    	    LagrLB_Local = ALPS_DBL_MAX;
	    	    ALVal_Local = ALPS_DBL_MAX;
cerr << "Subproblem " << tS << " infeasible on proc " << mpiRank << endl;
	    	    break;
		}
		subproblemSolvers[tS]->updateSolnInfo();
		//subproblemSolvers[tS]->compareLagrBds(omega_tilde[tS]);
		gapVal=subproblemSolvers[tS]->updateGapVal(omega_tilde[tS]);
		//updateTimer.addTime(updateTimeThisStep);
		scaleVec_[tS] = (gapVal <= 0.5*sqrDiscrNorm) ? 2.0:1.0;
		//Acumulate the values for the Lagrangian subproblems
	        recordKeeping[tS][0]=subproblemSolvers[tS]->getLagrBd();
	    	//cout << recordKeeping[tS][0] << " versus " << subproblemSolvers[tS]->evaluateVertexSolution(omega_tilde[tS]) << endl;
		LagrLB_Local += pr[tS]*subproblemSolvers[tS]->getLagrBd();//FWSPoptvalk;}
	 	if(useVertexHistory){
		    if(subproblemSolvers[tS]->getNumVertices() < nVerticesUsed || nVerticesUsed==0){subproblemSolvers[tS]->addVertex();}
		    else if(nVerticesUsed > 0){// && subproblemSolvers[tS]->getNumVertices() >= nVerticesUsed) {
			subproblemSolvers[tS]->replaceOldestVertex();
		    }
		}
//cerr << "Proc: " << mpiRank << ": Subproblem " << tS << " solve finished." << endl;
	}
	#ifdef USING_MPI
//if(mpiRank==0) cerr << "PerformColGenStep(): MPI before MPI_Barrier" << endl;;
//cerr << "*";
//MPI_Barrier(MPI_COMM_WORLD);
	if (mpiSize > 1) {
		//commsTimer.start();
		// Send decision variables to other processes.

		localReduceBuffer[0]=LagrLB_Local;
		localReduceBuffer[1]=ALVal_Local;
		localReduceBuffer[2]=localDiscrepNorm;

//if(mpiRank==0) cerr << "PerformColGenStep(): MPI before MPI_Allreduce" << endl;;
//cerr << "#";
		MPI_Allreduce(localReduceBuffer, reduceBuffer, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		//MPI::Comm::Allreduce(localReduceBuffer, reduceBuffer, 3, MPI::DOUBLE, MPI::SUM);
//if(mpiRank==0) cerr << "PerformColGenStep(): MPI after MPI_Allreduce" << endl;;

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
	for(int tS = 0; tS<nNodeSPs; tS++){
	    //cout << recordKeeping[tS][1] << " <=??? " << recordKeeping[tS][2] << "  (discr: " << recordKeeping[tS][3] << ")" << endl;
	    //cout << recordKeeping[tS][1] << " versus " << subproblemSolvers[tS]->evaluateVertexOptSolution(omega_current[tS]) << endl;;
	    //cout << " and " << recordKeeping[tS][2] << " versus " << subproblemSolvers[tS]->evaluateSolution(omega_current[tS])+0.5*recordKeeping[tS][3] << endl;
	    //if(recordKeeping[tS][1] > recordKeeping[tS][2]){
	    if(recordKeeping[tS][1] >  1e-2 + recordKeeping[tS][2] + 0.5*recordKeeping[tS][3]){
	        cout << recordKeeping[tS][1] << " <=??? " << recordKeeping[tS][2] << "  (discr: " << recordKeeping[tS][3] << ")" << endl;
	        cout << recordKeeping[tS][1] << " <=??? " << recordKeeping[tS][2]+0.5*recordKeeping[tS][3] << "  (discr: " << recordKeeping[tS][3] << ")" << endl;
		cout << "Current solution value: " << subproblemSolvers[tS]->evaluateSolution(omega_current[tS]) << endl;;
		subproblemSolvers[tS]->evaluateVertexHistory(omega_current[tS]);
	        subproblemSolvers[tS]->printWeights();	
	    }
	    //assert(recordKeeping[tS][1] <=  1e-2 + recordKeeping[tS][2] + 0.5*recordKeeping[tS][3]);
	}
	//updateModelStatusSP();
//if(mpiRank==0) cerr << endl;
	if(LagrLB > ALPS_INFINITY) modelStatus_[SP_STATUS]=SP_INFEAS; //for parallel
	else{modelStatus_[SP_STATUS]=SP_OPT;}
//if(mpiRank==0){cerr << "End performColGenStep(): "  << endl;}
	return modelStatus_[SP_STATUS];
}


AlpsTreeNode* PSCGModel::createRoot(){
    PSCGTreeNode *root = new PSCGTreeNode;
    PSCGNodeDesc *desc = new PSCGNodeDesc(this);
    desc->setZLBs(origVarLB_);
    desc->setZUBs(origVarUB_);
    root->setDesc(desc);
    root->setExplicit(1);
    return root;
}
//******************Command line paramater handler functions**********************



void PSCGModel::displayParameters(){
  if(mpiRank==0){
	#ifdef USING_MPI
	std::cout << "USINGMPI:TRUE" << endl;
	#else
	std::cout << "USINGMPI:FALSE" << endl;
	#endif

	std::cout << "ALGORITHM:PSCG with branch and bound" << std::endl;
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
	std::cout << "Using Algorithm variant to branch using z: " << AlgorithmZ << std::endl;

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

