/*
Finds the optimal solution to the Lagrangian dual of a two-stage stochastic
integer program using a Frank-Wolfe-based Method of Multipliers approach.
*/

#include "PSCG.h"
#include <utility>
//#include "CPLEXsolverSCG.h"
#include "ProblemDataBodur.h"
#include "PSCGScen.h"
#include "TssModel.h"



using namespace std;

void PSCG::initialiseParameters(){
	//par->readParameters(argc, argv);
	setMPIParams(par->mpiRank,par->mpiSize);

	std::string str_filename(par->filename);
	std::string str_outputFileID(par->outputFilename);
	strcpy(filepath,par->filename.c_str());
	initialiseFileName();
	nS = par->noScenarios;
	rho = par->penalty;
	algorithm = par->Algorithm;
	baselineRho=rho;
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
	//LBcalc = par->LBcalc;
	//AlgorithmZ = par->AlgorithmZ; 
	//disableHeuristic = par->disableHeuristic;
	ftype = par->filetype;

}

void PSCG::setupSolvers(){
	nNodeSPs = scenariosToThisModel.size();
	recordKeeping = new double[nNodeSPs][4];
	integrDiscr_ = new double[nNodeSPs];
	weights_ = new double[nNodeSPs];
	//scaleVec_ = new double[nNodeSPs];
	char logFileName[128];
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
#if 1
			subproblemSolvers.push_back( new PSCGScen_Bodur() );
		     	dynamic_cast<PSCGScen_Bodur*>(subproblemSolvers[tS])->initialiseBodur(par,pdBodur,scenariosToThisModel[tS]);
#endif
		      	break;
		    case 2: //SIPLIB problems
			subproblemSolvers.push_back( new PSCGScen_SMPS() );
		      	dynamic_cast<PSCGScen_SMPS*>(subproblemSolvers[tS])->initialiseSMPS(par,smpsModel,scenariosToThisModel[tS]); 
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
		scaling_matrix[tS][i] = rho ;
		//omega_saved[tS][i] = 0.0; //omega will be initialised from the node.
		omega_centre[tS][i] = 0.0; 
	    }
	    subproblemSolvers[tS]->setQuadraticTerm(scaling_matrix[tS]);
	}
	zeroOmega();
}

int PSCG::initialIteration(){
if(mpiRank==0){cerr << "Begin initialIteration()" << endl;}
	// This is part of the initialisation - initial Lagrangian subproblem computations
	double LagrLB_Local = 0.0;
        setPenalty(1.0);
	modelStatus_[SP_STATUS]=SP_ITER_LIM;
	for (int tS = 0; tS < nNodeSPs; tS++) {
 	    //updateTimer.start();
    	    //subproblemSolvers[tS]->setInitialSolution(NULL);

	   subproblemSolvers[tS]->resetDispersionsToZero();
    	   try{
	    spSolverStatuses_[tS] = subproblemSolvers[tS]->initialLPSolve(omega_centre[tS]);
            if(spSolverStatuses_[tS]==PSCG_PRIMAL_INF || spSolverStatuses_[tS]==PSCG_DUAL_INF){
	    	    LagrLB_Local = COIN_DBL_MAX;
cerr << "initialIteration(): Subproblem " << tS << " infeasible on proc " << mpiRank; 
cerr << " with CPLEX status: " << subproblemSolvers[tS]->getCPLEXErrorStatus() << endl;
	    	    continue;

	    }

	    spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_centre[tS]);
            if(spSolverStatuses_[tS]==PSCG_PRIMAL_INF || spSolverStatuses_[tS]==PSCG_DUAL_INF){
	    	    LagrLB_Local = COIN_DBL_MAX;
//cerr << "initialIteration(): Subproblem " << tS << " infeasible on proc " << mpiRank << endl;
cerr << "initialIteration(): Subproblem " << tS << " infeasible on proc " << mpiRank; 
cerr << " with CPLEX status: " << subproblemSolvers[tS]->getCPLEXErrorStatus() << endl;
//subproblemSolvers[tS]->printColBds();
	    	    continue;

	    }

	   }
	   catch(CoinError &e){
		cerr << "Exception thrown during MIP solve phase." << endl;
	   }

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
#if 0
	    for (int i = 0; i < n1; i++) {
	    	omega_current[tS][i] += scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
	    }
#endif
	}//for tS
	
//if(mpiRank==0){cout << "After z: " << endl;}
//printZ();
	//weightedAverage(x_current, pr, z_local, z_current, nNodeSPs, n1, mpiRank);

	#ifdef USING_MPI
	if (mpiSize > 1) {
		MPI_Allreduce(&LagrLB_Local, &trialLagrLB, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}
	#endif
	if(mpiSize == 1){
	    trialLagrLB = LagrLB_Local;
	}
	

	//currentLagrLB =-ALPS_DBL_MAX;
	currentLagrLB = trialLagrLB;
	centreLagrLB = trialLagrLB;

	//recordKeeping[0]=trialLagrLB;
	if(trialLagrLB > COIN_DBL_MAX/10.0){ modelStatus_[SP_STATUS]=SP_INFEAS;} //for parallel
	else{
	    //Update of z
	    //currentLagrLB=-ALPS_DBL_MAX;
	    modelStatus_[SP_STATUS]=SP_OPT;
	    updateZ();
#if 1
	    double sqrDiscrNorm = 0.0;
	    double sqrDiscrNormLocal = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++) {
	         subproblemSolvers[tS]->updateALValues(omega_centre[tS],z_current,scaling_matrix[tS]);
		 sqrDiscrNormLocal += 0.5*pr[tS]*subproblemSolvers[tS]->getSqrNormDiscr();
	    }
	    #ifdef USING_MPI
	    if (mpiSize > 1) {
		MPI_Allreduce(&sqrDiscrNormLocal, &sqrDiscrNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    }
	    #endif
	    if(mpiSize == 1){
	        sqrDiscrNorm = sqrDiscrNormLocal;
	    }
	    if(sqrDiscrNorm >= 1e-20){
		if(mpiRank==0) cout << "Initial penalty value: " << min(1e10,0.01*fabs(centreLagrLB)/sqrDiscrNorm) << endl;
		baselineRho = min(1e10,fabs(centreLagrLB)/sqrDiscrNorm);
		setPenalty( baselineRho );
	    }
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
		int solverStatus=subproblemSolvers[tS]->solveLagrangianWithXFixedToZ(z_rounded, NULL, currentVarLB_, currentVarUB_, colType_);

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
		MPI_Allreduce(&localObjVal, &objVal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}
	#endif

	if (mpiSize == 1) {
		objVal=localObjVal;
	}
if(mpiRank==0){cout << "End solveRecourseProblemGivenFixedZ() with value: " << objVal << endl;}
	return isFeas;
}

int PSCG::performColGenStep(){	
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
	        ALVal_tS = subproblemSolvers[tS]->updateALValues(omega_centre[tS],z_current,scaling_matrix[tS]);
		ALVal_Local += pr[tS]*ALVal_tS;
		sqrDiscrNorm_tS = subproblemSolvers[tS]->getSqrNormDiscr();
		localDiscrepNorm += pr[tS]*sqrDiscrNorm_tS;

		//****************** Compute Next Vertex **********************
	// Find next vertex and Lagrangian lower bound
//if(tS==0) subproblemSolvers[tS]->printC();
		recordKeeping[tS][2]=ALVal_tS;
		recordKeeping[tS][3]=sqrDiscrNorm_tS;

		//Solve Lagrangian MIP

#if 0

    double *xSoln = subproblemSolvers[tS]->getXVertex();
    memcpy(totalSoln_,xSoln,n1*sizeof(double));
    double *ySoln = subproblemSolvers[tS]->getYVertex();
    memcpy(totalSoln_+n1,ySoln,n2*sizeof(double));
    int *all_indices = new int[n1+n2];
    for(int jj=0; jj<n1+n2; jj++){ all_indices[jj]=jj;}
    //subproblemSolvers[tS]->setInitialSolution(all_indices,totalSoln_);
    delete [] all_indices;
#endif

		
		dynamic_cast<PSCGScen_SMPS*>(subproblemSolvers[tS])->getOSI()->resolve();

		//LagrLB_tS = subproblemSolvers[tS]->getLagrBd();	
		spSolverStatuses_[tS] = subproblemSolvers[tS]->solveLagrangianProblem(omega_tilde[tS]);
		if(spSolverStatuses_[tS]==PSCG_PRIMAL_INF || spSolverStatuses_[tS]==PSCG_DUAL_INF){
	    	    LagrLB_Local = COIN_DBL_MAX;
	    	    ALVal_Local = COIN_DBL_MAX;
cerr << "performColGenStep(): Subproblem " << tS << " infeasible on proc " << mpiRank << endl;
	    	    break;
		}

	double LagrLB_tS = subproblemSolvers[tS]->optimiseLagrOverVertexHistory(omega_tilde[tS]);
#ifdef KEEP_LOG
    	    *(logFiles[tS]) << "old LagrLB_ts: " << LagrLB_tS << "\tnew LagrLB_ts" << subproblemSolvers[tS]->getLagrBd() << "\tALVal_ts: " << ALVal_tS << "\tsqrDiscrNorm: " << sqrDiscrNorm_tS << endl;
    	    *(logFiles[tS]) << "Testing whether vertex is redundant: " << -(LagrLB_tS - subproblemSolvers[tS]->getLagrBd()) << endl;;
#endif
	if(LagrLB_tS <= subproblemSolvers[tS]->getLagrBd()+1e-10){
#ifdef KEEP_LOG
    	    *(logFiles[tS]) << "  FLAGGING: Redundant vertex found" << endl;;
	    *(logFiles[tS]) << " Is vertex really redundant? " << subproblemSolvers[tS]->checkWhetherVertexIsRedundant() << endl;
#endif
	}
	else{
	    updateVertexHistory(tS);
	}
#if 1
		LagrLB_tS = subproblemSolvers[tS]->getLagrBd();	
		lhsCritVal = ALVal_tS+0.5*sqrDiscrNorm_tS;
#ifdef KEEP_LOG
		*(logFiles[tS]) << "Printing penalities: ";
#endif
		for(int ii=0; ii<n1; ii++){
		    lhsCritVal += scaling_matrix[tS][ii]*z_current[ii]*(x_current[tS][ii]-z_current[ii]);
#ifdef KEEP_LOG
		    *(logFiles[tS]) << " " << scaling_matrix[tS][ii];
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
//MPI_Barrier(MPI_COMM_WORLD);
	if (mpiSize > 1) {
		localReduceBuffer[0]=LagrLB_Local;
		localReduceBuffer[1]=ALVal_Local;
		localReduceBuffer[2]=localDiscrepNorm;
		MPI_Allreduce(localReduceBuffer, reduceBuffer, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

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
	        cout << recordKeeping[tS][1] << " <=??? " << recordKeeping[tS][2] << "  (discr: " << recordKeeping[tS][3] << ")" << endl;
	        cout << recordKeeping[tS][1] << " <=??? " << recordKeeping[tS][2]+0.5*recordKeeping[tS][3] << "  (discr: " << recordKeeping[tS][3] << ")" << endl;
		cout << "Current solution value: " << subproblemSolvers[tS]->evaluateSolution(omega_centre[tS]) << endl;;

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

