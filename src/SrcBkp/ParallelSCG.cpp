/*
Finds the optimal solution to the Lagrangian dual of a two-stage stochastic
integer program using a Frank-Wolfe-based Method of Multipliers approach.
*/

#include "ParallelSCG.h"
//#include "CPLEXsolverSCG.h"
#include "ProblemDataBodur.h"
#include "CPLEXSCGSolver.h"
#include "Stopwatch.h"
#include "TssModel.h"

#define ROUNDING_TOLERANCE 1e-4
#define OUTER_LOOP_TERMINATION 1e-6
#define TIME_TYPES 2
#define DEFAULT_THREADS 1
#define SSC_PARAM 0.1
#define DEFAULT_PENALTY 100
#define DEFAULT_MAX_OUTER_LOOP 20 
#define KIWIEL_PENALTY 0 //set 1 to use Kiwel (2006) penalty update rule


using namespace std;

class myexception: public exception
{
  virtual const char* what() const throw()
  {
    std::cerr << "My exception happened" << endl;
	return "";
  }
} myex;

int main(int argc, char **argv) {

	//******************Wall Timing Setup****************

	//Start timing. More refined timing is possible for debugging.
	Stopwatch totalTimer;
   	Stopwatch interpTimer;
   	Stopwatch updateTimer;
   	double totalTimeQP [TIME_TYPES] = {0,0};
   	double totalTimeMIP [TIME_TYPES] = {0,0};

   	totalTimer.start();
   	double totalTimeAllStep [TIME_TYPES] = {0,0};

	//******************MPI Setup**********************
	int mpiRank;
	int mpiSize;
	bool mpiHead;
   	bool parallel;

   	//identifying whether the code is going to run in  parallel or not
	#ifdef USING_MPI
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
		parallel = (mpiSize > 1);
	#else
		mpiRank = 0;
		mpiSize = 1;
   		parallel = false;
	#endif

   	//Flag for the head node recognise itself as such
	mpiHead = (mpiRank == 0);

	if (mpiHead) {
		#ifdef USING_MPI
			std::cout << "USINGMPI:TRUE" << endl;
		#else
			std::cout << "USINGMPI:FALSE" << endl;
		#endif
	}

   	//******************Read Command Line Parameters**********************
	Params par;
	readParameters(&par, argc, argv);

	std::string str_filename(par.filename);
	std::string str_outputFileID(par.outputFilename);
	int nS = par.noScenarios;
	double penC = par.penalty;
	double penMult = par.penaltyMult;

	int maxStep = par.maxStep;
	int maxSeconds = par.maxSeconds;
	int fixInnerStep = par.fixInnerStep;
	int nVerticesUsed = par.UseVertexHistory;
	int nThreads = par.threads;
	if (nThreads < 0) { nThreads = DEFAULT_THREADS; }
	bool verbose = par.verbose;
	bool debug = par.debug;
	bool linRelaxFirst = par.linRelaxFirst;
	bool linRelaxSecond = par.linRelaxSecond;
	bool scaling = par.scaling;
	bool dataPathOverride = par.dataPathOverride;
	bool LBcalc = par.LBcalc;
	bool AlgorithmC = par.AlgorithmC; //obsolete, get rid of.
	bool disableHeuristic = par.disableHeuristic;
	int ftype = par.filetype;

	bool useVertexHistory;
	if (nVerticesUsed < 0) { useVertexHistory = false; }
	if (nVerticesUsed == 0) { useVertexHistory = true; }
	if (nVerticesUsed > 0) { useVertexHistory = true; }

   	//Hard coded standard input file location
   	//const std::string DATA_FILE_PATH("/short/ka3/Data");
   	const std::string DATA_FILE_PATH("/homes/bcdandurand/Data");
	std::string filepath(DATA_FILE_PATH);

   	//(FABS) This seems to be useless. No overriding seems to be performed.
   	if (dataPathOverride) {filepath = "";}

	//******************Display Command Line Parameters**********************
	if (mpiHead) {
		std::cout << "ALGORITHM:FWPH" << std::endl;
		std::cout << "Number of processors: " << mpiSize << std::endl;
		std::cout << "Data file path: " << filepath << std::endl;
		std::cout << "Problem filename: " << str_filename << std::endl;

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
	}

	//******************Penalty Parameter Setup*********************
	
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

	if (mpiHead) {
		std::cout << std::endl << std::endl;
	}

	//******************Termination, Bounds and Norms Setup**********************

	// Termination conditions
	bool terminate = false;
	bool outerConvCritReached = false;
	bool timeLimitReached = false;
	bool maxStepReached = false;
	int totalInnerSteps = 0;

	// Bounds
	double LagrLB_Local = 0.0;
	double LagrLB = 0.0;
	double currentLagrLB = 0.0;
	double ALVal_Local = 0.0;
	double ALVal = 0.0;
	
	// Norms
	double localDiscrepNorm = 1e9;
	double discrepNorm = 1e9;

	double localReduceBuffer[3]; //0-LagrLB_Local,  1-ALVal_Local,  2-localDiscrepNorm
	double reduceBuffer[3];	//0-LagrLB,  1-ALVal,  2-discrepNorm


	//******************Reading Problem Data**********************

	int n1;
	int n2;
	vector<double> pr;

    //fileRequest/ fileReply are structs with the info for instantiating the problem
	SMIP_fileRequest fileRequest(filepath.c_str(), str_filename.c_str(), ftype, nS, mpiHead, disableHeuristic, nThreads);
	SMIP_fileReply fileReply;
	std::string full_filename;

    //Flags for linear relaxations of 1st and 2nd stage variables
	fileRequest.lin1 = linRelaxFirst;
	fileRequest.lin2 = linRelaxSecond;

   	// We'll need this in case it is a CAP instance
   	ProblemDataBodur pdBodur(nS);
	// We'll need this in case it is read from SMPS
	TssModel smpsModel;

   	switch( ftype )	{
		case 1: //CAP Problem
      			pdBodur.initialise(&fileRequest, &fileReply);
			n1 = fileReply.n1;
			n2 = fileReply.n2;
			nS = fileReply.nS;
			break;
		case 2: //SMPS problem
			full_filename += filepath.c_str();
			full_filename += str_filename.c_str();
			smpsModel.readSmps(full_filename.c_str());
			n1 = smpsModel.getNumCols(0);
			n2 = smpsModel.getNumCols(1);
			nS = smpsModel.getNumScenarios();
			break;
		default:
			throw(-1);
			break;
	}
	

	if (mpiHead) {
		std::cout << "Problem data: " << std::endl;
		std::cout << "Number of first stage variables: " << n1 << std::endl;
		std::cout << "Number of second stage variables: " << n2 << std::endl;
		std::cout << "Number of scenarios: " << nS << std::endl;
    }

	//******************Assign Scenarios**********************

	int scenarioAssign[nS];
	for (int tS = 0; tS < nS; tS++) {
		scenarioAssign[tS] = tS % mpiSize;
	}
	
	vector<int> scenariosToThisNode;
	
	for (int sss = 0; sss < nS; sss++) {
	    if (scenarioAssign[sss] == mpiRank) {
		scenariosToThisNode.push_back(sss);
	    }
	}
	
	int nNodeSPs = scenariosToThisNode.size();
	
	//******************Decision Variable and Data Storage Setup**********************

	vector<double*> scaling_matrix; //glorified rho
	vector<double*> dvar_tilde; //dvar == omega
	vector<double*> dvar_current;

	vector<double*> x_current;
	vector<double*> y_current;

	double* z_current = new double[n1];
	double* z_local = new double[n1];

	//initialising zˆ0=0
	for (int i = 0; i < n1; i++) {
	    z_current[i] = 0;
	}

	//******************Decision Variable Initialisation**********************

	int step = 0; 

	#if 0
		forInitTimes.stop();
		forInitTimes.getTime(initializationTime);
		if(mpiRank==0) cout << "Initialization took: " << initializationTime[0] << " seconds." << endl;
		totalTimer.start();
	#endif

	vector<CPLEXSCGSolver*> subproblemSolver(nNodeSPs);
	
	double LagrLB_Init_Local = 0.0;

	for (int tS = 0; tS < nNodeSPs; tS++) {
	    //CPLEXsolverSCG &currentSPSolver = subproblemSolver[tS];
	    CPLEXSCGSolver *&currentSPSolver = subproblemSolver[tS];
	    switch( ftype ){
		    case 1: //CAP Problem
			currentSPSolver = new CPLEXSCGSolver_Bodur();
		     	dynamic_cast<CPLEXSCGSolver_Bodur*>(currentSPSolver)->initialiseBodur(pdBodur,&fileRequest,scenariosToThisNode[tS]);
		      	break;
		    case 2: //SIPLIB problems
			currentSPSolver = new CPLEXSCGSolver_SMPS();
		      	dynamic_cast<CPLEXSCGSolver_SMPS*>(currentSPSolver)->initialiseSMPS(&fileRequest,smpsModel,scenariosToThisNode[tS]); 
		      	break;
		    default:
				throw(-1);
				break;
	    }
	    currentSPSolver->finishInitialisation(); 
	    
	    x_current.push_back(currentSPSolver->getX());
	    y_current.push_back(currentSPSolver->getY());
	    pr.push_back(currentSPSolver->getProbabilities());
	    scaling_matrix.push_back(new double[n1]);
	    dvar_tilde.push_back(new double[n1]);
	    dvar_current.push_back(new double[n1]);
	    
	    for (int i = 0; i < n1; i++) {
			scaling_matrix[tS][i] = penC ;
			dvar_current[tS][i] = 0.0;
	    }
	    currentSPSolver->setQuadraticTerm(scaling_matrix[tS]);
	}
	#if 1 
		switch( ftype )	{
			case 1:
				pdBodur.cleanup(); //double-check whether this is being used
				break;
			case 2:
				break;
			default:
				throw(-1);
				break;
		}
	#endif
	
	// This is part of the initialisation - initial Lagrangian subproblem computations
	for (int tS = 0; tS < nNodeSPs; tS++) {
 	    //updateTimer.start();
	    subproblemSolver[tS]->solveLagrangianProblem(dvar_current[tS]);
	    subproblemSolver[tS]->setXToVertex();
	    subproblemSolver[tS]->setYToVertex();

	    //updateTimer.stop();
	    //updateTimer.addTime(updateTimeThisStep);

	    LagrLB_Init_Local += pr[tS]*(subproblemSolver[tS]->getLagrBd());//obj;

	    if (useVertexHistory) {
		subproblemSolver[tS]->updateVertexHistory();
	    }
	}
	
	//Update of z
	weightedAverage(x_current, pr, z_local, z_current, nNodeSPs, n1, mpiRank);

	#ifdef USING_MPI
	if (parallel) {
			// Send decision variables to other processes.
		MPI_Allreduce(&LagrLB_Init_Local, &LagrLB, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

			//commsTimer.stop();
			//commsTimer.addTime(commsTimeThisStep);
	}
	#endif
	
	LagrLB = LagrLB_Init_Local;
	currentLagrLB = LagrLB;

	//totalTimer.stop();
	//totalTimer.addTime(totalTimeThisStep);

	//******************Output Initialisation**********************

	if (mpiHead) {
		std::cout << "Initialisation step:" << std::endl;
		printf("Lagrangian Lower Bound: %0.9g\n", currentLagrLB);
		#if 0
			printf("QFLOB\tExplicit Lagrangian Lower Bound: %0.9g\n", LagrLB_Explicit);
			printf("QFTIM\tTime: %-7.3f %-7.3f\n", totalTimeThisStep[0], totalTimeThisStep[1]);
			printf("QFSPT\tTime spent on update subproblems by this process: %-7.3f %-7.3f\n", updateTimeThisStep[0], updateTimeThisStep[1]);
		#ifdef USING_MPI
			printf("QFCMT\tTime spent on communication: %-7.3f %-7.3f\n", commsTimeThisStep[0], commsTimeThisStep[1]);
			printf("QFBLK\tTime spent blocked by this process: %-7.3f %-7.3f\n", blockedTimeThisStep[0], blockedTimeThisStep[1]);
		#endif
			printf("QFLBT\tAdditional time spent on bound calculation: %-7.3f %-7.3f\n", boundTimeThisStep[0], boundTimeThisStep[1]);
		#endif
	
	#if 1
		if (verbose) {
			printf("\tConsensus values for first-stage decision variables:\n\t[");
			for (int i = 0; i < n1; i++) {
				printf("%0.6g ", z_current[i]);
			}
			printf("]\n");
		}
		printf("\n");
		cout << endl;
	#endif
	}

	//*******************************************************************
	//******************Outer FWPH Loop Begins Here**********************
	//*******************************************************************

	while (terminate == false) {
		#if 0
			for (int ttype = 0; ttype < TIME_TYPES; ttype++) {
				totalTimeThisStep[ttype] = 0;
				boundTimeThisStep[ttype] = 0;
				commsTimeThisStep[ttype] = 0;
				updateTimeThisStep[ttype] = 0;
				blockedTimeThisStep[ttype] = 0;
				interpTimeThisStep[ttype] = 0;
			}
		#endif
		//totalTimer.start();

		step++;

		//*******************************************************************
		//******************Inner FWPH Loop Begins Here**********************
		//*******************************************************************

		int innerStep = 0;

 		if(step>1){
			for(int itGS=0; itGS < fixInnerStep; itGS++) { //The inner loop has a fixed number of occurences
		    	for (int tS = 0; tS < nNodeSPs; tS++) {

					//*************************** Quadratic subproblem ***********************************
				    
				    if (useVertexHistory) { //Compute Next X, Y (With History)

						interpTimer.start();
						
						//Solve the quadratic master problem for x and y
						subproblemSolver[tS]->updatePrimalVariablesHistory_OneScenario(dvar_current[tS],z_current);
						interpTimer.stop();
						interpTimer.addTime(totalTimeQP);
						
						if(tS == nNodeSPs-1) {
							cout << "Walltime spent on QP (step: " << itGS+1 << "):" << totalTimeQP[0] << endl;
						}
						// note: the final weight corresponds to the existing x
				    }
				    else { //might not work correctly, untested! Compute Next X, Y (Without History)

						//SMIP_qu_singleVarUpdate question(tS, scaling_matrix[tS], x_current[tS], y_current[tS], NULL, NULL, z_current, dvar_current[tS]);
						//SMIP_ans_singleVarUpdate answer;
						//subproblemSolver[tS].updatePrimalVariables_OneScenario(&question, &answer);
						subproblemSolver[tS]->updatePrimalVariables_OneScenario(dvar_current[tS],z_current,scaling_matrix[tS]);
				    }
			    }
			    //clock_t end = clock();
			    //cout << "Processor " << mpiRank << " took " << (double)(end-start) / CLOCKS_PER_SEC  << " seconds for QSP." << endl;;
	    					
				// Update z_previous.
			    weightedAverage(x_current, pr, z_local, z_current, nNodeSPs, n1, mpiRank);
			}
		}

		//cout << "Done with step " << step << " QSP updates from node " << mpiRank << endl;
		//clock_t start2 = clock();
		
		LagrLB_Local = 0.0;
		ALVal_Local = 0.0;
		localDiscrepNorm = 0.0;
		
		for (int tS = 0; tS < nNodeSPs; tS++) {

			//****************** Compute Next Vertex **********************
			// Find next vertex and Lagrangian lower bound
			
			subproblemSolver[tS]->updateALValues(dvar_current[tS],z_current,scaling_matrix[tS]);
			ALVal_Local += pr[tS]*subproblemSolver[tS]->getALVal();
			localDiscrepNorm += pr[tS]*subproblemSolver[tS]->getSqrNormDiscr();
			
			for (int i = 0; i < n1; i++) {
			    dvar_tilde[tS][i] = dvar_current[tS][i] + scaling_matrix[tS][i] * (x_current[tS][i] - z_current[i]);
			}

			//Solve Lagrangian MIP
			updateTimer.start();
			subproblemSolver[tS]->solveLagrangianProblem(dvar_tilde[tS]);
			updateTimer.stop();
			updateTimer.addTime(totalTimeMIP);
			
			if(tS == nNodeSPs-1) {
				cout << "Walltime spent on MIP: " << totalTimeMIP[0] << endl;
			}
			//updateTimer.addTime(updateTimeThisStep);

			//Acumulate the values for the Lagrangian subproblems
			LagrLB_Local += pr[tS]*subproblemSolver[tS]->getLagrBd();//FWSPoptvalk;}

		 	if(useVertexHistory){
			    subproblemSolver[tS]->updateVertexHistory();

			    if (nVerticesUsed > 0 && subproblemSolver[tS]->getNumVertices() > nVerticesUsed) {
					subproblemSolver[tS]->removeBackVertex();
			    }
			}
		}

		//cout << "Done with step " << step << " MIP updates from node " << mpiRank << endl;
		// Iteration over scenarios ends here
		//clock_t end2 = clock();
		//cout << "Processor " << mpiRank << " took " << (double)(end2-start2) / CLOCKS_PER_SEC << " seconds for vertex-finding." << endl;;

		for (int tS = 0; tS < nNodeSPs; tS++) {
			for (int i = 0; i < n1; i++) {
				x_current[tS][i] = roundIfClose(x_current[tS][i]);
			}			
		}

		//Z UPDATE
		#ifdef USING_MPI
			if (parallel) {
				//commsTimer.start();
				// Send decision variables to other processes.

				localReduceBuffer[0]=LagrLB_Local;
				localReduceBuffer[1]=ALVal_Local;
				localReduceBuffer[2]=localDiscrepNorm;

				MPI_Allreduce(localReduceBuffer, reduceBuffer, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				LagrLB = reduceBuffer[0];
				ALVal = reduceBuffer[1];
				discrepNorm = reduceBuffer[2];
				if(mpiRank==0) cout << "Testing: " << ALVal+ discrepNorm  << " >= " << currentLagrLB << endl;

				//commsTimer.stop();
				//commsTimer.addTime(commsTimeThisStep);
			}
		#endif

		if (!parallel) {
			LagrLB = LagrLB_Local;
			ALVal = ALVal_Local;
			discrepNorm = localDiscrepNorm;
			cout << "Testing: " << ALVal+ discrepNorm  << " >= " << currentLagrLB << endl;
		}

		//********************** Serious Step Condition (SSC) **********************
		double SSCVal = (LagrLB-currentLagrLB)/(ALVal + 0.5*discrepNorm - currentLagrLB);
		
		if(SSCVal >= SSC_PARAM) {
		    for (int tS = 0; tS < nNodeSPs; tS++) {
				for (int i = 0; i < n1; i++) {
			    	dvar_current[tS][i] = dvar_tilde[tS][i];
				}
		    }
		    currentLagrLB = LagrLB;
		}
		else {
		    if(mpiRank==0) cout << "Null step taken." << endl;
		}

		
		//********************** Penalty update **********************
		//This should normally be turned off, this is a rule for updating the 
		//penalty from Kiwiel 2006 and Lubin et al.
	
		#if KIWIEL_PENALTY 
			penC = 1.0/min(max( (2.0/penC)*(1.0-SSCVal),  max(1.0/(10.0*penC),1e-4)    ), 10.0/penC);
			for (int tS = 0; tS < nNodeSPs; tS++) {
			    subproblemSolver[tS].setQuadraticTerm(penC);
		    	    for (int i = 0; i < n1; i++) {
				scaling_matrix[tS][i] = penC;
		    	    }
			}
		#endif

		//****************** Timing **********************

		//totalTimer.stop();
		//totalTimer.addTime(totalTimeThisStep);

		//****************** Share Timing Data ************************

		#if 0
			for (int ttype = 0; ttype < TIME_TYPES; ttype++) {
				totalTimeAllStep[ttype] += totalTimeThisStep[ttype];
				boundTimeAllStep[ttype] += boundTimeThisStep[ttype];
				commsTimeAllStep[ttype] += commsTimeThisStep[ttype];
				updateTimeAllStep[ttype] += updateTimeThisStep[ttype];
				blockedTimeAllStep[ttype] += blockedTimeThisStep[ttype];
			}
		#endif

		//****************** Check Termination Conditions **********************

		if (step >= maxStep) { terminate = true; maxStepReached = true; }

		if (ALVal + discrepNorm  - currentLagrLB < OUTER_LOOP_TERMINATION)
		{ terminate = true; outerConvCritReached = true; }

		#ifndef USING_MPI
			if (maxSeconds > 0 && totalTimeAllStep[0] > maxSeconds)
			{ terminate = true; timeLimitReached = true; }
		#endif

		//****************** Output **********************

		if (mpiHead) {
			printf("\nStep %d\n", step);
			//printf("\tNorm of z difference: %0.6g\n", zNorm);
			printf("Lagrangian Lower Bound: %0.9g\n", currentLagrLB);
			printf("Aug. Lagrangian value: %0.9g\n", ALVal);
			printf("Norm of primal discrepancy: %0.6g", discrepNorm);
			//printf("QFDES\tDescent on first inner step: %0.9g\n", descent);
			//printf("QFLOB\tExplicit Lagrangian Lower Bound: %0.9g\n", LagrLB_Explicit);
			#if 0
				printf("QFTIM\tTime: %-7.3f %-7.3f\n", totalTimeThisStep[0], totalTimeThisStep[1]);
				printf("QFSPT\tTime spent on update subproblems by this process: %-7.3f %-7.3f\n", updateTimeThisStep[0], updateTimeThisStep[1]);
				#ifdef USING_MPI
					printf("QFCMT\tTime spent on communication: %-7.3f %-7.3f\n", commsTimeThisStep[0], commsTimeThisStep[1]);
					printf("QFBLK\tTime spent blocked by this process: %-7.3f %-7.3f\n", blockedTimeThisStep[0], blockedTimeThisStep[1]);
				#endif
				printf("QFLBT\tAdditional time spent on bound calculation: %-7.3f %-7.3f\n", boundTimeThisStep[0], boundTimeThisStep[1]);
				printf("QFINT\tTime spent interpolating vertices: %-7.3f %-7.3f\n", interpTimeThisStep[0], interpTimeThisStep[1]);
			#endif


			if (verbose || terminate) {
				printf("\nConsensus values for first-stage decision variables:\n[");
				
				for (int i = 0; i < n1; i++) {
					printf("%0.2g ", z_current[i]);
				}
			
				printf("]\n");
			}	
		}
		
		printf("\n");

		if (terminate) {

			if (timeLimitReached)
			{ std::cout << "Time limit reached" << std::endl; }

			if (outerConvCritReached)
			{ std::cout << "Outer loop reached termination criterion" << std::endl; }

			if (maxStepReached)
			{ std::cout << "Outer loop terminated at maxStep" << std::endl; }

			std::cout << "Outer loop terminated at step: " << step << std::endl;
		}
		cout << endl;
	}
	
	//***********************************************************************
	//******************Outer FWPH Loop Terminates Here**********************
	//***********************************************************************

	#if 0 
		if (mpiHead) {
			double walltime = totalTimeAllStep[0];
			double updatetime = updateTimeAllStep[0];
			double commstime = commsTimeAllStep[0];
			double boundtime = boundTimeAllStep[0];
			double blocktime = blockedTimeAllStep[0];

			printf("Total inner steps required: %d\n", totalInnerSteps);

			printf("Total walltime required: %0.2f seconds\n", walltime);
			if (walltime > 600) {
				if (walltime > 7200)
					{ printf("which is: %0.2f hours\n", walltime/3600); }
				else
					{ printf("which is: %0.2f minutes\n", walltime/60); }
			}
			printf("Walltime spent on subproblems summed across processes: %0.2f seconds\n", updatetime);
			if (updatetime > 600) {
				if (updatetime > 7200)
					{ printf("which is: %0.2f hours\n", updatetime/3600); }
				else
					{ printf("which is: %0.2f minutes\n", updatetime/60); }
			}

			printf("Walltime spent on communication: %0.2f seconds\n", commstime);
			if (commstime > 600) {
				if (commstime > 7200)
					{ printf("which is: %0.2f hours\n", commstime/3600); }
				else
					{ printf("which is: %0.2f minutes\n", commstime/60); }
			}

			printf("Walltime spent blocked summed across processes: %0.2f seconds\n", blocktime);
			if (blocktime > 600) {
				if (blocktime > 7200)
					{ printf("which is: %0.2f hours\n", blocktime/3600); }
				else
					{ printf("which is: %0.2f minutes\n", blocktime/60); }
			}

			printf("Additional walltime spent on bounds: %0.2f seconds\n", boundtime);
			if (boundtime > 600) {
				if (boundtime > 7200)
					{ printf("which is: %0.2f hours\n", boundtime/3600); }
				else
					{ printf("which is: %0.2f minutes\n", boundtime/60); }
		}

		//printf("Explicit lower bound: %0.9f\n", LagrLB_Explicit);
		printf("Lower bound: %0.9f\n", LagrLB);
		printf("\n");

		std::cout << "RUN:SUCCESS" << std::endl;
		}
	#endif

	//outstrm.close();
	//clock_t startCleanup = clock();

	//Clean up structures
	for (int tS = 0; tS < nNodeSPs; tS++) {
		delete [] scaling_matrix[tS];
		delete [] dvar_tilde[tS];
		delete [] dvar_current[tS];
	}
	delete[] z_current;
	delete[] z_local;

	#if 0
		clock_t endCleanup = clock();
		if(mpiRank==0){
		    cout << "Processor " << mpiRank << " took " << (double)(endCleanup-startCleanup) / CLOCKS_PER_SEC << " seconds for final cleanup." << endl;;
		}
	#endif

	#ifdef USING_MPI
		MPI_Finalize();
	#endif
	
	totalTimer.stop();
	totalTimer.getTime(totalTimeAllStep);
	
	if (mpiHead) {
		double walltime = totalTimeAllStep[0];
		double cpuTime = totalTimeAllStep[1];
		printf("Total inner steps required: %d\n", totalInnerSteps);

		printf("Total walltime required: %0.2f seconds\n", walltime);
		if (walltime > 600) {
			if (walltime > 7200)
				{ printf("which is: %0.2f hours\n", walltime/3600); }
			else
				{ printf("which is: %0.2f minutes\n", walltime/60); }
		}
		printf("Total CPU time required: %0.2f seconds\n", cpuTime);
		if (cpuTime > 600) {
			if (cpuTime > 7200)
				{ printf("which is: %0.2f hours\n", cpuTime/3600); }
			else
				{ printf("which is: %0.2f minutes\n", cpuTime/60); }
		}
		
		printf("Lower bound: %0.9f\n", currentLagrLB);
		printf("\n");

		//Summary: walltime, no. iterations, ave. walltime, no. procs., no. GSIts, currentLagrLB
		cout << "Summary: " << walltime << "s,  " << step << " steps,  " << walltime/((double)step) << "s/step,  " 
		<< mpiSize << " node,  " << fixInnerStep << " innersteps,  Lag. LB: " << currentLagrLB << endl;
		std::cout << "RUN:SUCCESS" << std::endl;
	}
	return 0;
}



//******************FWPH functions**********************

void weightedAverage(const vector<double*> &x, const vector<double> &p, double *localZ, double* z, int nNodeSPs, int n1, int mpiRank) {
	for (int i = 0; i < n1; i++)
	{
	    localZ[i] = 0.0;
	    for (int tS = 0; tS < nNodeSPs; tS++)
	    {
		localZ[i] += x[tS][i] * p[tS];
	    }
	}
	#ifdef USING_MPI
		MPI_Allreduce(localZ, z, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#else
	    for (int i=0; i<n1; i++) z[i] = localZ[i]; //Only one node, trivially localZ == z
	#endif
}

//NOT CORRECTLY IMPLEMENTED
void computeFeasibleStartingPoint(int tS, double* x, double* yFeasible) {
	SMIP_qu_secondStage question(tS, x);
	SMIP_ans_secondStage answer(yFeasible);
	//CPLEXsolverSCG subproblemSolver;
	//subproblemSolver.solveForSecondStage(&question, &answer);
}

//******************Command line paramater handler functions**********************

Params* readParameters(Params* p, int argc, char** argv) {
	try {
		TCLAP::CmdLine cmdL1("Solves a SIP using Progressive Hedging", ' ', "0.1");
		CArgs args1(cmdL1);
		cmdL1.parse( argc, argv );

		// If there is a config file, read from that first
		// This could be done much more elegantly and rigorously
		if (args1.configFileArg.getValue().compare("") != 0) {
			std::string line;
			std::vector<std::string> arguments;

			ifstream confile (args1.configFileArg.getValue().c_str());
			if (!confile) {
				std::cerr << "ERROR: could not open config file '" << confile
				<< "' for reading" << std::endl;
				throw(-1);
			}
			confile.close();

			getline (confile, line);
			split(line, ' ', arguments);

			// Make file contents look like a command line
			int config_argc;
			const char** config_argv = 0;

			config_argc = arguments.size() + 1;

			config_argv = new const char*[config_argc];
			config_argv[0] = argv[0];
			for (int i = 1; i < config_argc; i++) {
				config_argv[i] = arguments.at(i-1).c_str();
			}

			TCLAP::CmdLine cmdL2("Internal config file reader - you should never see this!", ' ', "0.0");
			CArgs args2(cmdL2);
			cmdL2.parse( config_argc, config_argv );
			updateParams(p, &args2);
		}

		// Arguments from the command line now overwrite any from the config file
		updateParams(p, &args1);

	} catch (TCLAP::ArgException &e)
	{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

	validateParams(p);
	return p;
}

// ADDFLAG : put flag handling here
void updateParams(Params* p, CArgs* a) {
	if (a->filenameArg.getValue().compare("") != 0) {
		p->filename = a->filenameArg.getValue();
	}

	if (a->outputFileArg.getValue().compare("") != 0) {
		p->outputFilename = a->outputFileArg.getValue();
	}

	if (a->nsArg.getValue() != -1) {
		p->noScenarios = a->nsArg.getValue();
	}

	if (a->ppArg.getValue() != -1) {
		p->penalty = a->ppArg.getValue();
	}

	if (a->pmultArg.getValue() > 0) {
		p->penaltyMult = a->pmultArg.getValue();
	}

	if (a->stepArg.getValue() != -1) {
		p->maxStep = a->stepArg.getValue();
	}

	if (a->maxSecondsArg.getValue() != -1) {
		p->maxSeconds = a->maxSecondsArg.getValue();
	}

	if (a->innerStepArg.getValue() != -1) {
		p->fixInnerStep = a->innerStepArg.getValue();
	}

	if (a->vertexHistoryArg.getValue() != -1) {
		p->UseVertexHistory = a->vertexHistoryArg.getValue();
	}

	if (a->threadsArg.getValue() != -1) {
		p->threads = a->threadsArg.getValue();
	}

	if (a->verboseSwitch.getValue() == true) {
		p->verbose = true;
	}

	if (a->debugSwitch.getValue() == true) {
		p->debug = true;
	}

	if ((a->linSwitch.getValue() == true) || (a->linFirstSwitch.getValue() == true)) {
		p->linRelaxFirst = true;
	}

	if ((a->linSwitch.getValue() == true) || (a->linSecondSwitch.getValue() == true)) {
		p->linRelaxSecond = true;
	}

	if (a->scalingSwitch.getValue() == true) {
		p->scaling = true;
	}

	if (a->CAP_Switch.getValue() == true)
	{
		if (p->filetype == 0) {
			p->filetype = 1;
		}
		else
		{
			std::cerr << "ERROR: Tried to set file type more than once" << endl;
			throw(-1);
		}
	}

	if (a->SMPS_Switch.getValue() == true)
	{
		if (p->filetype == 0) {
			p->filetype = 2;
		}
		else
		{
			std::cerr << "ERROR: Tried to set file type more than once" << endl;
			throw(-1);
		}
	}

	if (a->LB_Switch.getValue() == true) {
		p->LBcalc = true;
	}

	if (a->AlgC_Switch.getValue() == true) {
		p->AlgorithmC = true;
	}

	if (a->Heur_Switch.getValue() == true) {
		p->disableHeuristic = true;
	}
}

// ADDFLAG: put flag checking here
int validateParams(Params* p) {

	if (p->filename.compare("") == 0) {
		std::cerr << "ERROR: No problem file specified" << endl;
		throw(-1);
	}

	if (p->outputFilename.compare("") == 0) {
		p->outputFilename = "output";
	}

	if (p->noScenarios <= 0) {
		//std::cerr << "ERROR: Number of scenarios not specified" << endl;
		//throw(-1);
		//This now signifies we should use all scenarios in the file
		//Will only work for SMPS files, though.
	}

	if (p->penalty <= 0) {
		//std::cerr << "ERROR: Starting penalty parameter not specified" << endl;
		//throw(-1);
	}

	if (p->maxStep <= 0) {
		std::cerr << "WARNING: maxStep not specified, defaults to " << DEFAULT_MAX_OUTER_LOOP << endl;
		p->maxStep = DEFAULT_MAX_OUTER_LOOP;
	}

	if (p->filetype == 0) {
		std::cerr << "ERROR: File type not specified, use e.g. --SMPS" << endl;
		throw(-1);
	}
	return 1;
}

//******************Helper functions**********************

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

double roundIfClose(double input) {
	double tmp = round(input);
	if (abs(tmp-input) < ROUNDING_TOLERANCE)
		{ return tmp; }
	else
		{ return input; }

}
