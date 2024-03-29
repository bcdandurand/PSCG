/*
Author: Brian C. Dandurand (2017-2018)
Based on and modified from code developed by Jeffrey Christiansen, Brian Dandurand, and Fabricio Oliveira
at RMIT in Melbourne Australia with funding under project ARC DP 140100985 during 2015-2017.
CIs of that projects were Prof. Andrew Eberhard, Prof. Natashia Boland, and PI Prof. Jeffrey Linderoth.
*/

/*
Finds the optimal solution to the Lagrangian dual of a two-stage stochastic
integer program using a Frank-Wolfe-based Method of Multipliers approach.
*/

#include "PSCGScen.h"

using namespace std;


int main(int argc, char **argv) {
#if 1
	//******************Wall Timing Setup****************

	//Start timing. More refined timing is possible for debugging.

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

	//PSCGParams par;
	//par.readParameters(argc, argv);
	//par.setMPIParams(mpiSize,mpiRank);
	
	smpsModel.readSmps(argv[1]);	

	//PSCG model(&par);
	#ifdef USING_MPI
	 PSCG model(smpsModel,MPI_COMM_WORLD);
	#else
	 PSCG model(smpsModel);
	#endif
	model.setPhase(2);
	model.setRho(50);
	model.setBaselineRho(50);
        model.setMaxNoSteps(400);
        model.setMaxNoGSSteps(1);
        model.setMaxNoInnerSteps(1);
	model.setInnerSSCParam(0.5);
	model.computeBound();
	//double *z = new double[model.n1];
	#if 0
	vector<double*> omega;
	for(int tS=0; tS<model.nNodeSPs; tS++){
	    omega.push_back(new double[model.n1]);
	    for(int ii=0;ii<model.n1; ii++) omega[tS][ii]=0.0;
	}
	#endif

	
	//******************Decision Variable and Data Storage Setup**********************


	//******************Decision Variable Initialisation**********************

	//int step = 0; 

	#if 0
		forInitTimes.stop();
		forInitTimes.getTime(initializationTime);
		if(mpiRank==0) cout << "Initialization took: " << initializationTime[0] << " seconds." << endl;
		totalTimer.start();
	#endif

	

//	model.setupSolvers();
	

	//totalTimer.stop();
	//totalTimer.addTime(totalTimeThisStep);
	//model.initialIteration();

	//*******************************************************************
	//******************Outer FWPH Loop Begins Here**********************
	//*******************************************************************

#if 0
	step++;
	bool terminate = false;
	bool timeLimitReached = false;
	bool maxStepReached = false;
	int maxStep = par.maxStep;
	double maxSeconds = par.maxSeconds;
	while (terminate == false) {
		//totalTimer.start();


		//int innerStep = 0;
		//model.regularIteration();
		//cout << "Proc " << mpiRank << " solution is feasible: " << model.checkZHasFullRecourse() << endl;


		//****************** Check Termination Conditions **********************

		if (step >= maxStep) { terminate = true; maxStepReached = true; }

		if(false)// (model.shouldTerminate())
		{ terminate = true; }

		#ifndef USING_MPI
			if (maxSeconds > 0 && totalTimeAllStep[0] > maxSeconds)
			{ terminate = true; timeLimitReached = true; }
		#endif

		//****************** Output **********************
#if 0
		if (mpiHead) {
			printf("\nStep %d\n", step);
			model.printStatus();


			if(par.verbose || terminate) {
				model.printZ();
			}	
		}
#endif
		

		if (terminate) {
		    if(mpiRank==0){
			std::cout << "Number of vertices: " << model.subproblemSolvers[0]->getNVertices() << endl;
			model.subproblemSolvers[0]->printColTypesFirstStage();
			if (timeLimitReached)
			{ std::cout << "Time limit reached" << std::endl; }

			if (model.shouldTerminate())
			{ std::cout << "Outer loop reached termination criterion" << std::endl; }

			if (maxStepReached)
			{ std::cout << "Outer loop terminated at maxStep" << std::endl; }

			std::cout << "Outer loop terminated at step: " << step << std::endl;
		    }
		    break;
		}
		step++;
	}//end while
	
#endif
	//***********************************************************************
	//******************Outer FWPH Loop Terminates Here**********************
	//***********************************************************************

	//if(mpiRank==0) std::cout << endl;

	#ifdef USING_MPI
		MPI_Finalize();
	#endif
	

#if 0
	totalTimer.stop();
	totalTimer.getTime(totalTimeAllStep);
	
	if (mpiHead) {
		double walltime = totalTimeAllStep[0];
		double cpuTime = totalTimeAllStep[1];
		printf("Total inner steps required: %d\n", model.totalNoGSSteps);

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
		
		//printf("Lower bound: %0.9f\n", model.currentLagrLB);
		printf("Lower bound: %0.9f\n", model.getIncumbentVal());
		printf("\n");

		//Summary: walltime, no. iterations, ave. walltime, no. procs., no. GSIts, currentLagrLB
		cout << "Summary: " << walltime << "s,  " << step << " steps,  " << walltime/((double)step) << "s/step,  " 
		<< mpiSize << " node,  " << model.totalNoGSSteps << " GS steps,  Lag. LB: " << model.currentLagrLB << endl;
		std::cout << "RUN:SUCCESS" << std::endl;
	}
#endif
	//delete [] z;
#endif
	return 0;
}

