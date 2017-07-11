/*Header for the main procedure ParallelSCG. */

#ifndef PARALLELSCG_H
#define PARALLELSCG_H

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <sys/time.h>
#include "tclap/CmdLine.h"
#include "StructureDefs.h"

#define smallNumber 0.000001

#ifdef USING_MPI
   #include <mpi.h>
#endif

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#define MAX_INNER_LOOP 200
#define DEFAULT_MAX_OUTER_LOOP 20 //TODO have this automatically update the help text below. Probably overkill though.

// Parameters for newConvergenceCriterion
#define BG_BETA 0.05

typedef double* ptrdouble;

/* TODO:
* plots
* Does scaling respect probabilities?
* Equality constraints
* Output precision
* dynamically allocating scenarios
* clustering
*/

// If you want to add a flag, search ADDFLAG in this file and the corresponding .cpp file and follow instructions.

// ADDFLAG : Add the variable here
typedef struct Params{
	std::string filename;
	std::string outputFilename;
	int noScenarios;
	int maxStep;
	int maxSeconds;
	int fixInnerStep;
	int UseVertexHistory;
	int threads;
	double penalty;
	double penaltyMult;

	int filetype;

	bool verbose;
	bool debug;
	bool linRelaxFirst;
	bool linRelaxSecond;
	bool scaling;
	bool dataPathOverride;
	bool LBcalc;
	bool AlgorithmC;
	bool disableHeuristic;
	Params() : filename(""), noScenarios(-1), maxStep(-1), maxSeconds(-1),
   fixInnerStep(-1), UseVertexHistory(-1), penalty(-1), penaltyMult(-1),
   filetype(0), verbose(false), debug(false), linRelaxFirst(false),
   linRelaxSecond(false), scaling(false), dataPathOverride(false),
   LBcalc(false), AlgorithmC(false), disableHeuristic(false) {}
} Params;

// ADDFLAG : Put the definition and declaration here.
typedef struct CArgs {
	TCLAP::ValueArg<std::string> filenameArg;
	TCLAP::ValueArg<std::string> configFileArg;
	TCLAP::ValueArg<std::string> outputFileArg;
	TCLAP::ValueArg<int> nsArg;
	TCLAP::ValueArg<int> innerStepArg;
	TCLAP::ValueArg<int> stepArg;
	TCLAP::ValueArg<int> maxSecondsArg;
	TCLAP::ValueArg<int> vertexHistoryArg;
	TCLAP::ValueArg<int> threadsArg;
	TCLAP::ValueArg<double> ppArg;
	TCLAP::ValueArg<double> pmultArg;

	TCLAP::SwitchArg verboseSwitch;
	TCLAP::SwitchArg debugSwitch;
	TCLAP::SwitchArg linSwitch;
	TCLAP::SwitchArg linFirstSwitch;
	TCLAP::SwitchArg linSecondSwitch;
	TCLAP::SwitchArg scalingSwitch;
	TCLAP::SwitchArg CAP_Switch;
	TCLAP::SwitchArg SMPS_Switch;
	TCLAP::SwitchArg LB_Switch;
	TCLAP::SwitchArg AlgC_Switch;
	TCLAP::SwitchArg Heur_Switch;

	CArgs(TCLAP::CmdLine &cmdL) :
		verboseSwitch("v","verbose","Print verbose output", false),
		debugSwitch("","debug","Debug output", false),
		linSwitch("l","linrelax","Toggle full linear relaxation in CAP problem", false),
		linFirstSwitch("","linfirst","Toggle first-stage linear relaxation in CAP problem", false),
		linSecondSwitch("","linsecond","Toggle first-stage linear relaxation in CAP problem", false),
		scalingSwitch("","scaling","Toggle penalty parameter scaling --DEPRECATED, DO NOT USE--", false),
		CAP_Switch("","CAP","Indicates the file to be read is a CAP file", false),
		SMPS_Switch("","SMPS","Indicates the files to be read comprise an SMPS file", false),
		LB_Switch("","lb","Toggle calculation of explicit lower bound", false),
		AlgC_Switch("","AlgC","Toggles use of Algorithm C variant", false),
		Heur_Switch("","disableHeur","Turns off use of CPLEX integer solution heuristic", false),

		ppArg("p", "penaltyParam", "Starting value for penalty --OVERRIDES pMultiplier--", false, -1, "float"),
		pmultArg("", "pMultiplier", "Multiplier for penalty chosen by heuristic", false, -1, "float"),
		nsArg("s", "noScenarios", "Number of scenarios to read", false, -1, "int"),
		stepArg("", "maxStep", "Maximum number of outer loop steps. Defaults to 20.", false, -1, "int"),
		maxSecondsArg("t", "maxTime", "Algorithm will terminate once it reaches this many seconds. Defaults to no limit. The job itself should have a time limit about 3 times this, just to be safe.", false, -1, "int"),
		innerStepArg("", "innerStep", "Fixes the number of inner-loop iterations for each outer-loop iteration. If not set, we default to 1.", false, -1, "int"),
		vertexHistoryArg("","vHistory", "Sets vertex history state (number of past vertices used, default is -1 for only most recent vertex. Zero means use all previous vertices.", false, -1, "int"),
		threadsArg("","threads", "Set number of threads CPLEX should use. Default is one.", false, -1, "int"),

		filenameArg("f","filename","Filename to read problem from",false,"","string"),
		outputFileArg("o","outputFile","Filename to write auxiliary output to --YOU SHOULD NOT NEED THIS--",false,"","string"),
		configFileArg("c","config","Filename to read additional flags from. Overridden by any command line flag.",false,"","string")
		{
		cmdL.add( verboseSwitch );
		cmdL.add( debugSwitch );
		cmdL.add( linSwitch );
		cmdL.add( linFirstSwitch );
		cmdL.add( linSecondSwitch );
		cmdL.add( scalingSwitch );
		cmdL.add( CAP_Switch );
		cmdL.add( SMPS_Switch );
		cmdL.add( LB_Switch );
		cmdL.add( AlgC_Switch );
		cmdL.add( Heur_Switch );
		cmdL.add( ppArg );
		cmdL.add( pmultArg );
		cmdL.add( nsArg );
		cmdL.add( stepArg );
		cmdL.add( maxSecondsArg );
		cmdL.add( innerStepArg );
		cmdL.add( vertexHistoryArg );
		cmdL.add( threadsArg );
		cmdL.add( configFileArg );
		cmdL.add( outputFileArg );
		cmdL.add( filenameArg );
	}
} CArgs;

//*** FWMM functions.
double computeNormedDiscrepancies(double** scaling_matrix, double** x, double* z, int nS, int n1);
void weightedAverage(const std::vector<double*> &x, const std::vector<double> &p, double *localZ, double* z, int nNodeSPs, int n1, int mpiRank);
void computeFeasibleStartingPoint(int tS, double* x, double* yFeasible);

//*** Command line parameter handlers.
Params* readParameters(Params* p, int argc, char** argv);
void updateParams(Params* p, CArgs* a);
int validateParams(Params* p);
//String manipulation helper
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

//*** Helper functions.
//Frobenius norm
double froNorm(double** array, int nS, int n1);
double froNorm(double* array, int n1);
//Deals with numerical imprecision from solver (MIPs)
double roundIfClose(double input);

#endif
