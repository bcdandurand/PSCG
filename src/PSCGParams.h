/*Header for the main procedure ParallelSCG. */

#ifndef PSCGPARAMS_H
#define PSCGPARAMS_H

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


#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
#define DEFAULT_MAX_OUTER_LOOP 20 //TODO have this automatically update the help text below. Probably overkill though.

enum Statuses{
    SP_STATUS=0,
    Z_STATUS
};

enum SPStatuses{
    SP_UNKNOWN=-1,
    SP_OPT=0, //PSCG termination criteria met
    SP_ITER_LIM, //otherwise feasible
    SP_INFEAS //at least one subproblem is infeasible
};

enum Z_Statuses{
    Z_OPT=0, //z is not only feasible, but optimal (PSCG meets termination criteria)
    Z_FEAS, //z is feasible (has both recourse and integer feas)
    Z_REC_INFEAS, //z has recourse (but its integer feas unknown)
    Z_INT_INFEAS, //z is integer feas (but its recourse is unknown)
    Z_INFEAS, //z is infeasible by both feasibility qualities
    Z_BOUNDED, //indicates status of z is irrelevant due to fathoming by bound
    Z_UNKNOWN
};




class PSCGParams{
public:
	std::string filename;
	std::string outputFilename;
	int noScenarios;
	int maxStep;
	int maxSeconds;
	int fixInnerStep;
	int UseVertexHistory;
	int threads;
	int mpiSize;
	int mpiRank;
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
	bool disableHeuristic;
	PSCGParams() : filename(""), noScenarios(-1), maxStep(-1), maxSeconds(-1),
   fixInnerStep(-1), UseVertexHistory(-1), penalty(-1), penaltyMult(-1),
   filetype(0), verbose(false), debug(false), linRelaxFirst(false),
   linRelaxSecond(false), scaling(false), dataPathOverride(false),
   LBcalc(false), disableHeuristic(false), mpiSize(1),mpiRank(0) {}

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
	TCLAP::ValueArg<int> algArg;
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
		Heur_Switch("","disableHeur","Turns off use of CPLEX integer solution heuristic", false),

		ppArg("p", "penaltyParam", "Starting value for penalty --OVERRIDES pMultiplier--", false, -1, "float"),
		pmultArg("", "pMultiplier", "Multiplier for penalty chosen by heuristic", false, -1, "float"),
		nsArg("s", "noScenarios", "Number of scenarios to read", false, -1, "int"),
		stepArg("", "maxStep", "Maximum number of outer loop steps. Defaults to 20.", false, -1, "int"),
		algArg("", "alg", "Algorithm to use. Defaults to 0 (baseline dual decomposition).", false, -1, "int"),
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
		cmdL.add( Heur_Switch );
		cmdL.add( ppArg );
		cmdL.add( pmultArg );
		cmdL.add( nsArg );
		cmdL.add( stepArg );
		cmdL.add( algArg );
		cmdL.add( maxSecondsArg );
		cmdL.add( innerStepArg );
		cmdL.add( vertexHistoryArg );
		cmdL.add( threadsArg );
		cmdL.add( configFileArg );
		cmdL.add( outputFileArg );
		cmdL.add( filenameArg );
	}
    } CArgs;


//*** Command line parameter handlers.
void readParameters(int argc, char** argv) {
	try {
		TCLAP::CmdLine cmdL1("Solves a SIP using Progressive Hedging", ' ', "0.1");
		CArgs args1(cmdL1);
		cmdL1.parse( argc, argv );

		// If there is a config file, read from that first
		// This could be done much more elegantly and rigorously
		if (args1.configFileArg.getValue().compare("") != 0) {
			std::string line;
			std::vector<std::string> arguments;

			std::ifstream confile (args1.configFileArg.getValue().c_str());
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
			updateParams(&args2);
		}

		// Arguments from the command line now overwrite any from the config file
		updateParams(&args1);

	} catch (TCLAP::ArgException &e)
	{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

	validateParams();
}

void setMPIParams(int size, int rank){
    mpiSize=size;
    mpiRank=rank;
}
//void updateParams(CArgs* a);
void updateParams(CArgs* a) {
	if (a->filenameArg.getValue().compare("") != 0) {
		filename = a->filenameArg.getValue();
	}

	if (a->outputFileArg.getValue().compare("") != 0) {
		outputFilename = a->outputFileArg.getValue();
	}

	if (a->nsArg.getValue() != -1) {
		noScenarios = a->nsArg.getValue();
	}

	if (a->ppArg.getValue() != -1) {
		penalty = a->ppArg.getValue();
	}

	if (a->pmultArg.getValue() > 0) {
		penaltyMult = a->pmultArg.getValue();
	}

	if (a->stepArg.getValue() != -1) {
		maxStep = a->stepArg.getValue();
	}

	if (a->maxSecondsArg.getValue() != -1) {
		maxSeconds = a->maxSecondsArg.getValue();
	}

	if (a->innerStepArg.getValue() != -1) {
		fixInnerStep = a->innerStepArg.getValue();
	}

	if (a->vertexHistoryArg.getValue() != -1) {
		UseVertexHistory = a->vertexHistoryArg.getValue();
	}

	if (a->threadsArg.getValue() != -1) {
		threads = a->threadsArg.getValue();
	}

	if (a->verboseSwitch.getValue() == true) {
		verbose = true;
	}

	if (a->debugSwitch.getValue() == true) {
		debug = true;
	}

	if ((a->linSwitch.getValue() == true) || (a->linFirstSwitch.getValue() == true)) {
		linRelaxFirst = true;
	}

	if ((a->linSwitch.getValue() == true) || (a->linSecondSwitch.getValue() == true)) {
		linRelaxSecond = true;
	}

	if (a->scalingSwitch.getValue() == true) {
		scaling = true;
	}

	if (a->CAP_Switch.getValue() == true)
	{
		if (filetype == 0) {
			filetype = 1;
		}
		else
		{
			std::cerr << "ERROR: Tried to set file type more than once" << std::endl;
			throw(-1);
		}
	}
	if (a->SMPS_Switch.getValue() == true)
	{
		if (filetype == 0) {
			filetype = 2;
		}
		else
		{
			std::cerr << "ERROR: Tried to set file type more than once" << std::endl;
			throw(-1);
		}
	}
	if (a->LB_Switch.getValue() == true) {
		LBcalc = true;
	}


	if (a->Heur_Switch.getValue() == true) {
		disableHeuristic = true;
	}
}
//String manipulation helper
int validateParams() {
	if (filename.compare("") == 0) {
		std::cerr << "ERROR: No problem file specified" << std::endl;
		throw(-1);
	}

	if (outputFilename.compare("") == 0) {
		outputFilename = "output";
	}

	if (noScenarios <= 0) {
		//std::cerr << "ERROR: Number of scenarios not specified" << endl;
		//throw(-1);
		//This now signifies we should use all scenarios in the file
		//Will only work for SMPS files, though.
	}

	if (penalty <= 0) {
		//std::cerr << "ERROR: Starting penalty parameter not specified" << endl;
		//throw(-1);
	}

	if (maxStep <= 0) {
		std::cerr << "WARNING: maxStep not specified, defaults to " << DEFAULT_MAX_OUTER_LOOP << std::endl;
		maxStep = DEFAULT_MAX_OUTER_LOOP;
	}

	if (filetype == 0) {
		std::cerr << "ERROR: File type not specified, use e.g. --SMPS" << std::endl;
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

};

#endif
