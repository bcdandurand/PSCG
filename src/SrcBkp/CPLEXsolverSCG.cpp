/**
*
*
*
*
*/

#include "CPLEXsolverSCG.h"

#include <cmath>
#include <sstream>
#include <string>
#include <stdio.h>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include <time.h>
#include <sys/time.h>

#include "ProblemDataBodur.h"

#include "CoinPragma.hpp"
#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"

// Initialises the problem data, reading it in from a file.
void CPLEXsolverSCG::initialiseBodur(ProblemDataBodur &pdBodur, SMIP_fileRequest *request, int scenario){
	n1 = pdBodur.get_n1();
	n2 = pdBodur.get_n2();
	nS = pdBodur.get_nS();

	filetype = request->ftype;
	tS = scenario;
	nThreads = request->nThreads;

	c = new double[n1];
	d = new double[n2];
	for (int i = 0; i < n1; i++) {
		c[i] = (pdBodur.get_c())[i];
	}
	for (int j = 0; j < n2; j++) {
		d[j] = (pdBodur.get_d(tS))[j];
	}
		
		pr = 1.0/nS;

		for(int i=0; i<n1; i++){ xVariables.add(IloNumVar(env,0, 1,ILOINT));}
		for(int i=0; i<n2; i++){ yVariables.add(IloNumVar(env,0.0, IloInfinity));}
		
		
		for (int i = 0; i < pdBodur.get_A().getSize(); i++) {
			slpModel.add(IloScalProd(xVariables, (pdBodur.get_A())[i]) >= (pdBodur.get_b())[i]);
		}
		
		for (int i = 0; i < (pdBodur.get_T(tS)).getSize(); i++) {
			slpModel.add(IloScalProd(xVariables, pdBodur.get_T(tS)[i]) + IloScalProd(yVariables, (pdBodur.get_W(tS))[i]) >= (pdBodur.get_h(tS))[i]) ;
		}
		
		for (int i = 0; i < n1; i++) {
		    p_times_c.add(c[i]);
		}
		for (int j = 0; j < n2; j++) {
		    p_times_d.add(d[j]);
		}
		slpObjective.setSense(IloObjective::Minimize);
		//IloNumExprArg origObj(IloScalProd(xVariables, p_times_c) + IloScalProd(yVariables,p_times_d)); 
		//slpObjective.setExpr( origObj );
		slpObjective.setExpr( IloScalProd(xVariables, p_times_c) + IloScalProd(yVariables,p_times_d) ); 
		slpModel.add(slpObjective);

		if (nThreads >= 0) { cplexMIP.setParam(IloCplex::Threads, nThreads); }
		cplexMIP.setParam(IloCplex::EpGap, 1e-6);
		cplexMIP.setOut(env.getNullStream());
		//cplexMIP.setParam();
		cplexMIP.extract(slpModel);
		//if(tS==0){cout << slpModel << endl;}

}

//Initialise SIPLIB
//int CPLEXsolverSCG::initialiseSMPS(SMIP_fileRequest *request, int scenario) {
int CPLEXsolverSCG::initialiseSMPS(SMIP_fileRequest *request, TssModel &smpsModel, int scenario) {
	tS = scenario;
	
	filetype = request->ftype;
	
	disableHeuristic = request->disableHeuristic;
	nThreads = request->nThreads;
	//ProblemDataSMPS *problemSMPS = ProblemDataSMPS::getSingleton();
	n1 = smpsModel.getNumCols(0);
	n2 = smpsModel.getNumCols(1);
	nS = smpsModel.getNumScenarios();

	c = new double[n1];
	d = new double[n2];
	CoinPackedMatrix * mat = new CoinPackedMatrix();
	double *clbd,*cubd,*obj,*rlbd,*rubd;
	char *ctype;
	
	int errcode = smpsModel.decompose(1, &tS, 0, NULL, NULL, NULL, 
		mat, clbd, cubd, ctype, obj, rlbd, rubd );
	LagrMIPInterface = new OsiCpxSolverInterface();
	pr = smpsModel.getProbability()[tS];
	for (int i = 0; i < n1; i++) {
		c[i] = obj[i];
	}
	for (int j = 0; j < n2; j++) {
		obj[j+n1] /= pr;  //More generally, divide by sum of probabilities. This is a hack; 
				//the TssModel::decompose may need to be modified
				//in the dual decomp. case to avoid the need for this.
		d[j] = obj[j+n1];
	}
	LagrMIPInterface->assignProblem(mat, clbd, cubd, obj, rlbd, rubd);
	for (int jj = 0; jj< n1+n2; jj++){
	    if(ctype[jj]=='I' || ctype[jj]=='B'){LagrMIPInterface->setInteger(jj);}
	    if(ctype[jj]=='B'){
		LagrMIPInterface->setColBounds(jj,0.0,1.0);	
	    }
	}
	delete [] ctype; //All other allocated data passed to assignProblem is owned by LagrMIPInterface

	//********Setting CPLEX Parameters*****************\\
	LagrMIPInterface->setIntParam(OsiParallelMode,0);
	LagrMIPInterface->setIntParam(OsiOutputControl,0);
	LagrMIPInterface->setIntParam(OsiMIPOutputControl,0);
	if (nThreads >= 0) { LagrMIPInterface->setIntParam(OsiParallelThreads, nThreads); }
	if (disableHeuristic) { LagrMIPInterface->setIntParam(OsiHeurFreq, -1); }
	LagrMIPInterface->setDblParam(OsiDualTolerance, 1e-6);
	LagrMIPInterface->setHintParam(OsiDoPresolveInInitial,true);
	LagrMIPInterface->setHintParam(OsiDoScale,true);
	LagrMIPInterface->setHintParam(OsiDoCrash,true);
	LagrMIPInterface->setHintParam(OsiDoReducePrint,true);
	LagrMIPInterface->setHintParam(OsiLastHintParam, true);
	//***********End setting parameters*********************//

	return errcode;
}

void CPLEXsolverSCG::finishInitialisation() {
	
	x = new double[n1];
	y = new double[n2];

	x_vertex.add(n1,0.0);
	y_vertex.add(n2,0.0);
		
	for(int i=0; i<n1; i++) xVertices.push_back(IloNumArray(env));// = new ptrIloNumArray[n1];
	for(int j=0; j<n2; j++) yVertices.push_back(vector<double>());// = new ptrDouble[n2];
		
		
	//add objective
	upvhObjective.setSense( IloObjective::Minimize );
		
	upvhModel.add(upvhObjective);
	upvhWeightConstraints.add(IloRange(env,1.0,upvhWeight0,1.0));
	upvhModel.add(upvhWeightConstraints);
		
	for(int i=0; i<n1; i++) {
	    upvhVertexConstraints.add(IloRange(env, 0.0, 0.0));
 	    upvhDummyVariables.add(IloNumVar(upvhVertexConstraints[i](-1.0)));
	    upvhDummyVariables[i].setLB(-IloInfinity);
	}
	upvhModel.add(upvhVertexConstraints);

	//We use dual simplex due to the nature of the problem
	cplexQP.setParam(IloCplex::RootAlg, IloCplex::Dual);
	//cplexQP.setParam(IloCplex::RootAlg, IloCplex::Primal);
	
	if (nThreads >= 0) { cplexQP.setParam(IloCplex::Threads, nThreads); }
	cplexQP.setOut(env.getNullStream());
	cplexQP.extract(upvhModel);
		
	cout << "Finish Initialisation: Total memory use for " << tS << ": " << env.getTotalMemoryUsage() << endl;
	initialised = true;
}

// Given a scenario index and a dual variable, find the anticipative solution for first and second stage variables.
int CPLEXsolverSCG::solveLagrangianProblem(const double *dual_var) {
	switch( filetype ) 
	{
	case 1:
		return solveLagrangianProblem_Bodur(dual_var);
		break;
	case 2:
		return solveLagrangianProblem_SMPS(dual_var);
		break;
	default:
		throw(-1);
		break;
	}
}

int CPLEXsolverSCG::solveLagrangianProblem_Bodur(const double *dual_var) {

	for (int i = 0; i < n1; i++) {
		//omega[i] = dual_var[i];
	        slpObjective.setLinearCoef(xVariables[i], p_times_c[i]+dual_var[i]);
	}
		
	if (!cplexMIP.solve()) {
		env.error() << "Failed to optimize in solveInitial" << endl;
		throw(-1);
	}

	cplexMIP.getValues(x_vertex, xVariables);
	cplexMIP.getValues(y_vertex, yVariables);
	
	LagrBd = cplexMIP.getObjValue();
	
	return 0;
}

int CPLEXsolverSCG::solveLagrangianProblem_SMPS(const double* dual_var) {
	
	OsiCpxSolverInterface* interface = LagrMIPInterface;
	
	for (int i = 0; i < n1; i++) {
		interface->setObjCoeff(i, c[i] + dual_var[i]);
	}

	interface->branchAndBound();
	
	const double* solution = interface->getColSolution();
	
	LagrBd = interface->getObjValue();
	
	for (int i = 0; i < n1; i++) {
		x_vertex[i] = solution[i];
	}
	
	for (int j = 0; j < n2; j++) {
		y_vertex[j] = solution[n1+j];
	}
	
	if (interface->isProvenOptimal() == false) {
		cerr << "Flagging: SMPS MIP solver indicated isProvenOptimal() == false." << endl;
	}
	
	return 0;
}


//Not tested!
void CPLEXsolverSCG::updatePrimalVariables_OneScenario(const double *dvar, const double *z, const double *scaling_vector) {
	
	double numerator = 0;
	double denominator = 0;
	
	for (int i = 0; i < n1; i++) {
		numerator -= (c[i] + dvar[i] + scaling_vector[i] * (x[i] - z[i])) * (x_vertex[i] - x[i]);
		denominator += (x_vertex[i] - x[i]) * scaling_vector[i] * (x_vertex[i] - x[i]);
	}
	
	for (int j = 0; j < n2; j++) {
		numerator -= d[j] * (y_vertex[j] - y[j]);
	}
	
	double a;
	if (denominator > 1e-9)	{
		a = numerator / denominator;
	}
	else {
		a = 1;
	}
	
	if (a > 1) {a = 1;}
	if (a < 0) {a = 0;}
	
	for (int i = 0; i < n1; i++) {
		x[i] = x[i] + a * (x_vertex[i] - x[i]);
	}

	for (int j = 0; j < n2; j++) {
		y[j] = y[j] + a * (y_vertex[j] - y[j]);
	}
}

void CPLEXsolverSCG::updatePrimalVariablesHistory_OneScenario(const double *dvar, const double *z) {

	IloNum weightObj0(0.0);
	for (int wI = 0; wI < nVertices; wI++) {
		weightObjective[wI] = baseWeightObj[wI];
	
		for (int i = 0; i < n1; i++) {
			weightObjective[wI] += xVertices[i][wI] * dvar[i];
		}
	}
	
	for (int i = 0; i < n1; i++) {
		weightObj0 += x[i] * (c[i] + dvar[i]);
	}

	for (int j = 0; j < n2; j++) {
		weightObj0 += y[j] * d[j];
	}

	upvhObjective.setExpr(IloScalProd(upvhWeightVariables, weightObjective) + weightObj0*upvhWeight0 + quadraticTerm);
		
	//modify vertex constraint
	for (int i = 0; i < n1; i++) {
		upvhVertexConstraints[i].setBounds(z[i], z[i]);
		upvhVertexConstraints[i].setLinearCoef(upvhWeight0, x[i]);
	}

	//cout << cplexQP.getModel() << endl;
	//cplexQP.exportModel("mymodel.lp");
	if (!cplexQP.solve()) {
		cout << "CPLEX status: " << cplexQP.getCplexStatus() << endl;
		//cout << "Num vars: " << cplexQP.getNcols() << endl;
		env.error() << "Failed to optimize in update step" << endl;
		throw(-1);
	}

	cplexQP.getValues(weightSoln, upvhWeightVariables);
	double weight0 = cplexQP.getValue(upvhWeight0);
	

	// note: the final weight corresponds to the existing x
	for (int i = 0; i < n1; i++) {
		x[i] = weight0 * x[i];
	}

	for (int j = 0; j < n2; j++) {
		y[j] = weight0 * y[j];
	}
			
	for(int wI=0; wI<nVertices; wI++) {
		
		for (int i = 0; i < n1; i++) {
			x[i] += weightSoln[wI] * xVertices[i][wI];
		}
	}
		
	for(int wI=0; wI<nVertices; wI++) {
		
		for (int j = 0; j < n2; j++) {
			y[j] += weightSoln[wI] * yVertices[j][wI];
		}
	}

	
}

//Not tested, not updated!
void CPLEXsolverSCG::getLagrangianGradient(SMIP_qu_getLagrangianGradient* question, SMIP_ans_getLagrangianGradient* answer) {
	
	//int tS = question->thisScenario;
	
	for (int i = 0; i < n1; i++)
	{
		//answer->x_gradient[i] = pr * c[i] + question->dvar[i];
		answer->x_gradient[i] = c[i] + question->dvar[i];
	}
	for (int j = 0; j < n2; j++)
	{
		//answer->y_gradient[j] = pr * d[j];
		answer->y_gradient[j] = d[j];
	}
	
}


double CPLEXsolverSCG::getDefaultPenaltyParameter() {

	double maxC = abs(c[0]);
	double meanC = abs(c[0]) / n1;
	
	for (int i = 1; i < n1; i++)
	{
		if (maxC < abs(c[i])) { maxC = abs(c[i]); }
		meanC = meanC + abs(c[i]) / n1;
	}
	
	double output = maxC * 1e-2;
	if (output < 0.5*meanC) {output = 0.5*meanC; }
	
	delete[] c;
	
	return output;
	
}

