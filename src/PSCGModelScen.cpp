/**
*
*
*
*
*/

#include "PSCGModelScen.h"

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
void PSCGModelScen_Bodur::initialiseBodur(PSCGParams *par, ProblemDataBodur &pdBodur, int scenario){
	n1 = pdBodur.get_n1();
	n2 = pdBodur.get_n2();
	nS = pdBodur.get_nS();

	tS = scenario;
	nThreads = par->threads;

	c = new double[n1];
	d = new double[n2];
	for (int i = 0; i < n1; i++) {
		c[i] = (pdBodur.get_c())[i];
	}
	for (int j = 0; j < n2; j++) {
		d[j] = (pdBodur.get_d(tS))[j];
	}
		
	pr = pdBodur.getProbabilities()[tS];

	for(int i=0; i<n1; i++){ xVariables.add(IloNumVar(env,0,1,ILOINT));}
	for(int i=0; i<n2; i++){ yVariables.add(IloNumVar(env,0.0,IloInfinity));}
	
	
	for (int i = 0; i < pdBodur.get_A().getSize(); i++) {
		slpModel.add(IloScalProd(xVariables, (pdBodur.get_A())[i]) >= (pdBodur.get_b())[i]);
	}
	
	for (int i = 0; i < (pdBodur.get_T(tS)).getSize(); i++) {
		slpModel.add(IloScalProd(xVariables, pdBodur.get_T(tS)[i]) + IloScalProd(yVariables, (pdBodur.get_W(tS))[i]) >= (pdBodur.get_h(tS))[i]) ;
	}
	
	for (int i = 0; i < n1; i++) {
	    c_vec.add(c[i]);
	}
	for (int j = 0; j < n2; j++) {
	    d_vec.add(d[j]);
	}
	slpObjective.setSense(IloObjective::Minimize);
	slpObjective.setExpr( IloScalProd(xVariables, c_vec) + IloScalProd(yVariables,d_vec) ); 
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
int PSCGModelScen_SMPS::initialiseSMPS(PSCGParams *par, TssModel &smpsModel, int scenario) {
	tS = scenario;
	
	disableHeuristic = par->disableHeuristic;
	nThreads = par->threads;
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
	LagrMIPInterface_ = new OsiCpxSolverInterface();
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
	LagrMIPInterface_->assignProblem(mat, clbd, cubd, obj, rlbd, rubd);
	for (int jj = 0; jj< n1+n2; jj++){
	    if(ctype[jj]=='I' || ctype[jj]=='B'){LagrMIPInterface_->setInteger(jj);}
	    if(ctype[jj]=='B'){
		LagrMIPInterface_->setColBounds(jj,0.0,1.0);	
	    }
	}
	delete [] ctype; //All other allocated data passed to assignProblem is owned by LagrMIPInterface_

	//********Setting CPLEX Parameters*****************\\
	LagrMIPInterface_->setIntParam(OsiParallelMode,0);
	LagrMIPInterface_->setIntParam(OsiOutputControl,0);
	LagrMIPInterface_->setIntParam(OsiMIPOutputControl,0);
	if (nThreads >= 0) { LagrMIPInterface_->setIntParam(OsiParallelThreads, nThreads); }
	if (disableHeuristic) { LagrMIPInterface_->setIntParam(OsiHeurFreq, -1); }
	LagrMIPInterface_->setDblParam(OsiDualTolerance, 1e-6);
	LagrMIPInterface_->setHintParam(OsiDoPresolveInInitial,true);
	LagrMIPInterface_->setHintParam(OsiDoScale,true);
	LagrMIPInterface_->setHintParam(OsiDoCrash,true);
	LagrMIPInterface_->setHintParam(OsiDoReducePrint,true);
	LagrMIPInterface_->setHintParam(OsiLastHintParam, true);
	//***********End setting parameters*********************//

	return errcode;
}

//n1 and n2 need to be set
void PSCGModelScen::finishInitialisation() {
	if (n1==0 && n2==0){
	     cerr << "PSCGModelScen::finishInitialisation(): n1==0 and n2==0, returning...." << endl;
	     return;
	}
	x_vertex = new double[n1];
	y_vertex = new double[n2];
	x = new double[n1];
	y = new double[n2];
	//x_vertex.add(n1,0.0);
	//y_vertex.add(n2,0.0);
		
	for(int i=0; i<n1; i++) xVertices.push_back(vector<double>());// = new ptrIloNumArray[n1];
	for(int j=0; j<n2; j++) yVertices.push_back(vector<double>());// = new ptrDouble[n2];
		
		
	//add objective
	mpObjective.setSense( IloObjective::Minimize );
		
	mpModel.add(mpObjective);
	mpWeightConstraints.add(IloRange(env,1.0,mpWeight0,1.0));
	mpModel.add(mpWeightConstraints);
		
	for(int i=0; i<n1; i++) {
	    mpVertexConstraints.add(IloRange(env, 0.0, 0.0));
 	    mpAuxVariables.add(IloNumVar(mpVertexConstraints[i](-1.0)));
	    mpAuxVariables[i].setLB(-IloInfinity);
	}
	mpModel.add(mpVertexConstraints);

	//We use dual simplex due to the nature of the problem
	cplexMP.setParam(IloCplex::RootAlg, IloCplex::Dual);
	//cplexQP.setParam(IloCplex::RootAlg, IloCplex::Primal);
	
	if (nThreads >= 0) { cplexMP.setParam(IloCplex::Threads, nThreads); }
	cplexMP.setOut(env.getNullStream());
	cplexMP.extract(mpModel);
		
	//cout << "Finish Initialisation: Total memory use for " << tS << ": " << env.getTotalMemoryUsage() << endl;
	initialised = true;
}

#if 1
int PSCGModelScen_Bodur::solveLagrangianProblem(const double *omega) {
	if(omega!=NULL){
	  for (int i = 0; i < n1; i++) {
	        slpObjective.setLinearCoef(xVariables[i], c_vec[i]+omega[i]);
	  }
	}
	else{
	  for (int i = 0; i < n1; i++) {
	        slpObjective.setLinearCoef(xVariables[i], c_vec[i]);
	  }
	}
		
	if (!cplexMIP.solve()) {
		env.error() << "Failed to optimize in solveInitial" << endl;
		throw(-1);
	}

	
	return 0;
}
#endif
int PSCGModelScen_Bodur::solveFeasibilityProblem(){
	
	for (int i = 0; i < n1; i++) {
		//omega[i] = omega[i];
	        slpObjective.setLinearCoef(xVariables[i], 0.0);
	}
	for (int i = 0; i < n2; i++) {
		//omega[i] = omega[i];
	        slpObjective.setLinearCoef(yVariables[i], 0.0);
	}
		
	if (!cplexMIP.solve()) {
		env.error() << "Failed to optimize in solveInitial" << endl;
		throw(-1);
	}

	//for(int ii=0; ii<n1; ii++) cplexMIP.getValue(xVariables[ii], x_vertex[ii]);
	//for(int jj=0; jj<n2; jj++) cplexMIP.getValue(yVariables[jj], y_vertex[jj]);
	
	//LagrBd = cplexMIP.getObjValue();
	for (int i = 0; i < n2; i++) {
		//omega[i] = omega[i];
	        slpObjective.setLinearCoef(yVariables[i], d_vec[i]);
	}
	
	return 0;
}

// Given a scenario index and a dual variable, find the anticipative solution for first and second stage variables.
int PSCGModelScen_SMPS::solveLagrangianProblem(const double* omega) {
	
	OsiCpxSolverInterface* osi = LagrMIPInterface_;
	if(omega!=NULL){
	  for (int i = 0; i < n1; i++) {
		osi->setObjCoeff(i, c[i] + omega[i]);
	  }
	}
	else{
	  for (int i = 0; i < n1; i++) {
		osi->setObjCoeff(i, c[i]);
	  }
	}

	osi->branchAndBound();
	LagrBd = osi->getObjValue()*osi->getObjSense();
if(omega==NULL){ 
 cout << "LagrBd " << LagrBd << endl;
 cout << "Printing solution: " << endl;
 for(int ii=0; ii<n1+n2; ii++){
  cout << " " << osi->getColSolution()[ii];
 }
 cout << endl;
}
	
	setSolverStatus();

	
	
	return solverStatus_;
}

int PSCGModelScen_SMPS::solveFeasibilityProblem(){

	OsiCpxSolverInterface* osi = LagrMIPInterface_;
	
#if 1
	for (int i = 0; i < n1; i++) {
		osi->setObjCoeff(i, 0.0);
	}
	for (int i = n1; i < n1+n2; i++) {
		osi->setObjCoeff(i, 0.0);
	}
#endif

	osi->branchAndBound();
	
	setSolverStatus();

	if(solverStatus_==PSCG_OPTIMAL || solverStatus_==PSCG_ITER_LIM){	
		if(solverStatus_==PSCG_ITER_LIM) cerr << "Flagging: SMPS MIP solver indicated iteration limit reached." << endl;
		const double* solution = osi->getColSolution();
for(int ii=0; ii<n1; ii++){
cout << " (" << solution[ii] << ")";
}
cout << endl;
#if 0
		const double* solution = osi->getColSolution();
		LagrBd = osi->getObjValue()*osi->getObjSense();
	
		for (int i = 0; i < n1; i++) {
			x_vertex[i] = solution[i];
		}
	
		for (int j = 0; j < n2; j++) {
			y_vertex[j] = solution[n1+j];
		}
#endif
	}
	else{
		cerr << "solveFeasibilityProblem(): Flagging: SMPS MIP solver indicated isProvenOptimal() == false." << endl;
	}
#if 1
	for (int i = n1; i < n1+n2; i++) {
		osi->setObjCoeff(i, d[i-n1]);
	}
#endif
	
	return solverStatus_;
}


//Not tested!
void PSCGModelScen::updatePrimalVariables_OneScenario(const double *omega, const double *z, const double *scaling_vector) {
	
	double numerator = 0.0;
	double denominator = 0.0;
	
	for (int i = 0; i < n1; i++) {
		numerator -= (c[i] + omega[i] + scaling_vector[i] * (x[i] - z[i])) * (x_vertex[i] - x[i]);
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

void PSCGModelScen::updatePrimalVariablesHistory_OneScenario(const double *omega, const double *z) {

	IloNum weightObj0(0.0);
	for (int wI = 0; wI < nVertices; wI++) {
		weightObjective[wI] = baseWeightObj[wI];
	
		for (int i = 0; i < n1; i++) {
			weightObjective[wI] += xVertices[i][wI] * omega[i];
		}
	}
	
	for (int i = 0; i < n1; i++) {
		weightObj0 += x[i] * (c[i] + omega[i]);
	}

	for (int j = 0; j < n2; j++) {
		weightObj0 += y[j] * d[j];
	}

	mpObjective.setExpr(IloScalProd(mpWeightVariables, weightObjective) + weightObj0*mpWeight0 + quadraticTerm);
		
	//modify vertex constraint
	for (int i = 0; i < n1; i++) {
		mpVertexConstraints[i].setBounds(z[i], z[i]);
		mpVertexConstraints[i].setLinearCoef(mpWeight0, x[i]);
	}

	//cout << cplexQP.getModel() << endl;
	//cplexQP.exportModel("mymodel.lp");
	if (!cplexMP.solve()) {
		//cout << "Num vars: " << cplexQP.getNcols() << endl;
		cout << "MP Subproblem CPLEX status: " << cplexMP.getCplexStatus() << endl;
		env.error() << "Failed to optimize in update step" << endl;
		throw(-1);
	}

	cplexMP.getValues(weightSoln, mpWeightVariables);
	double weight0 = cplexMP.getValue(mpWeight0);
	

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


double PSCGModelScen::getDefaultPenaltyParameter() {

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

