/*
Author: Brian C. Dandurand (2017-2018)
Based on and modified from code developed by Jeffrey Christiansen, Brian Dandurand, and Fabricio Oliveira
at RMIT in Melbourne Australia with funding under project ARC DP 140100985 during 2015-2017.
CIs of that projects were Prof. Andrew Eberhard, Prof. Natashia Boland, and PI Prof. Jeffrey Linderoth.
*/

/**
*
*
*
*
*/

#include "PSCGScen.h"

#include <cmath>
#include <sstream>
#include <string>
#include <stdio.h>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include <time.h>
#include <sys/time.h>


#include "CoinPragma.hpp"
#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"
//#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"

//Initialise SIPLIB
//int CPLEXsolverSCG::initialiseSMPS(SMIP_fileRequest *request, int scenario) {
int PSCGScen_SMPS::initialiseSMPS(DecTssModel &smpsModel, int scenario) {
	tS = scenario;
	
	//ProblemDataSMPS *problemSMPS = ProblemDataSMPS::getSingleton();
	n1 = smpsModel.getNumCols(0);
	n2 = smpsModel.getNumCols(1);
	nS = smpsModel.getNumScenarios();
	//nThreads = par->threads;
	origLBs_ = new double[n1+n2];
	origUBs_ = new double[n1+n2];
	aggrXY0 = new double[n1+n2];
	aggrXY1 = new double[n1+n2];

	c = new double[n1+n2];
	d = c+n1;
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
				//the DecTssModel::decompose may need to be modified
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
        for(int iii=0; iii<n1+n2; iii++){
	    origLBs_[iii]=LagrMIPInterface_->getColLower()[iii];
	    origUBs_[iii]=LagrMIPInterface_->getColUpper()[iii];
	}

	setCPXMIPParameters();
#if 0
	//********Setting CPLEX Parameters*****************\\
	LagrMIPInterface_->setIntParam(OsiParallelMode,0);
	LagrMIPInterface_->setIntParam(OsiOutputControl,0);
	LagrMIPInterface_->setIntParam(OsiMIPOutputControl,0);
	if (nThreads >= 0) { LagrMIPInterface_->setIntParam(OsiParallelThreads, nThreads); }
	LagrMIPInterface_->setDblParam(OsiDualTolerance, 1e-9);
	setGapTolerances(1e-9,1e-9);
	LagrMIPInterface_->setHintParam(OsiDoPresolveInInitial,true);
	LagrMIPInterface_->setHintParam(OsiDoScale,true);
	LagrMIPInterface_->setHintParam(OsiDoCrash,true);
	LagrMIPInterface_->setHintParam(OsiDoReducePrint,true);
	LagrMIPInterface_->setHintParam(OsiLastHintParam, true);
	//***********End setting parameters*********************//
#endif

	return errcode;
}

//n1 and n2 need to be set
void PSCGScen::finishInitialisation() {
	if (n1==0 && n2==0){
	     cerr << "PSCGScen::finishInitialisation(): n1==0 and n2==0, returning...." << endl;
	     return;
	}
	x_vertex = new double[n1+n2];
	//y_vertex = new double[n2];
	y_vertex = x_vertex+n1;
	x = new double[n1+n2];
	dispersions = new double[n1+n2];
	intDiscVec_ = new double[n1+n2];
	for(int ii=0; ii<n1+n2; ii++){ dispersions[ii]=0.0;}
	//y = new double[n2];
	y = x+n1;
	//x_vertex.add(n1,0.0);
	//y_vertex.add(n2,0.0);
		
	//for(int i=0; i<n1; i++) xVertices.push_back(vector<double>());// = new ptrIloNumArray[n1];
	//for(int j=0; j<n2; j++) yVertices.push_back(vector<double>());// = new ptrDouble[n2];
		
		
	//add objective
	//
#if 0
	try{

	mpObjective.setSense( IloObjective::Minimize );
	//mpObjective.setExpr(quadraticTerm);
		
	mpModel.add(mpObjective);
	mpWeightConstraints.add(IloRange(env,1.0,mpWeight0,1.0));
	mpModel.add(mpWeightConstraints);
		
	for(int i=0; i<n1; i++) {
	    mpVertexConstraints.add(IloRange(env, 0.0, 0.0));
 	    mpAuxVariables.add(IloNumVar(mpVertexConstraints[i](-1.0)));
	    mpAuxVariables[i].setLB(-1.0*IloInfinity);
	    mpAuxVariables[i].setUB(IloInfinity);
	}
	mpModel.add(mpVertexConstraints);

	cplexMP.setParam(IloCplex::RootAlg, IloCplex::Primal);
	//cplexMP.setParam(IloCplex::RootAlg, CPX_ALG_BARRIER);
        //cplexMP.setParam(IloCplex::Param::SolutionType, CPX_BASIC_SOLN);
	//cplexMP.setParam(IloCplex::RootAlg, 0);
	cplexMP.setParam(IloCplex::OptimalityTarget, 1);
	cplexMP.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-9);
        cplexMP.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
	
	if (nThreads >= 0) { cplexMP.setParam(IloCplex::Threads, nThreads); }
	cplexMP.setOut(env.getNullStream());
	cplexMP.setWarning(env.getNullStream());
	  cplexMP.extract(mpModel);
	}
   	catch(IloException& e){
		cout << "cplexMP.extract error: " << e.getMessage() << endl;
		cout << "Exception caught!" << endl;
	//refresh();
		e.end();
   	}
#endif
		
	//cout << "Finish Initialisation: Total memory use for " << tS << ": " << env.getTotalMemoryUsage() << endl;
	initialised = true;
}

#if 0
int PSCGScen_SMPS::initialLPSolve(const double* omega) {
	
	OsiCpxSolverInterface* osi = LagrMIPInterface_;
        //changeFromMIQPToMILP(); //This only does anything if the problem has not already been changed back
	//changeToMILP();
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
	osi->initialSolve();
	setSolverStatus();
	return solverStatus_;
}
#endif
	//setSolverStatus();
// Given a scenario index and a dual variable, find the anticipative solution for first and second stage variables.
int PSCGScen_SMPS::solveLagrangianProblem(const double* omega, bool doInitialSolve) {
	
	OsiCpxSolverInterface* osi = LagrMIPInterface_;
	double bestUB;

        //changeFromMIQPToMILP(); //This only does anything if the problem has not already been changed back
	//changeToMILP();
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
#if 0
	double startLagrBd = COIN_DBL_MAX;
	double *startSol=new double[n1+n2];
	if(nVertices>0){
	    //startSol = new double[n1+n2];
	    memcpy(startSol,x_vertex,n1*sizeof(double));
	    memcpy(startSol+n1,y_vertex,n2*sizeof(double));
	    osi->setColSolution(startSol);
	    startLagrBd = evaluateVertexSolution(omega);
	}
#endif

	if(doInitialSolve){
	    osi->initialSolve();
	}
	else{osi->resolve();}

#if 0
	setSolverStatus();
        if(solverStatus_==PSCG_PRIMAL_INF){ 
	    LagrBd =  COIN_DBL_MAX*osi->getObjSense();
	    cout << "Infeasible at initial LP solve for scenario " << tS << endl;
	}
    else
#endif
	{

	osi->branchAndBound();
#if 0
if(omega==NULL){ 
 cout << "LagrBd " << LagrBd << endl;
 cout << "Printing solution: " << endl;
 for(int ii=0; ii<n1+n2; ii++){
  cout << " " << osi->getColSolution()[ii];
 }
 cout << endl;
}
#endif
	//printLagrSoln();
	setSolverStatus();
        //if(tS==0) printXBounds();
        //if(tS==0) printYBounds();
        if(solverStatus_==PSCG_PRIMAL_INF){ 
	    LagrBd =  COIN_DBL_MAX*osi->getObjSense();
	    //printOrigBDs();
	    //printColBds();
	}
	else{
	    //osi->markHotStart();
	    //LagrBd = computeMIPVal(omega);

	    //LagrBd = osi->getObjValue()*osi->getObjSense();
	    bestUB = osi->getObjValue()*osi->getObjSense();
	    LagrBd = getMIPBestNodeVal();
	    if( bestUB - getMIPBestNodeVal() > 1e-6 ){
		    cout << "Gap in solver was not closed adequately: " << bestUB << " vs " 

			<<  getMIPBestNodeVal() << endl;
		    cout << "Solver status: " << getCPLEXErrorStatus() << endl;
	    }
	    updateSolnInfo();
	}
#if 0
	if(nVertices>0 && solverStatus_!=PSCG_PRIMAL_INF){
	  if(startLagrBd + 1e-5 < LagrBd && omega!=NULL){
	    cout << "Flagging " << LagrBd - startLagrBd << endl;
	    //assert(startLagrBd >= LagrBd);
	    //osi->setColSolution(startSol);
	    //LagrBd = startLagrBd;
	  }
	}
#endif
	
	//delete [] startSol;
    }//else initial LP solve was feasible
	
	return solverStatus_;
}
#if 0
int PSCGScen_SMPS::solveAugmentedLagrangianMIP(const double* omega, const double* z, const double rho, const double* scal){
	OsiCpxSolverInterface* osi = LagrMIPInterface_;
        changeFromMILPToMIQP(); //This only does anything if the problem has not already been changed back
	double *diag = new double[n1+n2];
	if(omega!=NULL){
	  for (int i = 0; i < n1; i++) {
		osi->setObjCoeff(i, c[i] + omega[i] - z[i]*rho*scal[i]);
		diag[i]=rho*scal[i];
	  }
	}
	else{
	  for (int i = 0; i < n1; i++) {
		osi->setObjCoeff(i, c[i] - z[i]*rho*scal[i]);
		diag[i]=rho*scal[i];
	  }
	}
	for (int j=0; j<n2; j++){
	    diag[n1+j]=0.0;
	}
	
	osi->setSeparableQuadraticObjectiveCoefficients(diag);
#if 1
    	osi->setIntParam(OsiSolutionTarget, 1);
     
    	osi->setIntParam(OsiParallelMode,1);
        osi->setIntParam(OsiOutputControl,0);
     	osi->setIntParam(OsiMIPOutputControl,0);
    	osi->setDblParam(OsiDualTolerance, 1e-3);
     
    	osi->setHintParam(OsiDoPresolveInInitial,true);
    	osi->setHintParam(OsiDoScale,true);
    	osi->setHintParam(OsiDoCrash,true);
    	osi->setHintParam(OsiDoReducePrint,true);
    	osi->setHintParam(OsiLastHintParam, true);
#endif
	osi->branchAndBound();
	setSolverStatus();
        if(solverStatus_==PSCG_PRIMAL_INF){ 
	    objVal =  COIN_DBL_MAX;
	    //osi->unmarkHotStart();
	}
	else{
	    memcpy(x_vertex,osi->getColSolution(),n1*sizeof(double));
	    memcpy(y_vertex,osi->getColSolution()+n1,n2*sizeof(double));
	    objVal = osi->getObjValue();
	    for(int ii=0; ii<n1; ii++){
		objVal += 0.5*rho*scal[ii]*z[ii]*z[ii];
	    }
	}
	for (int i = 0; i < n1; i++) {
	  osi->setObjCoeff(i, c[i]);
	  diag[i]=0.0;
	}
	osi->setSeparableQuadraticObjectiveCoefficients(diag);
	osi->setSolvingMIQP(false);
   	delete [] diag; 
}
#endif

int PSCGScen_SMPS::solveFeasibilityProblem(){

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

bool PSCGScen::addVertex(){
    if(checkWhetherVertexIsRedundant()){
	return false;
    }
    if(nVertices==maxNVertices){
	if(vecWeights.size()==0){vecWeights.push_back(1.0);}
	else{vecWeights.push_back(0.0);}
	xyVertices.push_back(vector<double>(n1+n2));
	//yVertices.push_back(vector<double>(n2));
	nVertices++;
	maxNVertices++;
	for(int i=0; i<n1+n2; i++) {
	   xyVertices[nVertices-1][i] = x_vertex[i];
	}
#if 0
	for(int i=0; i<n2; i++) { 
	   yVertices[nVertices-1][i] = y_vertex[i];
	}
#endif
	bestVertexIndex=nVertices-1;
#if 0
try{
	mpWeightVariables.add(IloNumVar(mpWeightConstraints[0](1.0)));
	//mpWeightVariables.add(IloNumVar());
	//mpWeightVariables[nVertices-1].setLB(0.0);
	
	for(int i=0; i<n1; i++) {
	    mpVertexConstraints[i].setLinearCoef(mpWeightVariables[nVertices-1], x_vertex[i]);
	}
	//mpWeightConstraints[0].setLinearCoef(mpWeightVariables[nVertices-1],1.0);

	//baseWeightObj.push_back(0.0);
	weightSoln.add(0.0);
	weightObjective.add(0.0);
}
catch(IloException &e){
    cout << "addVertex error: " << e.getMessage() << endl;
    cout << "Occurred while updating vertex constraint info" << endl;	
    exit(1);
}
#endif

#if 0
	for (int i = 0; i < n1; i++) {
		baseWeightObj[nVertices-1] += x_vertex[i] * c[i];
	}

	for (int j = 0; j < n2; j++) {
		baseWeightObj[nVertices-1] += y_vertex[j] * d[j];
	}
#endif

    }
    else{ //nVertices<maxNVertices
	assert(nVertices < maxNVertices);
	replaceVertexAtIndex(nVertices);
	nVertices++;
    }
    if(oldestVertexIndex==-1) oldestVertexIndex=0;
    return true;
}
bool PSCGScen::replaceVertexAtIndex(int iii){
    if(checkWhetherVertexIsRedundant()){
	return false;
    }

	//mpWeightVariables.add(IloNumVar(mpWeightConstraints[0](1.0)));
	//mpWeightVariables[nVertices].setBounds(0.0,1.0);
	assert(iii>=0);
	assert(iii<maxNVertices);
#if 0
try{
	mpWeightVariables[iii].setBounds(0.0,IloInfinity);
}
catch(IloException &e){
    cout << "replaceVertexAtIndex error: " << e.getMessage() << endl;
    cout << "Occurred while updating weight var bounds." << endl;	
    exit(1);
}
#endif
	for(int i=0; i<n1+n2; i++) {
	   xyVertices[iii][i] = x_vertex[i];
	}
#if 0
	for(int i=0; i<n2; i++) { 
	   yVertices[iii][i] = y_vertex[i];
	}
#endif
#if 0
try{
	for(int i=0; i<n1; i++) {
	    mpVertexConstraints[i].setLinearCoef(mpWeightVariables[iii], x_vertex[i]);
	}
}
catch(IloException &e){
    cout << "replaceVertexAtIndex error: " << e.getMessage() << endl;
    cout << "Occurred while updating vertex constraint column coefficients." << endl;	
    exit(1);
}

	weightObjective[iii]=0.0;//.add(0.0);
#endif

	bestVertexIndex=iii;
	vecWeights[iii]=0.0;
	return true;
}

//Not tested!
void PSCGScen::solveMPLineSearch(const double *omega, const double *z, const double rho, const double *scaling_vector, int vertexIndex, double *z_average) {
	
	double numerator = 0.0;
	double denominator = 0.0;
	double dir = 0.0;
	double *vertX, *vertY;
	if(vertexIndex==-1){
	    vertexIndex=bestVertexIndex;
	}
	vertX = &(xyVertices[vertexIndex][0]);
	//vertY = &(yVertices[vertexIndex][0]);
	vertY = vertX+n1;
	
	for (int i = 0; i < n1; i++) {
		if(omega==NULL){numerator -= (c[i]/rho + scaling_vector[i] * (x[i] - z[i])) * (vertX[i] - x[i]);}
		else{numerator -= ( (c[i] + omega[i])/rho + scaling_vector[i] * (x[i] - z[i])) * (vertX[i] - x[i]);}
		denominator += (vertX[i] - x[i]) * scaling_vector[i] * (vertX[i] - x[i]);
	}
	
	for (int j = 0; j < n2; j++) {
		numerator -= d[j] * (vertY[j] - y[j]) / rho;
	}
	
	double a;
	if (denominator > 1e-20)	{
		a = numerator / denominator;
	}
	else {
		a = 1;
	}
	
	if (a > 1) {a = 1;}
	if (a < 0) {a = 0;}
	double oldX;
	for (int i = 0; i < n1; i++) {
		dir = vertX[i] - x[i];
		oldX=x[i];
		x[i] = x[i] + a * (vertX[i] - x[i]);
		//if(updateDisp){dispersions[i] = max( (1.0-a)*(dispersions[i]+fabs(x[i]-oldX)) , a*fabs(x[i]-vertX[i]) );}
	}

	for (int j = 0; j < n2; j++) {
		y[j] = y[j] + a * (vertY[j] - y[j]);
	}

	for(int wI=0; wI<nVertices; wI++) {
	    vecWeights[wI] = (1.0-a)*vecWeights[wI];
	}
	vecWeights[vertexIndex]+=a;
}

void PSCGScen::solveMPVertices(const double *omega, const double *z, const double rho, const double *scaling_vector)
{
    //optimiseLagrOverVertexHistory(omega); //prepares next call of solveMPLineSearch(omega,z,rho,scaling_vector)
    solveMPLineSearch(omega,z,rho,scaling_vector);
#if 0
    for(int nn=0; nn<40; nn++){
      solveMPLineSearch(omega,z,rho,scaling_vector);
      for(int vv=0; vv<nVertices; vv++){
	solveMPLineSearch(omega,z,rho,scaling_vector,vv);
      }
      optimiseLagrOverVertexHistory(omega); //prepares next call of solveMPLineSearch(omega,z,rho,scaling_vector)
    }
#endif

}

#if 0
void PSCGScen::computeWeightsForCurrentSoln(const double *z) {
//double *oldDispersions = new double[n1];
//double absDiffSum=0.0;
//for(int ii=0; ii<n1; ii++) {oldDispersions[ii]=dispersions[ii];}
   if(z==NULL){z=x;}
   mpWeight0.setBounds(0.0,0.0);
   try{
	for (int wI = 0; wI < nVertices; wI++) {
		weightObjective[wI] = 1.0;
	}
	mpObjective.setExpr(IloScalProd(mpWeightVariables, weightObjective) + quadraticTerm);
	for (int i = 0; i < n1; i++) {
		//mpVertexConstraints[i].setBounds(x[i], x[i]);
		mpVertexConstraints[i].setBounds(z[i], z[i]);
	}
	if (!cplexMP.solve()) {
		cout << "MP Subproblem CPLEX status: " << cplexMP.getCplexStatus() << endl;
		env.error() << "Failed to optimize in update step" << endl;
		cout << "Number of vertices: " << getNVertices() << endl;
		cplexMP.exportModel("infeasModel.lp");
		//cout << "Refreshing solution..." << endl;
		//refresh(omega,z,scaling_vector);
	}
	if(cplexMP.getCplexStatus()!=IloCplex::Optimal) {
		cout << "MP Subproblem CPLEX status: " << cplexMP.getCplexStatus() << endl;
	}
	cplexMP.getValues(weightSoln, mpWeightVariables);
	polishWeights(); //fix any unfortunate numerical quirks
	for (int ii = 0; ii < n1; ii++) {
		dispersions[ii] = 0.0;
		for(int wI=0; wI<nVertices; wI++) {
		    dispersions[ii] += weightSoln[wI]*fabs(xVertices[wI][ii]-z[ii]) ; 
		}
	}
   }
   catch(IloException& e){
	//cout << "MPSolve error: " << e.getMessage() << endl;
	cout << "Exception caught...weights not computed accurately!" << endl;
	//refresh();
	e.end();
   }
   mpWeight0.setBounds(0.0,1.0);
//for(int ii=0; ii<n1; ii++) {absDiffSum+= fabs(oldDispersions[ii]-dispersions[ii]);}
//cout << "Difference between old and new dispersions " << absDiffSum << endl;
//delete [] oldDispersions;
}
#endif

void PSCGScen::solveMPHistory(const double *omega, const double *z, const double *zLBs, const double *zUBs, 
	const double rho, const double *scaling_vector, bool updateDisp) {

   IloEnv   env;
   try {
      IloModel model(env);
      IloNumVar a0(env);
      IloNumVarArray a(env);
      IloNumVarArray zeta(env);
      IloRangeArray con(env);
      IloNumArray linCoeff(env);
      IloNumArray ones(env);
      IloExpr objExpr(env);

      //Populate the model.
      IloNum weightObj0(0.0);
      if(omega!=NULL){	
	 for (int i = 0; i < n1; i++) {
	   weightObj0 += x[i] * (c[i] + omega[i]);
	 }
      }
      else{	
	 for (int i = 0; i < n1; i++) {
	   weightObj0 += x[i] * (c[i]);
         }
      }

      for (int j = 0; j < n2; j++) {
	weightObj0 += y[j] * d[j];
      }
      weightObj0 /= rho;
      objExpr += weightObj0*a0;
      for (int wI = 0; wI < nVertices; wI++) {
		//weightObjective[wI] = baseWeightObj[wI]/rho;
		linCoeff.add(0.0);
		ones.add(1.0);
		a.add(IloNumVar(env));
		if(omega!=NULL){	
		    for (int i = 0; i < n1; i++) {
			linCoeff[wI] += xyVertices[wI][i] * (c[i]+omega[i]);
		    }
		}
		else{
		    for (int i = 0; i < n1; i++) {
			linCoeff[wI] += xyVertices[wI][i] * (c[i]);
		    }
		}
		for (int j = n1; j < n1+n2; j++) {
		    linCoeff[wI] += xyVertices[wI][j] * c[j];
		}
		linCoeff[wI] /= rho;

		
      	objExpr += linCoeff[wI]*a[wI];
      }
   for (int ii = 0; ii < n1; ++ii) {
	 zeta.add(IloNumVar(env,-IloInfinity,IloInfinity));
         objExpr += 0.5*scaling_vector[ii]*zeta[ii]*zeta[ii];
   }
   IloObjective obj = IloMinimize(env, objExpr);
   model.add(obj);
   objExpr.end();
   con.add(1.0*a0 + IloScalProd(ones,a) == 1.0);
   for(int ii=0; ii<n1; ii++){
	IloExpr vertConstrExpr(env);
	vertConstrExpr += x[ii]*a0;
	for(int wI=0; wI<nVertices; wI++){
	    vertConstrExpr += xyVertices[wI][ii] * a[wI];
	}	
	vertConstrExpr += -1.0*zeta[ii];
	con.add(vertConstrExpr == z[ii]);
	vertConstrExpr.end();
   }
   model.add(con);

   
      IloCplex cplex(model);
	cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
	//cplexMP.setParam(IloCplex::RootAlg, CPX_ALG_BARRIER);
        //cplexMP.setParam(IloCplex::Param::SolutionType, CPX_BASIC_SOLN);
	//cplexMP.setParam(IloCplex::RootAlg, 0);
	cplex.setParam(IloCplex::OptimalityTarget, 1);
	cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-9);
        cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	
	cplex.setParam(IloCplex::Threads, nThreads); 

      // Optimize the problem and obtain solution.
      if ( !cplex.solve() ) {
         env.error() << "Failed to optimize" << endl;
         throw(-1);
      }
      IloNum weight0 = cplex.getValue(a0);
	for (int i = 0; i < n1+n2; i++) {
		x[i] = weight0 * x[i];
		//if(updateDisp){dispersions[i] = weight0*dispersions[i];}
	}

#if 0
	for (int j = 0; j < n2; j++) {
		y[j] = weight0 * y[j];
	}
#endif
	for(int wI=0; wI<nVertices; wI++) {
		vecWeights[wI] = weight0*vecWeights[wI] + cplex.getValue(a[wI]);
		for (int i = 0; i < n1+n2; i++) {
			x[i] += cplex.getValue(a[wI]) * xyVertices[wI][i];
		}
	#if 0	
		for (int j = 0; j < n2; j++) {
			y[j] += cplex.getValue(a[wI]) * yVertices[wI][j];
		}
	#endif
	}
//}

   }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }

   env.end();
}

void PSCGScen::solveMPHistory2Norm(const double *omega, const double *z, const double *zLBs, const double *zUBs, 
	const double rho, const double *scaling_vector, bool updateDisp) {

   IloEnv   env;
   try {
      IloModel model(env);
      IloNumVar a0(env);
      IloNumVarArray a(env);
      IloNumVar eta(env);
      IloNumVarArray zeta(env);
      IloRangeArray con(env);
      IloNumArray linCoeff(env);
      IloNumArray ones(env);
      IloExpr objExpr(env);

      //Populate the model.
      IloNum weightObj0(0.0);
	if(omega!=NULL){	
	    for (int i = 0; i < n1; i++) {
		weightObj0 += x[i] * (c[i] + omega[i]);
	    }
	}
	else{	
	    for (int i = 0; i < n1; i++) {
		weightObj0 += x[i] * (c[i]);
	    }
	}

	for (int j = 0; j < n2; j++) {
		weightObj0 += y[j] * d[j];
	}
      weightObj0 /= rho;
      objExpr += weightObj0*a0;
      for (int wI = 0; wI < nVertices; wI++) {
		//weightObjective[wI] = baseWeightObj[wI]/rho;
		linCoeff.add(0.0);
		ones.add(1.0);
		a.add(IloNumVar(env));
		if(omega!=NULL){	
		    for (int i = 0; i < n1; i++) {
			linCoeff[wI] += xyVertices[wI][i] * (c[i]+omega[i]);
		    }
		}
		else{
		    for (int i = 0; i < n1; i++) {
			linCoeff[wI] += xyVertices[wI][i] * (c[i]);
		    }
		}
		for (int j = n1; j < n1+n2; j++) {
		    linCoeff[wI] += xyVertices[wI][j] * c[j];
		}
		linCoeff[wI] /= rho;

		
      	objExpr += linCoeff[wI]*a[wI];
      }
      objExpr += eta;
      for (int ii = 0; ii < n1; ++ii) {
	 zeta.add(IloNumVar(env,-IloInfinity,IloInfinity));
         //objExpr += 0.5*scaling_vector[ii]*zeta[ii]*zeta[ii];
      }
   IloObjective obj = IloMinimize(env, objExpr);
   model.add(obj);
   objExpr.end();
   con.add(1.0*a0 + IloScalProd(ones,a) == 1.0);
   IloExpr etaBound(env);
   for(int ii=0; ii<n1; ii++){
	IloExpr vertConstrExpr(env);
	vertConstrExpr += x[ii]*a0;
	for(int wI=0; wI<nVertices; wI++){
	    vertConstrExpr += xyVertices[wI][ii] * a[wI];
	}	
	vertConstrExpr += -1.0*zeta[ii];
	con.add(vertConstrExpr == z[ii]);
	vertConstrExpr.end();
        etaBound += scaling_vector[ii]*zeta[ii]*zeta[ii];
   }
   con.add(etaBound - eta*eta <= 0.0);
   etaBound.end();
   model.add(con);

   IloCplex cplex(model);
   //cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
	cplex.setParam(IloCplex::RootAlg, CPX_ALG_BARRIER);
        //cplexMP.setParam(IloCplex::Param::SolutionType, CPX_BASIC_SOLN);
	//cplexMP.setParam(IloCplex::RootAlg, 0);
	//cplex.setParam(IloCplex::OptimalityTarget, 1);
	//cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-9);
        //cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	
	cplex.setParam(IloCplex::Threads, nThreads); 

      // Optimize the problem and obtain solution.
      if ( !cplex.solve() ) {
         env.error() << "Failed to optimize" << endl;
         throw(-1);
      }
      IloNum weight0 = cplex.getValue(a0);
	for (int i = 0; i < n1+n2; i++) {
		x[i] = weight0 * x[i];
		//if(updateDisp){dispersions[i] = weight0*dispersions[i];}
	}

#if 0
	for (int j = 0; j < n2; j++) {
		y[j] = weight0 * y[j];
	}
#endif
	for(int wI=0; wI<nVertices; wI++) {
		vecWeights[wI] = weight0*vecWeights[wI] + cplex.getValue(a[wI]);
		for (int i = 0; i < n1+n2; i++) {
			x[i] += cplex.getValue(a[wI]) * xyVertices[wI][i];
		}
		
#if 0
		for (int j = 0; j < n2; j++) {
			y[j] += cplex.getValue(a[wI]) * yVertices[wI][j];
		}
#endif
	}
//}


   }//try
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }

   env.end();
}

double PSCGScen::getDefaultPenaltyParameter() {

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

