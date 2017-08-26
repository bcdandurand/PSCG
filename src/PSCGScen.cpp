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

	setCPXMIPParameters();
#if 0
	//********Setting CPLEX Parameters*****************\\
	LagrMIPInterface_->setIntParam(OsiParallelMode,0);
	LagrMIPInterface_->setIntParam(OsiOutputControl,0);
	LagrMIPInterface_->setIntParam(OsiMIPOutputControl,0);
	if (nThreads >= 0) { LagrMIPInterface_->setIntParam(OsiParallelThreads, nThreads); }
	if (disableHeuristic) { LagrMIPInterface_->setIntParam(OsiHeurFreq, -1); }
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
void PSCGModelScen::finishInitialisation() {
	if (n1==0 && n2==0){
	     cerr << "PSCGModelScen::finishInitialisation(): n1==0 and n2==0, returning...." << endl;
	     return;
	}
	x_vertex = new double[n1];
	y_vertex = new double[n2];
	x = new double[n1];
	dispersions = new double[n1];
	dispersions2 = new double[n1];
	for(int ii=0; ii<n1; ii++){ dispersions[ii]=0.0;}
	for(int ii=0; ii<n1; ii++){ dispersions2[ii]=0.0;}
	y = new double[n2];
	//x_vertex.add(n1,0.0);
	//y_vertex.add(n2,0.0);
		
	//for(int i=0; i<n1; i++) xVertices.push_back(vector<double>());// = new ptrIloNumArray[n1];
	//for(int j=0; j<n2; j++) yVertices.push_back(vector<double>());// = new ptrDouble[n2];
		
		
	//add objective
	mpObjective.setSense( IloObjective::Minimize );
		
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

	//We use dual simplex due to the nature of the problem
	//cplexMP.setParam(IloCplex::RootAlg, IloCplex::Dual);
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
		
	//cout << "Finish Initialisation: Total memory use for " << tS << ": " << env.getTotalMemoryUsage() << endl;
	initialised = true;
}

#if 1

int PSCGModelScen_Bodur::solveLagrangianProblem(const double *omega, bool doInitialSolve) {

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
int PSCGModelScen_Bodur::solveAugmentedLagrangianMIP(const double* omega, const double* z, const double* scal){
//TODO
    return solveLagrangianProblem(omega);
}
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

int PSCGModelScen_SMPS::initialLPSolve(const double* omega) {
	
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
	return 0;
}
	//setSolverStatus();
// Given a scenario index and a dual variable, find the anticipative solution for first and second stage variables.
int PSCGModelScen_SMPS::solveLagrangianProblem(const double* omega, bool doInitialSolve) {
	
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
	    //osi->unmarkHotStart();
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
	
	return solverStatus_;
}
int PSCGModelScen_SMPS::solveAugmentedLagrangianMIP(const double* omega, const double* z, const double* scal){
	OsiCpxSolverInterface* osi = LagrMIPInterface_;
        changeFromMILPToMIQP(); //This only does anything if the problem has not already been changed back
	double *diag = new double[n1+n2];
	if(omega!=NULL){
	  for (int i = 0; i < n1; i++) {
		osi->setObjCoeff(i, c[i] + omega[i] - z[i]*scal[i]);
		diag[i]=scal[i];
	  }
	}
	else{
	  for (int i = 0; i < n1; i++) {
		osi->setObjCoeff(i, c[i] - z[i]*scal[i]);
		diag[i]=scal[i];
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
		objVal += 0.5*scal[ii]*z[ii]*z[ii];
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
void PSCGModelScen::updatePrimalVariables_OneScenario(const double *omega, const double *z, const double *scaling_vector, double *z_average, int vertexIndex) {
	
	double numerator = 0.0;
	double denominator = 0.0;
	double dir = 0.0;
	double *vertX, *vertY;
	if(vertexIndex==-1){
	    vertX = x_vertex;
	    vertY = y_vertex;
	}
	
	for (int i = 0; i < n1; i++) {
		if(omega==NULL){numerator -= (c[i] + scaling_vector[i] * (x[i] - z[i])) * (vertX[i] - x[i]);}
		else{numerator -= (c[i] + omega[i] + scaling_vector[i] * (x[i] - z[i])) * (vertX[i] - x[i]);}
		denominator += (vertX[i] - x[i]) * scaling_vector[i] * (vertX[i] - x[i]);
	}
	
	for (int j = 0; j < n2; j++) {
		numerator -= d[j] * (vertY[j] - y[j]);
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
	double oldX;
	for (int i = 0; i < n1; i++) {
		dir = vertX[i] - x[i];
		oldX=x[i];
		x[i] = x[i] + a * (vertX[i] - x[i]);
		if(z_average!=NULL){dispersions2[i] = (1.0-a)*dispersions2[i] + a*fabs(z_average[i]-vertX[i]) ;}
		//if(updateDisp){dispersions[i] = max( (1.0-a)*(dispersions[i]+fabs(x[i]-oldX)) , a*fabs(x[i]-x_vertex[i]) );}
	}

	for (int j = 0; j < n2; j++) {
		y[j] = y[j] + a * (vertY[j] - y[j]);
	}
}

void PSCGModelScen::computeWeightsForCurrentSoln(const double *z) {
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

void PSCGModelScen::updatePrimalVariablesHistory_OneScenario(const double *omega, const double *z, const double *zLBs, const double *zUBs, 
	const double *scaling_vector, bool updateDisp) {
   double direction=1.0;
   try{
#if 0
	for(int i=0; i<n1; i++){
	    mpAuxVariables[i].setLB(zLBs[i]-z[i]);
	    mpAuxVariables[i].setUB(zUBs[i]-z[i]);
	}
#endif
	if(updateDisp) mpWeight0.setBounds(0.0,0.0);
	IloNum weightObj0(0.0);
	for (int wI = 0; wI < nVertices; wI++) {
		weightObjective[wI] = baseWeightObj[wI];
		if(omega!=NULL){	
		    for (int i = 0; i < n1; i++) {
			weightObjective[wI] += xVertices[wI][i] * omega[i];
		    }
		}
	}
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

	mpObjective.setExpr(IloScalProd(mpWeightVariables, weightObjective) + weightObj0*mpWeight0 + quadraticTerm);
		
	//modify vertex constraint
	for (int i = 0; i < n1; i++) {
		mpVertexConstraints[i].setBounds(z[i], z[i]);
		mpVertexConstraints[i].setLinearCoef(mpWeight0, x[i]);
	}

	//cout << cplexQP.getModel() << endl;
	//cplexQP.exportModel("mymodel.lp");
	//cplexMP.presolve(IloCplex::Dual);
	//cout << "MP Subproblem CPLEX status: " << cplexMP.getCplexStatus() << endl;
		//cout << "MP Subproblem CPLEX status: " << cplexMP.getCplexStatus() << endl;
	//cout << mpModel << endl;
	//printVertices();
	if (!cplexMP.solve()) {
		//cout << "Num vars: " << cplexQP.getNcols() << endl;
		cout << "MP Subproblem CPLEX status: " << cplexMP.getCplexStatus() << endl;
		env.error() << "Failed to optimize in update step" << endl;
		//cout << mpModel << endl;
		//cplexMP.refineConflict();
		//cout << cplexMP.getConflict() << endl;
		cout << "Number of vertices: " << getNVertices() << endl;
		cplexMP.exportModel("infeasModel.lp");
		cout << "Refreshing solution..." << endl;
		refresh(omega,z,scaling_vector);
		//assert(false);
		//throw(-1);
	}
	if(cplexMP.getCplexStatus()!=IloCplex::Optimal) {
		cout << "MP Subproblem CPLEX status: " << cplexMP.getCplexStatus() << endl;
		//printVertices();
		//printWeights();
	}

	cplexMP.getValues(weightSoln, mpWeightVariables);
	weight0 = cplexMP.getValue(mpWeight0);
	
	polishWeights(); //fix any unfortunate numerical quirks
	// note: the final weight corresponds to the existing x
	for (int i = 0; i < n1; i++) {
		x[i] = weight0 * x[i];
		//if(updateDisp){dispersions[i] = weight0*dispersions[i];}
	}

	for (int j = 0; j < n2; j++) {
		y[j] = weight0 * y[j];
	}
	for(int wI=0; wI<nVertices; wI++) {
		vecWeights[wI] = weight0*vecWeights[wI] + weightSoln[wI];
		for (int i = 0; i < n1; i++) {
			x[i] += weightSoln[wI] * xVertices[wI][i];
		}
#if 0
	}


      if(updateDisp){
	for (int ii = 0; ii < n1; ii++) {
		dispersions[ii] = 0.0;
		for(int wI=0; wI<nVertices; wI++) {
		    dispersions[ii] += weightSoln[wI]*fabs(xVertices[ii][wI]-x[ii]) ; 
		}
	}
      }


	for(int wI=0; wI<nVertices; wI++) {
#endif		
		
		for (int j = 0; j < n2; j++) {
			y[j] += weightSoln[wI] * yVertices[wI][j];
		}
	}
	
	if(updateDisp) mpWeight0.setBounds(0.0,1.0);
   }
   catch(IloException& e){
	//cout << "MPSolve error: " << e.getMessage() << endl;
	cout << "Exception caught...Refreshing solution..." << endl;
	//refresh();
	refresh(omega,z,scaling_vector);
	e.end();
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

