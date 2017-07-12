#ifndef STRUCTUREDEFS_H
#define STRUCTUREDEFS_H

#include <list>

// WARNING WARNING
// Be VERY CAREFUL what you pass from exterior stuff compiled with MPI to interior compiled with GCC!
// Stick to things which you KNOW will work in two different compilers, in C and in C++ EXACTLY the same way!
// If you don't the struct will be filled with garbage and it will NOT be obvious why!

extern "C" {

typedef double* ptrdouble;
	
//enum FileType{ NONE = 0, CAP = 1, SMPS = 2};

//The canaries are there to make sure nothing has gone wrong 
//Unfortunately the file type cannot be an enum.
typedef struct SMIP_fileRequest{
	const char* filepath;
	const char* filename;
	int ftype;	
	int totalScenarios;
	bool lin1; // linear relaxation of first stage variables
	bool lin2; // linear relaxation of second stage variables
	bool mpiHead; // am I the head node which should do the talking?
	bool disableHeuristic; //should I disable the CPLEX integer solution heuristic?
	int nThreads; //number of CPLEX threads
	
	//SMIP_fileRequest() : ftype(0), totalScenarios(0), lin1(false), lin2(false), mpiHead(false), canary1(12345), canary2(56789) {;}
	SMIP_fileRequest(const char* fp,const char* fn, int ft, int nS, bool head, bool disHeur, int nT) : filepath(fp), filename(fn), ftype(ft), totalScenarios(nS), lin1(false), lin2(false), mpiHead(head), disableHeuristic(disHeur), nThreads(nT) {;}
} SMIP_fileRequest;

typedef struct SMIP_fileReply{
	bool fileLoaded;
	int n1;
	int n2;
	int nS;
	double* p;

	SMIP_fileReply() : fileLoaded(false), n1(0), n2(0), nS(0), p(0) {;}
	~SMIP_fileReply() {;}
} SMIP_fileReply;

typedef struct SMIP_qu_LagrProblem {
	unsigned int thisScenario;
	double* dvar;
	SMIP_qu_LagrProblem(unsigned int tS, double *dv) : thisScenario(tS), dvar(dv) {;}
	~SMIP_qu_LagrProblem() {;}
} SMIP_qu_LagrProblem;

typedef struct SMIP_qu_singleVarUpdate {
	int tS;
	double* scaling_vector;
	double* x;
	double* y;
	double* x_vertex;
	double* y_vertex;
	double* z;
	double* dvar;
	SMIP_qu_singleVarUpdate(int scen, double* scal, double* xx, double* yy, double* xv, double* yv, double* zz, double* d) : tS(scen), scaling_vector(scal), x(xx), y(yy), x_vertex(xv), y_vertex(yv), z(zz), dvar(d) {;}
	~SMIP_qu_singleVarUpdate() {;}
} SMPI_qu_singleVarUpdate;

typedef struct SMIP_ans_singleVarUpdate {
	double a;
	SMIP_ans_singleVarUpdate() : a(-1) {;}
	~SMIP_ans_singleVarUpdate() {;}
} SMIP_ans_singleVarUpdate;

typedef struct SMIP_qu_singleVarUpdateHistory {
	int tS;
	double* scaling_vector;
	double* x;
	double* y;
	std::list<double*> *x_vertices;
	std::list<double*> *y_vertices;
	double* z;
	double* dvar;
	SMIP_qu_singleVarUpdateHistory(int scen, double* scal, double* xx, double* yy, std::list<double*> * xv, std::list<double*> * yv, double* zz, double* d) : tS(scen), scaling_vector(scal), x(xx), y(yy), x_vertices(xv), y_vertices(yv), z(zz), dvar(d) {;}
	~SMIP_qu_singleVarUpdateHistory() {;}
} SMPI_qu_singleVarUpdateHistory;

typedef struct SMIP_ans_singleVarUpdateHistory {
	double* weights;
	SMIP_ans_singleVarUpdateHistory(double* w) : weights(w) {;}
	~SMIP_ans_singleVarUpdateHistory() {;}
} SMIP_ans_singleVarUpdateHistory;

typedef struct SMIP_qu_fwmmSolveTilde {
	int thisScenario;
	double* scaling_vector;
	double* x;
	double* z;
	double* dvar;
	
	SMIP_qu_fwmmSolveTilde(unsigned int tS, double* pC, double *xv, double *zv, double *dv) : thisScenario(tS), scaling_vector(pC), x(xv), z(zv), dvar(dv) {;}
	~SMIP_qu_fwmmSolveTilde() {;} // arguments are "loaned out" by the questioner and belong to it. They should not be deleted.
} SMIP_qu_fwmmSolveTilde;

typedef struct SMIP_qu_PHsolveUpdate {
	int thisScenario;
	double* scaling_vector;
	double* z;
	double* dvar;
	SMIP_qu_PHsolveUpdate(unsigned int s, double* scal, double* zz, double* dv) : thisScenario(s), scaling_vector(scal), z(zz), dvar(dv) {;}
	~SMIP_qu_PHsolveUpdate() {;}
} SMIP_qu_PHsolveUpdate;

typedef struct SMIP_qu_getLagrangianGradient {
	int thisScenario;
	double* dvar;
	SMIP_qu_getLagrangianGradient(int tS, double* d) : thisScenario(tS), dvar(d) {;}
	~SMIP_qu_getLagrangianGradient() {;}
} SMIP_qu_getLagrangianGradient;

typedef struct SMIP_ans_getLagrangianGradient {
	double* x_gradient;
	double* y_gradient;
	SMIP_ans_getLagrangianGradient(double* xg, double* yg) : x_gradient(xg), y_gradient(yg) {;}
	~SMIP_ans_getLagrangianGradient() {;}
} SMIP_ans_getLagrangianGradient;

typedef struct SMIP_qu_secondStage {
	int thisScenario;
	double* x;
	SMIP_qu_secondStage(int tS, double* xvar) : thisScenario(tS), x(xvar) {;}
	~SMIP_qu_secondStage() {;}
} SMIP_qu_secondStage;

typedef struct SMIP_ans_secondStage {
	double* y;
	SMIP_ans_secondStage(double* yver) : y(yver) {;}
	~SMIP_ans_secondStage() {;}
} SMIP_ans_secondStage;

typedef struct SMIP_ans {
	double* obj;
	double* x;
	double* y;
	
	SMIP_ans(double* o, double* xvar, double* yvar) : x(xvar), y(yvar), obj(o) {;}
	~SMIP_ans() {;} // obj, x and y are "loaned out" by the questioner and belong to it. They should not be deleted.
} SMIP_ans;

}
#endif
