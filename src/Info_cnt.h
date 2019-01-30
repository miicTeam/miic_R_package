#include <string>
#include <map>
#include <tuple>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////
//OPTIMIZATION FUNCTION MDL COMPLEXITY
double optfun_onerun_kmdl_coarse(int *sortidx_var, int *data, int nbrV, int **factors, int *r, double pxy,
								 double sc, double sc_levels1, double sc_levels2, int n, int nnr,int *cut,
								 int *r_opt, int maxbins, double* looklog, double* lookH,
								 double** cterms, int cplx, double* sample_weights, bool flag_effN);

//////////////////////////////////////////////////////////////////////////////////
double* compute_Ixy_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels, 
						 int n,int maxbins, int initbins, int **cut, int *r, double* c2terms,
						 double** cterms, double* looklog, double* lookH, int cplx, double* sample_weights, bool flag_effN);

double* compute_mi_cond_alg1(int** data, int** sortidx, int* AllLevels, int* ptr_cnt, int* ptrVarIdx,
							 int nbrUi,  int n, int maxbins, int initbins, double* c2terms, 
							 double** cterms, double* looklog, double* lookH, int cplx, double* sample_weights, bool flag_effN);

double* compute_Rscore_Ixyz_alg5(int** data, int** sortidx, int* AllLevels, int* ptr_cnt, int* ptrVarIdx,
								 int nbrUi, int ptrZiIdx, int n, int maxbins, int initbins, double* c2terms,
								 double** cterms, double* looklog, double*lookH, int cplx, double* sample_weights, bool flag_effN);

double* compute_Rscore_Ixyz_new_alg5(int** data, int** sortidx, int* AllLevels, int* ptr_cnt, int* ptrVarIdx,
									 int nbrUi, int ptrZiIdx, int n, int maxbins, int initbins, double* c2terms, 
									 double** cterms, double* looklog, double* lookH, int cplx,double* sample_weights, bool flag_effN);
