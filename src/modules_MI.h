//functions headers of the modules for dynamic programming optimization

/////////////////////////////////////////////////////////////////////////////
//INPUT:
//memory_cuts: vector length n with recorded recursvely best cuts,
//possible values 0...n :
//0->stops (one bins); -k -> stops (two bin ([0 k-1][k ..];k->continute [.. k-1][k ...])

//OUTPUT:
//r : number cuts
//cut: vector with cuts point-> [0 cut[0]][cut[0]+1 cut[1]]...[cut[r-2] cut[r-1]]
//int reconstruction_cut_coarse(int *memory_cuts, int *memory_cuts2,int np, int n, int *cut);
int reconstruction_cut_coarse(vector<int>& memory_cuts, vector<int>& memory_cuts2, int np, int n,int *cut);

// void reconstruction_cut_coarse(int *memory_cuts, int np, int d,int **cut,int *r);

//////////////////////////////////////////////////////////////////////////////////
//update datafactors
void update_datafactors(int **sortidx, int varidx, int **datafactors,int d, int n, int **cut);
// void update_datafactors(int **sortidx, int **repeated, int varidx, int **datafactors,int d, int n, int **cut);

/////////////////////////////////////////////////////////////////////////////////////////
// compute jointfactors functions

void jointfactors_uiyx(int **datafactors,int dui, int n, int Mui, int *r,int **,int *);

void jointfactors_uyx(int **datafactors,int* ptrVar, int n, int Mui, int *r, int **uyxfactors,int *ruyx);

void jointfactors_u(int **datafactors,int *ptrIdx,int n, int Mui, int *r, int *ufactors,int *ru);


////////////////////////////////////////////////////////////////////////////////


double* computeMIcond_knml(int **, int *, int*, int , double*, double* );
double* computeMI_knml(int* xfactors,int* ufactors,int* uxfactors,int* rux,int n,double* c2terms, double* looklog, int flag=0);


/////////////////////////////////////////////////////////////////////////////
//compute I with MDL COMPLEXITY

double* computeMIcond_kmdl(int **uiyxfactors, int *ruiyx, int *r,int n, double* looklog);

double* computeMI_kmdl(int* xfactors,int* ufactors,int* uxfactors,int* rux,int n,double* looklog, int flag=0);

//obsolete
// double* computeMI_kmdl(int **, int *, int*, int , int );
