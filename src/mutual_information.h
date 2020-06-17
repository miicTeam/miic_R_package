#ifndef MIIC_MUTUAL_INFORMATION_H_
#define MIIC_MUTUAL_INFORMATION_H_

#include <vector>

namespace miic {
namespace computation {

// functions headers of the modules for dynamic programming optimization
// INPUT:
// memory_cuts: vector length n with recorded recursvely best cuts,
// possible values 0...n :
// 0->stops (one bins); -k -> stops (two bin ([0 k-1][k ..];
// k -> continute [.. k-1][k ...])
// OUTPUT:
// r : number cuts
// cut: vector with cuts point-> [0 cut[0]][cut[0]+1 cut[1]]...[cut[r-2]
// cut[r-1]] int reconstruction_cut_coarse(int *memory_cuts, int
// *memory_cuts2,int np, int n, int *cut);
int reconstruction_cut_coarse(std::vector<int> &memory_cuts,
    std::vector<int> &memory_cuts2, int np, int n, int *cut);
void update_datafactors(
    std::vector<std::vector<int> > &sortidx, int varidx,
    int** datafactors, int d, int n, int **cut);
// compute jointfactors functions
void jointfactors_uiyx(
    int **datafactors, int dui, int n, int Mui, int *r, int **, int *);
void jointfactors_u(int **datafactors, int *ptrIdx, int n, int Mui, int *r,
    int *ufactors, int *ru);
double *computeMI_knml(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, double *c2terms, double *looklog, int flag = 0);
double *computeMI_knml(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, int n_eff, double *c2terms, double *looklog,
    std::vector<double> sample_weights, int flag = 0);
double *computeMI_kmdl(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, double *looklog, int flag = 0);
}  // namespace computation
}  // namespace miic

#endif  // MIIC_MUTUAL_INFORMATION_H_
