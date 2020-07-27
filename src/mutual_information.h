#ifndef MIIC_MUTUAL_INFORMATION_H_
#define MIIC_MUTUAL_INFORMATION_H_

#include <vector>

#include "computation_cache.h"

namespace miic {
namespace computation {

using std::vector;

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
int reconstruction_cut_coarse(vector<int> &memory_cuts,
    vector<int> &memory_cuts2, int np, int n, int *cut);
void update_datafactors(
    vector<vector<int> > &sortidx, int varidx,
    int** datafactors, int d, int n, int **cut);
// compute jointfactors functions
void jointfactors_uiyx(
    int **datafactors, int dui, int n, int Mui, int *r, int **, int *);
void jointfactors_u(int **datafactors, int *ptrIdx, int n, int Mui, int *r,
    int *ufactors, int *ru);
vector<double> computeMI_knml(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, std::shared_ptr<CtermCache> cache, int flag = 0);
vector<double> computeMI_knml(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, int n_eff, vector<double> sample_weights,
    std::shared_ptr<CtermCache> cache, int flag = 0);
vector<double> computeMI_kmdl(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, std::shared_ptr<CtermCache> cache, int flag = 0);
}  // namespace computation
}  // namespace miic

#endif  // MIIC_MUTUAL_INFORMATION_H_
