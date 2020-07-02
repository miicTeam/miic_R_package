// Module to compute the conditional mutual information for mixed continuous
// and discrete variables.
#ifndef MIIC_INFO_CNT_H_
#define MIIC_INFO_CNT_H_

#include <algorithm>
#include <functional>
#include <vector>

#include "environment.h"
#include "structure.h"

namespace miic {
namespace computation {

using std::vector;

double* compute_Ixy_alg1(vector<vector<int> > data, vector<vector<int> > sortidx,
    vector<int> ptr_cnt, vector<int> ptrVarIdx, vector<int> AllLevels, int n,
    int** cut, int* r, vector<double> sample_weights, bool flag_sample_weights,
    structure::Environment& environment, bool saveIterations = false);

double* compute_mi_cond_alg1(vector<vector<int> > data, vector<vector<int> > sortidx,
    vector<int> AllLevels, vector <int> ptr_cnt, vector<int> ptrVarIdx,
    int nbrUi, int n, vector<double> sample_weights,
    bool flag_sample_weights, structure::Environment& environment,
    bool saveIterations = false);

double* compute_Rscore_Ixyz_alg5(vector<vector<int> > data,
    vector<vector<int> > sortidx, vector<int> AllLevels, vector<int> ptr_cnt,
    vector<int> ptrVarIdx, int nbrUi, int ptrZiIdx, int n,
    vector<double> sample_weights, bool flag_sample_weights,
    structure::Environment& environment, bool saveIterations = false);

void optfun_onerun_kmdl_coarse(vector<int> sortidx_var, vector<int> data, int nbrV,
    int** factors, int* r, double sc, int sc_levels1, int previous_levels,
    int n, int nnr, int* cut, int* r_opt, vector<double> sample_weights,
    bool flag_sample_weights, structure::Environment& environment);

void reset_u_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx,
    int init_nbin, int maxbins, int lbin, int* r, int* AllLevels, int n);

void reset_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx,
    int init_nbin, int maxbins, int lbin, int* r, int* AllLevels, int n);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_INFO_CNT_H_
