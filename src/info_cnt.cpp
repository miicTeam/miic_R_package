#include "info_cnt.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <map>
#include <sstream>
#include <tuple>

#include "compute_info.h"
#include "mutual_information.h"
#include "structure.h"

#define STEPMAX 50
#define EPS 1e-5
#define FLAG_CPLX 1
//#define _MY_DEBUG_MInoU 1
//#define _MY_DEBUG_MI 1
//#define _MY_DEBUG_NEW_OPTFUN 1
//#define _MY_DEBUG_NEW 1

namespace miic {
namespace computation {

using uint = unsigned int;
using std::min;
using std::vector;
using namespace miic::structure;

void reset_u_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx,
    int init_nbin, int maxbins, int lbin, int* r, int* AllLevels, int n) {
  for (int l = 2; l < nbrUi + 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (int j = 0; j < init_nbin - 1; j++) {
        cut[l][j] = j * lbin + lbin - 1;
      }
      cut[l][init_nbin - 1] = n - 1;
      for (int j = init_nbin; j < maxbins; j++) {
        cut[l][j] = 0;
      }
      r[l] = init_nbin;
    } else {
      r[l] = AllLevels[ptrVarIdx[l]];
    }
  }
}

void reset_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx,
    int init_nbin, int maxbins, int lbin, int* r, int* AllLevels, int n) {
  for (int l = 0; l < nbrUi + 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (int j = 0; j < init_nbin - 1; j++) {
        cut[l][j] = j * lbin + lbin - 1;
      }
      cut[l][init_nbin - 1] = n - 1;
      for (int j = init_nbin; j < maxbins; j++) {
        cut[l][j] = 0;
      }
      r[l] = init_nbin;
    } else {
      r[l] = AllLevels[ptrVarIdx[l]];
    }
  }
}
// GENERAL VARIABLES NOTATION

// int n: number of sample
// int coarse : coarse graining step : minimum length of a bin : considering
// only n/coarse possible cuts

// double* looklog : lookup table for the logarithms of natural numbers up to n
// double* c2terms : precomputed vectors for the MDL complexity

/*  DYNAMIC PROGRAMMING OPTIMIZATION FUNCTION
 *
 *  Takes as input the factors and joint factors and optimizes a variable by
 * maximizing a given information. nbrV determines the number of factors to
 * consider.
 *
 * ////////////////////////////////////////////
 * INPUT
 *
 * The optimizing variable is defined by <sortidx_var> and <data>.
 * <sortidx> Contains the rank indices (ordering) of the optimized variable,
 * <data> contains the ranks.
 *
 * <nbrV> is the number of terms in the optimization functions, consisting in
 * variables and joint variables, which (discretized) values are stored in
 * <factors> : ( <nbrV> x n )  matrix : Example <nbrV>=2 : Optimizing X on
 * I(X;Y) xy factors : factors[0] x factors : factors[1]
 *
 * <r> is the number of unique values of each factor.
 *
 * <sc> is the stochastic complexity used for the simple MDL cost.
 * Its value is 0.5 * number_of_other_levels where number of other levels is the
 * product of number of unique observed levels of all variables except the
 * variable being optimized.
 *
 * <sc_levels> and <previous_levels> are respectively the number of levels of
 * the other variable in the mutual information term and the number of levels of
 * the otpimized variable in the previous iteration (or initial number for the
 * first iteration).
 *
 * <n> is the number of samples in total.
 *
 * <nnr> is the number of non repeated samples.
 *
 * <cut> is the vector containing cutpoints for discretizing the optimized
 * variables. Modified at the end wih the optimal solution with a call to
 * reconstruction_cut_coarse().
 *
 * <r_opt> is a pointer to the number of bins of the optimized variable, also
 * updated by the call to reconstruction_cut_coarse().
 *
 * <maxbins> is the maximum number of bins allowed. It controls the resolution
 * of possible cutpoints : For example, 50 maximum cutpoints on a 1000 samples
 * vector results in possible cutpoints every 1000/50 = 20 positions. This
 * speeds up significantly the dynamic programming at the cost of finding
 * approximated solutions.
 *
 * <looklog, lookH, cterms> are lookup tables for computing respectively log,
 * entropy and stochastic complexity (as in Kontkanen & Myllymäki, 2005)
 *
 * <cplx> is the choice of complexity : 0 for simple MDL (product of all
 * observed values, see <sc>) and 1 for NML with stochastic complexity and a
 * combinatorial term with the previous number of levels.
 *
 * <sample_weights> are unique weights for each sample to account for effective
 * N != N, if the <flag_efN> is set to true.
 *
 *
 * Use-cases :
 *  Optimizing Ik(X;Y) by choosing discretization on X:
 *      <sortidx_var>, <data> -> Contains the ranks of X
 *      <nbrV> -> 2
 *      <factors>
 *           [0] : The discretized Y
 *           [1] : A vector of one repeated level (starting from X with one
 * single bin) <sc_levels1>        -> the number of levels of discretized Y
 *      <previous_levels>   -> the number of levels of discretizedd X in the
 * previous iteration (for cost correction) <cut>               -> the vector
 * that will contain the cutpoints for X
 */
// inline __attribute__((always_inline))
std::vector<double> optfun_onerun_kmdl_coarse(int* sortidx_var, int* data, int nbrV,
    int** factors, int* r, double sc, int sc_levels1, int previous_levels,
    int n, int nnr, int* cut, int* r_opt, Environment& environment) {
  int coarse = ceil(1.0 * nnr / environment.maxbins);  // step coarse graining
  if (coarse < 1) coarse = 1;
  uint np = ceil(1.0 * nnr / coarse);  // number of possible cuts

  // temp variables to memorize optimization cuts
  vector<int> memory_cuts_idx(np);  // indexes of the cuts (1..np)
  vector<int> memory_cuts_pos(np);  // positions of the cuts (1..n)

  // dynamic programming optimize function and memorize of cuts
  double Imax;   // I
  double Ikmax;  // I-k_NML

  // function max value at each step
  // double* Ik=(double *)calloc(np,sizeof(double));
  vector<double> I(np);   // The optimal information value found at each idx.
  double I_kj;            // The information value for a bin fom idx k to j.
  vector<double> Ik(np);
  double Ik_kj;

  int njforward(0), nkforward(0);  // Indexes at current position of j and k

  // entropy in kj interval for the <nbrV>+1 terms
  vector<double> H_kj(nbrV);
  vector<double> Hk_kj(nbrV);

  vector<vector<uint> > counts(nbrV);
  for (int m = 0; m < nbrV; m++) {
    counts[m].resize(r[m]);
  }

  Grid2d<uint> coarse_counts_marginal(np, r[1], 0);
  Grid2d<uint> coarse_counts_joint(np, r[0], 0);

  vector<vector<uint> > counts_k(nbrV);
  for (int m = 0; m < nbrV; m++) {
    counts_k[m].resize(r[m]);
  }
  uint weighted_count;

  vector<bool> check_repet(n);
  for (int i = 0; i < (n - 1); i++) {
    check_repet[i] = (data[sortidx_var[i + 1]] != data[sortidx_var[i]]);
  }

  vector<uint> n_values(np);
  vector<double> sum_sample_weights(np);
  int ir(0);  // Iterator on non repeated values
  uint level_marginal(0), level_joint(0);
  for (uint i = 0; i < np; i++) {
    ir = 0;  // iterator on not repeated values
    if (environment.flag_sample_weights && i > 0)
      sum_sample_weights[i] = sum_sample_weights[i - 1];

    while ((ir < coarse) && (njforward < n)) {
      level_marginal = factors[1][sortidx_var[njforward]];
      level_joint = factors[0][sortidx_var[njforward]];
      coarse_counts_marginal(i, level_marginal)++;
      coarse_counts_joint(i, level_joint)++;

      if (environment.flag_sample_weights)
        sum_sample_weights[i] += environment.sampleWeights[njforward];
      if (njforward + 1 < n) {  // check no repetition
        ir += int(check_repet[njforward]);
      }
      njforward++;
    }
    n_values[i] = njforward;
  }

  double ef_nj = 0;
  double efN_factor;
  double ef_nk;

  double I_newbin(0);
  double Ik_newbin(0);

  njforward = 0;  // iterator on values

  // moving j over the np possible cuts
  for (uint j = 0; j < np; j++) {  // j=1...n-1

    njforward = n_values[j];
    ef_nj = sum_sample_weights[j];
    if (environment.flag_sample_weights) efN_factor = ef_nj / njforward;
    coarse_counts_joint.add_row(counts[0], j);
    coarse_counts_marginal.add_row(counts[1], j);

    Hk_kj[0] = 0;  // joint
    Hk_kj[1] = 0;  // marginal
    H_kj[0] = 0;  // joint
    H_kj[1] = 0;  // marginal
    for (int level = 0; level < r[0]; level++) {
      weighted_count = environment.flag_sample_weights
                           ? int(efN_factor * counts[0][level] + 0.5)
                           : counts[0][level];
      Hk_kj[0] += environment.lookH[weighted_count];
    }
    H_kj[0]   = Hk_kj[0];

    for (int level = 0; level < r[1]; level++) {
      weighted_count = environment.flag_sample_weights
                           ? int(efN_factor * counts[1][level] + 0.5)
                           : counts[1][level];
      Hk_kj[1] -= environment.lookH[weighted_count];
      H_kj[1]  -= environment.lookH[weighted_count];

      if (environment.cplx == 0 && counts[1][level] > 0)
        Hk_kj[1] -= sc * environment.looklog[n];
      else if (environment.cplx == 1) {
        Hk_kj[1] -= computeLogC(weighted_count, sc_levels1, environment.looklog,
            environment.cterms);
      }
    }

    I[j]  = 0;
    Ik[j] = 0;
    for (int m = 0; m < nbrV; m++) {
      I[j]  += H_kj[m];  // herve
      Ik[j] += Hk_kj[m];  // herve
    }

    Imax = -DBL_MAX;
    Ikmax = -DBL_MAX;

    for (int m = 0; m < nbrV; m++) {
      for (int level = 0; level < r[m]; level++)
        counts_k[m][level] = counts[m][level];
    }
    ef_nk = ef_nj;

    // moving k

    // k iterator on possible cuts
    nkforward = 0;  // iterator on values
    // it iterator on not repeated vales

    // Before trying to create bins : solution is one single bin from 0 to j.
    memory_cuts_idx[j] = 0;
    memory_cuts_pos[j] = 0;

    for (uint k = 0; k < j; k++) {  // k=1...n-2 possible cuts

      nkforward = n_values[k];
      ef_nk = sum_sample_weights[j] - sum_sample_weights[k];
      if (environment.flag_sample_weights)
        efN_factor = ef_nk / (njforward - nkforward);
      coarse_counts_marginal.subtract_row(counts_k[1], k);
      coarse_counts_joint.subtract_row(counts_k[0], k);

      H_kj[0]  = 0;
      Hk_kj[0] = 0;
      for (int level = 0; level < r[0]; level++) {
        weighted_count = environment.flag_sample_weights
                             ? int(efN_factor * counts_k[0][level] + 0.5)
                             : counts_k[0][level];
        // j_efN *
        // nxyu_k[m][xyu]*looklog[int(j_efN
        // * nxyu_k[m][xyu] + 0.5)];
        Hk_kj[0] += environment.lookH[weighted_count];
      }
      H_kj[0] = Hk_kj[0];

      H_kj[1]  = 0;
      Hk_kj[1] = 0;
      for (int level = 0; level < r[1]; level++) {
        weighted_count = environment.flag_sample_weights
                             ? int(efN_factor * counts_k[1][level] + 0.5)
                             : counts_k[1][level];
        // j_efN *
        // nxyu_k[m][xyu]*looklog[int(j_efN
        // * nxyu_k[m][xyu] + 0.5)];<Paste>
        Hk_kj[1] -= environment.lookH[weighted_count];
        H_kj[1]  -= environment.lookH[weighted_count];

        if (environment.cplx == 0 && counts_k[1][level] > 0)
          Hk_kj[1] -= sc * environment.looklog[n];
        else if (environment.cplx == 1) {
          Hk_kj[1] -= computeLogC(weighted_count, sc_levels1,
              environment.looklog, environment.cterms);
        }
      }

      I_kj = 0;
      Ik_kj = 0;
      for (int m = 0; m < nbrV; m++) { 
        I_kj  += H_kj[m];
        Ik_kj += Hk_kj[m];
      }
      if (environment.cplx == 1) {
        // Combinatorial approximation
        Ik_kj -= logchoose(np - 1, previous_levels - 1, environment.looklog,
                     environment.lookchoose) /
                 (previous_levels - 1);
      }

      if ((Ik[k] + Ik_kj) > Ik[j]) {
        I_newbin  = I[k]  + I_kj;  //[0.. cuts.. k-1][k j] //herve
        Ik_newbin = Ik[k] + Ik_kj;  //[0.. cuts.. k-1][k j] //herve
        if (Ikmax < Ik_newbin) {
          Imax  = I_newbin;
          Ikmax = Ik_newbin;
          I[j]  = Imax;  // optimized function for the interval [0 j]
          Ik[j] = Ikmax;  // optimized function for the interval [0 j]
          if (memory_cuts_idx[k + 1] == 0) {
            memory_cuts_idx[j] = -k - 1;  // index  of the (last) optimal cut
          } else {
            memory_cuts_idx[j] = k + 1;  // index  of the (last) optimal cut
          }
          memory_cuts_pos[j] =
              nkforward - 1;  // position  of the (last) optimal cut
        }
      }
    }
  }

  *r_opt =
      reconstruction_cut_coarse(memory_cuts_idx, memory_cuts_pos, np, n, cut);

  vector<double> res(2);
  res[0] = I[np-1];
  res[1] = Ik[np-1];
  return(res);
}

// compute I(x,y)
// dynamic programming for optimizing variables binning
// optimize on x I(x,y): Hx - Hxy - kmdl
// optimize on y I(x,y): Hy - Hxy - kmdl
// until convergence
double* compute_Ixy_alg1(int** data, int** sortidx, int* ptr_cnt,
    int* ptrVarIdx, int* AllLevels, int n, int** cut, int* r,
    Environment& environment, bool saveIterations) {
  int maxbins = environment.maxbins;
  int initbins = environment.initbins;
  double* c2terms = environment.c2terms;
  double** lookchoose = environment.lookchoose;
  double* looklog = environment.looklog;
  int cplx = environment.cplx;

  int j, l;

  // res_tempults res_temp[0]->I,res_temp[1]->I-k
  double* res_temp = (double*)calloc(2, sizeof(double));

  // allocation factors  x y
  int** datafactors;
  datafactors = (int**)calloc((2), sizeof(int*));

  for (l = 0; l < (2); l++) {
    datafactors[l] = (int*)calloc(n, sizeof(int));
  }

  int* r_temp = (int*)calloc(3, sizeof(int));
  int* ptr = (int*)calloc(2, sizeof(int));
  int* xy_factors = (int*)calloc(n, sizeof(int));
  int rxy = 0;
  ptr[0] = 0;
  ptr[1] = 1;
  // initialization of datafactors && sortidx
  for (l = 0; l < 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
    } else {
      for (j = 0; j <= n - 1; j++) {
        datafactors[l][j] = data[ptrVarIdx[l]][j];
      }
    }
  }
  jointfactors_u(datafactors, ptr, n, 2, r, xy_factors, &rxy);
  r_temp[0] = r[0];
  r_temp[1] = r[1];
  r_temp[2] = rxy;
  if (cplx == 1)
    res_temp = computeMI_knml(datafactors[0], datafactors[1], xy_factors,
        r_temp, n, c2terms, looklog, 0);
  else
    res_temp = computeMI_kmdl(
        datafactors[0], datafactors[1], xy_factors, r_temp, n, looklog, 0);

  // all discrete
  if (ptr_cnt[ptrVarIdx[0]] == 0 && ptr_cnt[ptrVarIdx[1]] == 0) {
    free(xy_factors);
    free(r_temp);
    free(ptr);

    for (l = 0; l < (2); l++) free(datafactors[l]);
    free(datafactors);

    return res_temp;
  }

  // Find the best initial conditions with the same number of bins (equalfreq)
  // on all continuous variables.
  double max_res = -DBL_MAX;
  int max_initbins;
  int min_unique_values = n;
  for (l = 0; l < 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1)
      min_unique_values = min(min_unique_values, AllLevels[ptrVarIdx[l]]);
  }
  for (int new_initbins = 2; (new_initbins < initbins) && (new_initbins < 20) &&
                             (new_initbins < min_unique_values);
       new_initbins++) {
    int lbin = floor(n / new_initbins);
    if (lbin < 1) {
      lbin = 1;
      new_initbins = n;
    }
    // Reinitialization cut and r
    for (l = 0; l < 2; l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        for (j = 0; j < new_initbins - 1; j++) {
          cut[l][j] = j * lbin + lbin - 1;
        }
        cut[l][new_initbins - 1] = n - 1;
        r[l] = new_initbins;
      } else {
        r[l] = AllLevels[ptrVarIdx[l]];
      }
    }
    // initialization of datafactors && sortidx
    for (l = 0; l < 2; l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
      } else {
        for (j = 0; j <= n - 1; j++) {
          datafactors[l][j] = data[ptrVarIdx[l]][j];
        }
      }
    }

    jointfactors_u(datafactors, ptr, n, 2, r, xy_factors, &rxy);
    r_temp[0] = r[0];
    r_temp[1] = r[1];
    r_temp[2] = rxy;
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors[0], datafactors[1], xy_factors,
          r_temp, n, c2terms, looklog, 0);
    else
      res_temp = computeMI_kmdl(
          datafactors[0], datafactors[1], xy_factors, r_temp, n, looklog, 0);

    if (res_temp[1] > max_res) {
      max_initbins = new_initbins;
      max_res = res_temp[1];
    }
    free(res_temp);
  }

  int lbin = floor(n / max_initbins);
  if (lbin < 1) {
    lbin = 1;
    max_initbins = n;
  }
  // Reinitialization cut and r
  for (l = 0; l < 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (j = 0; j < max_initbins - 1; j++) {
        cut[l][j] = j * lbin + lbin - 1;
      }
      cut[l][max_initbins - 1] = n - 1;
      r[l] = max_initbins;
    } else {
      r[l] = AllLevels[ptrVarIdx[l]];
    }
  }
  // initialization of datafactors && sortidx
  for (l = 0; l < 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
    } else {
      for (j = 0; j <= n - 1; j++) {
        datafactors[l][j] = data[ptrVarIdx[l]][j];
      }
    }
  }
  // Run dynamic optimization with the best initial conditions.
  double sc;
  double* MI = (double*)calloc(STEPMAX, sizeof(double));
  double* MIk = (double*)calloc(STEPMAX, sizeof(double));
  double I_av, Ik_av;

  int stop, i, flag;
  int** factors1 = (int**)calloc(2, sizeof(int*));
  int* singlefactor = (int*)calloc(n, sizeof(int));
  int* rt1 = (int*)calloc(2, sizeof(int));
  int sc_levels1, sc_levels2;

  MI[0] = 0;
  MIk[0] = -DBL_MAX;

  flag = 0;
  int sc_levels_x = r[0];  // Number of levels of the first variable
  int sc_levels_y = r[1];  // Number of levels of the second variable
  int rx = r[0];
  int ry = r[1];
  int np;  // number of possible cuts
  for (stop = 1; stop < STEPMAX; stop++) {
    if (ptr_cnt[ptrVarIdx[0]] == 1) {
      // optimize on x
      // I(x;y)
      factors1[0] = datafactors[1];  // y
      factors1[1] = singlefactor;    // One single bin at start
      rt1[0] = ry;  // y
      rt1[1] = 1;   // One single bin at start

      sc_levels_y = ry;
      sc_levels1 = sc_levels_y;
      sc_levels2 = sc_levels_x;
      sc = 0.5 * (sc_levels_y - 1);
      rx = r[0];
      // Optimization run on X.
      optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]], 2,
          factors1, rt1, sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[0]],
          cut[0], &(r[0]),
          environment);  // 2 factors
    }
    if (ptr_cnt[ptrVarIdx[1]] == 1) {
      // opt y
      // I(x;y)
      factors1[0] = datafactors[0];  // x before its optimization
      factors1[1] = singlefactor;    // One single bin at start

      rt1[0] = rx;  // x before its optimization
      rt1[1] = 1;   // One single bin at start

      sc_levels_x = rx;
      sc_levels1 = sc_levels_x;
      sc_levels2 = sc_levels_y;
      sc = 0.5 * (sc_levels_x - 1);

      ry = r[1];
      // Optimization run on Y.
      // 2 factors
      optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]], 2,
          factors1, rt1, sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[1]],
          cut[1], &(r[1]), environment);
    }

    // update both datafactors
    if (ptr_cnt[ptrVarIdx[0]] == 1) {
      update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);
      rx = r[0];
    }
    if (ptr_cnt[ptrVarIdx[1]] == 1) {
      update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut);
      ry = r[1];
    }

    if (saveIterations) {
      // Save cut points
      for (j = 0; j < maxbins; j++) {
        environment.iterative_cuts[stop - 1][j] = cut[0][j];
        environment.iterative_cuts[stop - 1][j + maxbins] = cut[1][j];
      }
    }
    jointfactors_u(datafactors, ptr, n, 2, r, xy_factors, &rxy);
    r_temp[0] = r[0];
    r_temp[1] = r[1];
    r_temp[2] = rxy;
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors[0], datafactors[1], xy_factors,
          r_temp, n, c2terms, looklog, 0);
    else
      res_temp = computeMI_kmdl(
          datafactors[0], datafactors[1], xy_factors, r_temp, n, looklog, 0);
    // Adding combinatorial term
    if (ptr_cnt[ptrVarIdx[0]] == 1 && rx > 1) {
      np = min(maxbins, AllLevels[ptrVarIdx[0]]);
      if (rx < np)
        res_temp[1] -= logchoose(np - 1, rx - 1, looklog, lookchoose) / n;
    }
    if (ptr_cnt[ptrVarIdx[1]] == 1 && ry > 1) {
      np = min(maxbins, AllLevels[ptrVarIdx[1]]);
      if (ry < np)
        res_temp[1] -= logchoose(np - 1, ry - 1, looklog, lookchoose) / n;
    }

    for (i = stop - 1; i > 0; i--) {
      if (fabs(res_temp[1] - MIk[i]) < EPS) {
        flag = 1;
        Ik_av = MIk[i];
        I_av = MI[i];

        for (j = i + 1; j < stop; j++) {
          Ik_av += MIk[j];
          I_av += MI[j];
        }
        Ik_av /= (stop - i);  // average over the periodic cycle
        I_av /= (stop - i);
        break;
      }
    }
    MIk[stop] = res_temp[1];
    MI[stop] = res_temp[0];
    free(res_temp);

    if (flag || (ptr_cnt[ptrVarIdx[0]] == 0) || (ptr_cnt[ptrVarIdx[1]] == 0)) {
      break;
    }
  }  // for

  double* return_res = (double*)calloc(2, sizeof(double));
  if (flag) {
    return_res[0] = I_av;
    return_res[1] = Ik_av;
  } else {
    return_res[0] = MI[stop];
    return_res[1] = MIk[stop];
  }
  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if ((return_res[1] < 0) &&
      ((ptr_cnt[ptrVarIdx[0]] == 1) && (ptr_cnt[ptrVarIdx[1]] == 1))) {
    return_res[0] = 0;
    return_res[1] = 0;
  }

  if (saveIterations) {
    environment.iterative_cuts[stop][0] =
        -1;  // mark where we stopped iterating
    environment.iterative_cuts[stop][1] =
        100000 * return_res[1];  // Pass Ik[X;Y]
    environment.iterative_cuts[stop][2] =
        100000 * return_res[0];  // Pass I[X;Y]
    environment.iterative_cuts[stop][3] =
        100000 * max_res;  // Pass max res before optimization with equal freq
    environment.iterative_cuts[stop][maxbins] = -1;
  }

#if _MY_DEBUG_MInoU
  cout << (ptr_cnt[ptrVarIdx[0]]) << " " << (ptr_cnt[ptrVarIdx[1]]) << endl;
  printf("    final: I_xy=%lf Ik_xy=%lf\n", return_res[0], return_res[1]);

  printf("    0 : r=%d ", r[0]);
  for (i = 0; i < r[0]; i++) {
    printf("    %d ", cut[0][i]);
  }
  printf("\n");
  printf("    1 : r=%d ", r[1]);
  for (i = 0; i < r[1]; i++) {
    printf("    %d ", cut[1][i]);
  }
  printf("\n");
#endif

  free(xy_factors);
  free(r_temp);
  free(ptr);

  for (l = 0; l < (2); l++) free(datafactors[l]);
  free(datafactors);

  free(factors1);
  free(rt1);
  free(MI);
  free(MIk);
  free(singlefactor);

  return return_res;
}

double* compute_Ixy_cond_u_new_alg1(int** data, int** sortidx, int* ptr_cnt,
    int* ptrVarIdx, int* AllLevels, int nbrUi, int n, int** cut, int* r,
    int lbin, Environment& environment, bool saveIterations) {
  int maxbins = environment.maxbins;
  int initbins = environment.initbins;
  double* c2terms = environment.c2terms;
  double* looklog = environment.looklog;
  int cplx = environment.cplx;

  int j, l;
  int STEPMAX1 = 50;
  double** lookchoose = environment.lookchoose;
  // res_tempults res_temp[0]->I,res_temp[1]->I-k
  double* res_temp = (double*)calloc(2, sizeof(double));
  // allocation factors  x y
  int** datafactors = (int**)calloc((nbrUi + 2), sizeof(int*));
  for (l = 0; l < (nbrUi + 2); l++) {
    datafactors[l] = (int*)calloc(n, sizeof(int));
  }

  int* r_temp = (int*)calloc(3, sizeof(int));

  // initialization of datafactors && sortidx
  for (l = 0; l < (nbrUi + 2); l++) {
    // compute datafactors based on the positions of cut points in vector <cut>
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
    } else {  // discrete case
      for (j = 0; j <= n - 1; j++) {
        datafactors[l][j] = data[ptrVarIdx[l]][j];
      }
    }
  }

  int** uiyxfactors;  //({ui},{uiy}),{uix},{uiyx})
  uiyxfactors = (int**)calloc(4, sizeof(int*));
  for (l = 0; l < 4; l++) {
    uiyxfactors[l] = (int*)calloc(n, sizeof(int));
  }
  int* ruiyx;
  ruiyx = (int*)calloc(4, sizeof(int));

  double sc;

  double* MI = (double*)calloc(STEPMAX, sizeof(double));
  double* MIk = (double*)calloc(STEPMAX, sizeof(double));
  double* MI1 = (double*)calloc(STEPMAX1, sizeof(double));
  double* MIk1 = (double*)calloc(STEPMAX1, sizeof(double));
  double I_av1, Ik_av1;
  MI[0] = 0;
  MIk[0] = -DBL_MAX;

  // Loop variables
  int stop1, i, flag1;
  // Mutual informations for computing I(x;y|u).
  double I_x_yu;
  double I_y_xu;
  double I_x_u;
  double I_y_u;
  double cond_I;
  double Ik_x_yu;
  double Ik_y_xu;
  double Ik_x_u;
  double Ik_y_u;
  double cond_Ik;
  std::vector<double> I_x_yu_part;
  std::vector<double> I_y_xu_part;
  std::vector<double> I_x_u_part;
  std::vector<double> I_y_u_part;
  double cond_I_bis;
  double cond_Ik_bis;

  // Complexity factor (# of levels) passed to each optimization
  int sc_levels1, sc_levels2;

  int** factors1 = (int**)calloc(2, sizeof(int*));

  // Find the best initial conditions with the same number of bins (equalfreq)
  // on all continuous variables.
  double max_res = -DBL_MAX;
  int max_initbins = initbins;
  int min_unique_values = n;
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1)
      min_unique_values = min(min_unique_values, AllLevels[ptrVarIdx[l]]);
  }
  for (int new_initbins = 2; (new_initbins < initbins) && (new_initbins < 20) &&
                             (new_initbins < min_unique_values);
       new_initbins++) {
    // initialization cut and r
    for (l = 0; l < (nbrUi + 2); l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        for (j = 0; j < new_initbins - 1; j++) {
          cut[l][j] = j * lbin + lbin - 1;
        }
        cut[l][new_initbins - 1] = n - 1;
        r[l] = new_initbins;
      } else {
        r[l] = AllLevels[ptrVarIdx[l]];
      }
    }

    // compute datafactors based on the positions of cut points in vector <cut>
    for (l = 0; l < (nbrUi + 2); l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
      } else {  // discrete case
        for (j = 0; j <= n - 1; j++) {
          datafactors[l][j] = data[ptrVarIdx[l]][j];
        }
      }
    }

    jointfactors_uiyx(datafactors, -1, n, nbrUi, r, uiyxfactors, ruiyx);
    r_temp[0] = r[1];
    r_temp[1] = ruiyx[2];
    r_temp[2] = ruiyx[3];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors[1], uiyxfactors[2], uiyxfactors[3],
          r_temp, n, c2terms, looklog, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors[1], uiyxfactors[2], uiyxfactors[3],
          r_temp, n, looklog, 0);
    I_y_xu = res_temp[0];  // Before optimization on X.
    Ik_y_xu = res_temp[1];
    free(res_temp);

    r_temp[0] = r[0];
    r_temp[1] = ruiyx[1];
    r_temp[2] = ruiyx[3];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors[0], uiyxfactors[1], uiyxfactors[3],
          r_temp, n, c2terms, looklog, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors[0], uiyxfactors[1], uiyxfactors[3],
          r_temp, n, looklog, 0);
    I_x_yu = res_temp[0];  // Before updating Y (and X).
    Ik_x_yu = res_temp[1];
    free(res_temp);

    if ((Ik_y_xu + Ik_x_yu) > max_res) {
      max_initbins = new_initbins;
      max_res = (Ik_y_xu + Ik_x_yu);
    }
  }

  lbin = floor(n / max_initbins);
  if (lbin < 1) {
    lbin = 1;
    max_initbins = n;
  }

  // initialization cut and r
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (j = 0; j < max_initbins - 1; j++) {
        cut[l][j] = j * lbin + lbin - 1;
      }
      cut[l][max_initbins - 1] = n - 1;
      r[l] = max_initbins;
    } else {
      r[l] = AllLevels[ptrVarIdx[l]];
    }
  }
  // compute datafactors based on the positions of cut points in vector <cut>
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
    } else {
      for (j = 0; j <= n - 1; j++) {  // discrete case
        datafactors[l][j] = data[ptrVarIdx[l]][j];
      }
    }
  }
  // Run optimization with best initial equal freq.
  int* rt1 = (int*)calloc(2, sizeof(int));
  int* r_old = (int*)calloc((nbrUi + 2), sizeof(int));
  for (int ll = 0; ll < (nbrUi + 2); ll++) {
    r_old[ll] = r[ll];
  }
  int U_counter;

  int sc_levels_x;  // Number of levels of the first variable
  int sc_levels_y;  // Number of levels of the second variable
  flag1 = 0;
  int max_U_counter = 3;
  int np;  // number of possible cuts for combinatorial term
  for (stop1 = 1; stop1 < STEPMAX1; stop1++) {
    // optimize I(y;xu) over x and u
    U_counter = 0;
    while (U_counter < max_U_counter) {
      for (l = 0; l < nbrUi; l++) {
        if (ptr_cnt[ptrVarIdx[l + 2]] == 1) {
          // opt u
          // I(y;xu)
          jointfactors_uiyx(
              datafactors, l + 2, n, nbrUi, r_old, uiyxfactors, ruiyx);
          // init variables for the optimization run
          factors1[0] = uiyxfactors[3];  // xyu
          factors1[1] = uiyxfactors[2];  // xu

          rt1[0] = ruiyx[3];  // xyu
          rt1[1] = ruiyx[2];  // xu
          sc_levels_x = r_old[0];
          sc_levels_y = r_old[1];

          sc = 0.5 * (sc_levels_y - 1) * ruiyx[2];

          sc_levels1 = sc_levels_y;
          sc_levels2 = r_old[l + 2];  // old nlevels for combinatorial term

          // Run optimization on U.
          I_y_xu_part = optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l + 2]],
              data[ptrVarIdx[l + 2]], 2, factors1, rt1, sc, sc_levels1,
              sc_levels2, n, AllLevels[ptrVarIdx[l + 2]], cut[l + 2],
              &(r[l + 2]), environment);  // 2 factors
        }
      }  // for all Uis
      for (int ll = 0; ll < nbrUi; ll++) {
        if (ptr_cnt[ptrVarIdx[ll + 2]] == 1)
          update_datafactors(
              sortidx, ptrVarIdx[ll + 2], datafactors, ll + 2, n, cut);
        r_old[ll + 2] = r[ll + 2];
      }
      U_counter++;
      if (nbrUi == 1) U_counter = max_U_counter;
    }  // U_counter loop

    jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
    r_temp[0] = r_old[1];
    r_temp[1] = ruiyx[2];
    r_temp[2] = ruiyx[3];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors[1], uiyxfactors[2], uiyxfactors[3],
          r_temp, n, c2terms, looklog, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors[1], uiyxfactors[2], uiyxfactors[3],
          r_temp, n, looklog, 0);
    I_y_xu = res_temp[0];  // Before optimization on X.
    Ik_y_xu = res_temp[1];
    free(res_temp);
    if ((ptr_cnt[ptrVarIdx[0]] == 1) && (r_old[0] > 1)) {
      np = min(AllLevels[ptrVarIdx[0]], maxbins);
      if (r_old[0] < np) {
        Ik_y_xu -= logchoose(np - 1, r_old[0] - 1, looklog, lookchoose) / n;
        I_y_xu_part[1] -= logchoose(np - 1, r_old[0] - 1, looklog, lookchoose) / n;
      }
    }

    if (ptr_cnt[ptrVarIdx[0]] == 1) {
      // opt x
      // I(y;xu)
      // compute joint factors u yu xu xyu
      jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
      // init variables for the optimization run
      factors1[0] = uiyxfactors[1];  // uy
      factors1[1] = uiyxfactors[0];  // u

      rt1[0] = ruiyx[1];  // uy
      rt1[1] = ruiyx[0];  // u
      sc_levels_x = r_old[0];
      sc_levels_y = r_old[1];
      sc = 0.5 * (sc_levels_y - 1) * ruiyx[0];

#if _MY_DEBUG_MInoU
      printf("start optfun\n ");
      fflush(stdout);
#endif
      sc_levels1 = sc_levels_y;  // herve
      sc_levels2 = sc_levels_x;  // herve
      // Run optimization on X.
      optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]], 2,
          factors1, rt1, sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[0]],
          cut[0], &(r[0]), environment);  // 2 factors
    }

    // Reset cutpoints on U
    reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin,
        r, AllLevels, n);
    for (l = 0; l < nbrUi; l++) {
      if (ptr_cnt[ptrVarIdx[l + 2]] == 1)
        update_datafactors(
            sortidx, ptrVarIdx[l + 2], datafactors, l + 2, n, cut);
      // r[l+2] is set to init_nbins during reset_u_cutpoints
      r_old[l + 2] = r[l + 2];
    }
    // optimize I(x;yu) over y and u
    U_counter = 0;
    while (U_counter < max_U_counter) {
      for (l = 0; l < nbrUi; l++) {
        if (ptr_cnt[ptrVarIdx[l + 2]] == 1) {
          // opt u
          // I(x;yu)
          jointfactors_uiyx(
              datafactors, l + 2, n, nbrUi, r_old, uiyxfactors, ruiyx);
          // init variables for the optimization run
          factors1[0] = uiyxfactors[3];  // xyu
          factors1[1] = uiyxfactors[1];  // yu

          rt1[0] = ruiyx[3];  // xyu
          rt1[1] = ruiyx[1];  // yu
          sc_levels_y = r_old[1];
          sc_levels_x = r_old[0];

          sc_levels1 = sc_levels_x;   // herve
          sc_levels2 = r_old[l + 2];  // herve
          sc = 0.5 * (sc_levels_x - 1) * ruiyx[1];

          // Run optimization on U.
          I_x_yu_part = optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l + 2]],
              data[ptrVarIdx[l + 2]], 2, factors1, rt1, sc, sc_levels1,
              sc_levels2, n, AllLevels[ptrVarIdx[l + 2]], cut[l + 2],
              &(r[l + 2]), environment);  // 2 factors
        }
      }  // for all Uis
      for (int ll = 0; ll < nbrUi; ll++) {
        if (ptr_cnt[ptrVarIdx[ll + 2]] == 1)
          update_datafactors(
              sortidx, ptrVarIdx[ll + 2], datafactors, ll + 2, n, cut);
        r_old[ll + 2] = r[ll + 2];
      }
      U_counter++;
      if (nbrUi == 1) U_counter = max_U_counter;
    }  // U_counter loop

    jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
    r_temp[0] = r_old[0];
    r_temp[1] = ruiyx[1];
    r_temp[2] = ruiyx[3];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors[0], uiyxfactors[1], uiyxfactors[3],
          r_temp, n, c2terms, looklog, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors[0], uiyxfactors[1], uiyxfactors[3],
          r_temp, n, looklog, 0);
    I_x_yu = res_temp[0];  // Before updating Y (and X).
    Ik_x_yu = res_temp[1];
    free(res_temp);
    if ((ptr_cnt[ptrVarIdx[1]] == 1) && (r_old[1] > 1)) {
      np = min(AllLevels[ptrVarIdx[1]], maxbins);
      if (r_old[1] < np) {
        Ik_x_yu -= logchoose(np - 1, r_old[1] - 1, looklog, lookchoose) / n;
        I_x_yu_part[1] -= logchoose(np - 1, r_old[1] - 1, looklog, lookchoose) / n;
      }
    }

    if (ptr_cnt[ptrVarIdx[1]] == 1) {
      // optimize on y
      // I(x;yu)
      jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
      // init variables for the optimization run
      factors1[0] = uiyxfactors[2];  // ux
      factors1[1] = uiyxfactors[0];  // u

      rt1[0] = ruiyx[2];  // ux
      rt1[1] = ruiyx[0];  // u

      sc_levels_x = r_old[0];
      sc_levels_y = r_old[1];

      sc_levels1 = sc_levels_x;  // herve
      sc_levels2 = sc_levels_y;  // herve
      sc = 0.5 * (sc_levels_x - 1) * ruiyx[0];
      // Run optimization on Y.
      optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]], 2,
          factors1, rt1, sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[1]],
          cut[1], &(r[1]), environment);  // 2 factors
    }
    // Reset cutpoints on U
    reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin,
        r, AllLevels, n);
    for (l = 0; l < nbrUi; l++) {
      if (ptr_cnt[ptrVarIdx[l + 2]] == 1)
        update_datafactors(
            sortidx, ptrVarIdx[l + 2], datafactors, l + 2, n, cut);
      r_old[l + 2] = r[l + 2];
    }
    // optimize I(x;u) over u
    U_counter = 0;
    while (U_counter < max_U_counter) {
      for (l = 0; l < nbrUi; l++) {
        if (ptr_cnt[ptrVarIdx[l + 2]] == 1) {
          jointfactors_uiyx(
              datafactors, l + 2, n, nbrUi, r_old, uiyxfactors, ruiyx);
          // init variables for the optimization run
          factors1[0] = uiyxfactors[2];  // xu
          factors1[1] = uiyxfactors[0];  // u

          rt1[0] = ruiyx[2];  // xu
          rt1[1] = ruiyx[0];  // u

          sc_levels1 = r_old[0];      // x
          sc_levels2 = r_old[l + 2];  // u
          sc = 0.5 * (sc_levels_x - 1) * ruiyx[0];
          // optimization run on ptrVarIdx[l+2]
          I_x_u_part = optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l + 2]],
              data[ptrVarIdx[l + 2]], 2, factors1, rt1, sc, sc_levels1,
              sc_levels2, n, AllLevels[ptrVarIdx[l + 2]], cut[l + 2],
              &(r[l + 2]), environment);  // 2 factors //herve
        }
      }  // for all Uis
      for (int ll = 0; ll < nbrUi; ll++) {
        if (ptr_cnt[ptrVarIdx[ll + 2]] == 1)
          update_datafactors(
              sortidx, ptrVarIdx[ll + 2], datafactors, ll + 2, n, cut);
        r_old[ll + 2] = r[ll + 2];
      }
      U_counter++;
      if (nbrUi == 1) U_counter = max_U_counter;
    }  // U_counter loop

    jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
    r_temp[0] = r_old[0];
    r_temp[1] = ruiyx[0];
    r_temp[2] = ruiyx[2];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors[0], uiyxfactors[0], uiyxfactors[2],
          r_temp, n, c2terms, looklog, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors[0], uiyxfactors[0], uiyxfactors[2],
          r_temp, n, looklog, 0);
    I_x_u = res_temp[0];  // After optimization on U.
    Ik_x_u = res_temp[1];
    free(res_temp);
    // Reset cutpoints on U
    reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin,
        r, AllLevels, n);
    for (l = 0; l < nbrUi; l++) {
      if (ptr_cnt[ptrVarIdx[l + 2]] == 1)
        update_datafactors(
            sortidx, ptrVarIdx[l + 2], datafactors, l + 2, n, cut);
      r_old[l + 2] = r[l + 2];
    }
    // optimize I(y;u) over u
    U_counter = 0;
    while (U_counter < max_U_counter) {
      for (l = 0; l < nbrUi; l++) {
        if (ptr_cnt[ptrVarIdx[l + 2]] == 1) {
          jointfactors_uiyx(
              datafactors, l + 2, n, nbrUi, r_old, uiyxfactors, ruiyx);

          factors1[0] = uiyxfactors[1];  // yu
          factors1[1] = uiyxfactors[0];  // u

          rt1[0] = ruiyx[1];  // yu
          rt1[1] = ruiyx[0];  // u

          sc_levels1 = r_old[1];      // y
          sc_levels2 = r_old[l + 2];  // u
          sc = 0.5 * (sc_levels1 - 1) * ruiyx[0];
          // optimization run on ptrVarIdx[l+2]
          I_y_u_part = optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l + 2]],
              data[ptrVarIdx[l + 2]], 2, factors1, rt1, sc, sc_levels1,
              sc_levels2, n, AllLevels[ptrVarIdx[l + 2]], cut[l + 2],
              &(r[l + 2]), environment);  // 2 factors

          // update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n,
          // cut);
        }
      }  // for all Uis
      for (int ll = 0; ll < nbrUi; ll++) {
        if (ptr_cnt[ptrVarIdx[ll + 2]] == 1)
          update_datafactors(
              sortidx, ptrVarIdx[ll + 2], datafactors, ll + 2, n, cut);
        r_old[ll + 2] = r[ll + 2];
      }
      U_counter++;
      if (nbrUi == 1) U_counter = max_U_counter;
    }  // U_counter loop

    jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
    r_temp[0] = r_old[1];
    r_temp[1] = ruiyx[0];
    r_temp[2] = ruiyx[1];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors[1], uiyxfactors[0], uiyxfactors[1],
          r_temp, n, c2terms, looklog, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors[1], uiyxfactors[0], uiyxfactors[1],
          r_temp, n, looklog, 0);
    I_y_u = res_temp[0];  // After optimization on U.
    Ik_y_u = res_temp[1];
    free(res_temp);
    // Reset cutpoints on U
    reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin,
        r, AllLevels, n);
    for (l = 0; l < nbrUi; l++) {
      if (ptr_cnt[ptrVarIdx[l + 2]] == 1)
        update_datafactors(
            sortidx, ptrVarIdx[l + 2], datafactors, l + 2, n, cut);
      r_old[l + 2] = r[l + 2];
    }
    // Update X and Y
    if (ptr_cnt[ptrVarIdx[0]] == 1) {
      update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);
      r_old[0] = r[0];
    }
    if (ptr_cnt[ptrVarIdx[1]] == 1) {
      update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut);
      r_old[1] = r[1];
    }

    // Compute I(X;Y|U)
    cond_I = 0.5 * (I_x_yu - I_x_u + I_y_xu - I_y_u);
    cond_Ik = 0.5 * (Ik_x_yu - Ik_x_u + Ik_y_xu - Ik_y_u);
    cond_I_bis  = 0.5 * (I_x_yu_part[0] - I_x_u_part[0] + I_y_xu_part[0] - I_y_u_part[0])/n;
    cond_Ik_bis = 0.5 * (I_x_yu_part[1] - I_x_u_part[1] + I_y_xu_part[1] - I_y_u_part[1])/n;

    if (saveIterations) {
      for (j = 0; j < maxbins; j++) {
        for (l = 0; l < (nbrUi + 2); l++) {
          environment.iterative_cuts[stop1 - 1][j + l * maxbins] = cut[l][j];
        }
      }
    }
    // Test stop condition on stop1
    for (i = stop1 - 1; i > 0; i--) {
      // If no real improvement over last information
      if (fabs(cond_Ik_bis - MIk1[i]) < EPS) {
        flag1 = 1;
        Ik_av1 = MIk1[i];
        I_av1 = MI1[i];

        for (j = i + 1; j < stop1; j++) {
          Ik_av1 += MIk1[j];
          I_av1 += MI1[j];
        }
        Ik_av1 /= (stop1 - i);  // average over the periodic cycle
        I_av1 /= (stop1 - i);
        break;
      }
    }
    if (flag1) break;
    MIk1[stop1] = cond_Ik_bis;
    MI1[stop1] = cond_I_bis;

  }  // end stop1
  double* return_res = (double*)calloc(2, sizeof(double));
  if (flag1) {
    return_res[0] = I_av1;
    return_res[1] = Ik_av1;
  } else {
    return_res[0] = cond_I_bis;
    return_res[1] = cond_Ik_bis;
  }
  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if ((return_res[1] < 0) &&
      ((ptr_cnt[ptrVarIdx[0]] == 1) && (ptr_cnt[ptrVarIdx[1]] == 1))) {
    return_res[0] = 0;
    return_res[1] = 0;
  }

  if (saveIterations) {
    for (l = 0; l < (nbrUi + 2); l++) {
      environment.iterative_cuts[stop1][l * maxbins] =
          -1;  // mark where we stopped iterating
      environment.iterative_cuts[stop1][l * maxbins + 1] =
          100000 * return_res[1];  // pass Ik[X;Y|U]
      environment.iterative_cuts[stop1][l * maxbins + 2] =
          100000 * return_res[0];  // pass I[X;Y|U]
    }
  }

  free(r_temp);
  for (l = 0; l < (nbrUi + 2); l++) free(datafactors[l]);
  free(datafactors);
  for (l = 0; l < 4; l++) free(uiyxfactors[l]);
  free(uiyxfactors);
  free(factors1);
  free(rt1);
  free(r_old);
  free(ruiyx);
  free(MI);
  free(MIk);
  free(MI1);
  free(MIk1);
    //cond_Ik_bis = 0.5 * (I_x_yu_part[1] - I_x_u_part[1] + I_y_xu_part[1] - I_y_u_part[1])/n;

  return return_res;
}

double* compute_mi_cond_alg1(int** data, int** sortidx, int* AllLevels,
    int* ptr_cnt, int* ptrVarIdx, int nbrUi, int n, Environment& environment,
    bool saveIterations) {
  int maxbins = environment.maxbins;
  int initbins = environment.initbins;

  double* res = new double[3]();  // results res[0]->I,res[1]->I-k
  double* res_temp;  //=(double *)calloc(2,sizeof(double));//results
                     //res[0]->I,res[1]->I-k

  int j, l;

  int* r = (int*)calloc((nbrUi + 2), sizeof(int));

  int** cut;
  cut = (int**)calloc(nbrUi + 2, sizeof(int*));
  for (l = 0; l < (nbrUi + 2); l++) {
    cut[l] = (int*)calloc(maxbins, sizeof(int));
  }

  int lbin = floor(n / initbins);
  if (lbin < 1) {
    lbin = 1;
    initbins = n;
  }

  // no conditioning, empty set of variables in u
  // calling funcion compute_Ixy_alg1
  // NO u
  if (nbrUi == 0) {
    // initialization cut and r
    for (l = 0; l < (nbrUi + 2); l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        for (j = 0; j < initbins - 1; j++) {
          cut[l][j] = j * lbin + lbin - 1;
        }
        cut[l][initbins - 1] = n - 1;
        r[l] = initbins;
      } else {
        r[l] = AllLevels[ptrVarIdx[l]];
      }
    }

    res_temp = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, n,
        cut, r, environment, saveIterations);

    res[0] = n;
    res[1] = res_temp[0];
    res[2] = res_temp[0] - res_temp[1];

    free(res_temp);
    free(r);
    for (l = 0; l < (nbrUi + 2); l++) {
      free(cut[l]);
    }
    free(cut);

    return res;
  } else {  // with U
    // initialization cut and r
    for (l = 0; l < (nbrUi + 2); l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        for (j = 0; j < initbins - 1; j++) {
          cut[l][j] = j * lbin + lbin - 1;
        }
        cut[l][initbins - 1] = n - 1;
        r[l] = initbins;
      } else {
        r[l] = AllLevels[ptrVarIdx[l]];
      }
    }

    res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx,
        AllLevels, nbrUi, n, cut, r, lbin, environment, saveIterations);

    res[0] = n;
    res[1] = res_temp[0];
    res[2] = res_temp[0] - res_temp[1];

    free(res_temp);
    free(r);
    for (l = 0; l < (nbrUi + 2); l++) {
      free(cut[l]);
    }
    free(cut);

    return res;
  }
}

// compute Rscore and three point mutual information I(x;y;z | u)
// input x y z u-> compute Rscore and Ixyz
// Returns:
// res[0]=Rscore
// res[1]=N*Ixyz
// res[2]=N*kxyz
double* compute_Rscore_Ixyz_new_alg5(int** data, int** sortidx, int* AllLevels,
    int* ptr_cnt, int* ptrVarIdx, int nbrUi, int ptrZiIdx, int n,
    Environment& environment, bool saveIterations) {
  int maxbins = environment.maxbins;
  int initbins = environment.initbins;

  double Rscore;
  double nv, dpi, first, second, xz, yz;

  int j, l, ll;

  double I_xy_u, I_xy_zu;
  double Ik_xy_u, Ik_xz_u, Ik_yz_u, Ik_xy_zu;
  double I_xyz_u, Ik_xyz_u;

  double* res_temp =
      (double*)calloc(2, sizeof(double));  // results res[0]->I,res[1]->I-k
  double* res = new double[3]();           // results res[0]->I,res[1]->I-k

  int* ptrVarIdx_t = (int*)calloc((nbrUi + 2), sizeof(int));

  int* r = (int*)calloc((nbrUi + 3), sizeof(int));
  int** cut;
  cut = (int**)calloc(nbrUi + 3, sizeof(int*));
  for (l = 0; l < (nbrUi + 3); l++) {
    cut[l] = (int*)calloc(maxbins, sizeof(int));
  }

  int lbin = floor(n / initbins);
  if (lbin < 1) {
    lbin = 1;
    initbins = n;
  }

  // initialitize cuts vectors
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (j = 0; j < initbins - 1; j++) {
        cut[l][j] = j * lbin + lbin - 1;
      }
      cut[l][initbins - 1] = n - 1;
      r[l] = initbins;
    } else {
      r[l] = AllLevels[ptrVarIdx[l]];
    }
  }
  // z
  l = nbrUi + 2;
  if (ptr_cnt[ptrZiIdx] == 1) {
    for (j = 0; j < initbins - 1; j++) {
      cut[l][j] = j * lbin + lbin - 1;
    }
    cut[l][initbins - 1] = n - 1;
    r[l] = initbins;
  } else {
    r[l] = AllLevels[ptrZiIdx];
  }

  int** datafactors = (int**)calloc((nbrUi + 3), sizeof(int*));
  for (l = 0; l < (nbrUi + 3); l++) {
    datafactors[l] = (int*)calloc(n, sizeof(int));
  }
  int* ptr_u2_t = (int*)calloc(nbrUi + 2, sizeof(int));
  int** factors4_t = (int**)calloc(4, sizeof(int*));
  for (l = 0; l < 4; l++) {
    factors4_t[l] = (int*)calloc(n, sizeof(int));
  }
  int* r4_t = (int*)calloc(4, sizeof(int));
  int* r_temp = (int*)calloc(2, sizeof(int));

  // if opt
  int** cut_t = (int**)calloc(nbrUi + 2, sizeof(int*));
  int* r_t = (int*)calloc((nbrUi + 2), sizeof(int));
  for (l = 0; l < (nbrUi + 2); l++) {
    cut_t[l] = (int*)calloc(maxbins, sizeof(int));
  }

  // Optimize variables for each MI estimation for the R score

  // I(x,y|u,z)
  res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx,
      AllLevels, nbrUi + 1, n, cut, r, lbin, environment, saveIterations);
  I_xy_zu = res_temp[0];
  Ik_xy_zu = res_temp[1];
  free(res_temp);

  // I(x,y|u)
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (j = 0; j < initbins - 1; j++) {
        cut_t[l][j] = j * lbin + lbin - 1;
      }
      cut_t[l][initbins - 1] = n - 1;
      for (int j = initbins; j < maxbins; j++) {
        cut_t[l][j] = 0;
      }
      r_t[l] = initbins;
    } else {
      r_t[l] = AllLevels[ptrVarIdx[l]];
    }
  }
  // Do opt run on I(X;Y|U)
  if (nbrUi > 0) {
    res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx,
        AllLevels, nbrUi, n, cut_t, r_t, lbin, environment, saveIterations);
  } else {
    res_temp = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, n,
        cut_t, r_t, environment, saveIterations);
  }
  I_xy_u = res_temp[0];
  Ik_xy_u = res_temp[1];
  free(res_temp);

  // I(z,x|u)
  ptrVarIdx_t[0] = ptrVarIdx[0];          // X
  ptrVarIdx_t[1] = ptrVarIdx[nbrUi + 2];  // Z
  for (ll = 0; ll < nbrUi; ll++) ptrVarIdx_t[ll + 2] = ll + 2;
  // Reset cut
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx_t[l]] == 1) {
      for (j = 0; j < initbins - 1; j++) {
        cut_t[l][j] = j * lbin + lbin - 1;
      }
      cut_t[l][initbins - 1] = n - 1;
      for (int j = initbins; j < maxbins; j++) {
        cut_t[l][j] = 0;
      }
      r_t[l] = initbins;
    } else {
      r_t[l] = AllLevels[ptrVarIdx_t[l]];
    }
  }
  // Do opt run on I(X;Z|U)
  if (nbrUi > 0) {
    res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t,
        AllLevels, nbrUi, n, cut_t, r_t, lbin, environment, saveIterations);
  } else {
    res_temp = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t, AllLevels,
        n, cut_t, r_t, environment, saveIterations);
  }
  Ik_xz_u = res_temp[1];
  free(res_temp);

  // I(z,y|u)
  ptrVarIdx_t[0] = ptrVarIdx[1];          // Y
  ptrVarIdx_t[1] = ptrVarIdx[nbrUi + 2];  // Z
  for (ll = 0; ll < nbrUi; ll++) ptrVarIdx_t[ll + 2] = ll + 2;
  // Reset cut
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx_t[l]] == 1) {
      for (j = 0; j < initbins - 1; j++) {
        cut_t[l][j] = j * lbin + lbin - 1;
      }
      cut_t[l][initbins - 1] = n - 1;
      for (int j = initbins; j < maxbins; j++) {
        cut_t[l][j] = 0;
      }
      r_t[l] = initbins;
    } else {
      r_t[l] = AllLevels[ptrVarIdx_t[l]];
    }
  }
  // Do opt run on I(Y;Z|U)
  if (nbrUi > 0) {
    res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t,
        AllLevels, nbrUi, n, cut_t, r_t, lbin, environment, saveIterations);
  } else {
    res_temp = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t, AllLevels,
        n, cut_t, r_t, environment, saveIterations);
  }
  Ik_yz_u = res_temp[1];
  free(res_temp);

  // compute conditional three point mutual information
  I_xyz_u = I_xy_u - I_xy_zu;
  Ik_xyz_u = Ik_xy_u - Ik_xy_zu;

  // compute Rscore
  // compute probability of
  // not v-structure: nv
  // not dpi inequality : dpi

  // nv
  nv = n * Ik_xyz_u;
  xz = n * (Ik_xz_u - Ik_xy_u);
  yz = n * (Ik_yz_u - Ik_xy_u);

  if (xz < yz) {
    first = xz;
    second = yz;
  } else {
    first = yz;
    second = xz;
  }

  // dpi
  dpi = first - log1p(exp(first - second));

  if (dpi < nv) {
    // Pdpi>Pnv => Rscore=Pdpi
    Rscore = dpi;
  } else {
    // Pdpi<Pnv => Rscore=Pnv
    Rscore = nv;
  }

  res[0] = Rscore;
  res[1] = n * I_xyz_u;
  res[2] = n * I_xyz_u - nv;

  for (l = 0; l < (nbrUi + 3); l++) free(cut[l]);
  for (l = 0; l < (nbrUi + 2); l++) free(cut_t[l]);
  free(cut);
  free(cut_t);

  free(r);
  free(ptrVarIdx_t);

  for (l = 0; l < (nbrUi + 3); l++) free(datafactors[l]);
  free(datafactors);
  free(ptr_u2_t);
  for (l = 0; l < 4; l++) free(factors4_t[l]);
  free(factors4_t);
  free(r4_t);
  free(r_temp);
  free(r_t);

  return res;
}

}  // namespace computation
}  // namespace miic
