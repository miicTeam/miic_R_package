#include "computation_continuous.h"

#include <numeric>  // std::accumulate
#include <tuple>    // std::tie

#include "linear_allocator.h"
#include "mutual_information.h"
#include "structure.h"

#define STEPMAX 50
#define EPS 1e-5
#define FLAG_CPLX 1
//#define _MY_DEBUG_MInoU 1

namespace miic {
namespace computation {

using std::begin;
using std::end;
using std::min;
using std::vector;
using namespace miic::structure;
using namespace miic::utility;

namespace {

void reset_u_cutpoints(TempGrid2d<int>& cut, int nbrUi,
    const TempVector<int>& ptr_cnt, const TempVector<int>& ptrVarIdx,
    int init_nbin, int maxbins, int lbin, TempVector<int>& r,
    const TempVector<int>& AllLevels, int n) {
  for (int l = 2; l < nbrUi + 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] != 1) {
      r[l] = AllLevels[ptrVarIdx[l]];
    } else {
      for (int j = 0; j < init_nbin - 1; j++) {
        cut(l, j) = j * lbin + lbin - 1;
      }
      cut(l, init_nbin - 1) = n - 1;
      for (int j = init_nbin; j < maxbins; j++) {
        cut(l, j) = 0;
      }
      r[l] = min(init_nbin, AllLevels[ptrVarIdx[l]]);
    }
  }
}

// INPUT:
// memory_cuts: vector length n with recorded recursvely best cuts,
// possible values 0...n :
// 0->stops (one bins); -k -> stops (two bin ([0 k-1][k ..];
// k->continute [.. k-1][k ...])
// OUTPUT:
// r : number cuts
// cut: vector with cuts point-> [0 cut[0]][cut[0]+1 cut[1]]...[cut[r-2]
// cut[r-1]]
template <typename Ccut, typename = IsIntContainer<Ccut>>
int reconstruction_cut_coarse(const TempVector<int>& memory_cuts,
    const TempVector<int>& memory_cuts2, int np, int n, Ccut& cut) {
  int ncuts = 0;
  int l, s;
  if (memory_cuts[np - 1] == 0) {
    cut[0] = n - 1;
    return 1;
  }

  l = memory_cuts[np - 1];
  while (l > 0) {
    ncuts++;
    l = memory_cuts[l - 1];
  }
  if (l < 0) ncuts++;

  cut[ncuts] = n - 1;  // conventional last cut
  l = ncuts - 1;
  s = memory_cuts[np - 1];
  cut[l] = memory_cuts2[np - 1];  // truly last cut
  l--;
  while (s > 0 && l >= 0) {
    cut[l] = memory_cuts2[s - 1];
    s = memory_cuts[s - 1];
    l--;
  }
  return ncuts + 1;  // number of levels (r)
}

// GENERAL VARIABLES NOTATION

// int n: number of sample
// int coarse : coarse graining step : minimum length of a bin : considering
// only n/coarse possible cuts

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
template <typename Cf0, typename Cf1, typename Ccut,
    typename =
        void_t<IsIntContainer<Cf0>, IsIntContainer<Cf1>, IsIntContainer<Ccut>>>
void optfun_onerun_kmdl_coarse(const TempGrid2d<int>::ConstRow& sortidx_var,
    const TempGrid2d<int>::ConstRow& data, int nbrV, const Cf0& factors0,
    const Cf1& factors1, int r0, int r1, double sc, int sc_levels1,
    int previous_levels, int n, int nnr, Ccut&& cut, int* r_opt,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    Environment& environment) {
  TempAllocatorScope scope;

  int coarse = ceil(1.0 * nnr / environment.maxbins);  // step coarse graining
  if (coarse < 1) coarse = 1;
  int np = ceil(1.0 * nnr / coarse);  // number of possible cuts <= max_level

  // temp variables to memorize optimization cuts
  TempVector<int> memory_cuts_idx(np);  // indexes of the cuts (1..np)
  TempVector<int> memory_cuts_pos(np);  // positions of the cuts (1..n)

  // dynamic programming optimize function and memorize of cuts
  double Imax;   // I
  double Ikmax;  // I-k_NML

  // function max value at each step
  TempVector<double> I(np);  // The optimal information value found at each idx.
  TempVector<double> Ik(np);
  double I_kj;  // The information value for a bin fom idx k to j.
  double Ik_kj;

  int njforward(0), nkforward(0);  // Indexes at current position of j and k

  // entropy in kj interval for the <nbrV>+1 terms
  TempVector<double> H_kj(nbrV);
  TempVector<double> Hk_kj(nbrV);

  TempVector<int> counts_0(r0);
  TempVector<int> counts_1(r1);
  TempVector<int> counts_k_0(r0);
  TempVector<int> counts_k_1(r1);

  TempGrid2d<int> coarse_counts_marginal(np, r1, 0);
  TempGrid2d<int> coarse_counts_joint(np, r0, 0);

  int weighted_count;

  TempVector<int> check_repet(n);
  for (int i = 0; i < (n - 1); i++) {
    check_repet[i] = (data[sortidx_var[i + 1]] != data[sortidx_var[i]]);
  }

  TempVector<int> n_values(np);
  TempVector<double> sum_sample_weights(np);
  int ir(0);  // Iterator on non repeated values
  int level_marginal(0), level_joint(0);
  for (int i = 0; i < np; i++) {
    ir = 0;  // iterator on not repeated values
    if (flag_sample_weights && i > 0)
      sum_sample_weights[i] = sum_sample_weights[i - 1];

    while ((ir < coarse) && (njforward < n)) {
      level_marginal = factors1[sortidx_var[njforward]];
      level_joint = factors0[sortidx_var[njforward]];
      coarse_counts_marginal(i, level_marginal)++;
      coarse_counts_joint(i, level_joint)++;

      if (flag_sample_weights)
        sum_sample_weights[i] += sample_weights[njforward];
      if (njforward + 1 < n) {  // check no repetition
        ir += check_repet[njforward];
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
  for (int j = 0; j < np; j++) {  // j=1...n-1

    njforward = n_values[j];
    ef_nj = sum_sample_weights[j];
    if (flag_sample_weights) efN_factor = ef_nj / njforward;
    std::transform(counts_0.begin(), counts_0.end(),
        coarse_counts_joint.row_begin(j), counts_0.begin(), std::plus<int>());
    std::transform(counts_1.begin(), counts_1.end(),
        coarse_counts_marginal.row_begin(j), counts_1.begin(),
        std::plus<int>());

    Hk_kj[0] = 0;  // joint
    Hk_kj[1] = 0;  // marginal
    H_kj[0] = 0;   // joint
    H_kj[1] = 0;   // marginal
    for (int level = 0; level < r0; level++) {
      weighted_count = flag_sample_weights
                           ? int(efN_factor * counts_0[level] + 0.5)
                           : counts_0[level];
      Hk_kj[0] += environment.cache.cterm->getH(weighted_count);
    }
    H_kj[0] = Hk_kj[0];

    for (int level = 0; level < r1; level++) {
      weighted_count = flag_sample_weights
                           ? int(efN_factor * counts_1[level] + 0.5)
                           : counts_1[level];
      Hk_kj[1] -= environment.cache.cterm->getH(weighted_count);
      H_kj[1] -= environment.cache.cterm->getH(weighted_count);

      if (environment.cplx == 0 && counts_1[level] > 0)
        Hk_kj[1] -= sc * environment.cache.cterm->getLog(n);
      else if (environment.cplx == 1) {
        Hk_kj[1] -=
            environment.cache.cterm->getLogC(weighted_count, sc_levels1);
      }
    }

    I[j] = 0;
    Ik[j] = 0;
    for (int m = 0; m < nbrV; m++) {
      I[j] += H_kj[m];    // herve
      Ik[j] += Hk_kj[m];  // herve
    }

    Imax = -DBL_MAX;
    Ikmax = -DBL_MAX;

    counts_k_0 = counts_0;
    counts_k_1 = counts_1;
    ef_nk = ef_nj;

    // moving k

    // k iterator on possible cuts
    nkforward = 0;  // iterator on values
    // it iterator on not repeated vales

    // Before trying to create bins : solution is one single bin from 0 to j.
    memory_cuts_idx[j] = 0;
    memory_cuts_pos[j] = 0;

    for (int k = 0; k < j; k++) {  // k=1...n-2 possible cuts

      nkforward = n_values[k];
      ef_nk = sum_sample_weights[j] - sum_sample_weights[k];
      if (flag_sample_weights) efN_factor = ef_nk / (njforward - nkforward);

      std::transform(counts_k_1.begin(), counts_k_1.end(),
          coarse_counts_marginal.row_begin(k), counts_k_1.begin(),
          std::minus<int>());
      std::transform(counts_k_0.begin(), counts_k_0.end(),
          coarse_counts_joint.row_begin(k), counts_k_0.begin(),
          std::minus<int>());

      H_kj[0] = 0;
      Hk_kj[0] = 0;
      for (int level = 0; level < r0; level++) {
        weighted_count = flag_sample_weights
                             ? int(efN_factor * counts_k_0[level] + 0.5)
                             : counts_k_0[level];
        Hk_kj[0] += environment.cache.cterm->getH(weighted_count);
      }
      H_kj[0] = Hk_kj[0];

      H_kj[1] = 0;
      Hk_kj[1] = 0;
      for (int level = 0; level < r1; level++) {
        weighted_count = flag_sample_weights
                             ? int(efN_factor * counts_k_1[level] + 0.5)
                             : counts_k_1[level];
        Hk_kj[1] -= environment.cache.cterm->getH(weighted_count);
        H_kj[1] -= environment.cache.cterm->getH(weighted_count);

        if (environment.cplx == 0 && counts_k_1[level] > 0)
          Hk_kj[1] -= sc * environment.cache.cterm->getLog(n);
        else if (environment.cplx == 1) {
          Hk_kj[1] -=
              environment.cache.cterm->getLogC(weighted_count, sc_levels1);
        }
      }

      I_kj = 0;
      Ik_kj = 0;
      for (int m = 0; m < nbrV; m++) {
        I_kj += H_kj[m];
        Ik_kj += Hk_kj[m];
      }
      if (environment.cplx == 1) {
        // Combinatorial approximation
        Ik_kj -=
            environment.cache.cterm->getLogChoose(np - 1, previous_levels - 1) /
            (previous_levels - 1);
      }

      if ((Ik[k] + Ik_kj) > Ik[j]) {
        I_newbin = I[k] + I_kj;     //[0.. cuts.. k-1][k j] //herve
        Ik_newbin = Ik[k] + Ik_kj;  //[0.. cuts.. k-1][k j] //herve
        if (Ikmax < Ik_newbin) {
          Imax = I_newbin;
          Ikmax = Ik_newbin;
          I[j] = Imax;    // optimized function for the interval [0 j]
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
}

// compute I(x,y)
// dynamic programming for optimizing variables binning
// optimize on x I(x,y): Hx - Hxy - kmdl
// optimize on y I(x,y): Hy - Hxy - kmdl
// until convergence
vector<double> compute_Ixy_alg1(const TempGrid2d<int>& data,
    const TempGrid2d<int>& sortidx, const TempVector<int>& ptr_cnt,
    const TempVector<int>& ptrVarIdx, const TempVector<int>& AllLevels, int n,
    TempGrid2d<int>& cut, TempVector<int>& r,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    Environment& environment, bool saveIterations) {
  TempAllocatorScope scope;

  int maxbins = environment.maxbins;
  int initbins = environment.initbins;
  int cplx = environment.cplx;
  double n_eff =
      std::accumulate(sample_weights.begin(), sample_weights.end(), 0.0);

  int j, l;

  // res_tempults res_temp[0]->I,res_temp[1]->I-k
  vector<double> res_temp;

  // allocation factors  x y
  TempGrid2d<int> datafactors(2, n);

  TempVector<int> r_temp(3);
  TempVector<int> xy_factors(n);
  int rxy = 0;
  // initialization of datafactors && sortidx
  for (size_t l = 0; l < 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      update_datafactors(sortidx.getConstRow(ptrVarIdx[l]),
          datafactors.getRow(l), cut.getRow(l));
    } else {
      for (j = 0; j <= n - 1; j++) {
        datafactors(l, j) = data(ptrVarIdx[l], j);
      }
    }
  }
  jointfactors_u(datafactors, r, xy_factors, rxy);
  r_temp[0] = r[0];
  r_temp[1] = r[1];
  r_temp[2] = rxy;
  if (cplx == 1)
    res_temp =
        computeMI_knml(datafactors.getRow(0), datafactors.getRow(1), xy_factors,
            r_temp, n, n_eff, sample_weights, environment.cache.cterm, 0);
  else
    res_temp = computeMI_kmdl(datafactors.getRow(0), datafactors.getRow(1),
        xy_factors, r_temp, n, environment.cache.cterm, 0);

  // all discrete
  if (ptr_cnt[ptrVarIdx[0]] == 0 && ptr_cnt[ptrVarIdx[1]] == 0) return res_temp;

  // Find the best initial conditions with the same number of bins (equalfreq)
  // on all continuous variables.
  double max_res = -DBL_MAX;
  int max_initbins = initbins;
  int min_unique_values = n;
  for (l = 0; l < 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1)
      min_unique_values = min(min_unique_values, AllLevels[ptrVarIdx[l]]);
  }
  for (int new_initbins = 2; (new_initbins < initbins) && (new_initbins < 20) &&
                             (new_initbins < min_unique_values);
       new_initbins++) {
    int lbin = n / new_initbins;
    if (lbin < 1) {
      lbin = 1;
      new_initbins = n;
    }
    // Reinitialization cut and r
    for (l = 0; l < 2; l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        for (j = 0; j < new_initbins - 1; j++) {
          cut(l, j) = j * lbin + lbin - 1;
        }
        cut(l, new_initbins - 1) = n - 1;
        r[l] = new_initbins;
      } else {
        r[l] = AllLevels[ptrVarIdx[l]];
      }
    }
    // initialization of datafactors && sortidx
    for (l = 0; l < 2; l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        update_datafactors(sortidx.getConstRow(ptrVarIdx[l]),
            datafactors.getRow(l), cut.getRow(l));
      } else {
        for (j = 0; j <= n - 1; j++) {
          datafactors(l, j) = data(ptrVarIdx[l], j);
        }
      }
    }

    jointfactors_u(datafactors, r, xy_factors, rxy);
    r_temp[0] = r[0];
    r_temp[1] = r[1];
    r_temp[2] = rxy;
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), datafactors.getRow(1),
          xy_factors, r_temp, n, n_eff, sample_weights, environment.cache.cterm,
          0);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), datafactors.getRow(1),
          xy_factors, r_temp, n, environment.cache.cterm, 0);

    if (res_temp[1] > max_res) {
      max_initbins = new_initbins;
      max_res = res_temp[1];
    }
  }

  int lbin = n / max_initbins;
  if (lbin < 1) {
    lbin = 1;
    max_initbins = n;
  }
  // Reinitialization cut and r
  for (l = 0; l < 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (j = 0; j < max_initbins - 1; j++) {
        cut(l, j) = j * lbin + lbin - 1;
      }
      cut(l, max_initbins - 1) = n - 1;
      r[l] = max_initbins;
    } else {
      r[l] = AllLevels[ptrVarIdx[l]];
    }
  }
  // initialization of datafactors && sortidx
  for (l = 0; l < 2; l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      update_datafactors(sortidx.getConstRow(ptrVarIdx[l]),
          datafactors.getRow(l), cut.getRow(l));
    } else {
      for (j = 0; j <= n - 1; j++) {
        datafactors(l, j) = data(ptrVarIdx[l], j);
      }
    }
  }
  // Run dynamic optimization with the best initial conditions.
  double sc;
  TempVector<double> MI(STEPMAX);
  TempVector<double> MIk(STEPMAX);
  double I_av, Ik_av;

  int stop, i, flag;
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
      // optimize I(x;y) on x
      int r0 = ry;  // y
      int r1 = 1;   // One single bin at start

      sc_levels_y = ry;
      sc_levels1 = sc_levels_y;
      sc_levels2 = sc_levels_x;
      sc = 0.5 * (sc_levels_y - 1);
      rx = r[0];
      // Optimization run on X. 2 factors
      optfun_onerun_kmdl_coarse(sortidx.getConstRow(ptrVarIdx[0]),
          data.getConstRow(ptrVarIdx[0]), 2, datafactors.getRow(1),
          TempVector<int>(n, 0), r0, r1, sc, sc_levels1, sc_levels2, n,
          AllLevels[ptrVarIdx[0]], cut.getRow(0), &(r[0]), sample_weights,
          flag_sample_weights, environment);
    }
    if (ptr_cnt[ptrVarIdx[1]] == 1) {
      // optimize I(x;y) on y
      int r0 = rx;  // x before its optimization
      int r1 = 1;   // One single bin at start

      sc_levels_x = rx;
      sc_levels1 = sc_levels_x;
      sc_levels2 = sc_levels_y;
      sc = 0.5 * (sc_levels_x - 1);

      ry = r[1];
      // Optimization run on Y. 2 factors
      optfun_onerun_kmdl_coarse(sortidx.getConstRow(ptrVarIdx[1]),
          data.getConstRow(ptrVarIdx[1]), 2, datafactors.getRow(0),
          TempVector<int>(n, 0), r0, r1, sc, sc_levels1, sc_levels2, n,
          AllLevels[ptrVarIdx[1]], cut.getRow(1), &(r[1]), sample_weights,
          flag_sample_weights, environment);
    }

    // update both datafactors
    if (ptr_cnt[ptrVarIdx[0]] == 1) {
      update_datafactors(sortidx.getConstRow(ptrVarIdx[0]),
          datafactors.getRow(0), cut.getRow(0));
      rx = r[0];
    }
    if (ptr_cnt[ptrVarIdx[1]] == 1) {
      update_datafactors(sortidx.getConstRow(ptrVarIdx[1]),
          datafactors.getRow(1), cut.getRow(1));
      ry = r[1];
    }

    if (saveIterations) {
      // Save cut points
      for (j = 0; j < maxbins; j++) {
        environment.iterative_cuts(stop - 1, j) = cut(0, j);
        environment.iterative_cuts(stop - 1, j + maxbins) = cut(1, j);
      }
    }
    jointfactors_u(datafactors, r, xy_factors, rxy);
    r_temp[0] = r[0];
    r_temp[1] = r[1];
    r_temp[2] = rxy;
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), datafactors.getRow(1),
          xy_factors, r_temp, n, n_eff, sample_weights, environment.cache.cterm,
          0);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), datafactors.getRow(1),
          xy_factors, r_temp, n, environment.cache.cterm, 0);
    // Adding combinatorial term
    if (ptr_cnt[ptrVarIdx[0]] == 1 && rx > 1) {
      np = min(maxbins, AllLevels[ptrVarIdx[0]]);
      if (rx < np)
        res_temp[1] -=
            environment.cache.cterm->getLogChoose(np - 1, rx - 1) / n;
    }
    if (ptr_cnt[ptrVarIdx[1]] == 1 && ry > 1) {
      np = min(maxbins, AllLevels[ptrVarIdx[1]]);
      if (ry < np)
        res_temp[1] -=
            environment.cache.cterm->getLogChoose(np - 1, ry - 1) / n;
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

    if (flag || (ptr_cnt[ptrVarIdx[0]] == 0) || (ptr_cnt[ptrVarIdx[1]] == 0)) {
      break;
    }
  }  // for

  vector<double> return_res(2);
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
    // mark where we stopped iterating
    environment.iterative_cuts(stop, 0) = -1;
    // Pass Ik[X;Y]
    environment.iterative_cuts(stop, 1) = 100000 * return_res[1];
    // Pass I[X;Y]
    environment.iterative_cuts(stop, 2) = 100000 * return_res[0];
    // Pass max res before optimization with equal freq
    environment.iterative_cuts(stop, 3) = 100000 * max_res;
    environment.iterative_cuts(stop, maxbins) = -1;
  }

#if _MY_DEBUG_MInoU
  Rcout << (ptr_cnt[ptrVarIdx[0]]) << " " << (ptr_cnt[ptrVarIdx[1]]) << endl;
  Rprintf("    final: I_xy=%lf Ik_xy=%lf\n", return_res[0], return_res[1]);

  Rprintf("    0 : r=%d ", r[0]);
  for (i = 0; i < r[0]; i++) {
    Rprintf("    %d ", cut(0, i));
  }
  Rprintf("\n");
  Rprintf("    1 : r=%d ", r[1]);
  for (i = 0; i < r[1]; i++) {
    Rprintf("    %d ", cut(1, i));
  }
  Rprintf("\n");
#endif
  return return_res;
}

vector<double> compute_Ixy_cond_u_new_alg1(const TempGrid2d<int>& data,
    const TempGrid2d<int>& sortidx, const TempVector<int>& ptr_cnt,
    const TempVector<int>& ptrVarIdx, const TempVector<int>& AllLevels,
    int nbrUi, int n, TempGrid2d<int>& cut, TempVector<int>& r, int lbin,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    Environment& environment, bool saveIterations) {
  TempAllocatorScope scope;

  int maxbins = environment.maxbins;
  int initbins = environment.initbins;
  int cplx = environment.cplx;
  double n_eff =
      std::accumulate(sample_weights.begin(), sample_weights.end(), 0.0);

  int j, l;
  int STEPMAX1 = 50;
  // res_tempults res_temp[0]->I,res_temp[1]->I-k
  vector<double> res_temp;
  // allocation factors  x y
  TempGrid2d<int> datafactors(nbrUi + 2, n);

  TempVector<int> r_temp(3);

  // initialization of datafactors && sortidx
  for (l = 0; l < (nbrUi + 2); l++) {
    // compute datafactors based on the positions of cut points in vector <cut>
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      update_datafactors(sortidx.getConstRow(ptrVarIdx[l]),
          datafactors.getRow(l), cut.getRow(l));
    } else {  // discrete case
      for (j = 0; j <= n - 1; j++) {
        datafactors(l, j) = data(ptrVarIdx[l], j);
      }
    }
  }

  TempGrid2d<int> uiyxfactors(4, n);  //({ui},{uiy}),{uix},{uiyx})
  TempVector<int> ruiyx(4);

  double sc;

  TempVector<double> MI1(STEPMAX1);
  TempVector<double> MIk1(STEPMAX1);
  double I_av1, Ik_av1;

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

  // Complexity factor (# of levels) passed to each optimization
  int sc_levels1, sc_levels2;

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
          cut(l, j) = j * lbin + lbin - 1;
        }
        cut(l, new_initbins - 1) = n - 1;
        r[l] = new_initbins;
      } else {
        r[l] = AllLevels[ptrVarIdx[l]];
      }
    }

    // compute datafactors based on the positions of cut points in vector <cut>
    for (l = 0; l < (nbrUi + 2); l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        update_datafactors(sortidx.getConstRow(ptrVarIdx[l]),
            datafactors.getRow(l), cut.getRow(l));
      } else {  // discrete case
        for (j = 0; j <= n - 1; j++) {
          datafactors(l, j) = data(ptrVarIdx[l], j);
        }
      }
    }

    jointfactors_uiyx(datafactors, r, -1, uiyxfactors, ruiyx);
    r_temp[0] = r[1];
    r_temp[1] = ruiyx[2];
    r_temp[2] = ruiyx[3];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(1), uiyxfactors.getRow(2),
          uiyxfactors.getRow(3), r_temp, n, n_eff, sample_weights,
          environment.cache.cterm, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(1), uiyxfactors.getRow(2),
          uiyxfactors.getRow(3), r_temp, n, environment.cache.cterm, 0);
    I_y_xu = res_temp[0];  // Before optimization on X.
    Ik_y_xu = res_temp[1];

    r_temp[0] = r[0];
    r_temp[1] = ruiyx[1];
    r_temp[2] = ruiyx[3];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), uiyxfactors.getRow(1),
          uiyxfactors.getRow(3), r_temp, n, n_eff, sample_weights,
          environment.cache.cterm, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), uiyxfactors.getRow(1),
          uiyxfactors.getRow(3), r_temp, n, environment.cache.cterm, 0);
    I_x_yu = res_temp[0];  // Before updating Y (and X).
    Ik_x_yu = res_temp[1];

    if ((Ik_y_xu + Ik_x_yu) > max_res) {
      max_initbins = new_initbins;
      max_res = (Ik_y_xu + Ik_x_yu);
    }
  }

  lbin = n / max_initbins;
  if (lbin < 1) {
    lbin = 1;
    max_initbins = n;
  }

  // initialization cut and r
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (j = 0; j < max_initbins - 1; j++) {
        cut(l, j) = j * lbin + lbin - 1;
      }
      cut(l, max_initbins - 1) = n - 1;
      r[l] = max_initbins;
    } else {
      r[l] = AllLevels[ptrVarIdx[l]];
    }
  }
  // compute datafactors based on the positions of cut points in vector <cut>
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      update_datafactors(sortidx.getConstRow(ptrVarIdx[l]),
          datafactors.getRow(l), cut.getRow(l));
    } else {
      for (j = 0; j <= n - 1; j++) {  // discrete case
        datafactors(l, j) = data(ptrVarIdx[l], j);
      }
    }
  }
  // Run optimization with best initial equal freq.
  TempVector<int> r_old(r);  // copy
  int U_counter;

  int sc_levels_x{0};  // Number of levels of the first variable
  int sc_levels_y{0};  // Number of levels of the second variable
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
          jointfactors_uiyx(datafactors, r_old, l + 2, uiyxfactors, ruiyx);
          // init variables for the optimization run
          int r0 = ruiyx[3];  // xyu
          int r1 = ruiyx[2];  // xu
          sc_levels_x = r_old[0];
          sc_levels_y = r_old[1];

          sc = 0.5 * (sc_levels_y - 1) * ruiyx[2];

          sc_levels1 = sc_levels_y;
          sc_levels2 = r_old[l + 2];  // old nlevels for combinatorial term

          // Run optimization on U. 2 factors xyu and xu
          optfun_onerun_kmdl_coarse(sortidx.getConstRow(ptrVarIdx[l + 2]),
              data.getConstRow(ptrVarIdx[l + 2]), 2, uiyxfactors.getRow(3),
              uiyxfactors.getRow(2), r0, r1, sc, sc_levels1, sc_levels2, n,
              AllLevels[ptrVarIdx[l + 2]], cut.getRow(l + 2), &(r[l + 2]),
              sample_weights, flag_sample_weights, environment);
        }
      }  // for all Uis
      for (int ll = 0; ll < nbrUi; ll++) {
        if (ptr_cnt[ptrVarIdx[ll + 2]] == 1) {
          update_datafactors(sortidx.getConstRow(ptrVarIdx[ll + 2]),
              datafactors.getRow(ll + 2), cut.getRow(ll + 2));
          r_old[ll + 2] = r[ll + 2];
        }
      }
      U_counter++;
      if (nbrUi == 1) U_counter = max_U_counter;
    }  // U_counter loop

    jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
    r_temp[0] = r_old[1];
    r_temp[1] = ruiyx[2];
    r_temp[2] = ruiyx[3];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(1), uiyxfactors.getRow(2),
          uiyxfactors.getRow(3), r_temp, n, n_eff, sample_weights,
          environment.cache.cterm, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(1), uiyxfactors.getRow(2),
          uiyxfactors.getRow(3), r_temp, n, environment.cache.cterm, 0);
    I_y_xu = res_temp[0];  // Before optimization on X.
    Ik_y_xu = res_temp[1];
    if ((ptr_cnt[ptrVarIdx[0]] == 1) && (r_old[0] > 1)) {
      np = min(AllLevels[ptrVarIdx[0]], maxbins);
      if (r_old[0] < np) {
        Ik_y_xu -=
            environment.cache.cterm->getLogChoose(np - 1, r_old[0] - 1) / n;
      }
    }

    if (ptr_cnt[ptrVarIdx[0]] == 1) {
      // opt x
      // I(y;xu)
      // compute joint factors u yu xu xyu
      jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
      // init variables for the optimization run
      int r0 = ruiyx[1];  // uy
      int r1 = ruiyx[0];  // u
      sc_levels_x = r_old[0];
      sc_levels_y = r_old[1];
      sc = 0.5 * (sc_levels_y - 1) * ruiyx[0];

#if _MY_DEBUG_MInoU
      Rprintf("start optfun\n ");
      R_FlushConsole();
#endif
      sc_levels1 = sc_levels_y;  // herve
      sc_levels2 = sc_levels_x;  // herve
      // Run optimization on X. 2 factors uy and u
      optfun_onerun_kmdl_coarse(sortidx.getConstRow(ptrVarIdx[0]),
          data.getConstRow(ptrVarIdx[0]), 2, uiyxfactors.getRow(1),
          uiyxfactors.getRow(0), r0, r1, sc, sc_levels1, sc_levels2, n,
          AllLevels[ptrVarIdx[0]], cut.getRow(0), &(r[0]), sample_weights,
          flag_sample_weights, environment);
    }

    // Reset cutpoints on U
    reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin,
        r, AllLevels, n);
    for (l = 0; l < nbrUi; l++) {
      if (ptr_cnt[ptrVarIdx[l + 2]] == 1)
        update_datafactors(sortidx.getConstRow(ptrVarIdx[l + 2]),
            datafactors.getRow(l + 2), cut.getRow(l + 2));
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
          jointfactors_uiyx(datafactors, r_old, l + 2, uiyxfactors, ruiyx);
          // init variables for the optimization run
          int r0 = ruiyx[3];  // xyu
          int r1 = ruiyx[1];  // yu
          sc_levels_y = r_old[1];
          sc_levels_x = r_old[0];

          sc_levels1 = sc_levels_x;   // herve
          sc_levels2 = r_old[l + 2];  // herve
          sc = 0.5 * (sc_levels_x - 1) * ruiyx[1];

          // Run optimization on U. 2 factors xyu and yu
          optfun_onerun_kmdl_coarse(sortidx.getConstRow(ptrVarIdx[l + 2]),
              data.getConstRow(ptrVarIdx[l + 2]), 2, uiyxfactors.getRow(3),
              uiyxfactors.getRow(1), r0, r1, sc, sc_levels1, sc_levels2, n,
              AllLevels[ptrVarIdx[l + 2]], cut.getRow(l + 2), &(r[l + 2]),
              sample_weights, flag_sample_weights, environment);
        }
      }  // for all Uis
      for (int ll = 0; ll < nbrUi; ll++) {
        if (ptr_cnt[ptrVarIdx[ll + 2]] == 1) {
          update_datafactors(sortidx.getConstRow(ptrVarIdx[ll + 2]),
              datafactors.getRow(ll + 2), cut.getRow(ll + 2));
          r_old[ll + 2] = r[ll + 2];
        }
      }
      U_counter++;
      if (nbrUi == 1) U_counter = max_U_counter;
    }  // U_counter loop

    jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
    r_temp[0] = r_old[0];
    r_temp[1] = ruiyx[1];
    r_temp[2] = ruiyx[3];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), uiyxfactors.getRow(1),
          uiyxfactors.getRow(3), r_temp, n, n_eff, sample_weights,
          environment.cache.cterm, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), uiyxfactors.getRow(1),
          uiyxfactors.getRow(3), r_temp, n, environment.cache.cterm, 0);
    I_x_yu = res_temp[0];  // Before updating Y (and X).
    Ik_x_yu = res_temp[1];
    if ((ptr_cnt[ptrVarIdx[1]] == 1) && (r_old[1] > 1)) {
      np = min(AllLevels[ptrVarIdx[1]], maxbins);
      if (r_old[1] < np) {
        Ik_x_yu -=
            environment.cache.cterm->getLogChoose(np - 1, r_old[1] - 1) / n;
      }
    }

    if (ptr_cnt[ptrVarIdx[1]] == 1) {
      // optimize on y
      // I(x;yu)
      jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
      // init variables for the optimization run
      int r0 = ruiyx[2];  // ux
      int r1 = ruiyx[0];  // u

      sc_levels_x = r_old[0];
      sc_levels_y = r_old[1];

      sc_levels1 = sc_levels_x;  // herve
      sc_levels2 = sc_levels_y;  // herve
      sc = 0.5 * (sc_levels_x - 1) * ruiyx[0];
      // Run optimization on Y. 2 factors ux and u
      optfun_onerun_kmdl_coarse(sortidx.getConstRow(ptrVarIdx[1]),
          data.getConstRow(ptrVarIdx[1]), 2, uiyxfactors.getRow(2),
          uiyxfactors.getRow(0), r0, r1, sc, sc_levels1, sc_levels2, n,
          AllLevels[ptrVarIdx[1]], cut.getRow(1), &(r[1]), sample_weights,
          flag_sample_weights, environment);
    }
    // Reset cutpoints on U
    reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin,
        r, AllLevels, n);
    for (l = 0; l < nbrUi; l++) {
      if (ptr_cnt[ptrVarIdx[l + 2]] == 1)
        update_datafactors(sortidx.getConstRow(ptrVarIdx[l + 2]),
            datafactors.getRow(l + 2), cut.getRow(l + 2));
      r_old[l + 2] = r[l + 2];
    }
    // optimize I(x;u) over u
    U_counter = 0;
    while (U_counter < max_U_counter) {
      for (l = 0; l < nbrUi; l++) {
        if (ptr_cnt[ptrVarIdx[l + 2]] == 1) {
          jointfactors_uiyx(datafactors, r_old, l + 2, uiyxfactors, ruiyx);
          // init variables for the optimization run
          int r0 = ruiyx[2];  // xu
          int r1 = ruiyx[0];  // u

          sc_levels1 = r_old[0];      // x
          sc_levels2 = r_old[l + 2];  // u
          sc = 0.5 * (sc_levels_x - 1) * ruiyx[0];
          // optimization run on ptrVarIdx[l+2], 2 factors xu and u
          optfun_onerun_kmdl_coarse(sortidx.getConstRow(ptrVarIdx[l + 2]),
              data.getConstRow(ptrVarIdx[l + 2]), 2, uiyxfactors.getRow(2),
              uiyxfactors.getRow(0), r0, r1, sc, sc_levels1, sc_levels2, n,
              AllLevels[ptrVarIdx[l + 2]], cut.getRow(l + 2), &(r[l + 2]),
              sample_weights, flag_sample_weights, environment);
        }
      }  // for all Uis
      for (int ll = 0; ll < nbrUi; ll++) {
        if (ptr_cnt[ptrVarIdx[ll + 2]] == 1) {
          update_datafactors(sortidx.getConstRow(ptrVarIdx[ll + 2]),
              datafactors.getRow(ll + 2), cut.getRow(ll + 2));
          r_old[ll + 2] = r[ll + 2];
        }
      }
      U_counter++;
      if (nbrUi == 1) U_counter = max_U_counter;
    }  // U_counter loop

    jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
    r_temp[0] = r_old[0];
    r_temp[1] = ruiyx[0];
    r_temp[2] = ruiyx[2];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), uiyxfactors.getRow(0),
          uiyxfactors.getRow(2), r_temp, n, n_eff, sample_weights,
          environment.cache.cterm, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), uiyxfactors.getRow(0),
          uiyxfactors.getRow(2), r_temp, n, environment.cache.cterm, 0);
    I_x_u = res_temp[0];  // After optimization on U.
    Ik_x_u = res_temp[1];
    // Reset cutpoints on U
    reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin,
        r, AllLevels, n);
    for (l = 0; l < nbrUi; l++) {
      if (ptr_cnt[ptrVarIdx[l + 2]] == 1)
        update_datafactors(sortidx.getConstRow(ptrVarIdx[l + 2]),
            datafactors.getRow(l + 2), cut.getRow(l + 2));
      r_old[l + 2] = r[l + 2];
    }
    // optimize I(y;u) over u
    U_counter = 0;
    while (U_counter < max_U_counter) {
      for (l = 0; l < nbrUi; l++) {
        if (ptr_cnt[ptrVarIdx[l + 2]] == 1) {
          jointfactors_uiyx(datafactors, r_old, l + 2, uiyxfactors, ruiyx);
          int r0 = ruiyx[1];  // yu
          int r1 = ruiyx[0];  // u

          sc_levels1 = r_old[1];      // y
          sc_levels2 = r_old[l + 2];  // u
          sc = 0.5 * (sc_levels1 - 1) * ruiyx[0];
          // optimization run on ptrVarIdx[l+2], 2 factors yu and u
          optfun_onerun_kmdl_coarse(sortidx.getConstRow(ptrVarIdx[l + 2]),
              data.getConstRow(ptrVarIdx[l + 2]), 2, uiyxfactors.getRow(1),
              uiyxfactors.getRow(0), r0, r1, sc, sc_levels1, sc_levels2, n,
              AllLevels[ptrVarIdx[l + 2]], cut.getRow(l + 2), &(r[l + 2]),
              sample_weights, flag_sample_weights, environment);
        }
      }  // for all Uis
      for (int ll = 0; ll < nbrUi; ll++) {
        if (ptr_cnt[ptrVarIdx[ll + 2]] == 1) {
          update_datafactors(sortidx.getConstRow(ptrVarIdx[ll + 2]),
              datafactors.getRow(ll + 2), cut.getRow(ll + 2));
          r_old[ll + 2] = r[ll + 2];
        }
      }
      U_counter++;
      if (nbrUi == 1) U_counter = max_U_counter;
    }  // U_counter loop

    jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
    r_temp[0] = r_old[1];
    r_temp[1] = ruiyx[0];
    r_temp[2] = ruiyx[1];
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(1), uiyxfactors.getRow(0),
          uiyxfactors.getRow(1), r_temp, n, n_eff, sample_weights,
          environment.cache.cterm, FLAG_CPLX);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(1), uiyxfactors.getRow(0),
          uiyxfactors.getRow(1), r_temp, n, environment.cache.cterm, 0);
    I_y_u = res_temp[0];  // After optimization on U.
    Ik_y_u = res_temp[1];
    // Reset cutpoints on U
    reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin,
        r, AllLevels, n);
    for (l = 0; l < nbrUi; l++) {
      if (ptr_cnt[ptrVarIdx[l + 2]] == 1)
        update_datafactors(sortidx.getConstRow(ptrVarIdx[l + 2]),
            datafactors.getRow(l + 2), cut.getRow(l + 2));
      r_old[l + 2] = r[l + 2];
    }
    // Update X and Y
    if (ptr_cnt[ptrVarIdx[0]] == 1) {
      update_datafactors(sortidx.getConstRow(ptrVarIdx[0]),
          datafactors.getRow(0), cut.getRow(0));
      r_old[0] = r[0];
    }
    if (ptr_cnt[ptrVarIdx[1]] == 1) {
      update_datafactors(sortidx.getConstRow(ptrVarIdx[1]),
          datafactors.getRow(1), cut.getRow(1));
      r_old[1] = r[1];
    }

    // Compute I(X;Y|U)
    cond_I = 0.5 * (I_x_yu - I_x_u + I_y_xu - I_y_u);
    cond_Ik = 0.5 * (Ik_x_yu - Ik_x_u + Ik_y_xu - Ik_y_u);

    if (saveIterations) {
      for (j = 0; j < maxbins; j++) {
        for (l = 0; l < (nbrUi + 2); l++) {
          environment.iterative_cuts(stop1 - 1, j + l * maxbins) = cut(l, j);
        }
      }
    }
    // Test stop condition on stop1
    for (i = stop1 - 1; i > 0; i--) {
      // If no real improvement over last information
      if (fabs(cond_Ik - MIk1[i]) < EPS) {
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
    MIk1[stop1] = cond_Ik;
    MI1[stop1] = cond_I;

  }  // end stop1
  vector<double> return_res(2);
  if (flag1) {
    return_res[0] = I_av1;
    return_res[1] = Ik_av1;
  } else {
    return_res[0] = cond_I;
    return_res[1] = cond_Ik;
  }
  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if ((return_res[1] < 0) &&
      ((ptr_cnt[ptrVarIdx[0]] == 1) && (ptr_cnt[ptrVarIdx[1]] == 1))) {
    return_res[0] = 0;
    return_res[1] = 0;
  }
  if (saveIterations) {
    for (l = 0; l < (nbrUi + 2); l++) {
      environment.iterative_cuts(stop1, l * maxbins) =
          -1;  // mark where we stopped iterating
      environment.iterative_cuts(stop1, l * maxbins + 1) =
          100000 * return_res[1];  // pass Ik[X;Y|U]
      environment.iterative_cuts(stop1, l * maxbins + 2) =
          100000 * return_res[0];  // pass I[X;Y|U]
    }
  }

  return return_res;
}

}  // anonymous namespace

InfoBlock computeCondMutualInfo(const TempGrid2d<int>& data,
    const TempGrid2d<int>& sortidx, const TempVector<int>& AllLevels,
    const TempVector<int>& ptr_cnt, const TempVector<int>& ptrVarIdx, int nbrUi,
    int n, const TempVector<double>& sample_weights, bool flag_sample_weights,
    Environment& environment, bool saveIterations) {
  TempAllocatorScope scope;

  int maxbins = environment.maxbins;
  int initbins = environment.initbins;

  vector<double> res_temp;      // res_temp[0]->I,res_temp[1]->I-k

  int j, l;

  TempVector<int> r(nbrUi + 2);
  TempGrid2d<int> cut(nbrUi + 2, maxbins);

  int lbin = n / initbins;
  if (lbin < 1) {
    lbin = 1;
    initbins = n;
  }
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1)
      initbins = min(initbins, AllLevels[ptrVarIdx[l]]);
  }

  // no conditioning, empty set of variables in u
  // calling funcion compute_Ixy_alg1
  // NO u
  if (nbrUi == 0) {
    // initialization cut and r
    for (l = 0; l < (nbrUi + 2); l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        for (j = 0; j < initbins - 1; j++) {
          cut(l, j) = j * lbin + lbin - 1;
        }
        cut(l, initbins - 1) = n - 1;
        r[l] = initbins;
      } else {
        r[l] = AllLevels[ptrVarIdx[l]];
      }
    }

    res_temp = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, n,
        cut, r, sample_weights, flag_sample_weights, environment,
        saveIterations);

    double Ixy_ui = n * res_temp[0];
    double kxy_ui = n * (res_temp[0] - res_temp[1]);
    return InfoBlock(n, Ixy_ui, kxy_ui);
  } else {  // with U
    // initialization cut and r
    for (l = 0; l < (nbrUi + 2); l++) {
      if (ptr_cnt[ptrVarIdx[l]] == 1) {
        for (j = 0; j < initbins - 1; j++) {
          cut(l, j) = j * lbin + lbin - 1;
        }
        cut(l, initbins - 1) = n - 1;
        r[l] = initbins;
      } else {
        r[l] = AllLevels[ptrVarIdx[l]];
      }
    }

    res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx,
        AllLevels, nbrUi, n, cut, r, lbin, sample_weights, flag_sample_weights,
        environment, saveIterations);

    double Ixy_ui = n * res_temp[0];
    double kxy_ui = n * (res_temp[0] - res_temp[1]);
    return InfoBlock(n, Ixy_ui, kxy_ui);
  }
}

// compute Rscore and three point mutual information I(x;y;z | u)
// return Info3PointBlock{score, N * Ixyz_ui, N * kxyz_ui}
Info3PointBlock compute_Rscore_Ixyz_alg5(const TempGrid2d<int>& data,
    const TempGrid2d<int>& sortidx, const TempVector<int>& AllLevels,
    const TempVector<int>& ptr_cnt, const TempVector<int>& ptrVarIdx, int nbrUi,
    int ptrZiIdx, int n, const TempVector<double>& sample_weights,
    bool flag_sample_weights, Environment& environment, bool saveIterations) {
  TempAllocatorScope scope;

  int maxbins = environment.maxbins;
  int initbins = environment.initbins;

  int j, l, ll;

  double I_xy_u, I_xy_zu;
  double Ik_xy_u, Ik_xz_u, Ik_yz_u, Ik_xy_zu;

  vector<double> res_temp;      // res_temp[0]->I,res_temp[1]->I-k

  TempVector<int> ptrVarIdx_t(nbrUi + 2);

  TempVector<int> r(nbrUi + 3);
  TempGrid2d<int> cut(nbrUi + 3, maxbins);

  int lbin = n / initbins;
  if (lbin < 1) {
    lbin = 1;
    initbins = n;
  }

  // initialitize cuts vectors
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (j = 0; j < initbins - 1; j++) {
        cut(l, j) = j * lbin + lbin - 1;
      }
      cut(l, initbins - 1) = n - 1;
      r[l] = initbins;
    } else {
      r[l] = AllLevels[ptrVarIdx[l]];
    }
  }
  // z
  l = nbrUi + 2;
  if (ptr_cnt[ptrZiIdx] == 1) {
    for (j = 0; j < initbins - 1; j++) {
      cut(l, j) = j * lbin + lbin - 1;
    }
    cut(l, initbins - 1) = n - 1;
    r[l] = initbins;
  } else {
    r[l] = AllLevels[ptrZiIdx];
  }

  // if opt
  TempVector<int> r_t(nbrUi + 2);
  TempGrid2d<int> cut_t(nbrUi + 2, maxbins);

  // Optimize variables for each MI estimation for the R score

  // I(x,y|u,z)
  res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx,
      AllLevels, nbrUi + 1, n, cut, r, lbin, sample_weights,
      flag_sample_weights, environment, saveIterations);
  I_xy_zu = res_temp[0];
  Ik_xy_zu = res_temp[1];

  // I(x,y|u)
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx[l]] == 1) {
      for (j = 0; j < initbins - 1; j++) {
        cut_t(l, j) = j * lbin + lbin - 1;
      }
      cut_t(l, initbins - 1) = n - 1;
      for (int j = initbins; j < maxbins; j++) {
        cut_t(l, j) = 0;
      }
      r_t[l] = initbins;
    } else {
      r_t[l] = AllLevels[ptrVarIdx[l]];
    }
  }
  // Do opt run on I(X;Y|U)
  if (nbrUi > 0) {
    res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx,
        AllLevels, nbrUi, n, cut_t, r_t, lbin, sample_weights,
        flag_sample_weights, environment, saveIterations);
  } else {
    res_temp = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, n,
        cut_t, r_t, sample_weights, flag_sample_weights, environment,
        saveIterations);
  }
  I_xy_u = res_temp[0];
  Ik_xy_u = res_temp[1];

  // I(z,x|u)
  ptrVarIdx_t[0] = ptrVarIdx[0];          // X
  ptrVarIdx_t[1] = ptrVarIdx[nbrUi + 2];  // Z
  for (ll = 0; ll < nbrUi; ll++)
    ptrVarIdx_t[ll + 2] = ll + 2;
  // Reset cut
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx_t[l]] == 1) {
      for (j = 0; j < initbins - 1; j++) {
        cut_t(l, j) = j * lbin + lbin - 1;
      }
      cut_t(l, initbins - 1) = n - 1;
      for (int j = initbins; j < maxbins; j++) {
        cut_t(l, j) = 0;
      }
      r_t[l] = initbins;
    } else {
      r_t[l] = AllLevels[ptrVarIdx_t[l]];
    }
  }
  // Do opt run on I(X;Z|U)
  if (nbrUi > 0) {
    res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t,
        AllLevels, nbrUi, n, cut_t, r_t, lbin, sample_weights,
        flag_sample_weights, environment, saveIterations);
  } else {
    res_temp = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t, AllLevels,
        n, cut_t, r_t, sample_weights, flag_sample_weights, environment,
        saveIterations);
  }
  Ik_xz_u = res_temp[1];

  // I(z,y|u)
  ptrVarIdx_t[0] = ptrVarIdx[1];          // Y
  ptrVarIdx_t[1] = ptrVarIdx[nbrUi + 2];  // Z
  for (ll = 0; ll < nbrUi; ll++)
    ptrVarIdx_t[ll + 2] = ll + 2;
  // Reset cut
  for (l = 0; l < (nbrUi + 2); l++) {
    if (ptr_cnt[ptrVarIdx_t[l]] == 1) {
      for (j = 0; j < initbins - 1; j++) {
        cut_t(l, j) = j * lbin + lbin - 1;
      }
      cut_t(l, initbins - 1) = n - 1;
      for (int j = initbins; j < maxbins; j++) {
        cut_t(l, j) = 0;
      }
      r_t[l] = initbins;
    } else {
      r_t[l] = AllLevels[ptrVarIdx_t[l]];
    }
  }
  // Do opt run on I(Y;Z|U)
  if (nbrUi > 0) {
    res_temp = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t,
        AllLevels, nbrUi, n, cut_t, r_t, lbin, sample_weights,
        flag_sample_weights, environment, saveIterations);
  } else {
    res_temp = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t, AllLevels,
        n, cut_t, r_t, sample_weights, flag_sample_weights, environment,
        saveIterations);
  }
  Ik_yz_u = res_temp[1];

  // compute conditional three point mutual information
  double I_xyz_u = I_xy_u - I_xy_zu;
  double Ik_xyz_u = Ik_xy_u - Ik_xy_zu;

  double xz = n * (Ik_xz_u - Ik_xy_u);
  double yz = n * (Ik_yz_u - Ik_xy_u);
  double lower{0}, higher{0};
  std::tie(lower, higher) = std::minmax(xz, yz);

  // Data processing inequality
  double dpi = lower - log1p(exp(lower - higher));
  // probability of not v-structure
  double nv = n * Ik_xyz_u;

  double Rscore = dpi < nv ? dpi : nv;
  return Info3PointBlock{Rscore, n * I_xyz_u, n * I_xyz_u - nv};
}

}  // namespace computation
}  // namespace miic
