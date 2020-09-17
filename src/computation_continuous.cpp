#include "computation_continuous.h"

#include <algorithm>  // std::minmax, std::transform
#include <limits>
#include <numeric>  // std::accumulate
#include <tuple>    // std::tie

#include "linear_allocator.h"
#include "mutual_information.h"
#include "structure.h"

constexpr double kEps = 1.e-5;
constexpr int kMaxIter = 50;
constexpr int kMaxIterOnU = 3;

namespace miic {
namespace computation {

using std::accumulate;
using std::min;
using namespace miic::structure;
using namespace miic::utility;

namespace {

void resetCutPoints(const TempVector<int>& levels,
    const TempVector<int>& is_continuous, const TempVector<int>& var_idx,
    int var_begin, int var_end, int init_nbin, int n, TempGrid2d<int>& cut) {
  for (int l = var_begin; l < var_end; ++l) {
    if (is_continuous[var_idx[l]] != 1) continue;

    int n_bins = min(init_nbin, levels[var_idx[l]]);
    int lbin = n / n_bins;
    if (lbin < 1) {
      lbin = 1;
      n_bins = n;
    }
    for (int j = 0; j < n_bins - 1; ++j) {
      cut(l, j) = j * lbin + lbin - 1;
    }
    cut(l, n_bins - 1) = n - 1;
    for (size_t j = n_bins; j < cut.n_cols(); ++j) {
      cut(l, j) = 0;
    }
  }
}

// update datafactors of a variable from the cut positions vector <cut>
void updateFactors(const TempGrid2d<int>& data_idx, const TempGrid2d<int>& cut,
    const TempVector<int>& is_continuous, const TempVector<int>& var_idx,
    int var_begin, int var_end, TempGrid2d<int>& factors, TempVector<int>& r) {
  int n_samples = data_idx.n_cols();
  for (int l = var_begin; l < var_end; ++l) {
    if (is_continuous[var_idx[l]] != 1) continue;

    int level = 0;
    for (int j = 0; j < n_samples; ++j) {
      if (j > cut(l, level)) ++level;
      int index = data_idx(var_idx[l], j);
      factors(l, index) = level;
    }
    r[l] = level + 1;
  }
  return;
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
void reconstructCutCoarse(const TempVector<int>& memory_cuts,
    const TempVector<int>& memory_cuts2, int np, int n, Ccut& cut) {
  if (memory_cuts[np - 1] == 0) {
    cut[0] = n - 1;
    return;
  }

  int ncuts = 0;
  int l = memory_cuts[np - 1];
  while (l > 0) {
    ncuts++;
    l = memory_cuts[l - 1];
  }
  if (l < 0) ncuts++;

  cut[ncuts] = n - 1;  // conventional last cut
  l = ncuts - 1;
  int s = memory_cuts[np - 1];
  cut[l] = memory_cuts2[np - 1];  // truly last cut
  l--;
  while (s > 0 && l >= 0) {
    cut[l] = memory_cuts2[s - 1];
    s = memory_cuts[s - 1];
    l--;
  }
}

// GENERAL VARIABLES NOTATION

// int coarse : coarse graining step : minimum length of a bin : considering
// only n_samples/coarse possible cuts

/*  DYNAMIC PROGRAMMING OPTIMIZATION FUNCTION
 *
 *  Takes as input the factors and joint factors and optimizes a variable by
 * maximizing a given information. n_vars determines the number of factors to
 * consider.
 *
 * ////////////////////////////////////////////
 * INPUT
 *
 * The optimizing variable is defined by <sortidx_var> and <data>.
 * <sortidx> Contains the rank indices (ordering) of the optimized variable,
 * <data> contains the ranks.
 *
 * <n_vars> is the number of terms in the optimization functions, consisting in
 * variables and joint variables, which (discretized) values are stored in
 * <factors> : ( <n_vars> x n_samples )  matrix : Example <n_vars>=2 :
 * Optimizing X on I(X;Y) xy factors : factors[0] x factors : factors[1]
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
 * <n_samples> is the number of samples in total.
 *
 * <nnr> is the number of non repeated samples.
 *
 * <cut> is the vector containing cutpoints for discretizing the optimized
 * variables. Modified at the end wih the optimal solution with a call to
 * reconstructCutCoarse().
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
 *      <n_vars> -> 2
 *      <factors>
 *           [0] : The discretized Y
 *           [1] : A vector of one repeated level (starting from X with one
 * single bin) <sc_levels1>        -> the number of levels of discretized Y
 *      <previous_levels>   -> the number of levels of discretizedd X in the
 * previous iteration (for cost correction) <cut>               -> the vector
 * that will contain the cutpoints for X
 */
template <typename Cf0, typename Cf1, typename Ccut,
    typename =
        void_t<IsIntContainer<Cf0>, IsIntContainer<Cf1>, IsIntContainer<Ccut>>>
void optimizeCutPoints(const TempGrid2d<int>::ConstRow& data,
    const TempGrid2d<int>::ConstRow& data_idx, const Cf0& factors0,
    const Cf1& factors1, int r0, int r1, int sc_levels1, int previous_levels,
    int nnr, const TempVector<double>& sample_weights, bool flag_sample_weights,
    int maxbins, std::shared_ptr<CtermCache> cache, int cplx, Ccut&& cut) {
  TempAllocatorScope scope;

  int coarse = ceil(1.0 * nnr / maxbins);  // step coarse graining
  if (coarse < 1) coarse = 1;
  int np = ceil(1.0 * nnr / coarse);  // number of possible cuts <= max_level
  int n_vars = 2;
  int n_samples = data.size();

  // temp variables to memorize optimization cuts
  TempVector<int> memory_cuts_idx(np);  // indexes of the cuts (1..np)
  TempVector<int> memory_cuts_pos(np);  // positions of the cuts (1..n_samples)

  // function max value at each step
  TempVector<double> I(np);  // The optimal information value found at each idx.
  TempVector<double> Ik(np);
  double I_kj;  // The information value for a bin fom idx k to j.
  double Ik_kj;

  int njforward(0), nkforward(0);  // Indexes at current position of j and k

  // entropy in kj interval for the <n_vars>+1 terms
  TempVector<double> H_kj(n_vars);
  TempVector<double> Hk_kj(n_vars);

  TempVector<int> counts_0(r0);
  TempVector<int> counts_1(r1);
  TempVector<int> counts_k_0(r0);
  TempVector<int> counts_k_1(r1);

  TempGrid2d<int> coarse_counts_marginal(np, r1, 0);
  TempGrid2d<int> coarse_counts_joint(np, r0, 0);

  int weighted_count;

  TempVector<int> check_repet(n_samples);
  for (int i = 0; i < (n_samples - 1); i++) {
    check_repet[i] = (data[data_idx[i + 1]] != data[data_idx[i]]);
  }

  TempVector<int> n_values(np);
  TempVector<double> sum_sample_weights(np);
  int ir(0);  // Iterator on non repeated values
  int level_marginal(0), level_joint(0);
  for (int i = 0; i < np; i++) {
    ir = 0;  // iterator on not repeated values
    if (flag_sample_weights && i > 0)
      sum_sample_weights[i] = sum_sample_weights[i - 1];

    while ((ir < coarse) && (njforward < n_samples)) {
      level_marginal = factors1[data_idx[njforward]];
      level_joint = factors0[data_idx[njforward]];
      coarse_counts_marginal(i, level_marginal)++;
      coarse_counts_joint(i, level_joint)++;

      if (flag_sample_weights)
        sum_sample_weights[i] += sample_weights[njforward];
      if (njforward + 1 < n_samples) {  // check no repetition
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

  double sc = 0.5 * (sc_levels1 - 1) * r1;
  // moving j over the np possible cuts
  for (int j = 0; j < np; j++) {  // j=1...n_samples-1

    njforward = n_values[j];
    ef_nj = sum_sample_weights[j];
    if (flag_sample_weights) efN_factor = ef_nj / njforward;
    std::transform(begin(counts_0), end(counts_0),
        coarse_counts_joint.row_begin(j), begin(counts_0), std::plus<int>());
    std::transform(begin(counts_1), end(counts_1),
        coarse_counts_marginal.row_begin(j), begin(counts_1), std::plus<int>());

    Hk_kj[0] = 0;  // joint
    Hk_kj[1] = 0;  // marginal
    H_kj[0] = 0;   // joint
    H_kj[1] = 0;   // marginal
    for (int level = 0; level < r0; level++) {
      weighted_count = flag_sample_weights
                           ? int(efN_factor * counts_0[level] + 0.5)
                           : counts_0[level];
      Hk_kj[0] += cache->getH(weighted_count);
    }
    H_kj[0] = Hk_kj[0];

    for (int level = 0; level < r1; level++) {
      weighted_count = flag_sample_weights
                           ? int(efN_factor * counts_1[level] + 0.5)
                           : counts_1[level];
      Hk_kj[1] -= cache->getH(weighted_count);
      H_kj[1] -= cache->getH(weighted_count);

      if (cplx == 0 && counts_1[level] > 0)
        Hk_kj[1] -= sc * cache->getLog(n_samples);
      else if (cplx == 1) {
        Hk_kj[1] -= cache->getLogC(weighted_count, sc_levels1);
      }
    }

    I[j] = 0;
    Ik[j] = 0;
    for (int m = 0; m < n_vars; m++) {
      I[j] += H_kj[m];    // herve
      Ik[j] += Hk_kj[m];  // herve
    }

    double Imax = std::numeric_limits<double>::lowest();
    double Ikmax = std::numeric_limits<double>::lowest();

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

    for (int k = 0; k < j; k++) {  // k=1...n_samples-2 possible cuts

      nkforward = n_values[k];
      ef_nk = sum_sample_weights[j] - sum_sample_weights[k];
      if (flag_sample_weights) efN_factor = ef_nk / (njforward - nkforward);

      std::transform(begin(counts_k_1), end(counts_k_1),
          coarse_counts_marginal.row_begin(k), begin(counts_k_1),
          std::minus<int>());
      std::transform(begin(counts_k_0), end(counts_k_0),
          coarse_counts_joint.row_begin(k), begin(counts_k_0),
          std::minus<int>());

      H_kj[0] = 0;
      Hk_kj[0] = 0;
      for (int level = 0; level < r0; level++) {
        weighted_count = flag_sample_weights
                             ? int(efN_factor * counts_k_0[level] + 0.5)
                             : counts_k_0[level];
        Hk_kj[0] += cache->getH(weighted_count);
      }
      H_kj[0] = Hk_kj[0];

      H_kj[1] = 0;
      Hk_kj[1] = 0;
      for (int level = 0; level < r1; level++) {
        weighted_count = flag_sample_weights
                             ? int(efN_factor * counts_k_1[level] + 0.5)
                             : counts_k_1[level];
        Hk_kj[1] -= cache->getH(weighted_count);
        H_kj[1] -= cache->getH(weighted_count);

        if (cplx == 0 && counts_k_1[level] > 0)
          Hk_kj[1] -= sc * cache->getLog(n_samples);
        else if (cplx == 1) {
          Hk_kj[1] -= cache->getLogC(weighted_count, sc_levels1);
        }
      }

      I_kj = 0;
      Ik_kj = 0;
      for (int m = 0; m < n_vars; m++) {
        I_kj += H_kj[m];
        Ik_kj += Hk_kj[m];
      }
      if (cplx == 1) {
        // Combinatorial approximation
        Ik_kj -= cache->getLogChoose(np - 1, previous_levels - 1) /
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

  reconstructCutCoarse(memory_cuts_idx, memory_cuts_pos, np, n_samples, cut);
}

// compute I(x,y)
// dynamic programming for optimizing variables binning
// optimize on x I(x,y): Hx - Hxy - kmdl
// optimize on y I(x,y): Hy - Hxy - kmdl
// until convergence
InfoBlock computeIxy(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& is_continuous,
    const TempVector<int>& var_idx, const TempVector<int>& levels,
    const TempVector<double>& weights, bool flag_sample_weights, int initbins,
    int maxbins, int cplx, std::shared_ptr<CtermCache> cache,
    std::shared_ptr<CutPointsInfo> cuts_info = nullptr) {
  TempAllocatorScope scope;

  bool save_cuts = cuts_info && !cuts_info->cutpoints.empty();
  int n_samples = data.n_cols();
  double n_eff = accumulate(begin(weights), end(weights), 0.0);

  TempGrid2d<int> datafactors(2, n_samples);
  TempGrid2d<int> cut(2, maxbins);
  TempVector<int> r(2);
  // Initialize discrete factors and n_levels
  for (int l = 0; l < 2; ++l) {
    if (is_continuous[var_idx[l]] != 1) {
      r[l] = levels[var_idx[l]];
      for (int j = 0; j < n_samples; ++j)
        datafactors(l, j) = data(var_idx[l], j);
    }
  }
  // Find the best initial conditions with the same number of bins (equalfreq)
  // on all continuous variables.
  double best_res{std::numeric_limits<double>::lowest()};
  int best_initbins{initbins};
  int n_levels_min{n_samples};
  for (int l = 0; l < 2; ++l) {
    if (is_continuous[var_idx[l]] == 1)
      n_levels_min = min(n_levels_min, levels[var_idx[l]]);
  }
  int n_test_max = min(min(initbins, 20), n_levels_min);
  for (int test_n_bins = 2; test_n_bins < n_test_max; ++test_n_bins) {
    TempAllocatorScope scope;
    // Reinitialization cut, datafactors and r
    resetCutPoints(
        levels, is_continuous, var_idx, 0, 2, test_n_bins, n_samples, cut);
    updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);

    int rxy{0};
    TempVector<int> xy_factors(n_samples);
    jointfactors_u(datafactors, r, xy_factors, rxy);

    TempVector<int> r_temp{r[0], r[1], rxy};
    InfoBlock res_temp;
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), datafactors.getRow(1),
          xy_factors, r_temp, n_eff, weights, cache, 0);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), datafactors.getRow(1),
          xy_factors, r_temp, cache, 0);
    // All set if both variables are discrete
    if (is_continuous[var_idx[0]] == 0 && is_continuous[var_idx[1]] == 0)
      return res_temp;

    double Ik = res_temp.Ixy_ui - res_temp.kxy_ui;
    if (Ik > best_res) {
      best_initbins = test_n_bins;
      best_res = Ik;
    }
  }
  // Reinitialization cut and r
  resetCutPoints(
      levels, is_continuous, var_idx, 0, 2, best_initbins, n_samples, cut);
  // initialization of datafactors && sortidx
  updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);
  // Run dynamic optimization with the best initial conditions.

  int iter_max = save_cuts ? cuts_info->cutpoints.n_rows() : kMaxIter;
  // Keep the result of each iteration
  TempVector<double> I_list(iter_max);
  TempVector<double> Ik_list(iter_max);
  double Ixy{0}, Ikxy{0};  // to be returned
  for (int step = 0; step < iter_max; ++step) {
    if (is_continuous[var_idx[0]]) {
      TempAllocatorScope scope;
      // max_x{ I(x;y) }
      int r0 = r[1];  // y
      int r1 = 1;   // One single bin at start
      int sc_levels1 = r0;
      int sc_levels2 = r[0];
      // Optimization run on X. 2 factors
      optimizeCutPoints(data.getConstRow(var_idx[0]),
          data_idx.getConstRow(var_idx[0]), datafactors.getRow(1),
          TempVector<int>(n_samples, 0), r0, r1, sc_levels1, sc_levels2,
          levels[var_idx[0]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(0));
    }
    if (is_continuous[var_idx[1]]) {
      TempAllocatorScope scope;
      // max_y{ I(x;y) }
      int r0 = r[0];  // x before its optimization
      int r1 = 1;   // One single bin at start
      int sc_levels1 = r0;
      int sc_levels2 = r[1];
      // Optimization run on Y. 2 factors
      optimizeCutPoints(data.getConstRow(var_idx[1]),
          data_idx.getConstRow(var_idx[1]), datafactors.getRow(0),
          TempVector<int>(n_samples, 0), r0, r1, sc_levels1, sc_levels2,
          levels[var_idx[1]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(1));
    }
    updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);

    if (save_cuts) {
      for (int j = 0; j < maxbins; ++j) {
        cuts_info->cutpoints(step, j) = cut(0, j);
        cuts_info->cutpoints(step, j + maxbins) = cut(1, j);
      }
    }
    TempAllocatorScope scope;

    int rxy = 0;
    TempVector<int> xy_factors(n_samples);
    jointfactors_u(datafactors, r, xy_factors, rxy);

    TempVector<int> r_temp{r[0], r[1], rxy};
    InfoBlock res_temp;
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), datafactors.getRow(1),
          xy_factors, r_temp, n_eff, weights, cache, 0);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), datafactors.getRow(1),
          xy_factors, r_temp, cache, 0);
    // Adding combinatorial term
    if (is_continuous[var_idx[0]] && r[0] > 1) {
      int np = min(maxbins, levels[var_idx[0]]);
      if (r[0] < np)
        res_temp.kxy_ui += cache->getLogChoose(np - 1, r[0] - 1);
    }
    if (is_continuous[var_idx[1]] && r[1] > 1) {
      int np = min(maxbins, levels[var_idx[1]]);
      if (r[1] < np)
        res_temp.kxy_ui += cache->getLogChoose(np - 1, r[1] - 1);
    }

    I_list[step] = res_temp.Ixy_ui;
    Ik_list[step] = res_temp.Ixy_ui - res_temp.kxy_ui;
    bool converged{false};
    for (int i = step - 1; i >= 0; --i) {
      if (fabs(Ik_list[step] - Ik_list[i]) < kEps) {
        converged = true;
        Ixy = accumulate(begin(I_list) + i, begin(I_list) + step, 0.0);
        Ikxy = accumulate(begin(Ik_list) + i, begin(Ik_list) + step, 0.0);
        Ikxy /= (step - i);  // average over the periodic cycle
        Ixy /= (step - i);
        break;
      }
    }

    if (converged) break;

    Ixy = res_temp.Ixy_ui;
    Ikxy = res_temp.Ixy_ui - res_temp.kxy_ui;
    // Already optimal if any one of them is discrete
    if (!is_continuous[var_idx[0]] || !is_continuous[var_idx[1]]) break;
  }  // for step

  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if (Ikxy < 0 && is_continuous[var_idx[0]] && is_continuous[var_idx[1]]) {
    Ixy = 0;
    Ikxy = 0;
  }
  if (save_cuts) {
    cuts_info->Ik = Ikxy / n_samples;
    cuts_info->I = Ixy / n_samples;
    cuts_info->I_equal_freq_max = best_res;
  }

  return InfoBlock(n_samples, Ixy, Ixy - Ikxy);
}

InfoBlock computeIxyui(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& is_continuous,
    const TempVector<int>& var_idx, const TempVector<int>& levels,
    const TempVector<double>& weights, bool flag_sample_weights, int initbins,
    int maxbins, int cplx, std::shared_ptr<CtermCache> cache,
    std::shared_ptr<CutPointsInfo> cuts_info = nullptr) {
  TempAllocatorScope scope;

  bool save_cuts = cuts_info && !cuts_info->cutpoints.empty();
  int n_samples = data.n_cols();
  int n_nodes = var_idx.size();
  int n_ui = n_nodes - 2;
  double n_eff = accumulate(begin(weights), end(weights), 0.0);
  // allocation factors  x y
  TempGrid2d<int> datafactors(n_nodes, n_samples);
  TempGrid2d<int> cut(n_nodes, maxbins);
  TempVector<int> r(n_nodes);  // n_levels of optimized variables

  // Initialize discrete factors and n_levels
  for (int l = 0; l < n_nodes; ++l) {
    if (is_continuous[var_idx[l]] != 1) {
      r[l] = levels[var_idx[l]];
      for (int j = 0; j < n_samples; ++j)
        datafactors(l, j) = data(var_idx[l], j);
    }
  }
  TempGrid2d<int> uiyxfactors(4, n_samples);  // ui, uiy, uix, uiyx
  TempVector<int> ruiyx(4);
  // Find the best initial conditions with the same number of bins (equalfreq)
  // on all continuous variables.
  double best_res{std::numeric_limits<double>::lowest()};
  int best_initbins{initbins};
  int n_levels_min{n_samples};
  for (int l = 0; l < n_nodes; ++l) {
    if (is_continuous[var_idx[l]] == 1)
      n_levels_min = min(n_levels_min, levels[var_idx[l]]);
  }
  int n_test_max = min(min(initbins, 20), n_levels_min);
  TempVector<int> r_temp(3);
  InfoBlock res_temp;
  for (int test_n_bins = 2; test_n_bins < n_test_max; ++test_n_bins) {
    // Initialize factors, cut and r
    resetCutPoints(levels, is_continuous, var_idx, 0, n_nodes, test_n_bins,
        n_samples, cut);
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 0, n_nodes, datafactors, r);

    jointfactors_uiyx(datafactors, r, -1, uiyxfactors, ruiyx);
    r_temp.assign({r[1], ruiyx[2], ruiyx[3]});
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(1), uiyxfactors.getRow(2),
          uiyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(1), uiyxfactors.getRow(2),
          uiyxfactors.getRow(3), r_temp, cache, 0);
    double Ik_y_xu = res_temp.Ixy_ui - res_temp.kxy_ui;

    r_temp.assign({r[0], ruiyx[1], ruiyx[3]});
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), uiyxfactors.getRow(1),
          uiyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), uiyxfactors.getRow(1),
          uiyxfactors.getRow(3), r_temp, cache, 0);
    double Ik_x_yu = res_temp.Ixy_ui - res_temp.kxy_ui;

    if ((Ik_y_xu + Ik_x_yu) > best_res) {
      best_initbins = test_n_bins;
      best_res = (Ik_y_xu + Ik_x_yu);
    }
  }
  // Initialize X and Y cuts with best_initbins
  resetCutPoints(
      levels, is_continuous, var_idx, 0, 2, best_initbins, n_samples, cut);
  updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);

  int iter_max = save_cuts ? cuts_info->cutpoints.n_rows() : kMaxIter;
  // Keep the result of each iteration
  TempVector<double> I_list(iter_max);
  TempVector<double> Ik_list(iter_max);
  TempVector<int> r_old(r);  // n_levels in the previous iteration
  double Ixy_ui{0}, Ikxy_ui{0};  // to be returned
  // Run optimization with best initial equal freq.
  for (int step = 0; step < iter_max; ++step) {
    // optimize I(y;xu) over x and u
    resetCutPoints(
        levels, is_continuous, var_idx, 2, n_nodes, initbins, n_samples, cut);
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
    std::copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; count < kMaxIterOnU; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] != 1) continue;
        // opt u, I(y;xu)
        jointfactors_uiyx(datafactors, r_old, l, uiyxfactors, ruiyx);
        // init variables for the optimization run
        int r0 = ruiyx[3];  // xyu
        int r1 = ruiyx[2];  // xu
        int sc_levels1 = r_old[1];
        int sc_levels2 = r_old[l];  // old nlevels for combinatorial term
        // Run optimization on U. 2 factors xyu and xu
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uiyxfactors.getRow(3),
            uiyxfactors.getRow(2), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors(
          data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
      std::copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
    r_temp.assign({r_old[1], ruiyx[2], ruiyx[3]});
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(1), uiyxfactors.getRow(2),
          uiyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(1), uiyxfactors.getRow(2),
          uiyxfactors.getRow(3), r_temp, cache, 0);
    double I_y_xu = res_temp.Ixy_ui;  // Before optimization on X.
    double Ik_y_xu = res_temp.Ixy_ui - res_temp.kxy_ui;
    if (is_continuous[var_idx[0]] && r_old[0] > 1) {
      int np = min(levels[var_idx[0]], maxbins);
      if (r_old[0] < np) Ik_y_xu -= cache->getLogChoose(np - 1, r_old[0] - 1);
    }

    if (is_continuous[var_idx[0]] == 1) {
      // I(y;xu), optimize on x
      int r0 = ruiyx[1];  // uy
      int r1 = ruiyx[0];  // u
      int sc_levels1 = r_old[1];
      int sc_levels2 = r_old[0];
      // Run optimization on X. 2 factors uy and u
      optimizeCutPoints(data.getConstRow(var_idx[0]),
          data_idx.getConstRow(var_idx[0]), uiyxfactors.getRow(1),
          uiyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
          levels[var_idx[0]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(0));
    }
    // optimize I(x;yu) over y and u
    // Reset cutpoints on U
    resetCutPoints(
        levels, is_continuous, var_idx, 2, n_nodes, initbins, n_samples, cut);
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
    std::copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; count < kMaxIterOnU; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] != 1) continue;

        jointfactors_uiyx(datafactors, r_old, l, uiyxfactors, ruiyx);
        // init variables for the optimization run
        int r0 = ruiyx[3];  // xyu
        int r1 = ruiyx[1];  // yu
        int sc_levels1 = r_old[0];
        int sc_levels2 = r_old[l];
        // Run optimization on U. 2 factors xyu and yu
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uiyxfactors.getRow(3),
            uiyxfactors.getRow(1), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors(
          data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
      std::copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
    r_temp.assign({r_old[0], ruiyx[1], ruiyx[3]});
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), uiyxfactors.getRow(1),
          uiyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), uiyxfactors.getRow(1),
          uiyxfactors.getRow(3), r_temp, cache, 0);
    double I_x_yu = res_temp.Ixy_ui;  // Before updating Y (and X).
    double Ik_x_yu = res_temp.Ixy_ui - res_temp.kxy_ui;
    if ((is_continuous[var_idx[1]] == 1) && (r_old[1] > 1)) {
      int np = min(levels[var_idx[1]], maxbins);
      if (r_old[1] < np) {
        Ik_x_yu -= cache->getLogChoose(np - 1, r_old[1] - 1);
      }
    }

    if (is_continuous[var_idx[1]] == 1) {
      // I(x;yu), optimize on y
      int r0 = ruiyx[2];  // ux
      int r1 = ruiyx[0];  // u
      int sc_levels1 = r_old[0];
      int sc_levels2 = r_old[1];
      // Run optimization on Y. 2 factors ux and u
      optimizeCutPoints(data.getConstRow(var_idx[1]),
          data_idx.getConstRow(var_idx[1]), uiyxfactors.getRow(2),
          uiyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
          levels[var_idx[1]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(1));
    }

    // Reset cutpoints on U
    resetCutPoints(
        levels, is_continuous, var_idx, 2, n_nodes, initbins, n_samples, cut);
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
    std::copy(begin(r) + 2, end(r), begin(r_old) + 2);
    // optimize I(x;u) over u
    for (int count = 0; count < kMaxIterOnU; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] != 1) continue;

        jointfactors_uiyx(datafactors, r_old, l, uiyxfactors, ruiyx);
        // init variables for the optimization run
        int r0 = ruiyx[2];  // xu
        int r1 = ruiyx[0];  // u
        int sc_levels1 = r_old[0];  // x
        int sc_levels2 = r_old[l];  // u
        // optimization run on var_idx[l], 2 factors xu and u
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uiyxfactors.getRow(2),
            uiyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors(
          data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
      std::copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
    r_temp.assign({r_old[0], ruiyx[0], ruiyx[2]});
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(0), uiyxfactors.getRow(0),
          uiyxfactors.getRow(2), r_temp, n_eff, weights, cache, cplx);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(0), uiyxfactors.getRow(0),
          uiyxfactors.getRow(2), r_temp, cache, 0);
    double I_x_u = res_temp.Ixy_ui;  // After optimization on U.
    double Ik_x_u = res_temp.Ixy_ui - res_temp.kxy_ui;

    // optimize I(y;u) over u
    resetCutPoints(
        levels, is_continuous, var_idx, 2, n_nodes, initbins, n_samples, cut);
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
    std::copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; count < kMaxIterOnU; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] != 1) continue;

        jointfactors_uiyx(datafactors, r_old, l, uiyxfactors, ruiyx);
        int r0 = ruiyx[1];  // yu
        int r1 = ruiyx[0];  // u
        int sc_levels1 = r_old[1];  // y
        int sc_levels2 = r_old[l];  // u
        // optimization run on var_idx[l], 2 factors yu and u
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uiyxfactors.getRow(1),
            uiyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors(
          data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
      std::copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    jointfactors_uiyx(datafactors, r_old, -1, uiyxfactors, ruiyx);
    r_temp.assign({r_old[1], ruiyx[0], ruiyx[1]});
    if (cplx == 1)
      res_temp = computeMI_knml(datafactors.getRow(1), uiyxfactors.getRow(0),
          uiyxfactors.getRow(1), r_temp, n_eff, weights, cache, cplx);
    else
      res_temp = computeMI_kmdl(datafactors.getRow(1), uiyxfactors.getRow(0),
          uiyxfactors.getRow(1), r_temp, cache, 0);
    double I_y_u = res_temp.Ixy_ui;  // After optimization on U.
    double Ik_y_u = res_temp.Ixy_ui - res_temp.kxy_ui;

    // End of iteration: update X and Y cuts
    updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);
    std::copy(begin(r), end(r), begin(r_old));

    if (save_cuts) {
      for (int j = 0; j < maxbins; ++j) {
        cuts_info->cutpoints(step, j) = cut(0, j);
        cuts_info->cutpoints(step, j + maxbins) = cut(1, j);
      }
    }

    // Compute I(X;Y|U)
    double cond_I = 0.5 * (I_x_yu - I_x_u + I_y_xu - I_y_u);
    double cond_Ik = 0.5 * (Ik_x_yu - Ik_x_u + Ik_y_xu - Ik_y_u);
    I_list[step] = cond_I;
    Ik_list[step] = cond_Ik;
    // Test stop condition on stop1
    bool converged{false};
    for (int i = step - 1; i >= 0; i--) {
      // If no real improvement over last information
      if (fabs(cond_Ik - Ik_list[i]) < kEps) {
        converged = true;
        Ixy_ui = accumulate(begin(I_list) + i, begin(I_list) + step, 0.0);
        Ikxy_ui = accumulate(begin(Ik_list) + i, begin(Ik_list) + step, 0.0);
        Ikxy_ui /= (step - i);  // average over the periodic cycle
        Ixy_ui /= (step - i);
        break;
      }
    }
    if (converged) break;

    Ixy_ui = cond_I;
    Ikxy_ui = cond_Ik;
  }  // for step
  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if (Ikxy_ui < 0 && is_continuous[var_idx[0]] && is_continuous[var_idx[1]]) {
    Ixy_ui = 0;
    Ikxy_ui = 0;
  }
  if (save_cuts) {
    cuts_info->Ik = Ikxy_ui / n_samples;
    cuts_info->I = Ixy_ui / n_samples;
  }

  return InfoBlock(n_samples, Ixy_ui, Ixy_ui - Ikxy_ui);
}

}  // anonymous namespace

InfoBlock computeCondMutualInfo(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& levels,
    const TempVector<int>& is_continuous, const TempVector<int>& var_idx,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    int initbins, int maxbins, int cplx, std::shared_ptr<CtermCache> cache,
    std::shared_ptr<CutPointsInfo> cuts_info) {

  if (data.n_rows() == 2) {
    return computeIxy(data, data_idx, is_continuous, var_idx, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx, cache,
        cuts_info);
  } else {
    return computeIxyui(data, data_idx, is_continuous, var_idx, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx, cache,
        cuts_info);
  }
}

// compute Rscore and three point mutual information I(x;y;z | u)
// return Info3PointBlock{score, N * Ixyz_ui, N * kxyz_ui}
Info3PointBlock computeInfo3PointAndScore(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& levels,
    const TempVector<int>& is_continuous, const TempVector<int>& var_idx,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    int initbins, int maxbins, int cplx, std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_ui = data.n_rows() - 3;
  // Optimize variables for each MI estimation for the R score
  // I(x,y|u,z)
  InfoBlock res_temp = computeIxyui(data, data_idx, is_continuous, var_idx,
      levels, sample_weights, flag_sample_weights, initbins, maxbins, cplx,
      cache);
  double I_xy_zu = res_temp.Ixy_ui;
  double Ik_xy_zu = res_temp.Ixy_ui - res_temp.kxy_ui;

  TempVector<int> var_idx_t(n_ui + 2);
  for (int l = 0; l < n_ui + 2; ++l)
    var_idx_t[l] = var_idx[l];
  // Do opt run on I(X;Y|U)
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx, cache);
  }
  double I_xy_u = res_temp.Ixy_ui;
  double Ik_xy_u = res_temp.Ixy_ui - res_temp.kxy_ui;

  // I(z,x|u)
  var_idx_t[0] = var_idx[0];         // X
  var_idx_t[1] = var_idx[n_ui + 2];  // Z
  for (int l = 2; l < n_ui + 2; ++l)
    var_idx_t[l] = l;
  // Do opt run on I(X;Z|U)
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx, cache);
  }
  double Ik_xz_u = res_temp.Ixy_ui - res_temp.kxy_ui;

  // I(z,y|u)
  var_idx_t[0] = var_idx[1];         // Y
  var_idx_t[1] = var_idx[n_ui + 2];  // Z
  for (int l = 2; l < n_ui + 2; ++l)
    var_idx_t[l] = l;
  // Do opt run on I(Y;Z|U)
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx, cache);
  }
  double Ik_yz_u = res_temp.Ixy_ui - res_temp.kxy_ui;

  // compute conditional three point mutual information
  double I_xyz_u = I_xy_u - I_xy_zu;
  double Ik_xyz_u = Ik_xy_u - Ik_xy_zu;

  double xz = Ik_xz_u - Ik_xy_u;
  double yz = Ik_yz_u - Ik_xy_u;
  double lower{0}, higher{0};
  std::tie(lower, higher) = std::minmax(xz, yz);

  // Data processing inequality
  double dpi = lower - log1p(exp(lower - higher));
  // Probability of not v-structure
  double nv = Ik_xyz_u;

  double Rscore = dpi < nv ? dpi : nv;
  return Info3PointBlock{Rscore, I_xyz_u, I_xyz_u - nv};
}

}  // namespace computation
}  // namespace miic
