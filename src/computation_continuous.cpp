#include "computation_continuous.h"

#include <algorithm>  // std::min, std::transform
#include <limits>
#include <numeric>  // std::accumulate

#include "linear_allocator.h"
#include "mutual_information.h"
#include "structure.h"

constexpr double kEpsilon = 1.e-5;
constexpr double kPrecision = 1.e-10;
constexpr int kMaxIter = 50;
constexpr int kMaxIterOnU = 3;

namespace miic {
namespace computation {

using std::accumulate;
using std::copy;
using std::min;
using std::max;
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
void reconstructCutCoarse(const TempVector<int>& memory_cuts_idx,
    const TempVector<int>& memory_cuts_pos, int n_samples, Ccut& cut) {
  if (memory_cuts_idx.back() == 0) {
    cut[0] = n_samples - 1;
    return;
  }

  int ncuts{0};
  int l{memory_cuts_idx.back()};
  while (l > 0) {
    ++ncuts;
    l = memory_cuts_idx[l - 1];
  }
  if (l < 0) ncuts++;

  cut[ncuts] = n_samples - 1;  // conventional last cut
  l = ncuts - 1;
  cut[l] = memory_cuts_pos.back();  // truly last cut
  --l;
  int s = memory_cuts_idx.back();
  while (s > 0 && l >= 0) {
    cut[l] = memory_cuts_pos[s - 1];
    s = memory_cuts_idx[s - 1];
    --l;
  }
}

// Given factors of target variable and joint variable, find a set
// of cut points for target variable <f_target> by maximizing mutual information
// max_cuts{ I(target|aux) }
//
// The optimizing variable is defined by <data_idx> and <data>.
// data_idx[i] is the [i]th smallest value of the optimized variable,
// <data> contains the ranks.
//
// <r> is the number of unique values of each factor.
//
// <sc> is the stochastic complexity used for the simple MDL cost.
// Its value is 0.5 * number_of_other_levels where number of other levels is the
// product of number of unique observed levels of all variables except the
// variable being optimized.
//
// <sc_levels> and <r1_prev> are respectively the number of levels of
// the other variable in the mutual information term and the number of levels of
// the otpimized variable in the previous iteration.
//
// <n_nr> is the number of non repeated values.
//
// <cut> is the vector containing cutpoints for discretizing the optimized
// variables. Modified at the end wih the optimal solution.
//
// <maxbins> is the maximum number of bins allowed. It controls the resolution
// of possible cutpoints : For example, 50 maximum cutpoints on a 1000 samples
// vector results in possible cutpoints every 1000/50 = 20 positions. This
// speeds up significantly the dynamic programming at the cost of finding
// approximated solutions.
//
// <cplx> is the choice of complexity : 0 for simple MDL (product of all
// observed values, see <sc>) and 1 for NML with stochastic complexity and a
// combinatorial term with the previous number of levels.
//
// <weights> are unique weights for each sample to account for effective
// N != N, if the <flag_efN> is set to true.
//
// Example: Discretize continuous variable X by maximizing Ik(X;Y).
//   <data>: Ranks of X without discretization
//   <data_idx>: Order of <data>, data[data_idx[0]] is the minimal value of X
//   <factors0>: Discret(e)ized variable Y
//   <factors1>: vector<int>(n_samples, 0) (X with one single bin)
//   <sc_levels1>: Number of levels of discret(e)ized variable Y
//   <r1_prev>: Number of levels of discretized X in the previous iteration,
//              used for complexity term
//   <cut>: Output, cut points for X
template <typename Cf0, typename Cf1, typename Ccut,
    typename =
        void_t<IsIntContainer<Cf0>, IsIntContainer<Cf1>, IsIntContainer<Ccut>>>
void optimizeCutPoints(const TempGrid2d<int>::ConstRow& data,
    const TempGrid2d<int>::ConstRow& data_idx, const Cf0& factors0,
    const Cf1& factors1, int r0, int r1, int sc_levels1, int r1_prev, int n_nr,
    const TempVector<double>& weights, bool flag_weights, int maxbins,
    std::shared_ptr<CtermCache> cache, int cplx, Ccut&& cut) {
  TempAllocatorScope scope;

  // coarse graining step, minimum length of a bin
  int coarse = std::max(ceil(1.0 * n_nr / maxbins), 1.0);
  // maximal number of possible cut points
  int n_cuts_max = ceil(1.0 * n_nr / coarse);
  int n_samples = data.size();

  // entropy in kj interval
  TempGrid2d<int> coarse_counts_target(n_cuts_max, r1, 0);
  TempGrid2d<int> coarse_counts_joint(n_cuts_max, r0, 0);

  TempVector<int> n_values(n_cuts_max);
  TempVector<double> sum_weights(n_cuts_max);
  for (int i = 0, it_j = 0; i < n_cuts_max; ++i) {
    if (flag_weights && i > 0) sum_weights[i] = sum_weights[i - 1];

    int ir{0};  // iterator on not repeated values
    while (ir < coarse && it_j < n_samples) {
      int level_target = factors1[data_idx[it_j]];
      int level_joint = factors0[data_idx[it_j]];
      ++coarse_counts_target(i, level_target);
      ++coarse_counts_joint(i, level_joint);

      if (flag_weights) sum_weights[i] += weights[it_j];
      if (it_j + 1 < n_samples) {  // check no repetition
        ir += data[data_idx[it_j + 1]] != data[data_idx[it_j]];
      }
      ++it_j;
    }
    n_values[i] = it_j;
  }

  // Keep the result of each iteration
  TempVector<double> I(n_cuts_max);
  TempVector<double> Ik(n_cuts_max);
  // indexes of the cuts (1..n_cuts_max)
  TempVector<int> memory_cuts_idx(n_cuts_max);
  // positions of the cuts (1..n_samples)
  TempVector<int> memory_cuts_pos(n_cuts_max);

  TempVector<int> counts_0(r0);
  TempVector<int> counts_1(r1);

  double k_mdl = 0.5 * (sc_levels1 - 1) * r1;
  for (int j = 0; j < n_cuts_max; ++j) {
    int it_j = n_values[j];
    double ef_nj = sum_weights[j];
    double eff_n_factors = ef_nj / it_j;
    std::transform(begin(counts_0), end(counts_0),
        coarse_counts_joint.row_begin(j), begin(counts_0), std::plus<int>());
    std::transform(begin(counts_1), end(counts_1),
        coarse_counts_target.row_begin(j), begin(counts_1), std::plus<int>());

    double Hk_j_joint{0}, H_j_joint{0};
    for (int level = 0; level < r0; level++) {
      int weighted_count = flag_weights
                               ? int(eff_n_factors * counts_0[level] + 0.5)
                               : counts_0[level];
      Hk_j_joint += cache->getH(weighted_count);
    }
    H_j_joint = Hk_j_joint;

    double Hk_j_target{0}, H_j_target{0};
    for (int level = 0; level < r1; level++) {
      int weighted_count = flag_weights
                               ? int(eff_n_factors * counts_1[level] + 0.5)
                               : counts_1[level];
      Hk_j_target -= cache->getH(weighted_count);
      H_j_target -= cache->getH(weighted_count);

      if (cplx == 0 && counts_1[level] > 0)
        Hk_j_target -= k_mdl * cache->getLog(n_samples);
      else if (cplx == 1)
        Hk_j_target -= cache->getLogC(weighted_count, sc_levels1);
    }

    I[j] = H_j_joint + H_j_target;
    Ik[j] = Hk_j_joint + Hk_j_target;

    TempAllocatorScope scope;
    // copy
    TempVector<int> counts_k_0(counts_0);
    TempVector<int> counts_k_1(counts_1);
    // Before trying to create bins : solution is one single bin from 0 to j.
    memory_cuts_idx[j] = 0;
    memory_cuts_pos[j] = 0;
    // moving k iterator on possible cuts
    for (int k = 0; k < j; ++k) {
      int it_k = n_values[k];
      if (flag_weights)
        eff_n_factors = (sum_weights[j] - sum_weights[k]) / (it_j - it_k);

      std::transform(begin(counts_k_1), end(counts_k_1),
          coarse_counts_target.row_begin(k), begin(counts_k_1),
          std::minus<int>());
      std::transform(begin(counts_k_0), end(counts_k_0),
          coarse_counts_joint.row_begin(k), begin(counts_k_0),
          std::minus<int>());

      double Hk_k_joint{0}, H_k_joint{0};
      for (int level = 0; level < r0; level++) {
        int weighted_count = flag_weights
                                 ? int(eff_n_factors * counts_k_0[level] + 0.5)
                                 : counts_k_0[level];
        Hk_k_joint += cache->getH(weighted_count);
      }
      H_k_joint = Hk_k_joint;

      double Hk_k_target{0}, H_k_target{0};
      for (int level = 0; level < r1; level++) {
        int weighted_count = flag_weights
                                 ? int(eff_n_factors * counts_k_1[level] + 0.5)
                                 : counts_k_1[level];
        Hk_k_target -= cache->getH(weighted_count);
        H_k_target -= cache->getH(weighted_count);

        if (cplx == 0 && counts_k_1[level] > 0)
          Hk_k_target -= k_mdl * cache->getLog(n_samples);
        else if (cplx == 1)
          Hk_k_target -= cache->getLogC(weighted_count, sc_levels1);
      }
      // The information value for a bin fom idx k to j
      double I_kj = H_k_joint + H_k_target;
      double Ik_kj = Hk_k_joint + Hk_k_target;
      if (cplx == 1) {
        // Combinatorial approximation
        Ik_kj -=
            cache->getLogChoose(n_cuts_max - 1, r1_prev - 1) / (r1_prev - 1);
      }

      double Ik_newbin = Ik[k] + Ik_kj;  // [0.. cuts.. k-1] + [k j]
      if ((Ik_newbin - Ik[j]) > kEpsilon) {
        I[j] = I[k] + I_kj;  // max info in the interval [0 j]
        Ik[j] = Ik_newbin;   // max info in the interval [0 j]
        if (memory_cuts_idx[k + 1] == 0) {
          memory_cuts_idx[j] = -k - 1;  // index  of the (last) optimal cut
        } else {
          memory_cuts_idx[j] = k + 1;  // index  of the (last) optimal cut
        }
        memory_cuts_pos[j] = it_k - 1;  // position  of the (last) optimal cut
      }
    }  // inner loop on k
  }    // outer loop on j
  reconstructCutCoarse(memory_cuts_idx, memory_cuts_pos, n_samples, cut);
}

// Initialize Ik(x,y) with equal bin discretization
// Repeat
//   optimize on x Ik(x,y): Hx - Hxy - kmdl
//   optimize on y Ik(x,y): Hy - Hxy - kmdl
// Until convergence
InfoBlock computeIxy(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& is_continuous,
    const TempVector<int>& var_idx, const TempVector<int>& levels,
    const TempVector<double>& weights, bool flag_sample_weights, int initbins,
    int maxbins, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache,
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
    auto var = var_idx[l];
    if (is_continuous[var] != 1) {
      r[l] = levels[var];
      copy(data.row_begin(var), data.row_end(var), datafactors.row_begin(l));
    }
  }
  TempVector<int> r_temp(3);
  int rxy{0};
  TempVector<int> xy_factors(n_samples);
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
    // Initialize cut, datafactors and r
    resetCutPoints(
        levels, is_continuous, var_idx, 0, 2, test_n_bins, n_samples, cut);
    updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);

    TempAllocatorScope scope;
    rxy = setJointFactors(datafactors, r, TempVector<int>{0, 1}, xy_factors);

    r_temp.assign({r[0], r[1], rxy});
    InfoBlock res_temp = computeMI(datafactors.getRow(0), datafactors.getRow(1),
        xy_factors, r_temp, n_eff, weights, cache, cplx, 0);
    // All set if both variables are discrete
    if (is_continuous[var_idx[0]] == 0 && is_continuous[var_idx[1]] == 0)
      return res_temp;

    double Ik = res_temp.I - res_temp.k;
    if (Ik > best_res) {
      best_initbins = test_n_bins;
      best_res = Ik;
    }
  }
  // Initialize cut, datafactors and r
  resetCutPoints(
      levels, is_continuous, var_idx, 0, 2, best_initbins, n_samples, cut);
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
      // max_X{ I(X;Y) }, X starts with a single bin
      optimizeCutPoints(data.getConstRow(var_idx[0]),
          data_idx.getConstRow(var_idx[0]), datafactors.getRow(1),
          TempVector<int>(n_samples, 0), r[1], 1, r[1], r[0],
          levels[var_idx[0]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(0));
    }
    if (is_continuous[var_idx[1]]) {
      TempAllocatorScope scope;
      // max_y{ I(x;y) }, Y starts with a single bin
      optimizeCutPoints(data.getConstRow(var_idx[1]),
          data_idx.getConstRow(var_idx[1]), datafactors.getRow(0),
          TempVector<int>(n_samples, 0), r[0], 1, r[0], r[1],
          levels[var_idx[1]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(1));
    }
    updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);

    if (save_cuts) {
      auto& cp = cuts_info->cutpoints;
      copy(cut.row_begin(0), cut.row_end(0), cp.row_begin(step));
      copy(cut.row_begin(1), cut.row_end(1), cp.row_begin(step) + maxbins);
    }

    TempAllocatorScope scope;
    rxy = setJointFactors(datafactors, r, TempVector<int>{0, 1}, xy_factors);
    r_temp.assign({r[0], r[1], rxy});
    InfoBlock res_temp = computeMI(datafactors.getRow(0), datafactors.getRow(1),
        xy_factors, r_temp, n_eff, weights, cache, cplx, 0);
    // Adding combinatorial term
    if (is_continuous[var_idx[0]] && r[0] > 1) {
      int n_cuts_max = min(maxbins, levels[var_idx[0]]);
      if (r[0] < n_cuts_max)
        res_temp.k += cache->getLogChoose(n_cuts_max - 1, r[0] - 1);
    }
    if (is_continuous[var_idx[1]] && r[1] > 1) {
      int n_cuts_max = min(maxbins, levels[var_idx[1]]);
      if (r[1] < n_cuts_max)
        res_temp.k += cache->getLogChoose(n_cuts_max - 1, r[1] - 1);
    }

    I_list[step] = res_temp.I;
    Ik_list[step] = res_temp.I - res_temp.k;
    bool converged{false};
    for (int i = step - 1; i >= 0; --i) {
      if (fabs(Ik_list[step] - Ik_list[i]) < kEpsilon) {
        converged = true;
        Ixy = accumulate(begin(I_list) + i, begin(I_list) + step, 0.0);
        Ikxy = accumulate(begin(Ik_list) + i, begin(Ik_list) + step, 0.0);
        Ikxy /= (step - i);  // average over the periodic cycle
        Ixy /= (step - i);
        break;
      }
    }

    if (save_cuts) cuts_info->n_iterations = step+1;
    if (converged) break;

    Ixy = res_temp.I;
    Ikxy = res_temp.I - res_temp.k;
    // Already optimal if any one of them is discrete
    if (!is_continuous[var_idx[0]] || !is_continuous[var_idx[1]]) break;
  }  // for step

  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if (!negative_info && Ikxy < 0) {
    Ixy = 0;
    Ikxy = 0;
  }
  if (save_cuts) {
    cuts_info->Ik = Ikxy / n_eff;
    cuts_info->I = Ixy / n_eff;
    cuts_info->I_equal_freq_max = best_res;
  }

  return InfoBlock(n_samples, Ixy, Ixy - Ikxy);
}

InfoBlock computeIxyui(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& is_continuous,
    const TempVector<int>& var_idx, const TempVector<int>& levels,
    const TempVector<double>& weights, bool flag_sample_weights, int initbins,
    int maxbins, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache,
    std::shared_ptr<CutPointsInfo> cuts_info = nullptr) {
  TempAllocatorScope scope;

  bool save_cuts = cuts_info && !cuts_info->cutpoints.empty();
  int n_samples = data.n_cols();
  int n_nodes = var_idx.size();
  int n_ui = n_nodes - 2;
  double n_eff = accumulate(begin(weights), end(weights), 0.0);
  int u_initbins = min(30, max(2, int(0.5 + pow(n_eff, 1.0/(n_nodes)))));
  // allocation factors  x y
  TempGrid2d<int> datafactors(n_nodes, n_samples);
  TempGrid2d<int> cut(n_nodes, maxbins);
  TempGrid2d<int> cut_y_u(n_nodes, maxbins);
  TempGrid2d<int> cut_x_u(n_nodes, maxbins);
  TempVector<int> r(n_nodes);  // n_levels of optimized variables

  // Initialize discrete factors and n_levels
  for (int l = 0; l < n_nodes; ++l) {
    auto var = var_idx[l];
    if (is_continuous[var] != 1) {
      r[l] = levels[var];
      copy(data.row_begin(var), data.row_end(var), datafactors.row_begin(l));
    }
  }
  TempGrid2d<int> uyxfactors(4, n_samples);  // u, uy, ux, uyx
  TempVector<int> ruyx(4);
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

    setUyxJointFactors(datafactors, r, -1, uyxfactors, ruyx);
    r_temp.assign({r[1], ruyx[2], ruyx[3]});
    res_temp = computeMI(datafactors.getRow(1), uyxfactors.getRow(2),
        uyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx, 1);
    double Ik_y_xu = res_temp.I - res_temp.k;

    r_temp.assign({r[0], ruyx[1], ruyx[3]});
    res_temp = computeMI(datafactors.getRow(0), uyxfactors.getRow(1),
        uyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx, 1);
    double Ik_x_yu = res_temp.I - res_temp.k;

    if ((Ik_y_xu + Ik_x_yu) > best_res) {
      best_initbins = test_n_bins;
      best_res = (Ik_y_xu + Ik_x_yu);
    }
  }
  // Initialize X and Y cuts with best_initbins
  resetCutPoints(
      levels, is_continuous, var_idx, 0, 2, best_initbins, n_samples, cut);
  updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);
  TempVector<int> r_old(r);  // n_levels in the previous iteration

  bool reuse_x_u_cuts = false;
  bool reuse_y_u_cuts = false;

  int iter_max = save_cuts ? cuts_info->cutpoints.n_rows() : kMaxIter;
  // Keep the result of each iteration
  TempVector<double> I_list(iter_max);
  TempVector<double> Ik_list(iter_max);
  double Ixy_ui{0}, Ikxy_ui{0};  // to be returned
  // Run optimization with best initial equal freq.
  for (int step = 0; step < iter_max; ++step) {
    // optimize I(y;xu) over x and u
    // Either reset cutpoints or re-use I(y;u) cutpoints
    if(reuse_y_u_cuts){
      for(int l = 2; l < n_nodes; ++l){
        copy(cut_y_u.row_begin(l), cut_y_u.row_end(l), cut.row_begin(l));
      }
    } else {
      resetCutPoints(
          levels, is_continuous, var_idx, 2, n_nodes, u_initbins,
          n_samples, cut);
    }
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
    copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; (count < kMaxIterOnU) && !reuse_y_u_cuts; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] == 0) continue;
        // opt u, I(y;xu)
        setUyxJointFactors(datafactors, r_old, l, uyxfactors, ruyx);
        // init variables for the optimization run
        int r0 = ruyx[3];  // xyu
        int r1 = ruyx[2];  // xu
        int sc_levels1 = r_old[1];
        int sc_levels2 = r_old[l];  // old nlevels for combinatorial term
        // Run optimization on U. 2 factors xyu and xu
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uyxfactors.getRow(3),
            uyxfactors.getRow(2), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors(
          data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
      copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    setUyxJointFactors(datafactors, r_old, -1, uyxfactors, ruyx);
    r_temp.assign({r_old[1], ruyx[2], ruyx[3]});
    res_temp = computeMI(datafactors.getRow(1), uyxfactors.getRow(2),
        uyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx, 1);
    double I_y_xu = res_temp.I;  // Before optimization on X.
    double Ik_y_xu = res_temp.I - res_temp.k;
    if (is_continuous[var_idx[0]] && r_old[0] > 1) {
      int n_cuts_max = min(levels[var_idx[0]], maxbins);
      if (r_old[0] < n_cuts_max)
        Ik_y_xu -= cache->getLogChoose(n_cuts_max - 1, r_old[0] - 1);
    }

    if (is_continuous[var_idx[0]] == 1) {
      // I(y;xu), optimize on x
      int r0 = ruyx[1];  // uy
      int r1 = ruyx[0];  // u
      int sc_levels1 = r_old[1];
      int sc_levels2 = r_old[0];
      // Run optimization on X. 2 factors uy and u
      optimizeCutPoints(data.getConstRow(var_idx[0]),
          data_idx.getConstRow(var_idx[0]), uyxfactors.getRow(1),
          uyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
          levels[var_idx[0]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(0));
    }

    // optimize I(x;yu) over y and u
    // Either reset cutpoints or re-use I(x;u) cutpoints
    if(reuse_x_u_cuts){ 
      for(int l = 2; l < n_nodes; ++l){
        copy(cut_x_u.row_begin(l), cut_x_u.row_end(l), cut.row_begin(l));
      }
    } else {
      resetCutPoints(
          levels, is_continuous, var_idx, 2, n_nodes, u_initbins,
          n_samples, cut);
    }
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
    copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; (count < kMaxIterOnU) && !reuse_x_u_cuts; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] == 0) continue;

        setUyxJointFactors(datafactors, r_old, l, uyxfactors, ruyx);
        // init variables for the optimization run
        int r0 = ruyx[3];  // xyu
        int r1 = ruyx[1];  // yu
        int sc_levels1 = r_old[0];
        int sc_levels2 = r_old[l];
        // Run optimization on U. 2 factors xyu and yu
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uyxfactors.getRow(3),
            uyxfactors.getRow(1), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors(
          data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
      copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    setUyxJointFactors(datafactors, r_old, -1, uyxfactors, ruyx);
    r_temp.assign({r_old[0], ruyx[1], ruyx[3]});
    res_temp = computeMI(datafactors.getRow(0), uyxfactors.getRow(1),
        uyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx, 1);
    double I_x_yu = res_temp.I;  // Before updating Y (and X).
    double Ik_x_yu = res_temp.I - res_temp.k;
    if ((is_continuous[var_idx[1]] == 1) && (r_old[1] > 1)) {
      int n_cuts_max = min(levels[var_idx[1]], maxbins);
      if (r_old[1] < n_cuts_max) {
        Ik_x_yu -= cache->getLogChoose(n_cuts_max - 1, r_old[1] - 1);
      }
    }

    if (is_continuous[var_idx[1]] == 1) {
      // I(x;yu), optimize on y
      int r0 = ruyx[2];  // ux
      int r1 = ruyx[0];  // u
      int sc_levels1 = r_old[0];
      int sc_levels2 = r_old[1];
      // Run optimization on Y. 2 factors ux and u
      optimizeCutPoints(data.getConstRow(var_idx[1]),
          data_idx.getConstRow(var_idx[1]), uyxfactors.getRow(2),
          uyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
          levels[var_idx[1]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(1));
    }

    // optimize I(x;u) over u
    // Either reset cutpoints or re-use I(y;u) cutpoints
    if(reuse_x_u_cuts){ 
      for(int l = 2; l < n_nodes; ++l){
        copy(cut_x_u.row_begin(l), cut_x_u.row_end(l), cut.row_begin(l));
      }
    } else {
      resetCutPoints(
          levels, is_continuous, var_idx, 2, n_nodes, u_initbins,
          n_samples, cut);
    }
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
    copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; (count < kMaxIterOnU) && !reuse_x_u_cuts; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] == 0) continue;

        setUyxJointFactors(datafactors, r_old, l, uyxfactors, ruyx);
        // init variables for the optimization run
        int r0 = ruyx[2];           // xu
        int r1 = ruyx[0];           // u
        int sc_levels1 = r_old[0];  // x
        int sc_levels2 = r_old[l];  // u
        // optimization run on var_idx[l], 2 factors xu and u
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uyxfactors.getRow(2),
            uyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors(
          data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
      copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    //save u cutpoints
    for(int l = 2; l < n_nodes; ++l){
      copy(cut.row_begin(l), cut.row_end(l), cut_x_u.row_begin(l));
    }
    setUyxJointFactors(datafactors, r_old, -1, uyxfactors, ruyx);
    r_temp.assign({r_old[0], ruyx[0], ruyx[2]});
    res_temp = computeMI(datafactors.getRow(0), uyxfactors.getRow(0),
        uyxfactors.getRow(2), r_temp, n_eff, weights, cache, cplx, 1);
    double I_x_u = res_temp.I;  // After optimization on U.
    double Ik_x_u = res_temp.I - res_temp.k;

    // optimize I(y;u) over u
    // Either reset cutpoints or re-use I(y;u) cutpoints
    if(reuse_y_u_cuts){
      for(int l = 2; l < n_nodes; ++l){
        copy(cut_y_u.row_begin(l), cut_y_u.row_end(l), cut.row_begin(l));
      }
    } else {
      resetCutPoints(
          levels, is_continuous, var_idx, 2, n_nodes, u_initbins,
          n_samples, cut);
    }
    updateFactors(
        data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
    copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; (count < kMaxIterOnU) && !reuse_y_u_cuts; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] == 0) continue;

        setUyxJointFactors(datafactors, r_old, l, uyxfactors, ruyx);
        int r0 = ruyx[1];           // yu
        int r1 = ruyx[0];           // u
        int sc_levels1 = r_old[1];  // y
        int sc_levels2 = r_old[l];  // u
        // optimization run on var_idx[l], 2 factors yu and u
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uyxfactors.getRow(1),
            uyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors(
          data_idx, cut, is_continuous, var_idx, 2, n_nodes, datafactors, r);
      copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    //save u cutpoints
    for(int l = 2; l < n_nodes; ++l){
      copy(cut.row_begin(l), cut.row_end(l), cut_y_u.row_begin(l));
    }
    setUyxJointFactors(datafactors, r_old, -1, uyxfactors, ruyx);
    r_temp.assign({r_old[1], ruyx[0], ruyx[1]});
    res_temp = computeMI(datafactors.getRow(1), uyxfactors.getRow(0),
        uyxfactors.getRow(1), r_temp, n_eff, weights, cache, cplx, 1);
    double I_y_u = res_temp.I;  // After optimization on U.
    double Ik_y_u = res_temp.I - res_temp.k;

    // End of iteration: update X and Y cuts
    updateFactors(data_idx, cut, is_continuous, var_idx, 0, 2, datafactors, r);
    copy(begin(r), end(r), begin(r_old));

    if (save_cuts) {
      auto& cp = cuts_info->cutpoints;
      copy(cut.row_begin(0), cut.row_end(0), cp.row_begin(step));
      copy(cut.row_begin(1), cut.row_end(1), cp.row_begin(step) + maxbins);
    }

    // Compute I(X;Y|U)
    double part_one = I_x_yu - I_x_u;
    double part_two = I_y_xu - I_y_u;
    // Either sum may be negative when we can't find U bins on I(X;YU) (or 
    //I(Y;XY)) that are as good as those found on I(X;U) (I(Y;U)). If that is
    //the case, we reuse the better U bins for all terms for the next iteration.
    reuse_x_u_cuts = part_one < kPrecision || (reuse_x_u_cuts && (part_two < kPrecision));
    reuse_y_u_cuts = part_two < kPrecision || (reuse_y_u_cuts && (part_one < kPrecision));

    if(reuse_x_u_cuts || reuse_y_u_cuts) {
      bool one_bin = r[0] == 1 || r[1] == 1;
      if(one_bin && (step != (iter_max-1))) {
        // Next iteration will have 0 information with no further optimization possible.
        iter_max = step+2;
        if (save_cuts) cuts_info->n_iterations = step+1;
      }
      if(!one_bin) continue;
      // With one bin, we can keep the result and check for convergence
    }

    double cond_I = 0.5 * (part_one + part_two);
    double cond_Ik = 0.5 * (Ik_x_yu - Ik_x_u + Ik_y_xu - Ik_y_u);
    I_list[step] = cond_I;
    Ik_list[step] = cond_Ik;

    // Test stop condition on stop1
    bool converged{false};
    for (int i = step - 1; i >= 0; i--) {
      // If no real improvement over last information
      if (fabs(cond_Ik - Ik_list[i]) < kEpsilon) {
        converged = true;
        Ixy_ui = accumulate(begin(I_list) + i, begin(I_list) + step, 0.0);
        Ikxy_ui = accumulate(begin(Ik_list) + i, begin(Ik_list) + step, 0.0);
        Ikxy_ui /= (step - i);  // average over the periodic cycle
        Ixy_ui /= (step - i);
        break;
      }
    }
    if (save_cuts) cuts_info->n_iterations = step;
    if (converged) break;

    Ixy_ui = cond_I;
    Ikxy_ui = cond_Ik;
  }  // for step
  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if (!negative_info && Ikxy_ui < 0) {
    Ixy_ui = 0;
    Ikxy_ui = 0;
  }
  if (save_cuts) {
    cuts_info->Ik = Ikxy_ui / n_eff;
    cuts_info->I = Ixy_ui / n_eff;
  }

  return InfoBlock(n_samples, Ixy_ui, Ixy_ui - Ikxy_ui);
}

}  // anonymous namespace

InfoBlock computeCondMutualInfo(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& levels,
    const TempVector<int>& is_continuous, const TempVector<int>& var_idx,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    int initbins, int maxbins, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache,
    std::shared_ptr<CutPointsInfo> cuts_info) {
  if (data.n_rows() == 2) {
    return computeIxy(data, data_idx, is_continuous, var_idx, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache, cuts_info);
  } else {
    return computeIxyui(data, data_idx, is_continuous, var_idx, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache, cuts_info);
  }
}

// compute Rscore and three point mutual information I(x;y;z | u)
// return Info3PointBlock{score, N * Ixyz_ui, N * kxyz_ui}
Info3PointBlock computeInfo3PointAndScore(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& levels,
    const TempVector<int>& is_continuous, const TempVector<int>& var_idx,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    int initbins, int maxbins, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_ui = data.n_rows() - 3;
  // Optimize variables for each MI estimation for the R score
  // I(x,y|u,z)
  InfoBlock res_temp = computeIxyui(data, data_idx, is_continuous, var_idx,
      levels, sample_weights, flag_sample_weights, initbins, maxbins, cplx,
      negative_info, cache);
  double I_xy_zu = res_temp.I;
  double Ik_xy_zu = res_temp.I - res_temp.k;

  TempVector<int> var_idx_t(begin(var_idx), begin(var_idx) + n_ui + 2);
  // Do opt run on I(X;Y|U)
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  }
  double I_xy_u = res_temp.I;
  double Ik_xy_u = res_temp.I - res_temp.k;

  // I(x,z|u)
  var_idx_t[1] = var_idx.back();  // Z
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  }
  double Ik_xz_u = res_temp.I - res_temp.k;

  // I(y,z|u)
  var_idx_t[0] = var_idx[1];  // Y
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  }
  double Ik_yz_u = res_temp.I - res_temp.k;

  double xz = Ik_xz_u - Ik_xy_u;
  double yz = Ik_yz_u - Ik_xy_u;
  // Data processing inequality
  double dpi = std::fmin(xz, yz) - std::log1p(exp(-std::fabs(xz - yz)));

  // Conditional three point information
  double I_xyz_u = I_xy_u - I_xy_zu;
  // Ik_xyz_u can be seen as the probability of non-v-structure
  double Ik_xyz_u = Ik_xy_u - Ik_xy_zu;

  double Rscore = std::fmin(dpi, Ik_xyz_u);
  return Info3PointBlock{Rscore, I_xyz_u, I_xyz_u - Ik_xyz_u};
}

}  // namespace computation
}  // namespace miic
