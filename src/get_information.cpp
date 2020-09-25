#include "get_information.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <tuple>

#include "computation_continuous.h"
#include "computation_discrete.h"
#include "structure.h"
#include "utilities.h"

namespace miic {
namespace computation {

using namespace miic::structure;
using namespace miic::utility;
using std::vector;

constexpr double kPrecision = 1.e-10;

// Return either (when bool get_info == true)
//   shifted 3-point information I(X;Y;Z|ui) - k(X;Y;Z|ui),
// or (when bool get_info == false)
//   score of Z as a candidate separating node w.r.t. (X,Y|ui),
// since the two quantities are both of type double and are always computed at
// the same time, but are never required at the same time.
double getInfo3PointOrScore(Environment& environment, int X, int Y, int Z,
    const vector<int>& ui_list, bool get_info) {
  TempAllocatorScope scope;
  auto& cache = environment.cache.info_score;

  if (get_info) {
    double info{0};
    bool found{false};
    std::tie(info, found) = cache->getInfo3Point(X, Y, Z, ui_list);
    if (found) return info;
  } else {
    double score{std::numeric_limits<double>::lowest()};
    bool found{false};
    std::tie(score, found) = cache->getScore(X, Y, Z, ui_list);
    if (found) return score;
  }

  bool any_na = environment.has_na[X] || environment.has_na[Y] ||
                environment.has_na[Z] ||
                std::any_of(begin(ui_list), end(ui_list),
                    [&environment](int u) { return environment.has_na[u]; });
  // Mark rows containing NAs and count the number of complete samples
  // 1: sample contains no NA, 0: sample contains NA
  TempVector<int> sample_is_not_NA(environment.n_samples, 1);
  // vector with the number of rows containing NAs seen at rank i
  TempVector<int> na_count(environment.n_samples, 0);
  int n_samples_non_na_z = environment.n_samples;
  if (any_na)
    n_samples_non_na_z = countNonNA(
        X, Y, Z, ui_list, environment.data_numeric, sample_is_not_NA, na_count);

  if (n_samples_non_na_z <= 2) {  // not sufficient statistics
    if (get_info) {
      cache->saveInfo3Point(X, Y, Z, ui_list, 0);
      return 0;
    } else {
      double score = std::numeric_limits<double>::lowest();
      cache->saveScore(X, Y, Z, ui_list, score);
      return score;
    }
  }

  int n_ui = ui_list.size();
  // Allocate data reducted *_red without rows containing NAs
  // All *_red variables are passed to the optimization routine
  TempGrid2d<int> data_red(n_ui + 3, n_samples_non_na_z);
  TempGrid2d<int> data_idx_red(n_ui + 3, n_samples_non_na_z);
  TempVector<int> levels_red(n_ui + 3);
  TempVector<int> is_continuous_red(n_ui + 3);
  TempVector<int> var_idx_red(n_ui + 3);
  TempVector<double> weights_red(n_samples_non_na_z);
  bool flag_sample_weights = filterNA(X, Y, Z, ui_list,
      environment.data_numeric, environment.data_numeric_idx,
      environment.levels, environment.is_continuous, environment.sample_weights,
      sample_is_not_NA, na_count, data_red, data_idx_red, levels_red,
      is_continuous_red, var_idx_red, weights_red, any_na);

  // If X or Y has 1 or less level, the 3-point information should be zero, and
  // the score should be low
  bool ok = levels_red[0] > 1 && levels_red[1] > 1;

  if (any_na) {
    int n_samples_nonNA =
        getNumSamplesNonNA(environment.data_numeric, X, Y, ui_list);
    if (ok && environment.test_mar && n_samples_nonNA != n_samples_non_na_z) {
      double kldiv = compute_kl_divergence(environment.data_numeric,
          environment.data_double, X, Y, ui_list, environment.levels,
          environment.is_continuous, n_samples_non_na_z, levels_red,
          sample_is_not_NA, environment.noise_vec);
      double cplxMdl = environment.cache.cterm->getLog(n_samples_non_na_z);

      if ((kldiv - cplxMdl) > 0) {
        // The sample is not representative of the population, hence for 3-point
        // information, we cannot draw conclusion the unshielded triple (X, Z,
        // Y), return 0; For contributing score, Z is not a good candidate.
        ok = false;
      }
    }
  }
  if (!ok) {
    if (get_info) {
      cache->saveInfo3Point(X, Y, Z, ui_list, 0);
      return 0;
    } else {
      double score = std::numeric_limits<double>::lowest();
      cache->saveScore(X, Y, Z, ui_list, score);
      return score;
    }
  }

  Info3PointBlock res{0, 0, 0};
  if (std::all_of(begin(is_continuous_red), end(is_continuous_red),
          [](int x) { return x == 0; })) {
    res = computeInfo3PointAndScoreDiscrete(data_red, levels_red, var_idx_red,
        weights_red, environment.cplx, environment.cache.cterm);
  } else {
    res = computeInfo3PointAndScore(data_red, data_idx_red, levels_red,
        is_continuous_red, var_idx_red, weights_red, flag_sample_weights,
        environment.initbins, environment.maxbins, environment.cplx,
        environment.cache.cterm);
  }
  double info = res.Ixyz_ui - res.kxyz_ui;  // I(x;y;z|u) - cplx I(x;y;z|u)
  double score = res.score;                 // R(X,Y;Z|ui)
  if (std::fabs(info) < kPrecision) info = 0;
  if (std::fabs(score) < kPrecision) score = 0;

  if (get_info) {
    cache->saveInfo3Point(X, Y, Z, ui_list, info);
    return info;
  } else {
    cache->saveScore(X, Y, Z, ui_list, score);
    return score;
  }
}

InfoBlock getCondMutualInfo(int X, int Y, const vector<int>& ui_list,
    const Grid2d<int>& data_numeric, const Grid2d<int>& data_numeric_idx,
    Environment& environment) {
  TempAllocatorScope scope;
  auto& cache = environment.cache.info_score;

  int n_ui = ui_list.size();
  if (n_ui != 0) {
    InfoBlock res{0, 0, 0};
    bool found{false};
    std::tie(res, found) = cache->getMutualInfo(X, Y, ui_list);
    if (found) return res;
  }

  bool any_na = environment.has_na[X] || environment.has_na[Y] ||
                std::any_of(begin(ui_list), end(ui_list),
                    [&environment](int u) { return environment.has_na[u]; });
  // Mark rows containing NAs and count the number of complete samples
  // 1: sample contains no NA, 0: sample contains NA
  TempVector<int> sample_is_not_NA(environment.n_samples, 1);
  // vector with the number of rows containing NAs seen at rank i
  TempVector<int> na_count(environment.n_samples, 0);
  int n_samples_non_na = environment.n_samples;
  if (any_na)
    n_samples_non_na = countNonNA(
        X, Y, /*Z*/ -1, ui_list, data_numeric, sample_is_not_NA, na_count);

  if (n_samples_non_na <= 2) {
    InfoBlock res{n_samples_non_na, 0, 0};
    if (n_ui != 0) cache->saveMutualInfo(X, Y, ui_list, res);
    return res;
  }

  // Allocate data reducted *_red without rows containing NAs
  // All *_red variables are passed to the optimization routine
  TempGrid2d<int> data_red(n_ui + 2, n_samples_non_na);
  TempGrid2d<int> data_idx_red(n_ui + 2, n_samples_non_na);
  TempVector<int> levels_red(n_ui + 2);
  TempVector<int> is_continuous_red(n_ui + 2);
  TempVector<int> var_idx_red(n_ui + 2);
  TempVector<double> weights_red(n_samples_non_na);

  bool flag_sample_weights = filterNA(X, Y, /*Z*/ -1, ui_list, data_numeric,
      data_numeric_idx, environment.levels, environment.is_continuous,
      environment.sample_weights, sample_is_not_NA, na_count, data_red,
      data_idx_red, levels_red, is_continuous_red, var_idx_red, weights_red,
      any_na);

  // If X or Y has only 1 level
  if (levels_red[0] <= 1 || levels_red[1] <= 1) {
    InfoBlock res{n_samples_non_na, 0, 0};
    if (n_ui != 0) cache->saveMutualInfo(X, Y, ui_list, res);
    return res;
  }

  InfoBlock res{0, 0, 0};
  if (std::all_of(begin(is_continuous_red), end(is_continuous_red),
          [](int x) { return x == 0; })) {
    res = computeCondMutualInfoDiscrete(data_red, levels_red, var_idx_red,
        weights_red, environment.cplx, environment.cache.cterm);
  } else {
    res = computeCondMutualInfo(data_red, data_idx_red, levels_red,
        is_continuous_red, var_idx_red, weights_red, flag_sample_weights,
        environment.initbins, environment.maxbins, environment.cplx,
        environment.cache.cterm);
  }
  if (std::fabs(res.I) < kPrecision) res.I = 0;
  if (std::fabs(res.k) < kPrecision) res.k = 0;

  if (n_ui != 0) cache->saveMutualInfo(X, Y, ui_list, res);
  return res;
}

// parallel: whether the search can be done in parallel
void searchForBestContributingNode(
    Environment& environment, int X, int Y, bool parallel) {
  auto info = environment.edges(X, Y).shared_info;
  auto& zi_list = info->zi_list;
  if (!environment.latent) {
    auto is_isolated = [&environment, X, Y](int Z) {
      return !environment.edges(X, Z).status && !environment.edges(Y, Z).status;
    };
    // remove zi that is not connected to neither x nor y
    zi_list.erase(
        remove_if(begin(zi_list), end(zi_list), is_isolated), end(zi_list));
  }

  int n_zi = zi_list.size();
  info->Rxyz_ui = 0;
#ifdef _OPENMP
#pragma omp parallel for if (parallel && n_zi > environment.n_threads)
#endif
  for (int i = 0; i < n_zi; i++) {
    int Z = zi_list[i];
    double score = getInfo3PointOrScore(
        environment, X, Y, Z, info->ui_list, /* get_info = */ false);
#ifdef _OPENMP
#pragma omp critical
#endif
    if (score > info->Rxyz_ui) {
      info->top_z = Z;
      info->Rxyz_ui = score;
    }
  }
}

}  // namespace computation
}  // namespace miic
