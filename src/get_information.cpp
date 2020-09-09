#include "get_information.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <tuple>

#include "computation_discrete.h"
#include "computation_continuous.h"
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

  // Mark rows containing NAs and count the number of complete samples
  // 1: sample contains no NA, 0: sample contains NA
  TempVector<int> sample_is_not_NA(environment.n_samples);
  // vector with the number of rows containing NAs seen at rank i
  TempVector<int> NAs_count(environment.n_samples);
  int n_samples_non_na_z = countNonNA(
      X, Y, Z, ui_list, environment.data_numeric, sample_is_not_NA, NAs_count);

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
  TempVector<double> sample_weights_red(n_samples_non_na_z);
  TempGrid2d<int> data_numeric_idx_red(n_ui + 3, n_samples_non_na_z);
  TempGrid2d<int> data_numeric_red(n_ui + 3, n_samples_non_na_z);
  TempVector<int> all_levels_red(n_ui + 3);
  TempVector<int> cnt_red(n_ui + 3);
  TempVector<int> posArray_red(n_ui + 3);

  bool flag_sample_weights = filterNA(X, Y, Z, ui_list,
      environment.data_numeric, environment.data_numeric_idx,
      environment.is_continuous, environment.sample_weights, sample_is_not_NA,
      NAs_count, data_numeric_red, data_numeric_idx_red, all_levels_red,
      cnt_red, posArray_red, sample_weights_red);

  // If X or Y has 1 or less level, the 3-point information should be zero, and
  // the score should be low
  bool ok = all_levels_red[0] > 1 && all_levels_red[1] > 1;

  int n_samples_nonNA = getNumSamplesNonNA(environment, X, Y);
  if (ok && environment.test_mar && n_samples_nonNA != n_samples_non_na_z) {
    double kldiv = compute_kl_divergence(X, Y, environment, n_samples_non_na_z,
        all_levels_red, sample_is_not_NA);
    double cplxMdl = environment.cache.cterm->getLog(n_samples_non_na_z);

    if ((kldiv - cplxMdl) > 0) {
      // The sample is not representative of the population, hence for 3-point
      // information, we cannot draw conclusion the unshielded triple (X, Z, Y),
      // return 0; For contributing score, Z is not a good candidate.
      ok = false;
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
  if (std::all_of(begin(cnt_red), end(cnt_red), [](int x) { return x == 0; })) {
    res = computeInfo3PointAndScoreDiscrete(data_numeric_red, all_levels_red,
        posArray_red, sample_weights_red, environment.cplx,
        environment.cache.cterm);
  } else {
    res = computeInfo3PointAndScore(data_numeric_red, data_numeric_idx_red,
        all_levels_red, cnt_red, posArray_red, n_ui, n_ui + 2,
        n_samples_non_na_z, sample_weights_red, flag_sample_weights,
        environment);
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
    const vector<vector<int>>& data_numeric,
    const vector<vector<int>>& data_numeric_idx, Environment& environment) {
  TempAllocatorScope scope;
  auto& cache = environment.cache.info_score;

  int n_ui = ui_list.size();
  if (n_ui != 0) {
    InfoBlock res{0, 0, 0};
    bool found{false};
    std::tie(res, found) = cache->getMutualInfo(X, Y, ui_list);
    if (found) return res;
  }

  // TODO : speedup by only removing NAs for marked columns
  // Mark rows containing NAs and count the number of complete samples
  TempVector<int> sample_is_not_NA(environment.n_samples);
  TempVector<int> NAs_count(environment.n_samples);
  int n_samples_non_na = countNonNA(
      X, Y, /*Z*/ -1, ui_list, data_numeric, sample_is_not_NA, NAs_count);

  if (n_samples_non_na <= 2) {
    InfoBlock res{n_samples_non_na, 0, 0};
    if (n_ui != 0) cache->saveMutualInfo(X, Y, ui_list, res);
    return res;
  }

  // Allocate data reducted *_red without rows containing NAs
  // All *_red variables are passed to the optimization routine
  TempVector<int> all_levels_red(n_ui + 2);
  TempVector<int> cnt_red(n_ui + 2);
  TempVector<int> posArray_red(n_ui + 2);
  TempVector<double> sample_weights_red(n_samples_non_na);
  TempGrid2d<int> data_numeric_red(n_ui + 2, n_samples_non_na);
  TempGrid2d<int> data_numeric_idx_red(n_ui + 2, n_samples_non_na);

  bool flag_sample_weights = filterNA(X, Y, /*Z*/ -1, ui_list, data_numeric,
      data_numeric_idx, environment.is_continuous, environment.sample_weights,
      sample_is_not_NA, NAs_count, data_numeric_red, data_numeric_idx_red,
      all_levels_red, cnt_red, posArray_red, sample_weights_red);

  // If X or Y has only 1 level
  if (all_levels_red[0] <= 1 || all_levels_red[1] <= 1) {
    InfoBlock res{n_samples_non_na, 0, 0};
    if (n_ui != 0) cache->saveMutualInfo(X, Y, ui_list, res);
    return res;
  }

  InfoBlock res{0, 0, 0};
  if (std::all_of(begin(cnt_red), end(cnt_red), [](int x) { return x == 0; })) {
    res = computeCondMutualInfoDiscrete(data_numeric_red, all_levels_red,
        posArray_red, sample_weights_red, environment.cplx,
        environment.cache.cterm);
  } else {
    res = computeCondMutualInfo(data_numeric_red, data_numeric_idx_red,
        all_levels_red, cnt_red, posArray_red, n_ui, n_samples_non_na,
        sample_weights_red, flag_sample_weights, environment);
  }
  if (std::fabs(res.Ixy_ui) < kPrecision) res.Ixy_ui = 0;
  if (std::fabs(res.kxy_ui) < kPrecision) res.kxy_ui = 0;

  if (n_ui != 0) cache->saveMutualInfo(X, Y, ui_list, res);
  return res;
}

// parallel: whether the search can be done in parallel
void searchForBestContributingNode(
    Environment& environment, int X, int Y, bool parallel) {
  auto info = environment.edges(X, Y).shared_info;
  auto& zi_list = info->zi_list;
  if (!environment.latent) {
    // remove zi that is not connected to neither x nor y
    zi_list.erase(std::remove_if(begin(zi_list), end(zi_list),
                      [&environment, X, Y](int Z) {
                        const auto& edges = environment.edges;
                        return !edges(X, Z).status && !edges(Y, Z).status;
                      }),
        zi_list.end());
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
