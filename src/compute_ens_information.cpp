#include "compute_ens_information.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <tuple>

#include "compute_info.h"
#include "info_cnt.h"
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

  double info{0}, score{std::numeric_limits<double>::lowest()};
  if (std::all_of(begin(cnt_red), end(cnt_red), [](int x) { return x == 0; })) {
    double* res = getAllInfoNEW(environment.data_numeric, environment.levels, X,
        Y, Z, ui_list, environment.n_eff, environment.cplx, environment.is_k23,
        environment.sample_weights, environment.test_mar,
        environment.cache.cterm);
    info = res[7] + res[8];  // I(X;Y;Z|ui) - k(X;Y;Z|ui), res[8] = -k
    score = res[6];          // R(X,Y;Z|ui)
    delete[] res;
  } else {
    double* res = compute_Rscore_Ixyz_alg5(data_numeric_red,
        data_numeric_idx_red, all_levels_red, cnt_red, posArray_red, n_ui,
        n_ui + 2, n_samples_non_na_z, sample_weights_red, flag_sample_weights,
        environment);
    info = res[1] - res[2];  // I(x;y;z|u) - cplx I(x;y;z|u)
    score = res[0];          // R(X,Y;Z|ui)
    delete[] res;
  }
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

double* getCondMutualInfo(
    Environment& environment, int X, int Y, const vector<int>& ui_list) {
  TempAllocatorScope scope;
  // TODO : speedup by only removing NAs for marked columns
  // Mark rows containing NAs and count the number of complete samples
  TempVector<int> sample_is_not_NA(environment.n_samples);
  TempVector<int> NAs_count(environment.n_samples);
  int n_samples_non_na = countNonNA(X, Y, /*Z*/ -1, ui_list,
      environment.data_numeric, sample_is_not_NA, NAs_count);

  if (n_samples_non_na <= 2) {
    double* res_new = new double[3];
    res_new[0] = n_samples_non_na;
    res_new[1] = 0;  // I(X;Y|ui)
    res_new[2] = 0;  // cplx
    return res_new;
  }

  int n_ui = ui_list.size();
  // Allocate data reducted *_red without rows containing NAs
  // All *_red variables are passed to the optimization routine
  TempVector<int> all_levels_red(n_ui + 2);
  TempVector<int> cnt_red(n_ui + 2);
  TempVector<int> posArray_red(n_ui + 2);
  TempVector<double> sample_weights_red(n_samples_non_na);
  TempGrid2d<int> data_numeric_red(n_ui + 2, n_samples_non_na);
  TempGrid2d<int> data_numeric_idx_red(n_ui + 2, n_samples_non_na);

  bool flag_sample_weights = filterNA(X, Y, /*Z*/ -1, ui_list,
      environment.data_numeric, environment.data_numeric_idx,
      environment.is_continuous, environment.sample_weights, sample_is_not_NA,
      NAs_count, data_numeric_red, data_numeric_idx_red, all_levels_red,
      cnt_red, posArray_red, sample_weights_red);

  // If X or Y has only 1 level
  if (all_levels_red[0] <= 1 || all_levels_red[1] <= 1) {
    double* res_new = new double[3];
    res_new[0] = n_samples_non_na;
    res_new[1] = 0;  // I(X;Y|ui)
    res_new[2] = 0;  // cplx
    return res_new;
  }

  double* res_new = NULL;
  if (std::all_of(begin(cnt_red), end(cnt_red), [](int x) { return x == 0; })) {
    res_new = getAllInfoNEW(environment.data_numeric,
        environment.levels, X, Y, -1, ui_list, environment.n_eff,
        environment.cplx, environment.is_k23, environment.sample_weights,
        environment.test_mar, environment.cache.cterm);
  } else {
    res_new = compute_mi_cond_alg1(data_numeric_red, data_numeric_idx_red,
        all_levels_red, cnt_red, posArray_red, n_ui, n_samples_non_na,
        sample_weights_red, flag_sample_weights, environment);
    res_new[1] = res_new[1] * res_new[0];  // Ixy|u
    res_new[2] = res_new[2] * res_new[0];  // cplx
  }
  for (int i = 0; i < 3; i++) {
    if (std::fabs(res_new[i]) < kPrecision) res_new[i] = 0.0;
  }

  return res_new;
}

// parallel: whether the search can be done in parallel
void searchForBestContributingNode(
    Environment& environment, int X, int Y, bool parallel) {
  auto info = environment.edges[X][Y].shared_info;
  auto& zi_list = info->zi_list;
  if (!environment.latent) {
    // remove zi that is not connected to neither x nor y
    zi_list.erase(std::remove_if(begin(zi_list), end(zi_list),
                      [&environment, X, Y](int Z) {
                        const auto& edges = environment.edges;
                        return !edges[X][Z].status && !edges[Y][Z].status;
                      }),
        zi_list.end());
  }
  int n_zi = zi_list.size();

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
      info->z_name_idx = i;
      info->Rxyz_ui = score;
    }
  }
}

}  // namespace computation
}  // namespace miic
