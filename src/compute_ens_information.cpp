#include "compute_ens_information.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <limits>

#include "compute_info.h"
#include "info_cnt.h"
#include "structure.h"
#include "utilities.h"

namespace miic {
namespace computation {

using namespace miic::structure;
using namespace miic::utility;
using std::vector;

namespace {

double computeContributingScores(Environment& environment, int X, int Y, int Z,
    const vector<int>& ui_list, int n_samples_nonNA) {
  TempAllocatorScope scope;

  double output_score = lookupScore(X, Y, ui_list, Z, environment);
  if (output_score != -1) return output_score;

  // Mark rows containing NAs and count the number of complete samples
  TempVector<int> sample_is_not_NA(environment.n_samples);
  TempVector<int> NAs_count(environment.n_samples);
  int n_samples_non_na_z =
      count_non_NAs(X, Y, ui_list, sample_is_not_NA, NAs_count, environment, Z);
  if (n_samples_non_na_z <= 2) {
    double output_score = std::numeric_limits<double>::lowest();
    saveScore(X, Y, ui_list, Z, output_score, environment);
    return output_score;
  }

  int n_ui = ui_list.size();
  // Allocate data reduced *_red without rows containing NAs
  // All *_red variables are passed to the optimization routine
  TempVector<double> sample_weights_red(n_samples_non_na_z);
  TempGrid2d<int> data_numeric_idx_red(n_ui + 3, n_samples_non_na_z);
  TempGrid2d<int> data_numeric_red(n_ui + 3, n_samples_non_na_z);
  TempVector<int> all_levels_red(n_ui + 3);
  TempVector<int> cnt_red(n_ui + 3);
  TempVector<int> posArray_red(n_ui + 3);

  bool flag_sample_weights = filter_NAs(X, Y, ui_list, all_levels_red, cnt_red,
      posArray_red, data_numeric_red, data_numeric_idx_red, sample_weights_red,
      sample_is_not_NA, NAs_count, environment, Z);

  // we do not want to add a z if x or y have only one bin
  bool ok = all_levels_red[0] > 1 && all_levels_red[1] > 1;

  if (ok && environment.test_mar && n_samples_nonNA != n_samples_non_na_z) {
    double kldiv = compute_kl_divergence(X, Y, environment, n_samples_non_na_z,
        all_levels_red, sample_is_not_NA);
    double cplxMdl = environment.cache.cterm->getLog(n_samples_non_na_z);

    if ((kldiv - cplxMdl) > 0) {
      // the sample is not representative of the population, hence we do
      // not want this z as possible z
      ok = false;
    }
  }
  if (!ok) {
    double output_score = std::numeric_limits<double>::lowest();
    saveScore(X, Y, ui_list, Z, output_score, environment);
    return output_score;
  }

  if (std::all_of(begin(cnt_red), end(cnt_red), [](int x) { return x == 0; })) {
    TempGrid2d<double> jointFreqs = getJointFreqs(environment, X, Y);
    // call discrete code
    double* res = getAllInfoNEW(environment.data_numeric, environment.levels, X,
        Y, ui_list, TempVector<int>{Z}, environment.n_eff, environment.cplx,
        environment.is_k23, environment.sample_weights, jointFreqs,
        environment.test_mar, environment.cache.cterm);

    output_score = res[6];
    delete[] res;
  } else {
    double* res = compute_Rscore_Ixyz_alg5(data_numeric_red,
        data_numeric_idx_red, all_levels_red, cnt_red, posArray_red, n_ui,
        n_ui + 2, n_samples_non_na_z, sample_weights_red, flag_sample_weights,
        environment);
    output_score = res[0];
    delete[] res;
  }

  saveScore(X, Y, ui_list, Z, output_score, environment);
  return output_score;
}

}  // anonymous namespace

// Computes the two point information X;Y|Ui and the three point information
// X;Y;Z|Ui
double get3PointInfo(
    Environment& environment, int X, int Y, int Z, const vector<int>& ui_list) {
  TempAllocatorScope scope;

  double* res_new = new double[3];
  res_new[0] = -1;

  lookupScore(X, Y, ui_list, Z, res_new, environment);
  if (res_new[0] != -1) {
    double I3 = res_new[1] - res_new[2];
    delete[] res_new;
    return I3;
  }

  // Mark rows containing NAs and count the number of complete samples
  // 1: sample contains NA, 0: sample contains no NA
  TempVector<int> sample_is_not_NA(environment.n_samples);
  // vector with the number of rows containing NAs seen at rank i
  TempVector<int> NAs_count(environment.n_samples);
  int n_samples_non_na_z =
      count_non_NAs(X, Y, ui_list, sample_is_not_NA, NAs_count, environment, Z);

  if (n_samples_non_na_z <= 2) {  // not sufficient statistics
    res_new[0] = n_samples_non_na_z;
    res_new[1] = 0;  // Ixyz
    res_new[2] = 0;  // cplx Ixyz
    saveScore(X, Y, ui_list, Z, res_new, environment);
    delete[] res_new;
    return 0;
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

  bool flag_sample_weights = filter_NAs(X, Y, ui_list, all_levels_red, cnt_red,
      posArray_red, data_numeric_red, data_numeric_idx_red, sample_weights_red,
      sample_is_not_NA, NAs_count, environment, Z);

  const auto& is_continuous = environment.is_continuous;
  if (!is_continuous[X] && !is_continuous[Z] && -!is_continuous[Y]) {
    TempGrid2d<double> jointFreqs = getJointFreqs(environment, X, Y);
    double* res = getAllInfoNEW(environment.data_numeric, environment.levels, X,
        Y, ui_list, TempVector<int>{Z}, environment.n_eff, environment.cplx,
        environment.is_k23, environment.sample_weights, jointFreqs,
        environment.test_mar, environment.cache.cterm);
    res_new[0] = res[0];
    res_new[1] = res[7];
    res_new[2] = -res[8];
    double I3 = res_new[1] - res_new[2];
    delete[] res_new;
    delete[] res;
    return I3;
  }

  // If X or Y has only 1 level
  bool ok = all_levels_red[0] > 1 && all_levels_red[1] > 1;

  //int n_samples_nonNA = getNumSamplesNonNA(environment, X, Y);
  //if (ok && environment.test_mar && n_samples_nonNA != n_samples_non_na_z) {
  //  double kldiv = compute_kl_divergence(X, Y, environment, n_samples_non_na_z,
  //      all_levels_red, sample_is_not_NA);
  //  double cplxMdl = environment.cache.cterm->getLog(n_samples_non_na_z);

  //  if ((kldiv - cplxMdl) > 0) {
  //    // the sample is not representative of the population, hence we do
  //    // not want this z as possible z
  //    ok = false;
  //  }
  //}
  if (!ok) {
    res_new[0] = n_samples_non_na_z;
    res_new[1] = 0;  // Ixyz
    res_new[2] = 0;  // cplx Ixyz
    saveScore(X, Y, ui_list, Z, res_new, environment);
    delete[] res_new;
    return 0;
  }
  //else {
  double* res = compute_Rscore_Ixyz_alg5(data_numeric_red, data_numeric_idx_red,
      all_levels_red, cnt_red, posArray_red, n_ui, n_ui + 2, n_samples_non_na_z,
      sample_weights_red, flag_sample_weights, environment);
  res_new[0] = n_samples_non_na_z;
  res_new[1] = res[1];  // I(x;y;z|u)
  res_new[2] = res[2];  // cplx I(x;y;z|u)
  //}
  delete[] res;

  for (int i = 0; i < 3; i++) {
    if (res_new[i] > -0.0000000001 && res_new[i] < 0.0000000001) {
      res_new[i] = 0.0;
    }
  }

  saveScore(X, Y, ui_list, Z, res_new, environment);
  double I3 = res_new[1] - res_new[2];
  delete[] res_new;
  return I3;
}

double* computeEnsInformationContinuous(Environment& environment, int X, int Y,
    const vector<int>& ui_list, const vector<int>& zi_list, int cplx) {
  TempAllocatorScope scope;

  double* res_new;
  int nbrRetValues = 3;
  int n_ui = ui_list.size();
  int n_zi = zi_list.size();

  // initialization part (no z)
  if (n_zi == 0) {
    // TODO : speedup by only removing NAs for marked columns
    // Mark rows containing NAs and count the number of complete samples
    TempVector<int> sample_is_not_NA(environment.n_samples);
    TempVector<int> NAs_count(environment.n_samples);
    int samplesNotNA =
        count_non_NAs(X, Y, ui_list, sample_is_not_NA, NAs_count, environment);

    if (samplesNotNA <= 2) {
      res_new = new double[3];
      res_new[0] = (double)samplesNotNA;  // N
      res_new[1] = 0;                     // Ixyu
      res_new[2] = 0;                     // cplx
      return res_new;
    } else {
      TempAllocatorScope scope;
      // Allocate data reducted *_red without rows containing NAs
      // All *_red variables are passed to the optimization routine
      TempVector<int> AllLevels_red(n_ui + 2);
      TempVector<int> cnt_red(n_ui + 2);
      TempVector<int> posArray_red(n_ui + 2);
      TempVector<double> sample_weights_red(samplesNotNA);
      TempGrid2d<int> dataNumeric_red(n_ui + 2, samplesNotNA);
      TempGrid2d<int> dataNumericIdx_red(n_ui + 2, samplesNotNA);

      bool flag_sample_weights = filter_NAs(X, Y, ui_list, AllLevels_red,
          cnt_red, posArray_red, dataNumeric_red, dataNumericIdx_red,
          sample_weights_red, sample_is_not_NA, NAs_count, environment);

      // If X or Y has only 1 level
      if (AllLevels_red[0] == 1 || AllLevels_red[1] == 1) {
        res_new = new double[3];
        res_new[0] = (double)samplesNotNA;
        res_new[1] = 0;  // Ixyz
        res_new[2] = 0;  // cplx Ixyz
        return res_new;
      } else {
        res_new = compute_mi_cond_alg1(dataNumeric_red, dataNumericIdx_red,
            AllLevels_red, cnt_red, posArray_red, n_ui, samplesNotNA,
            sample_weights_red, flag_sample_weights, environment);

        res_new[1] = res_new[1] * res_new[0];  // Ixy|u
        res_new[2] = res_new[2] * res_new[0];  // cplx
      }
    }
  } else {  // if nbrZi>0 : iteration part
    res_new = new double[3];
    res_new[0] = environment.n_samples;
    res_new[1] = -1;
    res_new[2] = std::numeric_limits<double>::lowest();

    int n_samples_nonNA = getNumSamplesNonNA(environment, X, Y);
#ifdef _OPENMP
    bool parallelizable =
        environment.first_iter_done && n_zi > environment.n_threads;
#pragma omp parallel for if (parallelizable)
#endif
    for (int iz = 0; iz < n_zi; iz++) {
      int Z = zi_list[iz];
      double score = computeContributingScores(
          environment, X, Y, Z, ui_list, n_samples_nonNA);
#ifdef _OPENMP
#pragma omp critical
#endif
      if (score > res_new[2]) {
        res_new[2] = score;
        res_new[1] = iz;
      }
    }
  }

  for (int i = 0; i < nbrRetValues; i++) {
    if (res_new[i] > -0.0000000001 && res_new[i] < 0.0000000001) {
      res_new[i] = 0.0;
    }
  }

  return res_new;
}

double* computeEnsInformationNew(Environment& environment, int X, int Y,
    const vector<int>& ui_list, const TempVector<int>& zi_list, int cplx) {
  TempAllocatorScope scope;

  TempGrid2d<double> jointFreqs = getJointFreqs(environment, X, Y);
  double* res_new = getAllInfoNEW(environment.data_numeric, environment.levels,
      X, Y, ui_list, zi_list, environment.n_eff, cplx, environment.is_k23,
      environment.sample_weights, jointFreqs, environment.test_mar,
      environment.cache.cterm);

  int nbrRetValues = zi_list.empty()? 3 : 9;
  // If !zi_list.empty(), return {nSample[z1]*I(..|{ui})[z1], NML(..|{ui})[z1],
  // nSample[z1],nSample[z2]*I(..|{ui})[z2], NML(..|{ui})[z2], nSample[z2], ...}
  for (int i = 0; i < nbrRetValues; i++) {
    if (res_new[i] > -0.0000000001 && res_new[i] < 0.0000000001) {
      res_new[i] = 0.0;
    }
  }

  return res_new;
}

void SearchForNewContributingNodeAndItsRank(
    Environment& environment, int X, int Y) {
  auto info = environment.edges[X][Y].shared_info;
  if (!environment.latent) {
    // remove zi that is not connected to neither x nor y
    info->zi_list.erase(
        std::remove_if(info->zi_list.begin(), info->zi_list.end(),
            [&environment, &X, &Y](int z) {
              return !environment.edges[X][z].status &&
                     !environment.edges[Y][z].status;
            }),
        info->zi_list.end());
  }
  if (info->zi_list.empty()) return;

  int cplx = -1;
  if (environment.is_k23 == true) cplx = environment.cplx;

  double* vect = computeEnsInformationContinuous(
      environment, X, Y, info->ui_list, info->zi_list, cplx);
  if (vect[2] > info->Rxyz_ui) {
    info->z_name_idx = vect[1];
    info->Rxyz_ui = vect[2];
  }
  delete[] vect;
}

}  // namespace computation
}  // namespace miic
