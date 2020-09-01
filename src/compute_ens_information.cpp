#include "compute_ens_information.h"

#include <float.h>
#include <math.h>

#include <set>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

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
  if (output_score != -1) {
    return output_score;
  }

  // Mark rows containing NAs and count the number of complete samples
  TempVector<int> sample_is_not_NA(environment.n_samples);
  TempVector<int> NAs_count(environment.n_samples);
  int samplesNotNA =
      count_non_NAs(X, Y, ui_list, sample_is_not_NA, NAs_count, environment, Z);

  if (samplesNotNA <= 2) {
    output_score = -DBL_MAX;
  } else {
    int n_ui = ui_list.size();
    TempAllocatorScope scope;
    // Allocate data reducted *_red without rows containing NAs
    // All *_red variables are passed to the optimization routine
    TempVector<double> sample_weights_red(samplesNotNA);
    TempGrid2d<int> dataNumericIdx_red(n_ui + 3, samplesNotNA);
    TempGrid2d<int> dataNumeric_red(n_ui + 3, samplesNotNA);
    TempVector<int> AllLevels_red(n_ui + 3);
    TempVector<int> cnt_red(n_ui + 3);
    TempVector<int> posArray_red(n_ui + 3);

    bool flag_sample_weights = filter_NAs(X, Y, ui_list, AllLevels_red, cnt_red,
        posArray_red, dataNumeric_red, dataNumericIdx_red, sample_weights_red,
        sample_is_not_NA, NAs_count, environment, Z);

    if (std::all_of(
            cnt_red.begin(), cnt_red.end(), [](int x) { return x == 0; })) {
      TempAllocatorScope scope;

      TempGrid2d<double> jointFreqs =
          getJointFreqs(environment, X, Y, sample_is_not_NA);

      // call discrete code
      double* res = getAllInfoNEW(environment.data_numeric, environment.levels,
          X, Y, ui_list, TempVector<int>{Z}, environment.n_eff,
          environment.cplx, environment.is_k23, environment.sample_weights,
          jointFreqs, environment.test_mar, environment.cache.cterm);

      output_score = res[6];
      delete[] res;
    } else {
      // we do not want to add a z if x or y have only one bin
      bool ok = true;  // ok : do we compute I or return 0?
      if (samplesNotNA < environment.n_samples) {
        std::set<int> s;
        for (int i = 0; i < 2 && ok; i++) {
          s.clear();
          for (int j = 0; j < samplesNotNA; j++) {
            s.insert(dataNumeric_red(i, j));
          }

          if (s.size() == 1) {
            ok = false;
            break;
          }
        }

        if (environment.test_mar && n_samples_nonNA != samplesNotNA) {
          double kldiv = compute_kl_divergence(
              X, Y, environment, samplesNotNA, AllLevels_red, sample_is_not_NA);
          double cplxMdl = environment.cache.cterm->getLog(samplesNotNA);

          if ((kldiv - cplxMdl) > 0) {
            // the sample is not representative of the population, hence we do
            // not want this z as possible z
            ok = false;
          }
        }
      }

      if (ok) {
        double* res = compute_Rscore_Ixyz_alg5(dataNumeric_red,
            dataNumericIdx_red, AllLevels_red, cnt_red, posArray_red, n_ui,
            n_ui + 2, samplesNotNA, sample_weights_red, flag_sample_weights,
            environment);
        output_score = res[0];
        delete[] res;
      } else {
        output_score = -DBL_MAX;
      }
    }
  }  // jump cond no statistics

  saveScore(X, Y, ui_list, Z, output_score, environment);
  return output_score;
}

}  // anonymous namespace

// Computes the two point information X;Y|Ui and the three point information
// X;Y;Z|Ui
double* computeEnsInformationContinuous_Orientation(Environment& environment,
    int X, int Y, const vector<int>& ui_list, int Z, const int cplx) {
  TempAllocatorScope scope;

  int nbrRetValues = 3;
  int n_ui = ui_list.size();

  double* res_new;
  res_new = new double[3];
  res_new[0] = -1;

  lookupScore(X, Y, ui_list, Z, res_new, environment);
  if (res_new[0] != -1) return res_new;

  // Mark rows containing NAs and count the number of complete samples
  // 1: sample contains NA, 0: sample contains no NA
  TempVector<int> sample_is_not_NA(environment.n_samples);
  // vector with the number of rows containing NAs seen at rank i
  TempVector<int> NAs_count(environment.n_samples);
  int samplesNotNA =
      count_non_NAs(X, Y, ui_list, sample_is_not_NA, NAs_count, environment, Z);

  if (samplesNotNA <= 2) {  // not sufficient statistics
    res_new[0] = samplesNotNA;
    res_new[1] = 0;  // Ixyz
    res_new[2] = 0;  // cplx Ixyz
  } else {
    TempAllocatorScope scope;
    // Allocate data reducted *_red without rows containing NAs
    // All *_red variables are passed to the optimization routine
    TempVector<double> sample_weights_red(samplesNotNA);
    TempGrid2d<int> dataNumericIdx_red(n_ui + 3, samplesNotNA);
    TempGrid2d<int> dataNumeric_red(n_ui + 3, samplesNotNA);
    TempVector<int> AllLevels_red(n_ui + 3);
    TempVector<int> cnt_red(n_ui + 3);
    TempVector<int> posArray_red(n_ui + 3);

    bool flag_sample_weights = filter_NAs(X, Y, ui_list, AllLevels_red, cnt_red,
        posArray_red, dataNumeric_red, dataNumericIdx_red, sample_weights_red,
        sample_is_not_NA, NAs_count, environment, Z);

    double* res;

    // If X or Y has only 1 level
    if (AllLevels_red[0] == 1 || AllLevels_red[1] == 1) {
      res_new[0] = (double)samplesNotNA;
      res_new[1] = 0;  // Ixyz
      res_new[2] = 0;  // cplx Ixyz

    } else {
      res = compute_Rscore_Ixyz_alg5(dataNumeric_red, dataNumericIdx_red,
          AllLevels_red, cnt_red, posArray_red, n_ui, n_ui + 2, samplesNotNA,
          sample_weights_red, flag_sample_weights, environment);

      res_new[0] = (double)samplesNotNA;
      res_new[1] = res[1];  // I(x;y;z|u)
      res_new[2] = res[2];  // cplx I(x;y;z|u)

      delete[] res;
    }

  }  // jump cond no sufficient statistics

  for (int i = 0; i < nbrRetValues; i++) {
    if (res_new[i] > -0.0000000001 && res_new[i] < 0.0000000001) {
      res_new[i] = 0.0;
    }
  }

  saveScore(X, Y, ui_list, Z, res_new, environment);
  return res_new;
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
    res_new[2] = -DBL_MAX;
    double* res;

    TempVector<int> zi_indices;

    // If x, y and uis are discrete we can put togheter all zi that are discrete
    // in a vector and evaluate them in one shot as discrete.
    if (!environment.is_continuous[X] && !environment.is_continuous[Y] &&
        std::all_of(begin(ui_list), end(ui_list),
            [&environment](int i) { return !environment.is_continuous[i]; })) {
      TempVector<int> zi_discrete;
      zi_discrete.reserve(n_zi);
      zi_indices.reserve(n_zi);
      for (int i = 0; i < n_zi; ++i) {
        if (environment.is_continuous[zi_list[i]])
          zi_indices.push_back(i);
        else
          zi_discrete.push_back(zi_list[i]);
      }
      zi_discrete.shrink_to_fit();
      zi_indices.shrink_to_fit();
      if (!zi_discrete.empty()) {
        TempAllocatorScope scope;

        TempVector<int> posZi(zi_discrete.size());
        int pos = 0;
        for (int i = 0; i < n_zi; i++) {
          if (!environment.is_continuous[zi_list[i]]) posZi[pos++] = i;
        }
        TempGrid2d<double> jointFreqs = getJointFreqs(environment, X, Y);
        res = getAllInfoNEW(environment.data_numeric, environment.levels, X, Y,
            ui_list, zi_discrete, environment.n_eff, cplx, environment.is_k23,
            environment.sample_weights, jointFreqs, environment.test_mar,
            environment.cache.cterm);

        // update res new, it will be compared to the continuous variables
        res_new[2] = res[6];
        if (res[3] >= 0) {
          res_new[1] = posZi[int(res[3])];
        }
        res_new[0] = res[0];
        delete[] res;
      }
      n_zi -= zi_discrete.size();
    }

    int n_samples_nonNA = getNumSamplesNonNA(environment, X, Y);
#ifdef _OPENMP
    bool parallelizable =
        environment.first_iter_done && n_zi > environment.n_threads;
#pragma omp parallel for if (parallelizable)
#endif
    for (int iz = 0; iz < n_zi; iz++) {
      int Z = zi_list[iz];
      if (!zi_indices.empty()) Z = zi_list[zi_indices[iz]];
      double score = computeContributingScores(
          environment, X, Y, Z, ui_list, n_samples_nonNA);
      if (score > res_new[2]) {
        res_new[2] = score;
        if (zi_indices.empty()) {
          res_new[1] = iz;
        } else {
          res_new[1] = zi_indices[iz];
        }
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
    Environment& environment, const int posX, const int posY) {
  auto info = environment.edges[posX][posY].shared_info;
  if (!environment.latent) {
    // remove zi that is not connected to neither x nor y
    info->zi_list.erase(
        std::remove_if(info->zi_list.begin(), info->zi_list.end(),
            [&environment, &posX, &posY](int z) {
              return !environment.edges[posX][z].status &&
                     !environment.edges[posY][z].status;
            }),
        info->zi_list.end());
  }
  if (info->zi_list.empty()) return;

  int argEnsInfo = -1;
  if (environment.is_k23 == true) argEnsInfo = environment.cplx;

  double* vect = NULL;

  if (std::all_of(environment.is_continuous.begin(),
          environment.is_continuous.end(), [](int i) { return i == 0; })) {
    TempAllocatorScope scope;
    TempVector<int> zi_list(begin(info->zi_list), end(info->zi_list));
    vect = computeEnsInformationNew(
        environment, posX, posY, info->ui_list, zi_list, argEnsInfo);
    if (vect[6] - info->Rxyz_ui > 0) {
      info->z_name_idx = vect[3];
      info->Rxyz_ui = vect[6];
    }
  } else {
    vect = computeEnsInformationContinuous(
        environment, posX, posY, info->ui_list, info->zi_list, argEnsInfo);
    if (vect[2] - info->Rxyz_ui > 0) {
      info->z_name_idx = vect[1];
      info->Rxyz_ui = vect[2];
    }
  }
  delete[] vect;
}

}  // namespace computation
}  // namespace miic
