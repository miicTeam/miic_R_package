#include "compute_ens_information.h"

#include <float.h>
#include <math.h>

#include <iostream>
#include <set>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "compute_info.h"
#include "info_cnt.h"
#include "structure.h"
#include "utilities.h"

// for memory space on continuous data
#define MAX_NBRUI 10
namespace miic {
namespace computation {

using namespace miic::structure;
using namespace miic::utility;
using std::cout;
using std::endl;
using std::vector;

// Computes the two point information X;Y|Ui and the three point information
// X;Y;Z|Ui
double* computeEnsInformationContinuous_Orientation(Environment& environment,
    int* myCond, int myNbrUi, int* myZi, const int myVarIdxX,
    const int myVarIdxY, const int cplx, MemorySpace& m) {
  vector<int> posArray(
      2 + environment.edges[myVarIdxX][myVarIdxY].shared_info->ui_list.size());
  posArray[0] = myVarIdxX;
  posArray[1] = myVarIdxY;

  int nbrRetValues = 3;

  if (myCond != NULL) {
    if (myNbrUi > 0) {
      // The index of the variables in this dataset
      for (int i = 0; i < myNbrUi; ++i) posArray[i + 2] = myCond[i];
    }
  }

  double* res_new;
  res_new = new double[3];
  res_new[0] = -1;
  int z = myZi[0];

  lookupScore(posArray, myNbrUi, z, res_new, environment);
  if(res_new[0] != -1) {
    return res_new;
  }

  // Mark rows containing NAs and count the number of complete samples
  //vector with zero or one according if the sample at position i contains NA or not
  vector <int> sample_is_not_NA(environment.n_samples);
  //vector with the number of rows containing NAs seen at rank i
  vector <int> NAs_count(environment.n_samples);
  int samplesNotNA = count_non_NAs(myNbrUi, sample_is_not_NA,
    NAs_count, posArray, environment, z);

  if (samplesNotNA <= 2) {  // not sufficient statistics
    res_new[0] = samplesNotNA;
    res_new[1] = 0;  // Ixyz
    res_new[2] = 1;  // cplx Ixyz
  } else {

    // Allocate data reducted *_red without rows containing NAs
    // All *_red variables are passed to the optimization routine
    vector<double> sample_weights_red(samplesNotNA);
    vector<vector<int> > dataNumericIdx_red(myNbrUi+3, vector<int>(samplesNotNA));
    vector<vector<int> > dataNumeric_red(myNbrUi+3, vector<int>(samplesNotNA));
    vector<int> AllLevels_red(myNbrUi+3);
    vector<int> cnt_red(myNbrUi+3);
    vector<int> posArray_red(myNbrUi+3);

    bool flag_sample_weights = filter_NAs(myNbrUi, AllLevels_red, cnt_red,
      posArray_red, posArray, dataNumeric_red, dataNumericIdx_red,
      sample_weights_red, sample_is_not_NA, NAs_count,
      environment, z);

    double* res;

    // If X or Y has only 1 level
    if (AllLevels_red[0] == 1 || AllLevels_red[1] == 1){
      res_new[0] = (double)samplesNotNA;
      res_new[1] = 0;  // Ixyz
      res_new[2] = 1;  // cplx Ixyz

    } else {
      res = compute_Rscore_Ixyz_alg5(dataNumeric_red, dataNumericIdx_red,
          AllLevels_red, cnt_red, posArray_red, myNbrUi, myNbrUi + 2,
          samplesNotNA, sample_weights_red, flag_sample_weights, environment);

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

#if _MY_DEBUG_NEW
  printf("\n# =====> after getAllInfoNEW \n");
  if (myNbrZi == 0) {
    printf("# N=res_new[%d]=%g ", 0, res_new[0]);
    for (int i = 1; i < 3; i++) printf("# res_new[%d]=%g ", i, res_new[i]);
  }
  if (myNbrZi > 0) {
    printf("# N=res_new[%d]=%g ", 0, res_new[0]);
    for (int i = 1; i < 3; i++) printf("# res_new[%d]=%g ", i, res_new[i]);
    printf("# z=res_new[%d]=%g ", 3, res_new[3]);
    for (int i = 4; i < 9; i++) printf(" res_new[%d]=%g \n", i, res_new[i]);
  }
  printf("\n");
#endif  // _MY_DEBUG_NEW

  saveScore(posArray, myNbrUi, z, res_new, environment);
  return res_new;
}

void computeContributingScores(Environment& environment, int* ziContPosIdx,
    int iz, int* myZi, int myNbrUi, int n_samples_nonNA,
    const vector<int>& posArray, double* scoresZ, MemorySpace m) {
  // progressive data rank with repetition for same values

  int cplx = environment.cplx;
  int z;
  if (ziContPosIdx == NULL)
    z = myZi[iz];
  else
    z = myZi[ziContPosIdx[iz]];

  double output_score = lookupScore(posArray, myNbrUi, z, environment);
  if(output_score != -1) {
    scoresZ[iz] = output_score;
    return;
  }

  // Mark rows containing NAs and count the number of complete samples
  vector <int> sample_is_not_NA(environment.n_samples);
  vector <int> NAs_count(environment.n_samples);
  int samplesNotNA = count_non_NAs(myNbrUi, sample_is_not_NA,
    NAs_count, posArray, environment, z);

  if (samplesNotNA <= 2) {
    output_score = -DBL_MAX;
  } else {
    // Allocate data reducted *_red without rows containing NAs
    // All *_red variables are passed to the optimization routine
    vector<double> sample_weights_red(samplesNotNA);
    vector<vector<int> > dataNumericIdx_red(myNbrUi+3, vector<int>(samplesNotNA));
    vector<vector<int> > dataNumeric_red(myNbrUi+3, vector<int>(samplesNotNA));
    vector<int> AllLevels_red(myNbrUi+3);
    vector<int> cnt_red(myNbrUi+3);
    vector<int> posArray_red(myNbrUi+3);

    bool flag_sample_weights = filter_NAs(myNbrUi, AllLevels_red, cnt_red,
      posArray_red, posArray, dataNumeric_red, dataNumericIdx_red,
      sample_weights_red, sample_is_not_NA, NAs_count,
      environment, z);

    if (std::all_of(cnt_red.begin(), cnt_red.end(), [](int x){ return x == 0; })) {
      // call discrete code
      int* zz = new int[1];
      zz[0] = z;

      double** jointFreqs = getJointFreqs(
          environment, posArray[0], posArray[1], sample_is_not_NA);

      double* res = getAllInfoNEW(environment.oneLineMatrix,
          environment.levels, posArray, myNbrUi, zz, 1, -1,
          environment.n_samples, environment.n_eff, cplx, environment.is_k23,
          environment.looklog, environment.c2terms, &m,
          environment.sample_weights, jointFreqs, environment.test_mar);

      output_score = res[6];
      delete[] res;

      for (int level0 = 0; level0 < environment.levels[posArray[0]];
           level0++)
        delete[] jointFreqs[level0];
      delete[] jointFreqs;
      delete[] zz;

    } else {
      // we do not want to add a z if x or y have only one bin
      bool ok = true;  // ok : do we compute I or return 0?
      if (samplesNotNA < environment.n_samples) {
        std::set<int> s;
        for (int i = 0; i < 2 && ok; i++) {
          s.clear();
          for (int j = 0; j < samplesNotNA; j++) {
            s.insert(dataNumeric_red[i][j]);
          }

          if (s.size() == 1) {
            ok = false;
            break;
          }
        }

        if (environment.test_mar && n_samples_nonNA != samplesNotNA) {
          double kldiv = compute_kl_divergence(posArray, environment,
              samplesNotNA, AllLevels_red, sample_is_not_NA);
          double cplxMdl = log(samplesNotNA);

          if ((kldiv - cplxMdl) > 0) {
            // the sample is not representative of the population, hence we do
            // not want this z as possible z
            ok = false;
          }
        }
      }

      if (ok) {
        double* res = compute_Rscore_Ixyz_alg5(dataNumeric_red,
            dataNumericIdx_red, AllLevels_red, cnt_red, posArray_red, myNbrUi,
            myNbrUi + 2, samplesNotNA, sample_weights_red, flag_sample_weights,
            environment);
        output_score = res[0];
        delete[] res;
      } else {
        output_score = -DBL_MAX;
      }
    }
  }  // jump cond no statistics

  saveScore(posArray, myNbrUi, z, output_score, environment);
  scoresZ[iz] = output_score;
}

double* computeEnsInformationContinuous(Environment& environment, int* myCond,
    int myNbrUi, int* myZi, int myNbrZi, int myZiPos, const int myVarIdxX,
    const int myVarIdxY, const int cplx, MemorySpace& m) {
  vector<int> posArray(
      2 + environment.edges[myVarIdxX][myVarIdxY].shared_info->ui_list.size());
  posArray[0] = myVarIdxX;
  posArray[1] = myVarIdxY;

  int nbrRetValues = 3;

  if (myCond != NULL) {
    if (myNbrUi > 0) {
      // The index of the variables in this dataset
      for (int i = 0; i < myNbrUi; ++i) posArray[i + 2] = myCond[i];
    }
  }

  double* res_new;

  // initialization part (no z)
  if (myNbrZi == 0) {

    // TODO : speedup by only removing NAs for marked columns
    // Mark rows containing NAs and count the number of complete samples
    vector <int> sample_is_not_NA(environment.n_samples);
    vector <int> NAs_count(environment.n_samples);
    int samplesNotNA = count_non_NAs(myNbrUi, sample_is_not_NA,
      NAs_count, posArray, environment);

    if (samplesNotNA <= 2) {
      res_new = new double[3];
      res_new[0] = (double)samplesNotNA;  // N
      res_new[1] = 0;             // Ixyu
      res_new[2] = 1;             // cplx
    } else {

      // Allocate data reducted *_red without rows containing NAs
      // All *_red variables are passed to the optimization routine
      vector<int> AllLevels_red(myNbrUi+2);
      vector<int> cnt_red(myNbrUi+2);
      vector<int> posArray_red(myNbrUi+2);
      vector<double> sample_weights_red(samplesNotNA);
      vector<vector<int> > dataNumeric_red((myNbrUi+2), vector<int> (samplesNotNA));
      vector<vector<int> > dataNumericIdx_red((myNbrUi+2), vector<int> (samplesNotNA));

      bool flag_sample_weights = filter_NAs(myNbrUi, AllLevels_red, cnt_red,
        posArray_red, posArray, dataNumeric_red, dataNumericIdx_red,
        sample_weights_red, sample_is_not_NA, NAs_count,
        environment);

      // If X or Y has only 1 level
      if (AllLevels_red[0] == 1 || AllLevels_red[1] == 1){
        res_new = new double[3];
        res_new[0] = (double)samplesNotNA;
        res_new[1] = 0;  // Ixyz
        res_new[2] = 1;  // cplx Ixyz
      } else{
        res_new = compute_mi_cond_alg1(dataNumeric_red, dataNumericIdx_red,
          AllLevels_red, cnt_red, posArray_red, myNbrUi, samplesNotNA,
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

    int z;
    int* ziContPosIdx = NULL;

    // If x, y and uis are discrete we can put togheter all zi that are discrete
    // in a vector and evaluate them in one shot as discrete.
    if (std::all_of(posArray.cbegin(), posArray.cend(),
            [&environment](int i) { return !environment.is_continuous[i]; })) {
      // search for z that are discrete
      int countZDiscrete = 0;
      for (int iz = 0; iz < myNbrZi; iz++) {
        z = myZi[iz];
        if (!environment.is_continuous[z]) countZDiscrete++;
      }

      if (countZDiscrete > 0) {
        int* posZi = new int[countZDiscrete];
        int* zz = new int[countZDiscrete];
        int pos = 0;
        for (int iz = 0; iz < myNbrZi; iz++) {
          z = myZi[iz];
          if (!environment.is_continuous[z]) {
            zz[pos] = z;
            posZi[pos] = iz;
            pos++;
          }
        }
        double** jointFreqs = getJointFreqs(
            environment, posArray[0], posArray[1]);
        res = getAllInfoNEW(environment.oneLineMatrix, environment.levels,
            posArray, myNbrUi, zz, countZDiscrete, -1, environment.n_samples,
            environment.n_eff, cplx, environment.is_k23, environment.looklog,
            environment.c2terms, &m, environment.sample_weights, jointFreqs,
            environment.test_mar);

        for (int level0 = 0; level0 < environment.levels[posArray[0]];
             level0++)
          delete[] jointFreqs[level0];
        delete[] jointFreqs;
        delete[] zz;

        // keep in ziContPos only the position of the continuous variables
        ziContPosIdx = new int[myNbrZi - countZDiscrete];
        pos = 0;
        for (int iz = 0; iz < myNbrZi; iz++) {
          z = myZi[iz];
          if (environment.is_continuous[z]) {
            ziContPosIdx[pos] = iz;
            pos++;
          }
        }
        myNbrZi = pos;

        // update res new, it will be compared to the continuous variables
        res_new[2] = res[6];
        res_new[1] = posZi[int(res[3])];
        res_new[0] = res[0];
        free(res);
        delete[] posZi;
      }
    }

    double* scoresZ = new double[myNbrZi];
#ifdef _OPENMP
    bool parallelizable =
        environment.first_iter_done && myNbrZi > environment.n_threads;
#pragma omp parallel for if (parallelizable)
#endif
    for (int iz = 0; iz < myNbrZi; iz++) {
      MemorySpace privateM = m;
#ifdef _OPENMP
      if (parallelizable)
        privateM = environment.memoryThreads[omp_get_thread_num()];
#endif
      int n_samples_nonNA =
          getNumSamplesNonNA(environment, posArray[0], posArray[1]);

      computeContributingScores(environment, ziContPosIdx, iz, myZi, myNbrUi,
          n_samples_nonNA, posArray, scoresZ, privateM);
    }  // parallel for on z

    for (int iz = 0; iz < myNbrZi; iz++) {  // find optimal z
      if (scoresZ[iz] > res_new[2]) {
        res_new[2] = scoresZ[iz];
        if (ziContPosIdx == NULL) {
          res_new[1] = iz;
        } else {
          res_new[1] = ziContPosIdx[iz];
        }
        // TO DEFINE
        // res_new[0]=(double) samplesNotNA;
      }
    }  // optimal z search

    if (ziContPosIdx != NULL) delete[] ziContPosIdx;
  }

  for (int i = 0; i < nbrRetValues; i++) {
    if (res_new[i] > -0.0000000001 && res_new[i] < 0.0000000001) {
      res_new[i] = 0.0;
    }
  }

#if _MY_DEBUG_NEW
  printf("\n# =====> after getAllInfoNEW \n");
  if (myNbrZi == 0) {
    printf("# N=res_new[%d]=%g ", 0, res_new[0]);
    for (int i = 1; i < 3; i++) printf("# res_new[%d]=%g ", i, res_new[i]);
  }
  if (myNbrZi > 0) {
    printf("# N=res_new[%d]=%g ", 0, res_new[0]);
    for (int i = 1; i < 3; i++) printf("# res_new[%d]=%g ", i, res_new[i]);
    printf("# z=res_new[%d]=%g ", 3, res_new[3]);
    for (int i = 4; i < 9; i++) printf(" res_new[%d]=%g \n", i, res_new[i]);
  }
  printf("\n");
#endif  // _MY_DEBUG_NEW

  return res_new;
}

double* computeEnsInformationNew(Environment& environment, int* myCond,
    int myNbrUi, int* myZi, int myNbrZi, int myZiPos, const int myVarIdxX,
    const int myVarIdxY, const int cplx, MemorySpace& m) {
  vector<int> posArray(
      2 + environment.edges[myVarIdxX][myVarIdxY].shared_info->ui_list.size());
  posArray[0] = myVarIdxX;
  posArray[1] = myVarIdxY;

  if (myCond != NULL) {
    if (myNbrUi > 0) {
      // The index of the variables in this dataset
      for (int i = 0; i < myNbrUi; ++i) posArray[i + 2] = myCond[i];
    }
  }

// Compute the mutual information
#if _MY_DEBUG_NEW
  printf("\n# =====> before getAllInfoNEW \n");
#endif  // _MY_DEBUG_NEW

  double** jointFreqs =
      getJointFreqs(environment, posArray[0], posArray[1]);

  double* res_new = getAllInfoNEW(environment.oneLineMatrix,
      environment.levels, posArray, myNbrUi, myZi, myNbrZi, myZiPos,
      environment.n_samples, environment.n_eff, cplx, environment.is_k23,
      environment.looklog, environment.c2terms, &m, environment.sample_weights,
      jointFreqs, environment.test_mar);

  for (int level0 = 0; level0 < environment.levels[posArray[0]]; level0++)
    delete[] jointFreqs[level0];
  delete[] jointFreqs;

  int nbrRetValues = 3;

  // If nbrZi > 0, return {nSample[z1]*I(..|{ui})[z1], NML(..|{ui})[z1],
  // nSample[z1],nSample[z2]*I(..|{ui})[z2], NML(..|{ui})[z2], nSample[z2], ...}
  if (myNbrZi > 0) nbrRetValues = 9;

  for (int i = 0; i < nbrRetValues; i++) {
    if (res_new[i] > -0.0000000001 && res_new[i] < 0.0000000001) {
      res_new[i] = 0.0;
    }
  }

#if _MY_DEBUG_NEW
  printf("\n# =====> after getAllInfoNEW \n");
  if (myNbrZi == 0) {
    printf("# N=res_new[%d]=%g ", 0, res_new[0]);
    for (int i = 1; i < 3; i++) printf("# res_new[%d]=%g ", i, res_new[i]);
  }
  if (myNbrZi > 0) {
    printf("# N=res_new[%d]=%g ", 0, res_new[0]);
    for (int i = 1; i < 3; i++) printf("# res_new[%d]=%g ", i, res_new[i]);
    printf("# z=res_new[%d]=%g ", 3, res_new[3]);
    for (int i = 4; i < 9; i++) printf(" res_new[%d]=%g \n", i, res_new[i]);
  }
  printf("\n");
#endif  // _MY_DEBUG_NEW

  return res_new;
}

void SearchForNewContributingNodeAndItsRank(
    Environment& environment, const int posX, const int posY, MemorySpace& m) {
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

  int* ui;
  int* zi;

  if (info->ui_list.empty())
    ui = NULL;
  else
    ui = &info->ui_list[0];

  if (info->zi_list.empty())
    zi = NULL;
  else
    zi = &info->zi_list[0];

  int argEnsInfo = -1;
  if (environment.is_k23 == true) argEnsInfo = environment.cplx;

  double* vect = NULL;

  if (std::all_of(environment.is_continuous.begin(),
          environment.is_continuous.end(), [](int i) { return i == 0; })) {
    vect = computeEnsInformationNew(environment, ui, info->ui_list.size(), zi,
        info->zi_list.size(), info->ui_list.size() + 2, posX, posY, argEnsInfo,
        m);
    if (vect[6] - info->Rxyz_ui > 0) {
      if (environment.verbose) {
        cout << "\n"
             << posX << "    " << posY << "# -----> possible zi: "
             << environment
                    .nodes[environment.edges[posX][posY]
                               .shared_info->zi_list[vect[3]]]
                    .name
             << "(" << vect[6] << " > " << info->Rxyz_ui << ")\n";
      }

      // The order matters: set first the z.name.idx, than get the corresponding
      // zi from the original vect / Doing this way, we make sure that the
      // z.name has the right bin xyzi key
      info->z_name_idx = vect[3];
      info->Rxyz_ui = vect[6];
      free(vect);

    } else if (environment.verbose) {
      cout << "# --!!--> Rxyz_ui.tmp = " << vect[6]
           << " < Rxyz_ui = " << info->Rxyz_ui << "\n";
    }

  } else {
    vect = computeEnsInformationContinuous(environment, ui,
        info->ui_list.size(), zi, info->zi_list.size(),
        info->ui_list.size() + 2, posX, posY, argEnsInfo, m);
    if (vect[2] - info->Rxyz_ui > 0) {
      if (environment.verbose) {
        cout << "\n"
             << posX << " " << posY << " # -----> possible zi: "
             << environment.nodes[info->zi_list[vect[1]]].name << "(" << vect[2]
             << " > " << info->Rxyz_ui << ")\n";
      }

      // The order matters: set first the z.name.idx, than get the corresponding
      // zi from the original vect / Doing this way, we make sure
      // that the z.name has the right bin xyzi key
      info->z_name_idx = vect[1];
      info->Rxyz_ui = vect[2];

    } else if (environment.verbose) {
      cout << "# --!!--> Rxyz_ui.tmp = " << vect[2]
           << " < Rxyz_ui = " << info->Rxyz_ui << "\n";
    }
    delete[] vect;
  }
}

}  // namespace computation
}  // namespace miic
