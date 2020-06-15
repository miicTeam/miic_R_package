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
using uint = unsigned int;
using std::cout;
using std::endl;
using std::vector;

// for gaussian variables
double correlations(Environment& environment, const int a, const int b,
    const int* subset, const int l, double** p) {
  int dim = l + 2;
  // initialization of p (looks like rho)
  for (int i = 0; i < dim - 1; i++) {
    for (int j = i + 1; j < dim; j++) {
      int first, second;

      if (i == 0)
        first = a;
      else if (i == 1)
        first = b;
      else
        first = subset[i - 2];

      if (j == 1)
        second = b;
      else
        second = subset[j - 2];

      p[i][j] = p[j][i] = environment.rho[first][second];
    }
  }

  for (int k = 1; k <= l; k++) {
    for (int i = 0; i <= (l - k); i++) {
      for (int j = i + 1; j < (dim - k); j++) {
        p[i][j] = p[j][i] =
            (p[i][j] - p[i][dim - k] * p[j][dim - k]) /
            (sqrt((1 - pow(p[i][dim - k], 2)) * (1 - pow(p[j][dim - k], 2))));
      }
    }
  }

  return p[0][1];
}

double getCorrelationCoefficient(Environment& environment, const int r,
    const int c, const int* subset, const int l, double** p) {
  double correlationCoefficient;

  if (l == 2) {
    double rhoRC, rhoRH1, rhoCH1;
    int h1 = subset[0];
    int h2 = subset[1];

    rhoRC = (environment.rho[r][c] -
                (environment.rho[r][h2] * environment.rho[c][h2])) /
            sqrt((1 - pow(environment.rho[r][h2], 2)) *
                 (1 - pow(environment.rho[c][h2], 2)));
    rhoRH1 = (environment.rho[r][h1] -
                 (environment.rho[r][h2] * environment.rho[h1][h2])) /
             sqrt((1 - pow(environment.rho[r][h2], 2)) *
                  (1 - pow(environment.rho[h1][h2], 2)));
    rhoCH1 = (environment.rho[c][h1] -
                 (environment.rho[c][h2] * environment.rho[h1][h2])) /
             sqrt((1 - pow(environment.rho[c][h2], 2)) *
                  (1 - pow(environment.rho[h1][h2], 2)));

    correlationCoefficient = (rhoRC - (rhoRH1 * rhoCH1)) /
                             sqrt((1 - pow(rhoRH1, 2)) * (1 - pow(rhoCH1, 2)));
  } else if (l == 3) {
    double rhoRC_H2H3, rhoRH1_H2H3, rhoCH1_H2H3, rhoRC_H3, rhoRH1_H3, rhoRH2_H3,
        rhoCH1_H3, rhoCH2_H3, rhoH1H2_H3;
    int h1 = subset[0];
    int h2 = subset[1];
    int h3 = subset[2];

    rhoRC_H3 = (environment.rho[r][c] -
                   (environment.rho[r][h3] * environment.rho[c][h3])) /
               sqrt((1 - pow(environment.rho[r][h3], 2)) *
                    (1 - pow(environment.rho[c][h3], 2)));
    rhoRH1_H3 = (environment.rho[r][h1] -
                    (environment.rho[r][h3] * environment.rho[h1][h3])) /
                sqrt((1 - pow(environment.rho[r][h3], 2)) *
                     (1 - pow(environment.rho[h1][h3], 2)));
    rhoRH2_H3 = (environment.rho[r][h2] -
                    (environment.rho[r][h3] * environment.rho[h2][h3])) /
                sqrt((1 - pow(environment.rho[r][h3], 2)) *
                     (1 - pow(environment.rho[h2][h3], 2)));
    rhoCH1_H3 = (environment.rho[c][h1] -
                    (environment.rho[c][h3] * environment.rho[h1][h3])) /
                sqrt((1 - pow(environment.rho[c][h3], 2)) *
                     (1 - pow(environment.rho[h1][h3], 2)));
    rhoCH2_H3 = (environment.rho[c][h2] -
                    (environment.rho[c][h3] * environment.rho[h2][h3])) /
                sqrt((1 - pow(environment.rho[c][h3], 2)) *
                     (1 - pow(environment.rho[h2][h3], 2)));
    rhoH1H2_H3 = (environment.rho[h1][h2] -
                     (environment.rho[h1][h3] * environment.rho[h2][h3])) /
                 sqrt((1 - pow(environment.rho[h1][h3], 2)) *
                      (1 - pow(environment.rho[h2][h3], 2)));

    rhoRC_H2H3 = (rhoRC_H3 - (rhoRH2_H3 * rhoCH2_H3)) /
                 sqrt((1 - pow(rhoRH2_H3, 2)) * (1 - pow(rhoCH2_H3, 2)));
    rhoRH1_H2H3 = (rhoRH1_H3 - (rhoRH2_H3 * rhoH1H2_H3)) /
                  sqrt((1 - pow(rhoRH2_H3, 2)) * (1 - pow(rhoH1H2_H3, 2)));
    rhoCH1_H2H3 = (rhoCH1_H3 - (rhoCH2_H3 * rhoH1H2_H3)) /
                  sqrt((1 - pow(rhoCH2_H3, 2)) * (1 - pow(rhoH1H2_H3, 2)));

    correlationCoefficient =
        (rhoRC_H2H3 - (rhoRH1_H2H3 * rhoCH1_H2H3)) /
        sqrt((1 - pow(rhoRH1_H2H3, 2)) * (1 - pow(rhoCH1_H2H3, 2)));
  } else {
    correlationCoefficient = correlations(environment, r, c, subset, l, p);
  }

  return correlationCoefficient;
}

// if whatToDo == 0 -> I(x;y)    or I(x;y|ui)
double* corrMutInfo(Environment& environment, double** dataset, int* uis,
    int numUis, int* zis, int numZis, int i, int j, int whatToDo) {
  double* r;
  if (numZis > 0)
    r = new double[numZis];
  else
    r = new double[1];

  if (whatToDo == -2) {
    if (numUis == 0 && numZis == 0) {
      r[0] = environment.rho[i][j];
      // cout << "corr:" << r[0] << endl;
      r[0] = (-log(1 - pow(r[0], 2)) / 2) * environment.numSamples;

    } else if (numZis == 0) {
      r[0] = getCorrelationCoefficient(
          environment, i, j, uis, numUis, environment.pMatrix);
      r[0] = (-log(1 - pow(r[0], 2)) / 2) * environment.numSamples;
    }
  }

  if (whatToDo == -1) {
    // z is not used
    r[0] = getCorrelationCoefficient(
        environment, i, j, uis, numUis, environment.pMatrix);

    r[0] = (-log(1 - pow(r[0], 2)) / 2) * environment.numSamples;
    for (int z = 1; z < numZis; z++) {
      r[z] = r[0];
    }
  } else if (whatToDo == 0) {
    for (int z = 0; z < numZis; z++) {
      r[z] = getCorrelationCoefficient(
          environment, zis[z], j, uis, numUis, environment.pMatrix);
      r[z] = (-log(1 - pow(r[z], 2)) / 2) * environment.numSamples;
    }
  } else if (whatToDo == 1) {
    for (int z = 0; z < numZis; z++) {
      r[z] = getCorrelationCoefficient(
          environment, i, zis[z], uis, numUis, environment.pMatrix);
      r[z] = (-log(1 - pow(r[z], 2)) / 2) * environment.numSamples;
    }
  } else if (whatToDo == 2) {
    int* uis1 = new int[numUis + 1];
    for (int k = 0; k < numUis; k++) {
      uis1[k] = uis[k];
    }

    for (int z = 0; z < numZis; z++) {
      // add z to uis
      uis1[numUis] = zis[z];
      r[z] = getCorrelationCoefficient(
          environment, i, j, uis1, numUis + 1, environment.pMatrix);
      r[z] = (-log(1 - pow(r[z], 2)) / 2) * environment.numSamples;
    }
  }

  return r;
}

// Compute the rank and find the minimum of each z
double* computeRankAndFindMin(double* xz, double* zy, double* xyz, int length) {
  double* res = new double[length];
  // score increasing
  double first;
  double second;
  for (int i = 0; i < length; ++i) {
    if (xz[i] < zy[i]) {
      first = xz[i];
      second = zy[i];
    } else {
      first = zy[i];
      second = xz[i];
    }

    double myPb = first - log1p(exp(-(second - first)));
    res[i] = std::min(xyz[i], myPb);
  }
  return res;
}

// Compute the difference of two arrays every 3 cells (only on information), and
// return the array of differences
double* computeDifferenceByStep(
    double* first, double* second, int position, int length) {
  double* res = new double[length / 3];
  int pos = 0;
  int step = position;
  while (step < length) {
    res[pos] = first[step] - second[step];
    step += 3;
    pos++;
  }
  return res;
}

// Compute the sum of two arrays cell by cell, but mulltiply the first by -1
double* computeMinus1TimeAndSum(double* first, double* second, int length) {
  double* res = new double[length];
  for (int pos = 0; pos < length; pos++)
    res[pos] = (-1 * first[pos]) + second[pos];

  return res;
}

double* computeSum(double* first, double* second, int length) {
  double* res = new double[length];
  for (int pos = 0; pos < length; pos++) {
    res[pos] = first[pos] + second[pos];
  }
  return res;
}

double* computeDifference(double* first, double* second, int length) {
  double* res = new double[length];
  for (int pos = 0; pos < length; pos++) res[pos] = first[pos] - second[pos];

  return res;
}

// Computes the two point information X;Y|Ui and the three point information
// X;Y;Z|Ui
double* computeEnsInformationContinuous_Orientation(Environment& environment,
    int* myCond, int myNbrUi, int* myZi, const int myVarIdxX,
    const int myVarIdxY, const int cplx, MemorySpace& m) {
  int* posArray = new int[2 + environment.edges[myVarIdxX][myVarIdxY]
                                  .shared_info->ui_vect_idx.size()];
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
  int z = myZi[0];

  // Mark rows containing NAs and count the number of complete samples
  //vector with zero or one according if the sample at position i contains NA or not
  vector <int> sample_is_not_NA(environment.numSamples);
  //vector with the number of rows containing NAs seen at rank i
  vector <int> NAs_count(environment.numSamples);
  uint samplesNotNA = count_non_NAs(myNbrUi, sample_is_not_NA,
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

      free(res);
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

  return res_new;
}

void computeContributingScores(Environment& environment, int* ziContPosIdx,
    int iz, int* myZi, int myNbrUi, uint numSamples_nonNA, int* posArray,
    double* res, double* scoresZ, MemorySpace m) {
  // progressive data rank with repetition for same values

  int cplx = environment.cplx;
  int z;
  if (ziContPosIdx == NULL)
    z = myZi[iz];
  else
    z = myZi[ziContPosIdx[iz]];

  // Mark rows containing NAs and count the number of complete samples
  vector <int> sample_is_not_NA(environment.numSamples);
  vector <int> NAs_count(environment.numSamples);
  uint samplesNotNA = count_non_NAs(myNbrUi, sample_is_not_NA,
    NAs_count, posArray, environment, z);

  if (samplesNotNA <= 2) {
    res = new double[3];
    res[0] = -DBL_MAX;  // Rscore
    res[1] = 0;         // I
    res[2] = 1;         // k
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

      res = getAllInfoNEW(environment.oneLineMatrix, environment.allLevels,
          posArray, myNbrUi, zz, 1, -1, environment.numSamples,
          environment.effN, cplx, environment.isK23, environment.looklog,
          environment.c2terms, &m, environment.sampleWeights, jointFreqs,
          environment.testDistribution);

      res[0] = res[6];

      for (uint level0 = 0; level0 < environment.allLevels[posArray[0]];
           level0++)
        delete[] jointFreqs[level0];
      delete[] jointFreqs;
      delete[] zz;

    } else {
      // res[0]=Rscore
      // res[1]=N*Ixyz
      // res[2]=N*kxyz
      // we do not want to add a z if x or y have only one bin
      bool ok = true;  // ok : do we compute I or return 0?
      if (samplesNotNA < environment.numSamples) {
        std::set<int> s;
        for (uint i = 0; i < 2 && ok; i++) {
          s.clear();
          for (uint j = 0; j < samplesNotNA; j++) {
            s.insert(dataNumeric_red[i][j]);
          }

          if (s.size() == 1) {
            ok = false;
            break;
          }
        }

        if (environment.testDistribution && numSamples_nonNA != samplesNotNA) {
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
        res = compute_Rscore_Ixyz_alg5(dataNumeric_red, dataNumericIdx_red,
            AllLevels_red, cnt_red, posArray_red, myNbrUi, myNbrUi + 2,
            samplesNotNA, sample_weights_red, flag_sample_weights, environment);
      } 
      else{
        res = new double[3];
        res[0] = -DBL_MAX;  // Rscore
        res[1] = 0;         // I
        res[2] = 1;         // k
      }
    }
  }  // jump cond no statistics

  scoresZ[iz] = res[0];
}

double* computeEnsInformationContinuous(Environment& environment, int* myCond,
    int myNbrUi, int* myZi, uint myNbrZi, int myZiPos, const int myVarIdxX,
    const int myVarIdxY, const int cplx, MemorySpace& m) {
  int* posArray = new int[2 + environment.edges[myVarIdxX][myVarIdxY]
                                  .shared_info->ui_vect_idx.size()];
  posArray[0] = myVarIdxX;
  posArray[1] = myVarIdxY;

  int nbrRetValues = 3;

  if (myCond != NULL) {
    if (myNbrUi > 0) {
      // The index of the variables in this dataset
      for (int i = 0; i < myNbrUi; ++i) posArray[i + 2] = myCond[i];
    }
  }

  int numSamples_nonNA =
      getNumSamples_nonNA(environment, posArray[0], posArray[1]);
  double* res_new;

  // initialization part (no z)
  if (myNbrZi == 0) {

    // TODO : speedup by only removing NAs for marked columns
    // Mark rows containing NAs and count the number of complete samples
    vector <int> sample_is_not_NA(environment.numSamples);
    vector <int> NAs_count(environment.numSamples);
    uint samplesNotNA = count_non_NAs(myNbrUi, sample_is_not_NA,
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
    res_new[2] = -DBL_MAX;
    res_new[1] = -1;
    res_new[0] = environment.numSamples;
    double* res;

    int z;
    int* ziContPosIdx = NULL;

    // If x, y and uis are discrete we can put togheter all zi that are discrete
    // in a vector and evaluate them in one shot as discrete.
    if (allVariablesDiscrete(
            environment.columnAsContinuous, posArray, (myNbrUi + 2))) {
      // search for z that are discrete
      int countZDiscrete = 0;
      for (uint iz = 0; iz < myNbrZi; iz++) {
        z = myZi[iz];
        if (environment.columnAsContinuous[z] == 0) countZDiscrete++;
      }

      if (countZDiscrete > 0) {
        int* posZi = new int[countZDiscrete];
        int* zz = new int[countZDiscrete];
        int pos = 0;
        for (uint iz = 0; iz < myNbrZi; iz++) {
          z = myZi[iz];
          if (environment.columnAsContinuous[z] == 0) {
            zz[pos] = z;
            posZi[pos] = iz;
            pos++;
          }
        }
        double** jointFreqs = getJointFreqs(
            environment, posArray[0], posArray[1]);
        res = getAllInfoNEW(environment.oneLineMatrix, environment.allLevels,
            posArray, myNbrUi, zz, countZDiscrete, -1, environment.numSamples,
            environment.effN, cplx, environment.isK23, environment.looklog,
            environment.c2terms, &m, environment.sampleWeights, jointFreqs,
            environment.testDistribution);

        for (uint level0 = 0; level0 < environment.allLevels[posArray[0]];
             level0++)
          delete[] jointFreqs[level0];
        delete[] jointFreqs;
        delete[] zz;

        // keep in ziContPos only the position of the continuous variables
        ziContPosIdx = new int[myNbrZi - countZDiscrete];
        pos = 0;
        for (uint iz = 0; iz < myNbrZi; iz++) {
          z = myZi[iz];
          if (environment.columnAsContinuous[z] == 1) {
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
        environment.firstIterationDone && myNbrZi > environment.nThreads;
#pragma omp parallel for if (parallelizable)
#endif
    for (uint iz = 0; iz < myNbrZi; iz++) {
      MemorySpace privateM = m;
#ifdef _OPENMP
      if (parallelizable)
        privateM = environment.memoryThreads[omp_get_thread_num()];
#endif
      int numSamples_nonNA =
          getNumSamples_nonNA(environment, posArray[0], posArray[1]);
      double* res;
      computeContributingScores(environment, ziContPosIdx, iz, myZi, myNbrUi,
          numSamples_nonNA, posArray, res, scoresZ, privateM);
    }  // parallel for on z

    for (uint iz = 0; iz < myNbrZi; iz++) {  // find optimal z
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

  delete[] posArray;

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
  bool test = false;

  int* posArray = new int[2 + environment.edges[myVarIdxX][myVarIdxY]
                                  .shared_info->ui_vect_idx.size()];
  posArray[0] = myVarIdxX;
  posArray[1] = myVarIdxY;

  if (myCond != NULL) {
    if (myNbrUi > 0) {
      // The index of the variables in this dataset
      for (int i = 0; i < myNbrUi; ++i) posArray[i + 2] = myCond[i];
    }
  }

  if (test) {
    if (myZi != NULL) {
      // The index of the variables in this dataset
      cout << "Number of zi :" << myNbrZi << endl;
      for (int i = 0; i < myNbrZi; i++) {
        cout << myZi[i] << " ";
      }
      cout << endl;
    }

    if (myCond != NULL) {
      cout << "Ui passed: ";

      for (int i = 0; i < myNbrUi; i++)
        cout << environment.nodes[myCond[i]].name << " ";
      cout << endl;

      // The index of the variables in this dataset
      cout << "Number of indexes with ui :" << myNbrUi << endl;
      for (int i = 0; i < myNbrUi + 2; i++) {
        cout << posArray[i] << " ";
      }
      cout << endl;
    }

    cout << "Levels\n";

    for (uint i = 0; i < environment.numNodes; i++) {
      cout << environment.allLevels[i] << " ";
    }
    cout << endl;

    cout << "Number of samples " << environment.numSamples << endl;
    cout << "Position x and y: " << posArray[0] << " " << posArray[1] << endl;
  }

// Compute the mutual information
#if _MY_DEBUG_NEW
  printf("\n# =====> before getAllInfoNEW \n");
#endif  // _MY_DEBUG_NEW

  int numSamples_nonNA =
      getNumSamples_nonNA(environment, posArray[0], posArray[1]);
  double** jointFreqs =
      getJointFreqs(environment, posArray[0], posArray[1]);

  double* res_new = getAllInfoNEW(environment.oneLineMatrix,
      environment.allLevels, posArray, myNbrUi, myZi, myNbrZi, myZiPos,
      environment.numSamples, environment.effN, cplx, environment.isK23,
      environment.looklog, environment.c2terms, &m, environment.sampleWeights,
      jointFreqs, environment.testDistribution);

  for (uint level0 = 0; level0 < environment.allLevels[posArray[0]]; level0++)
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

  delete[] posArray;

  return res_new;
}

void SearchForNewContributingNodeAndItsRank(
    Environment& environment, const int posX, const int posY, MemorySpace& m) {
  if (environment.edges[posX][posY].shared_info->zi_vect_idx.empty())
    return;

  int* ui;
  int* zi;

  if (environment.edges[posX][posY].shared_info->ui_vect_idx.empty())
    ui = NULL;
  else
    ui = &environment.edges[posX][posY].shared_info->ui_vect_idx[0];

  if (environment.edges[posX][posY].shared_info->zi_vect_idx.empty())
    zi = NULL;
  else
    zi = &environment.edges[posX][posY].shared_info->zi_vect_idx[0];

  int argEnsInfo = -1;
  if (environment.isK23 == true) argEnsInfo = environment.cplx;

  double* vect = NULL;

  if (environment.typeOfData == 0) {
    vect = computeEnsInformationNew(environment, ui,
        environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi,
        environment.edges[posX][posY].shared_info->zi_vect_idx.size(),
        environment.edges[posX][posY].shared_info->ui_vect_idx.size() + 2, posX,
        posY, argEnsInfo, m);
    if (vect[6] - environment.edges[posX][posY].shared_info->Rxyz_ui > 0) {
      if (environment.isVerbose) {
        cout << "\n"
             << posX << "    " << posY << "# -----> possible zi: "
             << environment
                    .nodes[environment.edges[posX][posY]
                               .shared_info->zi_vect_idx[vect[3]]]
                    .name
             << "(" << vect[6] << " > "
             << environment.edges[posX][posY].shared_info->Rxyz_ui << ")\n";
      }

      // The order matters: set first the z.name.idx, than get the corresponding
      // zi from the original vect / Doing this way, we make sure that the
      // z.name has the right bin xyzi key
      environment.edges[posX][posY].shared_info->z_name_idx = vect[3];
      environment.edges[posX][posY].shared_info->Rxyz_ui = vect[6];
      free(vect);

    } else if (environment.isVerbose) {
      cout << "# --!!--> Rxyz_ui.tmp = " << vect[6] << " < Rxyz_ui = "
           << environment.edges[posX][posY].shared_info->Rxyz_ui << "\n";
    }

  } else if (environment.typeOfData == 2 ||
             (environment.typeOfData == 1 && environment.isAllGaussian == 0)) {
    vect = computeEnsInformationContinuous(environment, ui,
        environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi,
        environment.edges[posX][posY].shared_info->zi_vect_idx.size(),
        environment.edges[posX][posY].shared_info->ui_vect_idx.size() + 2, posX,
        posY, argEnsInfo, m);
    if (vect[2] - environment.edges[posX][posY].shared_info->Rxyz_ui > 0) {
      if (environment.isVerbose) {
        cout << "\n"
             << posX << " " << posY << " # -----> possible zi: "
             << environment
                    .nodes[environment.edges[posX][posY]
                               .shared_info->zi_vect_idx[vect[1]]]
                    .name
             << "(" << vect[2] << " > "
             << environment.edges[posX][posY].shared_info->Rxyz_ui << ")\n";
      }

      // The order matters: set first the z.name.idx, than get the corresponding
      // zi from the original vect / Doing this way, we make sure
      // that the z.name has the right bin xyzi key
      environment.edges[posX][posY].shared_info->z_name_idx = vect[1];
      environment.edges[posX][posY].shared_info->Rxyz_ui = vect[2];

    } else if (environment.isVerbose) {
      cout << "# --!!--> Rxyz_ui.tmp = " << vect[2] << " < Rxyz_ui = "
           << environment.edges[posX][posY].shared_info->Rxyz_ui << "\n";
    }
    delete[] vect;
  }
}

void SearchForNewContributingNodeAndItsRankGaussian(
    Environment& environment, const int posX, const int posY, MemorySpace& m) {
  if (environment.edges[posX][posY].shared_info->zi_vect_idx.size() == 0)
    return;

  int nbrZi = environment.edges[posX][posY].shared_info->zi_vect_idx.size();

  if (nbrZi == 0) return;

  int* ui;
  int* zi;

  if (environment.edges[posX][posY].shared_info->ui_vect_idx.empty())
    ui = NULL;
  else
    ui = &environment.edges[posX][posY].shared_info->ui_vect_idx[0];

  if (environment.edges[posX][posY].shared_info->zi_vect_idx.empty())
    zi = NULL;
  else
    zi = &environment.edges[posX][posY].shared_info->zi_vect_idx[0];

  double* Ixy_ui_z = corrMutInfo(environment, environment.dataDouble, ui,
      environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi,
      environment.edges[posX][posY].shared_info->zi_vect_idx.size(), posX, posY,
      -1);

#if _MY_DEBUG_
  if (test) {
    cout << "Ixy_ui_z: ";
    for (int i = 0; i < 3 * nbrZi; i += 3) {
      cout << Ixy_ui_z[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
        // Get all I(zy|ui)[xyuiz]
  double* Izy_ui = corrMutInfo(environment, environment.dataDouble, ui,
      environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi,
      environment.edges[posX][posY].shared_info->zi_vect_idx.size(), posX, posY,
      0);
#if _MY_DEBUG_
  if (test) {
    cout << "Izy_ui: ";
    for (int i = 0; i < 3 * nbrZi; i += 3) {
      cout << Izy_ui[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_

  // Get all I(xz|ui)[xyuiz]
  double* Ixz_ui = corrMutInfo(environment, environment.dataDouble, ui,
      environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi,
      environment.edges[posX][posY].shared_info->zi_vect_idx.size(), posX, posY,
      1);
#if _MY_DEBUG_
  if (test) {
    cout << "Ixz_ui: ";
    for (int i = 0; i < 3 * nbrZi; i += 3) {
      cout << Ixz_ui[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
        // Get all I(xy|ui,z)[xyuiz]
  double* Ixy_uiz = corrMutInfo(environment, environment.dataDouble, ui,
      environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi,
      environment.edges[posX][posY].shared_info->zi_vect_idx.size(), posX, posY,
      2);
#if _MY_DEBUG_
  if (test) {
    cout << "Ixy_uiz: ";
    for (int i = 0; i < 3 * nbrZi; i += 3) {
      cout << Ixy_uiz[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
  // Compute I(xyz|ui)[xyuiz] = I(xy|ui)[xyuiz] - I(xy|ui,z)[xyuiz]
  double* Ixyz_ui_vect = computeDifference(Ixy_ui_z, Ixy_uiz, nbrZi);

#if _MY_DEBUG_
  if (test) {
    cout << "Ixyz_ui_vect: ";
    for (int i = 0; i < nbrZi; i++) {
      cout << Ixyz_ui_vect[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
  double* Ixz_ui_vect = Ixz_ui;
#if _MY_DEBUG_
  if (test) {
    cout << "Ixz_ui_vect: ";
    for (int i = 0; i < nbrZi; i++) {
      cout << Ixz_ui_vect[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
  double* Izy_ui_vect = Izy_ui;
#if _MY_DEBUG_
  if (test) {
    cout << "Izy_ui_vect: ";
    for (int i = 0; i < nbrZi; i++) {
      cout << Izy_ui_vect[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
  double* Ixy_ui_vect = Ixy_ui_z;
#if _MY_DEBUG_
  if (test) {
    cout << "Ixy_ui_vect: ";
    for (int i = 0; i < nbrZi; i++) {
      cout << Ixy_ui_vect[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
  double* xz = computeDifference(Ixz_ui_vect, Ixy_ui_vect, nbrZi);
#if _MY_DEBUG_
  if (test) {
    cout << "xz: ";
    for (int i = 0; i < nbrZi; i++) {
      cout << xz[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
  double* zy = computeDifference(Izy_ui_vect, Ixy_ui_vect, nbrZi);
#if _MY_DEBUG_
  if (test) {
    cout << "zy: ";
    for (int i = 0; i < nbrZi; i++) {
      cout << zy[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_
  double* xyz = Ixyz_ui_vect;
#if _MY_DEBUG_
  if (test) {
    cout << "xyz: ";
    for (int i = 0; i < nbrZi; i++) {
      cout << xyz[i] << ",";
    }
    cout << endl;
  }
#endif  // _MY_DEBUG_

  // For each zi, compute the rank
  double* score_vect = computeRankAndFindMin(xz, zy, xyz, nbrZi);
  // Display all the scores, the score max and the mutual information
  // corresponding to the two non-base edges
  if (environment.isVerbose) {
    cout << "# ------> scores: ";
    for (int i = 0; i < nbrZi; ++i) {
      cout << environment
                  .nodes[environment.edges[posX][posY]
                             .shared_info->zi_vect_idx[i]]
                  .name
           << ": " << score_vect[i];
      if (i + 1 < nbrZi) cout << ", ";
    }
    cout << endl;
  }
  // Get the index corresponding to the best rank
  int myZiBest_idx = 0;
  for (int i = 0; i < nbrZi; i++) {
    if (score_vect[i] > score_vect[myZiBest_idx]) myZiBest_idx = i;
  }
  // There can be more than one zi with the same rank... so, arbitrarly take
  // the first one
  if (score_vect[myZiBest_idx] > 0) {
    if (environment.isVerbose) {
      cout << "\n# -----> possible zi: "
           << environment
                  .nodes[environment.edges[posX][posY]
                             .shared_info->zi_vect_idx[myZiBest_idx]]
                  .name
           << "(" << score_vect[myZiBest_idx] << " > "
           << environment.edges[posX][posY].shared_info->Rxyz_ui << ")\n";
    }

    // The order matters: set first the z.name.idx, than get the corresponding
    // zi from the original vect / Doing this way, we make sure that the z.name
    // has the right bin xyzi key
    environment.edges[posX][posY].shared_info->z_name_idx = myZiBest_idx;
    environment.edges[posX][posY].shared_info->Rxyz_ui =
        score_vect[myZiBest_idx];

  } else if (environment.isVerbose) {
    cout << "# --!!--> Rxyz_ui.tmp = " << score_vect[myZiBest_idx]
         << " < Rxyz_ui = "
         << environment.edges[posX][posY].shared_info->Rxyz_ui << "\n";
  }

  delete[] Ixy_ui_z;
  delete[] Izy_ui;
  delete[] Ixz_ui;
  delete[] Ixy_uiz;
  delete[] xz;
  delete[] zy;
  delete[] score_vect;
  delete[] Ixyz_ui_vect;
  delete[] Ixz_ui_vect;
  delete[] Izy_ui_vect;
  delete[] Ixy_ui_vect;
  delete[] xyz;
}

double computeEnsInformationContinuous_Gaussian(
    Environment& environment, const int posX, const int posY, const int posZ) {
  if (environment.edges[posX][posY].shared_info->zi_vect_idx.size() == 0)
    return true;

  int* ui;
  if (environment.edges[posX][posY].shared_info->ui_vect_idx.empty())
    ui = NULL;
  else
    ui = &environment.edges[posX][posY].shared_info->ui_vect_idx[0];

  int* zi = new int[1];
  zi[0] = posZ;

  int nbrZi = 1;

  double* Ixy_ui = corrMutInfo(environment, environment.dataDouble, ui,
      environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi, 1,
      posX, posY, -1);

  // Get all I(zy|ui)[xyuiz]
  double* Izy_ui = corrMutInfo(environment, environment.dataDouble, ui,
      environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi, 1,
      posX, posY, 0);

  // Get all I(xz|ui)[xyuiz]
  double* Ixz_ui = corrMutInfo(environment, environment.dataDouble, ui,
      environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi, 1,
      posX, posY, 1);

  // Get all I(xy|ui,z)[xyuiz]
  double* Ixy_uiz = corrMutInfo(environment, environment.dataDouble, ui,
      environment.edges[posX][posY].shared_info->ui_vect_idx.size(), zi, 1,
      posX, posY, 2);

  // Compute I(xyz|ui)[xyuiz] = I(xy|ui)[xyuiz] - I(xy|ui,z)[xyuiz]
  double* Ixyz_ui_vect = computeDifference(Ixy_ui, Ixy_uiz, nbrZi);

  delete[] zi;
  delete[] Ixy_ui;
  delete[] Izy_ui;
  delete[] Ixz_ui;
  delete[] Ixy_uiz;

  delete[] Ixyz_ui_vect;

  return Ixyz_ui_vect[0];
}

}  // namespace computation
}  // namespace miic
