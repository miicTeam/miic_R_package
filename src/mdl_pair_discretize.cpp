#include <Rcpp.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>

#include "compute_info.h"
#include "info_cnt.h"
#include "mutual_information.h"
#include "utilities.h"

#define STEPMAX 50
#define STEPMAXOUTER 50
#define EPS 1e-5
#define FLAG_CPLX 1

using std::map;
using std::vector;
using namespace Rcpp;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;
using uint = unsigned int;

// Dealing with input variables
void transformToFactorsContinuous(
    double** data, int** dataNumeric, int i, int n) {
  std::multimap<double, int> myMap;

  // clean the dictionary since it is used column by column
  myMap.clear();

  vector<double> clmn;
  for (int j = 0; j < n; j++) {
    double entry = data[j][i];
    clmn.push_back(entry);
  }

  sort(clmn.begin(), clmn.end());

  for (size_t j = 0; j < clmn.size(); j++) {
    myMap.insert(std::pair<double, int>(clmn[j], j));
  }

  for (int j = 0; j < n; j++) {
    double entry = data[j][i];
    dataNumeric[j][i] = myMap.find(entry)->second;

    typedef std::multimap<double, int>::iterator iterator;
    std::pair<iterator, iterator> iterpair = myMap.equal_range(entry);
    iterator it = iterpair.first;
    for (; it != iterpair.second; ++it) {
      if (it->second == dataNumeric[j][i]) {
        myMap.erase(it);
        break;
      }
    }
  }
}

void transformToFactors(double** data, int** dataNumeric, int n, int i) {
  // create a dictionary to store the factors of the strings
  map<double, int> myMap;
  // clean the dictionary since it is used column by column
  myMap.clear();
  int factor = 0;

  for (int j = 0; j < n; j++) {
    map<double, int>::iterator it = myMap.find(data[j][i]);
    if (it != myMap.end()) {
      dataNumeric[j][i] = it->second;
    } else {
      myMap[data[j][i]] = factor;
      dataNumeric[j][i] = factor;
      factor++;
    }
  }
}

void transformToFactorsContinuousIdx(
    int** dataNumeric, int** dataNumericIdx, int n, int i) {
  map<int, int> myMap;
  // clean the dictionary since it is used column by column
  myMap.clear();
  // vector <int> clmn;
  for (int j = 0; j < n; j++) {
    int entry = dataNumeric[j][i];
    if (entry != -1) myMap[entry] = j;
  }

  int j = 0;
  for (map<int, int>::iterator it = myMap.begin(); it != myMap.end(); ++it) {
    dataNumericIdx[i][j] = it->second;
    j++;
  }
}

// [[Rcpp::export]]
List mydiscretizeMutual(SEXP RmyDist1, SEXP RmyDist2, SEXP RflatU,
    SEXP RnbrU, SEXP RmaxBins, SEXP Rinitbin, SEXP Rcplx, SEXP Rcnt_vec,
    SEXP Rnlevels, SEXP ReffN, SEXP RsampleWeights) {
  vector<double> myDist1Vec = Rcpp::as<vector<double> >(RmyDist1);
  vector<double> myDist2Vec = Rcpp::as<vector<double> >(RmyDist2);
  vector<double> cnt_vec = Rcpp::as<vector<double> >(Rcnt_vec);
  vector<double> nlevels = Rcpp::as<vector<double> >(Rnlevels);
  int maxbins = Rcpp::as<int>(RmaxBins);
  int init_bin = Rcpp::as<int>(Rinitbin);
  int cplx = Rcpp::as<int>(Rcplx);
  int nbrU = Rcpp::as<int>(RnbrU);
  int n = myDist1Vec.size();
  int effN = Rcpp::as<int>(ReffN);

  double* myDist1 = &myDist1Vec[0];
  double* myDist2 = &myDist2Vec[0];

  std::vector<double> vectorflatU = Rcpp::as<vector<double> >(RflatU);
  double** matrixU = new double*[nbrU];
  for (int l = 0; l < nbrU; l++) {
    matrixU[l] = new double[n];
  }
  if (nbrU > 0) {
    for (int l = 0; l < nbrU; l++) {
      matrixU[l] = new double[n];
      copy(vectorflatU.begin() + (l * n), vectorflatU.begin() + ((l + 1) * n),
          matrixU[l]);
    }
  }
  // create the data matrix for factors
  int i;
  int j;
  int l;
  double** data = new double*[n];
  for (i = 0; i < n; i++) {
    data[i] = new double[nbrU + 2];
    data[i][0] = myDist1[i];
    data[i][1] = myDist2[i];
    for (l = 0; l < nbrU; l++) {
      data[i][l + 2] = matrixU[l][i];
    }
  }

  int** dataNumeric = new int*[n];
  for (j = 0; j < n; j++) {
    dataNumeric[j] = new int[nbrU + 2];
  }
  for (i = 0; i < (nbrU + 2); i++) {
    if (cnt_vec[i] == 1)
      // update environment.dataNumeric not taking into account repetition
      transformToFactorsContinuous(data, dataNumeric, i, n);
    else
      transformToFactors(data, dataNumeric, n, i);
  }

  // sortidx for continuous
  // create the data matrix for factors indexes
  int** dataNumericIdx = new int*[nbrU + 2];
  for (i = 0; i < (nbrU + 2); i++) {
    dataNumericIdx[i] = new int[n];
    for (j = 0; j < n; j++) {
      dataNumericIdx[i][j] = -1;
    }
  }

  for (i = 0; i < (nbrU + 2); i++) {
    if (cnt_vec[i] == 1) {
      transformToFactorsContinuousIdx(dataNumeric, dataNumericIdx, n, i);
      // update environment.dataNumeric taking into account repetition
      transformToFactors(data, dataNumeric, n, i);
    }
  }
  // AllLevels
  uint* AllLevels = new uint[nbrU + 2];
  for (i = 0; i < (nbrU + 2); i++) {
    AllLevels[i] = nlevels[i];
  }
  // nrbUi
  int nbrUi = nbrU;
  // c2terms
  double* c2terms = new double[n + 1];
  for (int i = 0; i < n + 1; i++) {
    c2terms[i] = -1;
  }
  // lookup tables
  double* looklog = new double[n + 2];
  looklog[0] = 0.0;
  for (int i = 1; i < n + 2; i++) {
    looklog[i] = log(1.0 * i);
  }

  double* lookH = new double[n + 2];
  lookH[0] = 0.0;
  for (int i = 1; i < n + 2; i++) {
    lookH[i] = i * looklog[i];  //-(i+1)*looklog[(i+1)];
  }
  double** looklbc = new double*[n + 1];
  for (int k = 0; k < n + 1; k++) {
    looklbc[k] = new double[maxbins + 1];
    for (int z = 0; z < maxbins + 1; z++) {
      looklbc[k][z] = 0;
    }
  }

  // Declare the lookup table of the parametric complexity
  int ncol = std::max(maxbins, init_bin) + 1;
  for (int j = 0; j < (nbrU + 2); j++) {
    if (cnt_vec[j] == 0) ncol = std::max(ncol, int(AllLevels[j] + 1));
  }
  ncol = 1000;  // combinations of Us can exceed ncol
  double** sc_look = new double*[ncol];
  for (int K = 0; K < (ncol); K++) {
    sc_look[K] = new double[n + 1];
    for (int i = 0; i < (n + 1); i++) {
      if (K == 1)
        sc_look[K][i] = 0;
      else if (i == 0)
        sc_look[K][i] = 0;
      else
        sc_look[K][i] = -1;
    }
  }
  for (int i = 0; i < n; i++) {
    computeLogC(i, 2, looklog, sc_look);  // Initialize the c2 terms
  }

  double** lookchoose = new double*[ncol];
  for (int K = 0; K < (ncol); K++) {
    lookchoose[K] = new double[n + 1];
    for (int i = 0; i < (n + 1); i++) {
      lookchoose[K][i] = -1;
    }
  }

  miic::structure::Environment environment;


  environment.sample_weights = Rcpp::as<vector<double>>(RsampleWeights);
  if (environment.sample_weights.empty()) {
    double uniform_weight(1);
    if (effN != n)
      uniform_weight = (effN * 1.0) / n;
    environment.sample_weights.resize(n, uniform_weight);
  }

  environment.maxbins = maxbins;
  environment.initbins = init_bin;
  environment.lookH = lookH;
  environment.looklog = looklog;
  environment.cterms = sc_look;
  environment.c2terms = c2terms;
  environment.lookchoose = lookchoose;
  environment.cplx = cplx;
  environment.n_eff = effN;
  environment.n_samples = n;
  environment.data_numeric = dataNumeric;
  environment.data_numeric_idx = dataNumericIdx;
  environment.levels = AllLevels;
  environment.is_continuous.resize(nbrU+2);

  // Declare tables_red
  vector<int> posArray(nbrU + 2);
  for (i = 0; i < (nbrU + 2); i++) {
    posArray[i] = i;
    environment.is_continuous[i] = cnt_vec[i];
  }

  // Mark rows containing NAs and count the number of complete samples
  vector<int> sample_nonNA(environment.n_samples);
  vector<int> NAs_count(environment.n_samples);
  uint samplesNotNA = count_non_NAs(nbrU, sample_nonNA,
    NAs_count, posArray, environment);

  // Allocate data reducted *_red without rows containing NAs
  // All *_red variables are passed to the optimization routine
  vector<int> AllLevels_red(nbrU+2);
  vector<int> cnt_red(nbrU+2);
  vector<int> posArray_red(nbrU+2);
  vector<double> sample_weights_red(samplesNotNA);
  vector<vector<int> > dataNumeric_red(nbrU+2, vector<int> (samplesNotNA));
  vector<vector<int> > dataNumericIdx_red(nbrU+2, vector<int> (samplesNotNA));

  bool flag_sample_weights = filter_NAs(nbrU, AllLevels_red, cnt_red,
    posArray_red, posArray, dataNumeric_red, dataNumericIdx_red,
    sample_weights_red, sample_nonNA, NAs_count,
    environment);

  int** iterative_cuts = (int**)calloc(STEPMAX + 1, sizeof(int*));
  for (int i = 0; i < STEPMAX + 1; i++) {
    iterative_cuts[i] = (int*)calloc(maxbins * (2 + nbrU), sizeof(int));
  }
  environment.iterative_cuts = iterative_cuts;

  double* res = compute_mi_cond_alg1(dataNumeric_red, dataNumericIdx_red,
      AllLevels_red, cnt_red, posArray_red, nbrUi, n, sample_weights_red,
      flag_sample_weights, environment, true);

  for (int i = 0; i < n; i++) {
    delete[] dataNumeric[i];
  }
  delete[] dataNumeric;

  for (int i = 0; i < (nbrU + 2); i++) {
    delete[] dataNumericIdx[i];
  }
  delete[] dataNumericIdx;

  delete[] AllLevels;

  int niterations = 0;
  double max_res_ef;
  int** iterative_cutpoints = new int*[STEPMAX * maxbins];
  for (int i = 0; i < (STEPMAX * maxbins); i++) {
    iterative_cutpoints[i] = new int[nbrUi + 2];
  }
  for (int l = 0; l < STEPMAX + 1; l++) {
    if (iterative_cuts[l][0] == -1) {
      niterations = l;
      res[1] = iterative_cuts[l][1] / 100000.0;
      res[0] = iterative_cuts[l][2] / 100000.0;
      max_res_ef = iterative_cuts[l][3] / 100000.0;
      break;
    }
    for (int k = 0; k < (nbrUi + 2); k++) {
      i = 0;
      while (iterative_cuts[l][i + maxbins * k] <
             iterative_cuts[l][i + maxbins * k + 1]) {
        iterative_cutpoints[maxbins * l + i][k] =
            iterative_cuts[l][i + maxbins * k];
        i++;
      }
      for (int j = i; j < maxbins; j++) {
        iterative_cutpoints[maxbins * l + j][k] = -1;
      }
    }
  }

  NumericMatrix cutpoints(niterations * maxbins, nbrUi + 2);
  for (i = 0; i < cutpoints.nrow(); i++) {
    for (j = 0; j < (nbrUi + 2); j++) {
      cutpoints[i + j * cutpoints.nrow()] = iterative_cutpoints[i][j];  //
    }
  }

  List result = List::create(
      _["cutpointsmatrix"] = cutpoints,
      _["info"]            = res[0],
      _["infok"]           = res[1],
      _["efinfo"]          = max_res_ef);

  for (int i = 0; i < (STEPMAX + 1); i++) {
    free(iterative_cuts[i]);
  }
  free(iterative_cuts);

  for (int i = 0; i < ncol; i++) {
    delete[] sc_look[i];
    delete[] lookchoose[i];
  }
  delete[] sc_look;
  delete[] lookchoose;
  delete[] c2terms;

  for (int i = 0; i < n; i++) {
    delete[] data[i];
  }
  delete[] data;

  delete[] looklog;
  delete[] lookH;

  for (int l = 0; l < nbrU; l++) {
    delete[] matrixU[l];
  }
  delete[] matrixU;

  for (int i = 0; i < n + 1; i++) {
    delete[] looklbc[i];
  }
  delete[] looklbc;
  for (int i = 0; i < STEPMAX * maxbins; i++) {
    delete[] iterative_cutpoints[i];
  }
  delete[] iterative_cutpoints;

  return result;
}
