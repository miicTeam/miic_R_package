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

using Rcpp::_;
using Rcpp::as;
using Rcpp::DataFrame;
using Rcpp::List;
using std::string;
using std::vector;
using namespace Rcpp;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;

// [[Rcpp::export]]
List mydiscretizeMutual(List input_data, List arg_list){

  //List mydiscretizeMutual(SEXP RmyDist1, SEXP RmyDist2, SEXP RflatU,
  //  SEXP RnbrU, SEXP RmaxBins, SEXP Rinitbin, SEXP Rcplx, SEXP Rcnt_vec,
  //  SEXP Rnlevels, SEXP ReffN, SEXP RsampleWeights) {

  miic::structure::Environment environment(input_data, arg_list);
  int maxbins = environment.maxbins;
  int nbrU = environment.n_nodes-2;

  vector<int> posArray(nbrU+2);
  for(int i=0; i<(nbrU+2); i++) posArray[i] = i;
  vector<int> temp_ui_list(nbrU);
  std::iota(begin(temp_ui_list), end(temp_ui_list), 2);

  // Mark rows containing NAs and count the number of complete samples
  vector<int> sample_nonNA(environment.n_samples);
  vector<int> NAs_count(environment.n_samples);
  int samplesNotNA =
      count_non_NAs(0, 1, temp_ui_list, sample_nonNA, NAs_count, environment);

  // Allocate data reducted *_red without rows containing NAs
  // All *_red variables are passed to the optimization routine
  vector<int> AllLevels_red(nbrU+2);
  vector<int> cnt_red(nbrU+2);
  vector<int> posArray_red(nbrU+2);
  vector<double> sample_weights_red(samplesNotNA);
  vector<vector<int> > dataNumeric_red(nbrU+2, vector<int> (samplesNotNA));
  vector<vector<int> > dataNumericIdx_red(nbrU+2, vector<int> (samplesNotNA));

  bool flag_sample_weights = filter_NAs(0, 1, temp_ui_list, AllLevels_red,
      cnt_red, posArray_red, dataNumeric_red, dataNumericIdx_red,
      sample_weights_red, sample_nonNA, NAs_count, environment);

  int** iterative_cuts = (int**)calloc(STEPMAX + 1, sizeof(int*));
  for (int i = 0; i < STEPMAX + 1; i++) {
    iterative_cuts[i] = (int*)calloc(maxbins * (2 + nbrU), sizeof(int));
  }
  environment.iterative_cuts = iterative_cuts;

  double* res = compute_mi_cond_alg1(dataNumeric_red, dataNumericIdx_red,
      AllLevels_red, cnt_red, posArray_red, nbrU, environment.n_samples,
      sample_weights_red, flag_sample_weights, environment, true);

  int niterations = 0;
  int i = 0;
  double max_res_ef;
  vector<vector<int> > iterative_cutpoints(
    STEPMAX * maxbins, vector<int>(nbrU + 2));
  for (int l = 0; l < STEPMAX + 1; l++) {
    if (iterative_cuts[l][0] == -1) {
      niterations = l;
      res[1] = iterative_cuts[l][1] / 100000.0;
      res[0] = iterative_cuts[l][2] / 100000.0;
      max_res_ef = iterative_cuts[l][3] / 100000.0;
      break;
    }
    for (int k = 0; k < (nbrU + 2); k++) {
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

  NumericMatrix cutpoints(niterations * maxbins, nbrU + 2);
  for (int i = 0; i < cutpoints.nrow(); i++) {
    for (int j = 0; j < (nbrU + 2); j++) {
      cutpoints[i + j * cutpoints.nrow()] = iterative_cutpoints[i][j];
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

  return result;
}
