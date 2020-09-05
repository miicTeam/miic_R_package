#include <Rcpp.h>

#include <algorithm>  // std::max_element
#include <numeric>  // std::iota

#include "info_cnt.h"
#include "linear_allocator.h"
#include "utilities.h"

#define STEPMAX 50

using Rcpp::_;
using Rcpp::List;
using Rcpp::NumericMatrix;
using std::vector;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;

// [[Rcpp::export]]
List mydiscretizeMutual(List input_data, List arg_list) {
  Environment environment(input_data, arg_list);
  int maxbins = environment.maxbins;
  int nbrU = environment.n_nodes - 2;

  int max_level =
      *std::max_element(begin(environment.levels), end(environment.levels));
  size_t li_alloc_size = 8192;  // in bytes, extra space for fragmentation
  size_t n_integers = 4 * (environment.n_samples + 2) + 6 * (max_level + 1) +
                      3 * ((environment.n_samples + 1) * 7) +
                      4 * environment.n_nodes +
                      2 * (max_level + 1) * (max_level + 1) +
                      2 * max_level * max_level * max_level;  // coarse_counts
  size_t n_doubles = 3 * (max_level + 1) + 2 * environment.n_nodes;
  li_alloc_size += sizeof(int) * n_integers;
  li_alloc_size += sizeof(double) * n_doubles;

  li_alloc_ptr = std::make_unique<LinearAllocator>(li_alloc_size);

  vector<int> ui_list(nbrU);
  std::iota(begin(ui_list), end(ui_list), 2);

  // Mark rows containing NAs and count the number of complete samples
  TempVector<int> sample_nonNA(environment.n_samples);
  TempVector<int> NAs_count(environment.n_samples);
  int samplesNotNA = countNonNA(0, 1, /*Z*/ -1, ui_list,
      environment.data_numeric, sample_nonNA, NAs_count);

  // Allocate data reducted *_red without rows containing NAs
  // All *_red variables are passed to the optimization routine
  TempVector<int> AllLevels_red(nbrU + 2);
  TempVector<int> cnt_red(nbrU + 2);
  TempVector<int> posArray_red(nbrU + 2);
  TempVector<double> sample_weights_red(samplesNotNA);
  TempGrid2d<int> dataNumeric_red(nbrU + 2, samplesNotNA);
  TempGrid2d<int> dataNumericIdx_red(nbrU + 2, samplesNotNA);

  bool flag_sample_weights = filterNA(/*X*/ 0, /*Y*/ 1, /*Z*/ -1, ui_list,
      environment.data_numeric, environment.data_numeric_idx,
      environment.is_continuous, environment.sample_weights, sample_nonNA,
      NAs_count, dataNumeric_red, dataNumericIdx_red, AllLevels_red, cnt_red,
      posArray_red, sample_weights_red);

  Grid2d<int> iterative_cuts(STEPMAX + 1, maxbins * (2 + nbrU));
  environment.iterative_cuts = iterative_cuts;

  double* res = compute_mi_cond_alg1(dataNumeric_red, dataNumericIdx_red,
      AllLevels_red, cnt_red, posArray_red, nbrU, environment.n_samples,
      sample_weights_red, flag_sample_weights, environment, true);

  int niterations = 0;
  int i = 0;
  double max_res_ef = -1;
  TempGrid2d<int> iterative_cutpoints(STEPMAX * maxbins, nbrU + 2);
  for (int l = 0; l < STEPMAX + 1; l++) {
    if (iterative_cuts(l, 0) == -1) {
      niterations = l;
      res[1] = iterative_cuts(l, 1) / 100000.0;
      res[0] = iterative_cuts(l, 2) / 100000.0;
      max_res_ef = iterative_cuts(l, 3) / 100000.0;
      break;
    }
    for (int k = 0; k < (nbrU + 2); k++) {
      i = 0;
      while (iterative_cuts(l, i + maxbins * k) <
             iterative_cuts(l, i + maxbins * k + 1)) {
        iterative_cutpoints(maxbins * l + i, k) =
            iterative_cuts(l, i + maxbins * k);
        i++;
      }
      for (int j = i; j < maxbins; j++) {
        iterative_cutpoints(maxbins * l + j, k) = -1;
      }
    }
  }

  NumericMatrix cutpoints(niterations * maxbins, nbrU + 2);
  for (int i = 0; i < cutpoints.nrow(); i++) {
    for (int j = 0; j < (nbrU + 2); j++) {
      cutpoints[i + j * cutpoints.nrow()] = iterative_cutpoints(i, j);
    }
  }

  List result = List::create(
      _["cutpointsmatrix"] = cutpoints,
      _["info"]            = res[0],
      _["infok"]           = res[1],
      _["efinfo"]          = max_res_ef);

  return result;
}
