#include <Rcpp.h>

#include <algorithm>  // std::max_element
#include <array>
#include <numeric>  // std::iota

#include "computation_continuous.h"
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
  int nbrU = environment.n_nodes - 2;
  int maxbins = environment.maxbins;

  size_t li_alloc_size = getLinearAllocatorSize(environment.n_samples,
      environment.n_nodes, maxbins, environment.initbins,
      environment.is_continuous, environment.levels);
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

  environment.iterative_cuts = Grid2d<int>(STEPMAX + 1, maxbins * (2 + nbrU));

  computeCondMutualInfo(dataNumeric_red, dataNumericIdx_red, AllLevels_red,
      cnt_red, posArray_red, nbrU, environment.n_samples, sample_weights_red,
      flag_sample_weights, environment, true);

  int niterations = 0;
  int i = 0;
  double max_res_ef = -1;
  TempGrid2d<int> iterative_cutpoints(STEPMAX * maxbins, nbrU + 2);
  std::array<double, 2> res{0, 0};
  for (int l = 0; l < STEPMAX + 1; l++) {
    if (environment.iterative_cuts(l, 0) == -1) {
      niterations = l;
      res[1] = environment.iterative_cuts(l, 1) / 100000.0;
      res[0] = environment.iterative_cuts(l, 2) / 100000.0;
      max_res_ef = environment.iterative_cuts(l, 3) / 100000.0;
      break;
    }
    for (int k = 0; k < (nbrU + 2); k++) {
      i = 0;
      while (environment.iterative_cuts(l, i + maxbins * k) <
             environment.iterative_cuts(l, i + maxbins * k + 1)) {
        iterative_cutpoints(maxbins * l + i, k) =
            environment.iterative_cuts(l, i + maxbins * k);
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
