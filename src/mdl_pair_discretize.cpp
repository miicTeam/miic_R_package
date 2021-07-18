#include <Rcpp.h>

#include <algorithm>  // std::max_element
#include <array>
#include <numeric>  // std::iota

#include "computation_continuous.h"
#include "environment.h"
#include "linear_allocator.h"
#include "r_cpp_interface.h"
#include "utilities.h"

constexpr int kStepMax = 50;

using Rcpp::_;
using Rcpp::as;
using Rcpp::List;
using Rcpp::NumericMatrix;
using std::vector;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;

// [[Rcpp::export]]
List mydiscretizeMutual(List input_data, List arg_list) {
  // Initialize Environment with mandatory inputs
  Environment environment(as<int>(arg_list["n_samples"]),
      as<int>(arg_list["n_nodes"]), as<vector<int>>(input_data["factor"]),
      as<vector<int>>(input_data["order"]),
      as<vector<int>>(arg_list["is_continuous"]),
      as<vector<int>>(arg_list["levels"]));

  // Set optional parameters
  setEnvironmentFromR(input_data, arg_list, environment);
  int nbrU = environment.n_nodes - 2;
  int maxbins = environment.maxbins;

  size_t li_alloc_size = getLinearAllocatorSize(environment.n_samples,
      environment.n_nodes, maxbins, environment.initbins,
      environment.is_continuous, environment.levels);
  li_alloc_ptr = new LinearAllocator(li_alloc_size);

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

  bool any_na = environment.has_na[0] || environment.has_na[1];
  bool flag_sample_weights = filterNA(/*X*/ 0, /*Y*/ 1, /*Z*/ -1, ui_list,
      environment.data_numeric, environment.data_numeric_idx,
      environment.levels, environment.is_continuous, environment.sample_weights,
      sample_nonNA, NAs_count, dataNumeric_red, dataNumericIdx_red,
      AllLevels_red, cnt_red, posArray_red, sample_weights_red, any_na);

  auto cuts_ptr = std::make_shared<CutPointsInfo>(kStepMax, maxbins * 2);

  computeCondMutualInfo(dataNumeric_red, dataNumericIdx_red, AllLevels_red,
      cnt_red, posArray_red, sample_weights_red, flag_sample_weights,
      environment.initbins, environment.maxbins, environment.cplx,
      environment.cache.cterm, cuts_ptr);

  int niterations = cuts_ptr->n_iterations;
  TempGrid2d<int> iterative_cuts(kStepMax * maxbins, 2);
  const auto& cuts = cuts_ptr->cutpoints;
  for (int l = 0; l < kStepMax; ++l) {
    for (int k = 0; k < 2; ++k) {
      int i = 0;
      while (cuts(l, i + maxbins * k) < cuts(l, i + maxbins * k + 1)) {
        iterative_cuts(maxbins * l + i, k) = cuts(l, i + maxbins * k);
        ++i;
      }
      for (int j = i; j < maxbins; j++) {
        iterative_cuts(maxbins * l + j, k) = -1;
      }
    }
  }

  NumericMatrix cutpoints(niterations * maxbins, 2);
  for (int i = 0; i < cutpoints.nrow(); ++i) {
    for (int j = 0; j < 2; ++j) {
      cutpoints[i + j * cutpoints.nrow()] = iterative_cuts(i, j);
    }
  }

  List result = List::create(
      _["cutpointsmatrix"] = cutpoints,
      _["info"]            = cuts_ptr->I,
      _["infok"]           = cuts_ptr->Ik,
      _["efinfo"]          = cuts_ptr->I_equal_freq_max);

  delete li_alloc_ptr;
  return result;
}
