#include <Rcpp.h>

#include <algorithm>  // std::max_element
#include <array>
#include <numeric>  // std::iota

#include "get_information.h"
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
  int maxbins = environment.maxbins;

  size_t li_alloc_size = getLinearAllocatorSize(environment.n_samples,
      environment.n_nodes, maxbins, environment.initbins,
      environment.is_continuous, environment.levels);
  li_alloc_ptr = std::make_unique<LinearAllocator>(li_alloc_size);

  vector<int> ui_list(environment.n_nodes - 2);
  std::iota(begin(ui_list), end(ui_list), 2);

  std::shared_ptr<CutPointsInfo> cuts_ptr = nullptr;
  // X or Y is continuous, keep info on the discretization.
  if (environment.is_continuous[0] || environment.is_continuous[1]) {
    cuts_ptr = std::make_shared<CutPointsInfo>(kStepMax, maxbins * 2);
  }
  auto res = getCondMutualInfo(0, 1, ui_list, environment.data_numeric,
      environment.data_numeric_idx, environment, cuts_ptr);

  List result = List::create(
      _["info"]            = res.I,
      _["infok"]           = res.I - res.k);

  if (cuts_ptr != nullptr) {
    // Prepare cut points matrix
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

    result.push_back(cutpoints, "cutpointsmatrix");
    result.push_back(cuts_ptr->I_equal_freq_max, "efinfo");
  }

  delete li_alloc_ptr;
  return result;
}

// [[Rcpp::export]]
List miicRGetInfo3Point(List input_data, List arg_list) {
  // Initialize Environment with mandatory inputs
  Environment environment(as<int>(arg_list["n_samples"]),
      as<int>(arg_list["n_nodes"]), as<vector<int>>(input_data["factor"]),
      as<vector<int>>(input_data["order"]),
      as<vector<int>>(arg_list["is_continuous"]),
      as<vector<int>>(arg_list["levels"]));

  // Set optional parameters
  setEnvironmentFromR(input_data, arg_list, environment);
  int maxbins = environment.maxbins;

  size_t li_alloc_size = getLinearAllocatorSize(environment.n_samples,
      environment.n_nodes, maxbins, environment.initbins,
      environment.is_continuous, environment.levels);
  li_alloc_ptr = new LinearAllocator(li_alloc_size);

  vector<int> ui_list(environment.n_nodes - 3);
  std::iota(begin(ui_list), end(ui_list), 3);

  auto res = getInfo3Point(environment, 0, 1, 2, ui_list);

  List result = List::create(
      _["I3"]            = res.Ixyz_ui,
      _["I3k"]           = res.Ixyz_ui - res.kxyz_ui,
      _["I2"]            = res.Ixy_ui,
      _["I2k"]           = res.Ixy_ui - res.kxy_ui);

  return result;
}
