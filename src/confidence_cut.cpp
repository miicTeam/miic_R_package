#include "confidence_cut.h"

#include <Rcpp.h>

#include <set>

#include "compute_ens_information.h"
#include "environment.h"

namespace miic {
namespace reconstruction {

using std::string;
using std::vector;
using Rcpp::Rcout;
using namespace miic::computation;
using namespace miic::structure;

vector<vector<string>> confidenceCut(Environment& environment) {
  // Create a back up of the data, for later randomization
  vector<vector<int>> original_data(environment.data_numeric);
  vector<vector<int>> original_data_idx(environment.data_numeric_idx);

  // [0, 1, ..., n_samples - 1]
  Rcpp::IntegerVector indices(Rcpp::seq(0, environment.n_samples - 1));
  Rcpp::IntegerVector shuffled(indices);

  vector<EdgeID>& edge_list = environment.connected_list;
  int n_connected = edge_list.size();
  std::set<int> columns_to_shuffle;
  for (const auto& edge : edge_list) {
    columns_to_shuffle.insert(edge.i);
  }

  vector<double> confVect(n_connected, 0);
  for (int nb = 0; nb < environment.n_shuffles; nb++) {
    // Shuffle the dataset for selected columns
    for (auto col : columns_to_shuffle) {
      // random permutation
      shuffled = Rcpp::sample(indices, indices.size());
      for (int row = 0; row < environment.n_samples; row++) {
        environment.data_numeric[shuffled[row]][col] = original_data[row][col];
        if (environment.is_continuous[col]) {
          if (original_data_idx[col][row] == -1)
            environment.data_numeric_idx[col][row] = -1;
          else
            environment.data_numeric_idx[col][row] =
                shuffled[original_data_idx[col][row]];
        }
      }
    }

    for (int i = 0; i < environment.n_samples; i++) {
      for (int j = 0; j < environment.n_nodes; j++) {
        environment.oneLineMatrix[j * environment.n_samples + i] =
            environment.data_numeric[i][j];
      }
    }

    double NIxy_ui, k_xy_ui;
    // evaluate the mutual information for every edge
    for (int i = 0; i < n_connected; i++) {
      int X = edge_list[i].i;
      int Y = edge_list[i].j;
      if (!environment.is_continuous[X] && !environment.is_continuous[Y]) {
        // discrete case
        double* res = computeEnsInformationNew(environment, NULL, 0, NULL, 0,
            -1, X, Y, environment.cplx, environment.m);
        NIxy_ui = res[1];
        k_xy_ui = res[2];
        delete res;
      } else {
        // mixed case
        double* res = computeEnsInformationContinuous(environment, NULL, 0,
            NULL, 0, -1, X, Y, environment.cplx, environment.m);
        NIxy_ui = res[1];
        k_xy_ui = res[2];
        delete res;
      }

      double I_prime_shuffle = NIxy_ui - k_xy_ui;
      if (I_prime_shuffle < 0) {
        I_prime_shuffle = 0;
      }
      confVect[i] += exp(-I_prime_shuffle);  // n_shuffles times
    }
  }
  // evaluate the average confidence
  for (auto& c : confVect) {
    c /= environment.n_shuffles;
  }
  // remove edges based on confidence cut
  auto to_delete = [&environment, &confVect, &edge_list](EdgeID& id) {
    int X = id.i, Y = id.j;
    auto info = environment.edges[X][Y].shared_info;
    double I_prime_original = info->Ixy_ui - info->cplx;
    auto index = &id - &*begin(edge_list);
    // exp(I_shuffle - I_original)
    confVect[index] = exp(-I_prime_original) / confVect[index];
    if (confVect[index] > environment.conf_threshold) {
      info->connected = 0;
      environment.edges[X][Y].status = 0;
      environment.edges[Y][X].status = 0;
      return true;
    } else {
      return false;
    }
  };
  edge_list.erase(
      remove_if(begin(edge_list), end(edge_list), to_delete), end(edge_list));
  Rcout << "# -- number of edges cut: " << n_connected - edge_list.size()
            << "\n";

  // Copy data back
  environment.data_numeric = std::move(original_data);
  environment.data_numeric_idx = std::move(original_data_idx);

  for (int i = 0; i < environment.n_samples; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      environment.oneLineMatrix[j * environment.n_samples + i] =
          environment.data_numeric[i][j];
    }
  }

  vector<vector<string>> res;
  res.emplace_back(std::initializer_list<string>{"x", "y", "confidence_ratio"});
  for (int i = 0; i < n_connected; i++) {
    res.emplace_back(std::initializer_list<string>{
        environment.nodes[edge_list[i].i].name,
        environment.nodes[edge_list[i].j].name,
        std::to_string(confVect[i])});
  }
  std::sort(edge_list.begin(), edge_list.end());
  environment.numNoMore = edge_list.size();

  return res;
}
}  // namespace reconstruction
}  // namespace miic
