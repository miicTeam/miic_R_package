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

void setConfidence(Environment& environment) {
  vector<vector<int>> original_data(environment.data_numeric);
  vector<vector<int>> original_data_idx(environment.data_numeric_idx);

  vector<EdgeID> edge_list;
  for (int i = 1; i < environment.n_nodes; ++i) {
    for (int j = 0; j < i; ++j) {
      if (!environment.edges[i][j].status) continue;
      edge_list.emplace_back(j, i, environment.edges[i][j]);
    }
  }

  std::set<int> columns_to_shuffle;
  for (const auto& edge : edge_list) {
    columns_to_shuffle.insert(edge.i);
  }
  // [0, 1, ..., n_samples - 1]
  Rcpp::IntegerVector indices(Rcpp::seq(0, environment.n_samples - 1));
  Rcpp::IntegerVector shuffled(indices);
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
    for (const auto& edge : edge_list) {
      int X = edge.i, Y = edge.j;
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
      // n_shuffles times
      environment.edges[X][Y].shared_info->exp_shuffle += exp(-I_prime_shuffle);
    }
  }
  // evaluate average
  for (const auto& edge : edge_list) {
    environment.edges[edge.i][edge.j].shared_info->exp_shuffle /=
        environment.n_shuffles;
  }
  // Copy data back
  environment.data_numeric = std::move(original_data);
  environment.data_numeric_idx = std::move(original_data_idx);

  for (int i = 0; i < environment.n_samples; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      environment.oneLineMatrix[j * environment.n_samples + i] =
          environment.data_numeric[i][j];
    }
  }
}

vector<vector<string>> confidenceCut(Environment& environment) {
  vector<EdgeID>& edge_list = environment.connected_list;
  size_t n_connected = edge_list.size();
  // remove edges based on confidence
  auto to_delete = [&environment](EdgeID& id) {
    int X = id.i, Y = id.j;
    auto info = environment.edges[X][Y].shared_info;
    double I_prime_original = info->Ixy_ui - info->cplx;
    // exp(I_shuffle - I_original)
    auto confidence = exp(-I_prime_original) / info->exp_shuffle;
    if (confidence > environment.conf_threshold) {
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

  vector<vector<string>> res;
  res.emplace_back(std::initializer_list<string>{"x", "y", "confidence_ratio"});
  for (const auto& edge : edge_list) {
    auto info = environment.edges[edge.i][edge.j].shared_info;
    auto confidence = exp(info->cplx - info->Ixy_ui) / info->exp_shuffle;
    res.emplace_back(std::initializer_list<string>{
        environment.nodes[edge.i].name,
        environment.nodes[edge.j].name,
        std::to_string(confidence)});
  }
  std::sort(edge_list.begin(), edge_list.end());
  environment.numNoMore = edge_list.size();

  return res;
}

}  // namespace reconstruction
}  // namespace miic
