#include "confidence_cut.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>

#include <algorithm>  // std::sort
#include <set>

#include "get_information.h"
#include "environment.h"

namespace miic {
namespace reconstruction {

using std::string;
using std::vector;
using namespace miic::computation;
using namespace miic::structure;

// shuffle with R rng (std::random_shuffle is deprecated)
// https://en.cppreference.com/w/cpp/algorithm/random_shuffle
template <class RandomIt>
void rShuffle(RandomIt first, RandomIt last) {
  typename std::iterator_traits<RandomIt>::difference_type i, n;
  n = last - first;
  for (i = n - 1; i > 0; --i) {
    std::swap(first[i], first[floor(R::unif_rand() * (i + 1))]);
  }
}

void setConfidence(Environment& environment) {
  vector<EdgeID> edge_list;
  std::set<int> columns_to_shuffle;
  for (int i = 1; i < environment.n_nodes; ++i) {
    for (int j = 0; j < i; ++j) {
      auto& edge = environment.edges(i, j);
      if (!edge.status || edge.shared_info->exp_shuffle != -1) continue;

      edge.shared_info->exp_shuffle = 0;
      edge_list.emplace_back(i, j, edge);
      columns_to_shuffle.insert(j);
    }
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    Grid2d<int> shuffled_data(environment.data_numeric);
    Grid2d<int> shuffled_data_idx(environment.data_numeric_idx);
    vector<int> indices(environment.n_samples);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int nb = 0; nb < environment.n_shuffles; nb++) {
      // Shuffle the dataset for selected columns
      for (auto col : columns_to_shuffle) {
        // [0, 1, ..., n_samples - 1]
        std::iota(begin(indices), end(indices), 0);
        // random permutation
        rShuffle(begin(indices), end(indices));
        for (int row = 0; row < environment.n_samples; row++) {
          shuffled_data(col, indices[row]) = environment.data_numeric(col, row);
          if (environment.is_continuous[col]) {
            if (environment.data_numeric_idx(col, row) == -1)
              shuffled_data_idx(col, row) = -1;
            else
              shuffled_data_idx(col, row) =
                  indices[environment.data_numeric_idx(col, row)];
          }
        }
      }
      double NIxy_ui, k_xy_ui;
      // evaluate the mutual information for every edge
      for (const auto& edge : edge_list) {
        int X = edge.X, Y = edge.Y;
        auto xy_ui = getCondMutualInfo(
            X, Y, vector<int>(), shuffled_data, shuffled_data_idx, environment);
        NIxy_ui = xy_ui.I;
        k_xy_ui = xy_ui.k;

        double I_prime_shuffle = NIxy_ui - k_xy_ui;
        if (I_prime_shuffle < 0) {
          I_prime_shuffle = 0;
        }
        // n_shuffles times
#ifdef _OPENMP
#pragma omp atomic
#endif
        edge.getEdge().shared_info->exp_shuffle += exp(-I_prime_shuffle);
      }
    }
  }  // omp parallel
  // evaluate average
  for (const auto& edge : edge_list)
    edge.getEdge().shared_info->exp_shuffle /= environment.n_shuffles;
}

void confidenceCut(Environment& environment) {
  vector<EdgeID>& edge_list = environment.connected_list;
  // remove edges based on confidence
  auto to_delete = [&environment](EdgeID& id) {
    int X = id.X, Y = id.Y;
    auto info = id.getEdge().shared_info;
    double I_prime_original = info->Ixy_ui - info->kxy_ui;
    // exp(I_shuffle - I_original)
    auto confidence = exp(-I_prime_original) / info->exp_shuffle;
    if (confidence > environment.conf_threshold) {
      info->connected = 0;
      environment.edges(X, Y).status = 0;
      environment.edges(Y, X).status = 0;
      environment.edges(X, Y).proba_head = -1;
      environment.edges(Y, X).proba_head = -1;
      return true;
    } else {
      return false;
    }
  };
  edge_list.erase(
      remove_if(begin(edge_list), end(edge_list), to_delete), end(edge_list));
  std::sort(edge_list.begin(), edge_list.end());
}

}  // namespace reconstruction
}  // namespace miic
