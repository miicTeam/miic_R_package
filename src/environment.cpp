#include "environment.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>  // std::generate, std::any_of
#include <utility>    // std::move

namespace miic {
namespace structure {

using std::vector;
using namespace miic::utility;
using namespace miic::computation;

Environment::Environment(int n_samples, int n_nodes, vector<int> vec_numeric,
    vector<int> vec_index, vector<int> is_continuous_, vector<int> levels_)
    : n_samples(n_samples),
      n_nodes(n_nodes),
      data_numeric(n_nodes, n_samples, std::move(vec_numeric)),
      data_numeric_idx(n_nodes, n_samples, std::move(vec_index)),
      is_continuous(std::move(is_continuous_)),
      levels(std::move(levels_)),
      has_na(n_nodes, 0),
      n_eff(n_samples),
      edges(n_nodes, n_nodes),
      noise_vec(2 * n_samples),
      initbins(std::min(30, int(0.5 + std::cbrt(n_samples)))),
      cache(n_samples) {
  for (int i = 0; i < n_nodes; ++i) {
    has_na[i] = std::any_of(data_numeric.row_begin(i), data_numeric.row_end(i),
        [](int value) { return value == -1; });
  }
  for (int i = 0; i < n_nodes; ++i) {
    for (int j = 0; j < n_nodes; ++j) {
      if ((!is_continuous[i] && levels[i] == n_samples) ||
          (!is_continuous[j] && levels[j] == n_samples)) {
        // If a node is discrete with as many levels as there are samples, its
        // information with other nodes is null.
        edges(i, j).status = 0;
        edges(i, j).status_prev = 0;
      } else {
        // Initialise all other edges.
        edges(i, j).status = 1;
        edges(i, j).status_prev = 1;
      }
    }
  }
  for (int i = 0; i < n_nodes; ++i) {
    edges(i, i).status = 0;
    edges(i, i).status_prev = 0;
  }
}

void Environment::readBlackbox(const Grid2d<int>& node_list) {
  int n_pairs = node_list.n_rows();
  for (int i = 0; i < n_pairs; ++i) {
    const auto pair = node_list.getConstRow(i);
    edges(pair[0], pair[1]).status = 0;
    edges(pair[0], pair[1]).status_prev = 0;
    edges(pair[1], pair[0]).status = 0;
    edges(pair[1], pair[0]).status_prev = 0;
  }
}

}  // namespace structure
}  // namespace miic
