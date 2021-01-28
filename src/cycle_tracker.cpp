#include "cycle_tracker.h"

#include <algorithm>  // std::sort, std::none_of
#include <set>
#include <tuple>  // std::tie

#include "linear_allocator.h"
#include "structure.h"

namespace miic {
namespace reconstruction {
namespace detail {

using structure::TempVector;
using utility::TempAllocatorScope;

bool CycleTracker::hasCycle() {
  TempAllocatorScope scope;

  int n_edge = edge_list_.size();
  int n_nodes = edges_.n_rows();
  // Before saving the current iteration, search among previous iterations
  // those with the same number of edges
  auto range = edge_index_map_.equal_range(n_edge);
  bool no_cycle_found = range.first == range.second;
  // Indices of iteration that is possibly the end point of a cycle in the
  // backtracking sense. Example: suppose that iteration #1, #3, #6 have the
  // same number of edges as the current iteration #8, and each of them is a
  // possible start point of a cycle, therefore #2, #4, #7 are the corresponding
  // possible end points of the cycle (#8 #7 ... #2), (#8 #7 ... #4), (#8 #7).
  std::deque<int> iter_indices;
  for (auto it = range.first; it != range.second; ++it)
    iter_indices.push_back(it->second + 1);
  saveIteration();
  if (no_cycle_found) return false;
  // Backtracking requires starting from the largest index first
  std::sort(begin(iter_indices), end(iter_indices), std::greater<int>());
  // Whether the status of an edge has changed for at least once among all
  // iterations, if true the edge will be marked as connected and undirected
  TempGrid2d<int> has_changed(n_nodes, n_nodes, 0);
  // Whether the status of an edge is different between the current iteration
  // and the iteration before the one under consideration.
  TempVector<int> changed(n_nodes * n_nodes, 0);
  // Backtracking over iteration to get changed_edges
  cycle_size = 0;
  for (const auto& iter : iterations_) {
    ++cycle_size;
    for (const auto& k : iter.changed_edges) {
      int i{0}, j{0};
      std::tie(i, j) = getIndex2D(k.first, n_nodes);
      has_changed(i, j) = 1;
      // compare the status in the previous iteration against the latest status
      changed[k.first] = (k.second != edges_(i, j).status);
    }
    if (iter.index != iter_indices.front()) continue;
    iter_indices.pop_front();
    using std::none_of;
    if (none_of(begin(changed), end(changed), [](int j) { return j != 0; })) {
      for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
          if (has_changed(i, j) == 0) continue;

          edges_(i, j).status = 1;
          edges_(i, j).shared_info->setUndirected();
        }
      }
      return true;
    }
    if (iter_indices.empty()) return false;  // no cycle
  }
  // Never reached since iter_indices.size() < iterations_.size() by definition
  return false;
}

vector<vector<int>> CycleTracker::getAdjMatrices(int size) {
  vector<vector<int>> adj_matrices;
  for (auto i = iterations_.begin(); i != iterations_.begin() + size; ++i) {
    adj_matrices.push_back(i->adj_matrix_1d);
  }
  return adj_matrices;
}

vector<vector<double>> CycleTracker::getProbaAdjMatrices(int size) {
  vector<vector<double>> adj_matrices;
  for (auto i = iterations_.begin(); i != iterations_.begin() + size; ++i) {
    adj_matrices.push_back(i->proba_adj_matrix_1d);
  }
  return adj_matrices;
}

}  // namespace detail
}  // namespace reconstruction
}  // namespace miic
