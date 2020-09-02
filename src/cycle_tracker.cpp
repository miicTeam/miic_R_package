#include "cycle_tracker.h"

#include <algorithm>  // std::sort, std::none_of

namespace miic {
namespace reconstruction {
namespace detail {

bool CycleTracker::hasCycle() {
  int n_edge = env_.numNoMore;
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
  std::sort(iter_indices.begin(), iter_indices.end(), std::greater<int>());
  // Set of edges that are to be marked as connected and undirected
  std::set<int> edges_union;
  // Check if an edge is changed. vector is chosen over map for quicker access
  // and simpler syntax, at the cost of extra memory trace and (possible) extra
  // time complexity (in practice there are very few changes between each pair
  // of iterations).
  vector<int> changed(env_.n_nodes * (env_.n_nodes - 1) / 2, 0);
  // backtracking over iteration to get changed_edges
  cycle_size = 0;
  for (const auto& iter : iterations_) {
    ++cycle_size;
    for (const auto& k : iter.changed_edges) {
      edges_union.insert(k.first);
      // compare the status in the previous iteration against the latest status
      std::pair<int, int> p = getEdgeIndex2D(k.first);
      changed[k.first] = (k.second != env_.edges[p.first][p.second].status);
    }
    if (iter.index != iter_indices.front()) continue;
    iter_indices.pop_front();
    using std::none_of;
    if (none_of(begin(changed), end(changed), [](int j) { return j != 0; })) {
      for (auto& k : edges_union) {
        std::pair<int, int> p = getEdgeIndex2D(k);
        env_.edges[p.first][p.second].status = 1;
        env_.edges[p.second][p.first].status = 1;
        env_.edges[p.first][p.second].shared_info->setUndirected();
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
  for (auto i = iterations_.begin(), e = iterations_.begin() + size; i != e;
       ++i) {
    adj_matrices.push_back(i->adj_matrix_1d);
  }
  return adj_matrices;
}

}  // namespace detail
}  // namespace reconstruction
}  // namespace miic
