#ifndef MIIC_CYCLE_TRACKER_H_
#define MIIC_CYCLE_TRACKER_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <deque>
#include <map>
#include <utility>  // std::pair
#include <vector>

#include "structure.h"

namespace miic {
namespace reconstruction {

namespace detail {
using std::vector;
using namespace structure;

// During each consistent iteration of the network reconstruction, keep track of
// number of edges in the graph, edges with modified status with respect to the
// previous iteration, and the corresponding edge status in the previous
// iteration.
class CycleTracker {
 public:
  CycleTracker(Grid2d<Edge>& edges, const vector<EdgeID>& edge_list)
      : edges_(edges), edge_list_(edge_list) {}
  // convert lower triangular indices to 1d index
  static int getEdgeIndex1D(int i, int j) {
    return (j < i ? j + i * (i - 1) / 2 : i + j * (j - 1) / 2);
  }
  // convert 1d index to lower triangular indices
  std::pair<int, int> getEdgeIndex2D(int k) {
    int i = std::floor(0.5 + std::sqrt(0.25 + 2 * k));
    int j = k - i * (i - 1) / 2;
    return std::make_pair(i, j);
  }
  vector<vector<int>> getAdjMatrices(int size);
  int getCycleSize() { return cycle_size; }
  // check if a cycle exists between the current and the past iterations
  bool hasCycle();

 private:
  struct Iteration {
    int index;
    // key: index of edge
    // value: status of edge in the previous iteration
    std::map<int, int> changed_edges;
    vector<int> adj_matrix_1d;

    Iteration(const Grid2d<Edge>& edges, int i)
        : index(i), adj_matrix_1d(edges.size(), 0) {
      int n_nodes = edges.n_rows();
      for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
          adj_matrix_1d[j + i * n_nodes] = edges(i, j).status;
        }
      }
      // Keep track of the lower triangular part
      for (int i = 1; i < n_nodes; ++i) {
        for (int j = 0; j < i; ++j) {
          const auto& edge = edges(i, j);
          if (edge.status_prev == edge.status) continue;

          auto index_1d = getEdgeIndex1D(i, j);
          changed_edges.insert(std::make_pair(index_1d, edge.status_prev));
        }
      }
    }
  };

  // Internal class keeping track of iterations.
  class IterationList {
   public:
    template <class... Args>
    void add(Args&&... args) {
      iteration_list_.emplace_front(std::forward<Args>(args)...);
    }
    Iteration& get(int i) { return iteration_list_[i]; }

    size_t size() { return iteration_list_.size(); }

    auto begin() { return iteration_list_.begin(); }
    auto cbegin() const { return iteration_list_.cbegin(); }
    auto end() { return iteration_list_.end(); }
    auto cend() const { return iteration_list_.cend(); }

   private:
    std::deque<Iteration> iteration_list_;
  };

  Grid2d<Edge>& edges_;
  const vector<EdgeID>& edge_list_;
  IterationList iterations_;
  int n_saved{0};  // Number of saving operations performed
  int cycle_size{0};
  // key: number of edges in the graph
  // value: index of iteration
  std::multimap<int, int> edge_index_map_;

  void saveIteration() {
    int n_edge = edge_list_.size();
    // Index of the iteration starting from 0
    int index = n_saved++;
    edge_index_map_.insert(std::make_pair(n_edge, index));
    // skip the first iteration as its previous step is the initial graph
    if (index != 0) iterations_.add(edges_, index);
  }
};

}  // namespace detail
using detail::CycleTracker;
}  // namespace reconstruction
}  // namespace miic
#endif
