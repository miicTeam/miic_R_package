#ifndef MIIC_CYCLE_TRACKER_H_
#define MIIC_CYCLE_TRACKER_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <deque>
#include <map>
#include <utility>  // std::pair
#include <vector>

#include "structure.h"
#include "utilities.h"

namespace miic {
namespace reconstruction {

namespace detail {
using std::vector;
using namespace structure;
using namespace utility;

// convert 2d index in a n_nodes * n_nodes grid to 1d index
inline int getIndex1D(int i, int j, int n_nodes) { return j + i * n_nodes; }

// convert 1d index to 2d index in a n_nodes * n_nodes grid
inline std::pair<int, int> getIndex2D(int k, int n_nodes) {
  return std::make_pair(k / n_nodes, k % n_nodes);
}

// During each consistent iteration of the network reconstruction, keep track of
// number of edges in the graph, edges with modified status with respect to the
// previous iteration, and the corresponding edge status in the previous
// iteration.
class CycleTracker {
 public:
  CycleTracker(Grid2d<Edge>& edges, const vector<EdgeID>& edge_list)
      : edges_(edges), edge_list_(edge_list) {}
  vector<vector<int>> getAdjMatrices(int size);
  vector<vector<double>> getProbaAdjMatrices(int size);
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
    vector<double> proba_adj_matrix_1d;

    Iteration(const Grid2d<Edge>& edges, int i)
        : index(i),
          adj_matrix_1d(getAdjMatrix(edges)),
          proba_adj_matrix_1d(getProbaAdjMatrix(edges)) {
      int n_nodes = edges.n_rows();
      for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
          const auto& edge = edges(i, j);
          if (edge.status_prev == edge.status) continue;

          auto index_1d = getIndex1D(i, j, n_nodes);
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
    auto begin() const { return iteration_list_.cbegin(); }
    auto end() { return iteration_list_.end(); }
    auto end() const { return iteration_list_.cend(); }

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
