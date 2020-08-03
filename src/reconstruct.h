#ifndef MIIC_RECONSTRUCT_H_
#define MIIC_RECONSTRUCT_H_

#include <cmath>
#include <deque>
#include <map>
#include <set>
#include <stack>

#include "environment.h"
#include "structure.h"

namespace miic {
namespace reconstruction {

namespace reconstruction_impl {
using std::set;
using std::vector;

using structure::Environment;

// Class for biconnected component analysis
class BCC {
  Environment& environment;
  vector<int> is_cp;
  vector<int> degree_of;
  vector<int> bc_tree_rep;
  vector<int> bc_tree_inverse_index;
  vector<int> bc_tree_node_is_cp;
  vector<set<int> > bcc_list;
  vector<set<int> > bcc_set_indices;
  vector<set<int> > bc_tree_adj_list;

  void bcc_aux(int u, int& time, vector<int>& parent, vector<int>& lowest,
      vector<int>& depth, std::stack<std::pair<int, int> >& st);
  void bcc();
  vector<int> bc_tree_bfs(int start, int end) const;

 public:
  BCC(Environment& env)
      : environment(env),
        is_cp(env.n_nodes, 0),
        degree_of(env.n_nodes, 0),
        bc_tree_rep(env.n_nodes, -1),
        bcc_set_indices(env.n_nodes, set<int>()) {}

  void analyse() {
    // Reset members
    std::fill(is_cp.begin(), is_cp.end(), 0);
    std::fill(degree_of.begin(), degree_of.end(), 0);
    std::fill(bc_tree_rep.begin(), bc_tree_rep.end(), -1);
    bcc_list.clear();
    for (auto& s : bcc_set_indices) s.clear();
    // Set biconnected components
    bcc();
  }

  std::set<int> get_candidate_z(int x, int y) const;
  void set_candidate_z(int x, int y);
  bool is_consistent(int x, int y, const vector<int>& vect_z) const;
};

// During each consistent iteration of the network reconstruction, keep track of
// number of edges in the graph, edges with modified status with respect to the
// previous iteration, and the corresponding edge status in the previous
// iteration.
class CycleTracker {
 private:
  struct Iteration {
    int index;
    // key: index of edge
    // value: status of edge in the previous iteration
    std::map<int, int> changed_edges;
    vector<int> adj_matrix_1d;

    Iteration(const Environment& env, int i)
        : index(i), adj_matrix_1d(env.n_nodes * env.n_nodes, 0) {
      int n_node(env.n_nodes);
      for (int i = 0; i < n_node; ++i) {
        for (int j = 0; j < n_node; ++j) {
          adj_matrix_1d[j + i * n_node] = env.edges[i][j].status;
        }
      }
      // Keep track of the lower triangular part
      for (int i = 1; i < n_node; ++i) {
        for (int j = 0; j < i; ++j) {
          const auto& edge = env.edges[i][j];
          if (edge.status_prev == edge.status) continue;

          auto index_1d = getEdgeIndex1D(i, j);
          changed_edges.insert(std::make_pair(index_1d, edge.status_prev));
        }
      }
    }
  };

  // Internal class keeping track of iterations.
  class IterationList {
    std::deque<Iteration> iteration_list_;

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
  };

  Environment& env_;
  IterationList iterations_;
  int n_saved = 0;  // Number of saving operations performed
  // key: number of edges in the graph
  // value: index of iteration
  std::multimap<int, int> edge_index_map_;

  void saveIteration() {
    int n_edge = env_.numNoMore;
    // Index of the iteration starting from 0
    int index = n_saved++;
    edge_index_map_.insert(std::make_pair(n_edge, index));
    // skip the first iteration as its previous step is the initial graph
    if (index != 0) iterations_.add(env_, index);
  }

 public:
  CycleTracker(Environment& env) : env_(env) {}
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
  // check if a cycle exists between the current and the past iterations
  bool hasCycle();
  // when a cycle is found, hold the adj_matrix of all iterations in the cycle
  vector<vector<int>> adj_matrices;
};
}  // namespace reconstruction_impl
using reconstruction_impl::BCC;
using reconstruction_impl::CycleTracker;

bool skeletonInitialization(structure::Environment&);
bool firstStepIteration(structure::Environment&, BCC&);
bool skeletonIteration(structure::Environment&);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_RECONSTRUCT_H_
