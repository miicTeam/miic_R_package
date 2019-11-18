#ifndef MIIC_RECONSTRUCT_H_
#define MIIC_RECONSTRUCT_H_

#include <cmath>
#include <deque>
#include <map>
#include <set>
#include <stack>

#include "structure.h"

namespace miic {
namespace reconstruction {

namespace reconstruction_impl {
using std::set;
using std::vector;
using uint = unsigned int;
using edge_index_1d = uint;

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
        is_cp(env.numNodes, 0),
        degree_of(env.numNodes, 0),
        bc_tree_rep(env.numNodes, -1),
        bcc_set_indices(env.numNodes, set<int>()) {}

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
  // Max number of iterations to track in a cycle, cycle of larger size won't
  // be recognized.
  static constexpr uint max_cycle_size = 100;

  struct Iteration {
    uint index;
    // key: index of edge
    // value: status of edge in the previous iteration
    std::map<uint, int> changed_edges;

    Iteration(const Environment& env, uint i) : index(i) {
      // Keep track of the lower triangular part
      for (uint i = 1; i < env.numNodes; ++i) {
        for (uint j = 0; j < i; ++j) {
          const auto& edge = env.edges[i][j];
          if (edge.status_prev == edge.status) continue;

          auto index_1d = getEdgeIndex1D(i, j);
          changed_edges.insert(std::make_pair(index_1d, edge.status_prev));
        }
      }
    }
  };

  // Internal class keeping track of iterations, max size is given
  // by max_cycle_size, new iteration is inserted at the front
  class IterationList {
    std::deque<Iteration> iteration_list_;

   public:
    void add(Iteration&& i) {
      iteration_list_.push_front(std::move(i));
      if (iteration_list_.size() > max_cycle_size) iteration_list_.pop_back();
    }
    Iteration& get(uint i) { return iteration_list_[i]; }

    size_t size() { return iteration_list_.size(); }

    auto begin() { return iteration_list_.begin(); }
    auto cbegin() const { return iteration_list_.cbegin(); }
    auto end() { return iteration_list_.end(); }
    auto cend() const { return iteration_list_.cend(); }
  };

  Environment& env_;
  IterationList iterations_;
  uint n_saved = 0;  // Number of saving operations performed
  // key: number of edges in the graph
  // value: index of iteration
  std::multimap<uint, uint> edge_index_map_;

  void saveIteration() {
    uint n_edge = env_.numNoMore;
    // Index of the iteration starting from 0
    uint index = n_saved++;
    edge_index_map_.insert(std::make_pair(n_edge, index));
    // skip the first iteration as its previous step is the initial graph
    if (index != 0) iterations_.add(Iteration(env_, index));
  }

 public:
  CycleTracker(Environment& env) : env_(env) {}
  // convert lower triangular indices to 1d index
  static edge_index_1d getEdgeIndex1D(uint i, uint j) {
    return (j < i ? j + i * (i - 1) / 2 : i + j * (j - 1) / 2);
  }
  // convert 1d index to lower triangular indices
  std::pair<uint, uint> getEdgeIndex2D(uint k) {
    // floor is equivalent to a int cast for positive number
    uint i = (uint)(0.5 + std::sqrt(0.25 + 2 * k));
    uint j = k - i * (i - 1) / 2;
    return std::make_pair(i, j);
  }
  // check if a cycle exists between the current and the past iterations
  bool hasCycle();
};
}  // namespace reconstruction_impl
using reconstruction_impl::BCC;
using reconstruction_impl::CycleTracker;

bool skeletonInitialization(structure::Environment&);
bool firstStepIteration(structure::Environment&, BCC&);
bool skeletonIteration(structure::Environment&);
bool reconstruct(structure::Environment&, std::string, double, int, char*[]);
std::vector<std::vector<std::string> > confidenceCut(structure::Environment&);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_RECONSTRUCT_H_
