#ifndef MIIC_BICONNECTED_COMPONENT_H_
#define MIIC_BICONNECTED_COMPONENT_H_

#include <set>
#include <stack>
#include <utility>  // std::pair
#include <vector>

#include "environment.h"
#include "structure.h"

namespace miic {
namespace reconstruction {

namespace detail {
using std::pair;
using std::set;
using std::stack;
using std::vector;
using namespace miic::structure;
using namespace miic::utility;

class BiconnectedComponent {
 public:
  BiconnectedComponent(Environment& env)
      : environment(env),
        consistent(env.consistent),
        latent(env.latent),
        is_cut_point(env.n_nodes, 0),
        degree_of(env.n_nodes, 0),
        bc_tree_rep(env.n_nodes, -1),
        bcc_set_indices(env.n_nodes, set<int>()) {}

  void analyse() {
    // Reset members
    std::fill(begin(is_cut_point), end(is_cut_point), 0);
    std::fill(begin(degree_of), end(degree_of), 0);
    std::fill(begin(bc_tree_rep), end(bc_tree_rep), -1);
    bcc_list.clear();
    for (auto& s : bcc_set_indices)
      s.clear();
    // Set biconnected components
    bcc();
  }

  // Find and set all candidate Z for a given pair of vertices
  // using biconnected components and block-cut tree.
  set<int> getCandidateZ(int x, int y) const;
  void setCandidateZ(int x, int y, vector<int>& zi_list);
  // For each node z in vect_z, check if
  // (1) z lies on a path between node x and node y, and
  // (2) z is a non-child of either x or y.
  bool isConsistent(int x, int y, const vector<int>& vect_z) const;

 private:
  Environment& environment;
  const bool consistent;
  const bool latent;
  vector<int> is_cut_point;
  vector<int> degree_of;
  vector<int> bc_tree_rep;
  vector<int> bc_tree_inverse_index;
  vector<int> bc_tree_node_is_cp;
  vector<set<int>> bcc_list;
  vector<set<int>> bcc_set_indices;
  vector<set<int>> bc_tree_adj_list;

  void bccAux(int u, int& time, TempVector<int>& parent,
      TempVector<int>& lowest, TempVector<int>& depth,
      stack<pair<int, int>>& st);
  void bcc();
  TempVector<int> bcTreeBfs(int start, int end) const;
};
}  // namespace detail
using detail::BiconnectedComponent;
}  // namespace reconstruction
}  // namespace miic

#endif
