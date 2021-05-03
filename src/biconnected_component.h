#ifndef MIIC_BICONNECTED_COMPONENT_H_
#define MIIC_BICONNECTED_COMPONENT_H_

#include <set>
#include <stack>
#include <utility>  // std::pair
#include <vector>

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
  BiconnectedComponent(const Grid2d<Edge>& edges, int consistent, bool latent)
      : edges_(edges),
        n_nodes_(edges.n_rows()),
        consistent_(consistent),
        latent_(latent),
        is_cut_point_(n_nodes_, 0),
        degree_of_(n_nodes_, 0),
        bc_tree_rep_(n_nodes_, -1),
        bcc_set_indices_(n_nodes_, set<int>()) {}

  void analyse() {
    // Reset members
    std::fill(begin(is_cut_point_), end(is_cut_point_), 0);
    std::fill(begin(degree_of_), end(degree_of_), 0);
    std::fill(begin(bc_tree_rep_), end(bc_tree_rep_), -1);
    bcc_list_.clear();
    for (auto& s : bcc_set_indices_)
      s.clear();
    // Set biconnected components
    bcc();
  }

  // Find and set all candidate Z for a given pair of vertices,
  // use biconnected components and block-cut tree if consistent is true
  void setCandidateZ(int x, int y, vector<int>& zi_list);
  // For each node z in vect_z, check if
  // (1) z lies on a path between node x and node y, and
  // (2) z is a non-child of either x or y.
  bool isConsistent(int x, int y, const vector<int>& vect_z) const;

 private:
  const Grid2d<Edge>& edges_;
  const int n_nodes_;
  const int consistent_;
  const bool latent_;
  vector<int> is_cut_point_;
  vector<int> degree_of_;
  vector<int> bc_tree_rep_;
  vector<int> bc_tree_inverse_index_;
  vector<int> bc_tree_node_is_cp_;
  vector<set<int>> bcc_list_;
  vector<set<int>> bcc_set_indices_;
  vector<set<int>> bc_tree_adj_list_;

  void bccAux(int u, int& time, TempVector<int>& parent,
      TempVector<int>& lowest, TempVector<int>& depth,
      stack<pair<int, int>>& st);
  void bcc();
  TempVector<int> bcTreeBfs(int start, int end) const;
  set<int> getCandidateZ(int x, int y) const;
};

}  // namespace detail
using detail::BiconnectedComponent;
}  // namespace reconstruction
}  // namespace miic

#endif
