#include "biconnected_component.h"

#include <algorithm>  // std::copy_if, std::count
#include <queue>

#include "linear_allocator.h"

namespace miic {
namespace reconstruction {
namespace detail {

using std::copy_if;
using std::make_pair;
using std::min;

set<int> BiconnectedComponent::getCandidateZ(int x, int y) const {
  TempAllocatorScope scope;

  set<int> set_z;
  auto insert_it = inserter(set_z, set_z.begin());

  if (degree_of_[x] < 1 || degree_of_[y] < 1) return set_z;

  TempVector<int> common_bcc;
  set_intersection(bcc_set_indices_[x].begin(), bcc_set_indices_[x].end(),
      bcc_set_indices_[y].begin(), bcc_set_indices_[y].end(),
      inserter(common_bcc, common_bcc.begin()));
  if (common_bcc.empty()) {
    int start = bc_tree_rep_[x];
    int end = bc_tree_rep_[y];
    TempVector<int> bc_tree_path = bcTreeBfs(start, end);
    for (auto& node : bc_tree_path) {
      if (bc_tree_node_is_cp_[node]) continue;

      auto& bcc_set = bcc_list_[bc_tree_inverse_index_[node]];
      copy_if(bcc_set.cbegin(), bcc_set.cend(), insert_it,
          [x, y](int i) { return i != x && i != y; });
    }
    return set_z;
  } else {
    auto& bcc_set = bcc_list_[common_bcc[0]];
    copy_if(bcc_set.cbegin(), bcc_set.cend(), insert_it,
        [x, y](int i) { return i != x && i != y; });
    return set_z;
  }
}

void BiconnectedComponent::setCandidateZ(int x, int y, vector<int>& zi_list) {
  zi_list.clear();
  if (consistent_ != 0) {
    auto set_z = getCandidateZ(x, y);
    auto is_consistent = [this, x, y](int z){
      // consistent_ -> 1: orientation consistency; 2: skeleton consistency
      if (latent_ || consistent_ == 2) return true;
      // For double arrow headed edge (x <-> z), z is considered consistent
      if ((edges_(x, z).status_prev == 2 && edges_(z, x).status_prev == 2) ||
          (edges_(y, z).status_prev == 2 && edges_(z, y).status_prev == 2))
        return true;
      // status is either 0 (not connected) or 2 (z is the child of x or y)
      if (edges_(x, z).status_prev != 1 && edges_(y, z).status_prev != 1)
        return false;
      return true;
    };
    copy_if(begin(set_z), end(set_z), back_inserter(zi_list), is_consistent);
  } else {
    for (int z = 0; z < n_nodes_; ++z) {
      if (z == x || z == y) continue;
      if (latent_ || edges_(x, z).status_prev ||
          edges_(y, z).status_prev) {
        zi_list.push_back(z);
      }
    }
  }
}

bool BiconnectedComponent::isConsistent(
    int x, int y, const vector<int>& vect_z) const {
  if (vect_z.empty()) return true;
  set<int> set_z = getCandidateZ(x, y);
  for (auto& z : vect_z) {
    // Not in the consistent set
    if (set_z.find(z) == set_z.end())
      return false;
    // consistent_ -> 1: orientation consistency; 2: skeleton consistency
    if (consistent_ == 2)
      continue;
    // For double arrow headed edge (x <-> z), z is considered consistent
    if ((edges_(x, z).status == 2 && edges_(z, x).status == 2) ||
        (edges_(y, z).status == 2 && edges_(z, y).status == 2))
      continue;
    // status is either 0 (not connected) or 2 (z is the child of x or y)
    if (edges_(x, z).status != 1 && edges_(y, z).status != 1)
      return false;
  }
  return true;
}

// Biconnected components decomposition of the graph contained in the
// environment, allowing for search of candidate vertices for separation.
void BiconnectedComponent::bcc() {
  TempAllocatorScope scope;

  int time = 0;
  int n_nodes = n_nodes_;
  TempVector<int> depth(n_nodes, -1), lowest(n_nodes, -1), parent(n_nodes, -1);
  stack<pair<int, int>> st;

  for (int u = 0; u < n_nodes; u++) {
    if (depth[u] == -1) bccAux(u, time, parent, lowest, depth, st);

    if (!st.empty()) {
      set<int> s;
      int i = -1, j = -1;
      do {
        i = st.top().first;
        j = st.top().second;
        s.insert(i);
        s.insert(j);
        st.pop();
      } while (!st.empty());
      bcc_list_.push_back(s);
    }
  }

  int bc_tree_size = std::count(begin(is_cut_point_), end(is_cut_point_), 1) +
                     bcc_list_.size();
  bc_tree_adj_list_.assign(bc_tree_size, set<int>());
  bc_tree_inverse_index_.assign(bc_tree_size, -1);
  bc_tree_node_is_cp_.assign(bc_tree_size, 0);

  int bc_tree_index = 0;
  for (size_t index = 0; index < bcc_list_.size(); index++) {
    int rep = bc_tree_index++;
    bc_tree_inverse_index_[rep] = index;

    for (auto& node : bcc_list_[index]) {
      bcc_set_indices_[node].insert(index);

      if (is_cut_point_[node]) {
        if (bc_tree_rep_[node] == -1) {
          bc_tree_rep_[node] = bc_tree_index;
          bc_tree_node_is_cp_[bc_tree_index] = 1;
          bc_tree_inverse_index_[bc_tree_index] = node;
          bc_tree_index++;
        }

        bc_tree_adj_list_[bc_tree_rep_[node]].insert(rep);
        bc_tree_adj_list_[rep].insert(bc_tree_rep_[node]);
      } else
        bc_tree_rep_[node] = rep;
    }
  }

  for (int i = 0; i < n_nodes; i++) {
    for (int j = 0; j < n_nodes; j++) {
      if (i == j) continue;
      degree_of_[i] += edges_(i, j).shared_info->connected;
    }
  }
}

// Auxiliary recurrent method for biconnected component decomposition.
//
// @param int u Current node under consideration.
// @param int& time global time used to set the depth of each vertex.
// @param vector<int> parent Parent vertex of each vertex in the dfs search.
// @param vector<int> lowest Lowest point of each vertex.
// @param vector<int> depth Time when each vertex is visited in the dfs
// search.
// @param stack<pair<int>> st Stack for the dfs search.
void BiconnectedComponent::bccAux(int u, int& time, TempVector<int>& parent,
    TempVector<int>& lowest, TempVector<int>& depth,
    stack<pair<int, int>>& st) {
  int children = 0;
  depth[u] = lowest[u] = ++time;

  for (int v = 0; v < n_nodes_; v++) {
    // graph maybe (partially) directed, whereas biconnected component
    // concerns only the skeleton
    if (!edges_(u, v).status && !edges_(v, u).status)
      continue;

    if (depth[v] == -1) {
      parent[v] = u;
      children++;
      st.push(make_pair(u, v));

      bccAux(v, time, parent, lowest, depth, st);

      lowest[u] = min(lowest[u], lowest[v]);
      if ((parent[u] == -1 && children > 1) ||
          (parent[u] != -1 && lowest[v] >= depth[u])) {
        is_cut_point_[u] = 1;
        set<int> s;
        int i = -1, j = -1;
        do {
          i = st.top().first;
          j = st.top().second;
          s.insert(i);
          s.insert(j);
          st.pop();
        } while (i != u || j != v);
        bcc_list_.push_back(s);
      }
    } else if (v != parent[u] && depth[u] > depth[v]) {
      lowest[u] = min(lowest[u], depth[v]);
      st.push(make_pair(u, v));
    }
  }
}

// Find the shortest path between two nodes in the block-cut tree.
TempVector<int> BiconnectedComponent::bcTreeBfs(int start, int end) const {
  TempAllocatorScope scope;

  int n_nodes = bc_tree_adj_list_.size();
  TempVector<int> visited(n_nodes, 0);

  std::queue<pair<int, TempVector<int>>> bfs_queue;
  bfs_queue.push(make_pair(start, TempVector<int>{start}));

  while (!bfs_queue.empty()) {
    auto& p = bfs_queue.front();
    visited[p.first] = 1;
    for (auto& i : bc_tree_adj_list_[p.first]) {
      if (!visited[i]) {
        TempVector<int> new_path(p.second);
        new_path.push_back(i);
        if (i == end)
          return new_path;
        else
          bfs_queue.push(make_pair(i, new_path));
      }
    }
    bfs_queue.pop();
  }
  return TempVector<int>();
}

}  // namespace detail
}  // namespace reconstruction
}  // namespace miic
