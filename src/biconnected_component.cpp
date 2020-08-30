#include "biconnected_component.h"

#include <queue>

namespace miic {
namespace reconstruction {
namespace detail {

std::set<int> BiconnectedComponent::getCandidateZ(int x, int y) const {
  set<int> set_z;
  auto insert_it = inserter(set_z, set_z.begin());

  if (degree_of[x] < 1 || degree_of[y] < 1) return set_z;

  vector<int> common_bcc;
  set_intersection(bcc_set_indices[x].begin(), bcc_set_indices[x].end(),
      bcc_set_indices[y].begin(), bcc_set_indices[y].end(),
      inserter(common_bcc, common_bcc.begin()));
  if (common_bcc.empty()) {
    int start = bc_tree_rep[x];
    int end = bc_tree_rep[y];
    vector<int> bc_tree_path = bcTreeBfs(start, end);
    for (auto& node : bc_tree_path) {
      if (bc_tree_node_is_cp[node]) continue;

      auto& bcc_set = bcc_list[bc_tree_inverse_index[node]];
      copy_if(bcc_set.cbegin(), bcc_set.cend(), insert_it,
          [x, y](int i) { return i != x && i != y; });
    }
    return set_z;
  } else {
    auto& bcc_set = bcc_list[common_bcc[0]];
    copy_if(bcc_set.cbegin(), bcc_set.cend(), insert_it,
        [x, y](int i) { return i != x && i != y; });
    return set_z;
  }
}

void BiconnectedComponent::setCandidateZ(int x, int y) {
  auto& vect_z = environment.edges[x][y].shared_info->zi_list;
  auto set_z = getCandidateZ(x, y);
  copy_if(
      set_z.begin(), set_z.end(), back_inserter(vect_z), [this, x, y](int i) {
        return (environment.latent || environment.edges[i][x].status_prev > 0 ||
                environment.edges[i][y].status_prev > 0);
      });
}

bool BiconnectedComponent::isConsistent(
    int x, int y, const vector<int>& vect_z) const {
  if (vect_z.empty()) return true;
  std::set<int> set_z = getCandidateZ(x, y);
  for (auto& z : vect_z) {
    if (set_z.find(z) == set_z.end() ||
        (environment.edges[z][x].status <= 0 &&
            environment.edges[z][y].status <= 0))
      return false;
  }
  return true;
}

// Biconnected components decomposition of the graph contained in the
// environment, allowing for search of candidate vertices for separation.
void BiconnectedComponent::bcc() {
  int time = 0;
  int n_nodes = environment.n_nodes;
  vector<int> depth(n_nodes, -1), lowest(n_nodes, -1), parent(n_nodes, -1);
  std::stack<std::pair<int, int> > st;

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
      bcc_list.push_back(s);
    }
  }

  int bc_tree_size =
      count(is_cut_point.begin(), is_cut_point.end(), 1) + bcc_list.size();
  bc_tree_adj_list.assign(bc_tree_size, set<int>());
  bc_tree_inverse_index.assign(bc_tree_size, -1);
  bc_tree_node_is_cp.assign(bc_tree_size, 0);

  int bc_tree_index = 0;
  for (size_t index = 0; index < bcc_list.size(); index++) {
    int rep = bc_tree_index++;
    bc_tree_inverse_index[rep] = index;

    for (auto& node : bcc_list[index]) {
      bcc_set_indices[node].insert(index);

      if (is_cut_point[node]) {
        if (bc_tree_rep[node] == -1) {
          bc_tree_rep[node] = bc_tree_index;
          bc_tree_node_is_cp[bc_tree_index] = 1;
          bc_tree_inverse_index[bc_tree_index] = node;
          bc_tree_index++;
        }

        bc_tree_adj_list[bc_tree_rep[node]].insert(rep);
        bc_tree_adj_list[rep].insert(bc_tree_rep[node]);
      } else
        bc_tree_rep[node] = rep;
    }
  }

  for (int i = 0; i < n_nodes; i++) {
    for (int j = 0; j < n_nodes; j++) {
      if (i == j) continue;
      degree_of[i] += environment.edges[i][j].shared_info->connected;
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
void BiconnectedComponent::bccAux(int u, int& time, vector<int>& parent,
    vector<int>& lowest, vector<int>& depth,
    std::stack<std::pair<int, int> >& st) {
  int n_nodes = environment.n_nodes;
  int children = 0;
  depth[u] = lowest[u] = ++time;

  for (int v = 0; v < n_nodes; v++) {
    // graph maybe (partially) directed, whereas biconnected component
    // concerns only the skeleton
    if (!environment.edges[u][v].status && !environment.edges[v][u].status)
      continue;

    if (depth[v] == -1) {
      parent[v] = u;
      children++;
      st.push(std::make_pair(u, v));

      bccAux(v, time, parent, lowest, depth, st);

      lowest[u] = std::min(lowest[u], lowest[v]);
      if ((parent[u] == -1 && children > 1) ||
          (parent[u] != -1 && lowest[v] >= depth[u])) {
        is_cut_point[u] = 1;
        set<int> s;
        int i = -1, j = -1;
        do {
          i = st.top().first;
          j = st.top().second;
          s.insert(i);
          s.insert(j);
          st.pop();
        } while (i != u || j != v);
        bcc_list.push_back(s);
      }
    } else if (v != parent[u] && depth[u] > depth[v]) {
      lowest[u] = std::min(lowest[u], depth[v]);
      st.push(std::make_pair(u, v));
    }
  }
}

// Find the shortest path between two nodes in the block-cut tree.
vector<int> BiconnectedComponent::bcTreeBfs(int start, int end) const {
  int n_nodes = bc_tree_adj_list.size();
  vector<int> visited(n_nodes, 0);

  std::queue<std::pair<int, vector<int> > > bfs_queue;
  bfs_queue.push(std::make_pair(start, vector<int>{start}));

  while (!bfs_queue.empty()) {
    auto& p = bfs_queue.front();
    visited[p.first] = 1;
    for (auto& i : bc_tree_adj_list[p.first]) {
      if (!visited[i]) {
        vector<int> new_path(p.second);
        new_path.push_back(i);
        if (i == end)
          return new_path;
        else
          bfs_queue.push(make_pair(i, new_path));
      }
    }
    bfs_queue.pop();
  }
  return vector<int>();
}

}  // namespace detail
}  // namespace reconstruction
}  // namespace miic
