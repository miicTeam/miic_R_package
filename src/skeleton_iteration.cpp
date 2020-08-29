#include <Rcpp.h>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "reconstruct.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "compute_ens_information.h"
#include "structure.h"
#include "utilities.h"

using Rcpp::Rcout;
using Rcpp::Rcerr;
using std::endl;
using std::vector;
using namespace miic::computation;
using namespace miic::reconstruction;
using namespace miic::structure;
using namespace miic::utility;

template <bool latent = false>
void searchAndSetZi(
    Environment& environment, const int posX, const int posY) {
  // Search candidate nodes for the separation set of X and Y
  int numZiPos = 0;
  for (int c = 0; c < environment.n_nodes; c++) {
    if (c == posX || c == posY) continue;
    if (!latent && !environment.edges[posX][c].status_prev &&
        !environment.edges[posY][c].status_prev)
      continue;
    environment.edges[posX][posY].shared_info->zi_list.push_back(c);
    numZiPos++;
  }

  if (environment.verbose)
    Rcout << "The number of neighbours is: " << numZiPos << endl;
}

void searchAndSetZi(
    Environment& environment, const int posX, const int posY) {
  if (environment.latent)
    return searchAndSetZi<true>(environment, posX, posY);
  else
    return searchAndSetZi<false>(environment, posX, posY);
}

namespace miic {
namespace reconstruction {

bool firstStepIteration(Environment& environment, BCC& bcc) {
  // During first step iteration, search for U contributors is not
  // parallelizable see flag "parallelizable" in
  // computeEnsInformationContinuous() in computeEnsInformation.cpp
  environment.first_iter_done = false;

  environment.connected_list.clear();
  environment.unsettled_list.clear();

  // create and fill the unsettled_list struct, that keep track of i and j
  // positions of searchMore Edges
  environment.numSearchMore = 0;
  environment.numNoMore = 0;
  for (int i = 0; i < environment.n_nodes - 1; i++) {
    for (int j = i + 1; j < environment.n_nodes; j++) {
      // Do dot consider edges removed with unconditional independence
      if (!environment.edges[i][j].status) continue;
      environment.edges[i][j].shared_info->reset();
      environment.unsettled_list.emplace_back(i, j, environment.edges[i][j]);
      environment.numSearchMore++;
    }
  }

  if (environment.numSearchMore > 0) {
    if (environment.verbose)
      Rcout << "\n# -> searchMore edges, to get zi and noMore...\n";

    for (int i = 0; i < environment.numSearchMore; i++) {
      int posX = environment.unsettled_list[i].i;
      int posY = environment.unsettled_list[i].j;
      if (environment.verbose) {
        Rcout << "\n# --------------------\n# ----> EDGE: "
             << environment.nodes[posX].name << "--"
             << environment.nodes[posY].name << "\n# --------------------";
      }
      if (environment.consistent > 0)
        bcc.set_candidate_z(posX, posY);
      else
        searchAndSetZi(environment, posX, posY);
    }

    if (environment.verbose) {
      Rcout << "SEARCH OF BEST Z: ";
    }

    bool interrupt = false;
    int progress_percentile = -1;
    int n_jobs_done{0};
    auto loop_start_time = getLapStartTime();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < environment.numSearchMore; i++) {
      if (interrupt) {
        continue;  // will continue until out of for loop
      }
      int threadnum = 0;
#ifdef _OPENMP
      threadnum = omp_get_thread_num();
#endif
      if (threadnum == 0) {
        if (checkInterrupt()) interrupt = true;
      }
      const auto& edgeid = environment.unsettled_list[i];
      if (edgeid.getEdge().shared_info->zi_list.size() > 0)
        SearchForNewContributingNodeAndItsRank(environment, edgeid.i, edgeid.j);

#ifdef _OPENMP
#pragma omp atomic
      ++n_jobs_done;
#endif
      if (threadnum == 0)
        printProgress(
            static_cast<double>(n_jobs_done) / environment.numSearchMore,
            loop_start_time, progress_percentile);
    }
    // Print finished progress bar
    printProgress(1.0, loop_start_time, progress_percentile);
    Rcerr << '\n';

    if (interrupt) return false;

    for (int i = 0; i < environment.numSearchMore; i++) {
      int posX = environment.unsettled_list[i].i;
      int posY = environment.unsettled_list[i].j;
      if (environment.edges[posX][posY].shared_info->z_name_idx != -1) {
        if (environment.verbose) {
          Rcout << "## ------!!--> Update the edge element in 'searchMore': "
               << environment
                      .nodes[environment.edges[posX][posY].shared_info->zi_list
                                 [environment.edges[posX][posY]
                                         .shared_info->z_name_idx]]
                      .name
               << " is a good zi candidate\n";
        }
      } else {
        if (environment.verbose) {
          Rcout << "## ------!!--> Remove the edge element from searchMore.\n## "
                  "------!!--> Add edge to 'noMore' (no good zi candidate)\n";
        }
        // Move this edge element to "noMore"
        environment.connected_list.push_back(environment.unsettled_list[i]);
        environment.numNoMore++;
        environment.unsettled_list.erase(
            environment.unsettled_list.begin() + i);
        environment.numSearchMore--;
        i--;
        // Update the status
        environment.edges[posX][posY].shared_info->connected = 1;
      }
      if (environment.verbose) Rcout << "\n";
    }

    if (checkInterrupt()) {
      return false;
    }
    // sort the ranks
    std::sort(
        environment.unsettled_list.begin(), environment.unsettled_list.end());
  }
  environment.first_iter_done = true;
  return true;
}

bool skeletonIteration(Environment& environment) {
  int iIteration_count = 0;
  int max = 0;

  if (environment.verbose)
    Rcout << "Number of numSearchMore: " << environment.numSearchMore << endl;

  auto loop_start_time = getLapStartTime();
  int start_numSearchMore = environment.numSearchMore;

  int progress_percentile = -1;

  while (environment.numSearchMore > 0) {
    if (checkInterrupt()) {
      return false;
    }
    iIteration_count++;
    if (environment.verbose) {
      Rcout << "\n# Iteration " << iIteration_count << "\n";
    }
    // Get the first edge
    int posX = environment.unsettled_list[max].i;
    int posY = environment.unsettled_list[max].j;

    if (environment.verbose)
      Rcout << "Pos x : " << posX << " , pos y: " << posY << endl;

    auto topEdgeElt = environment.edges[posX][posY].shared_info;

    if (environment.verbose) Rcout << "# Before adding new zi to {ui}: ";

    // Keep the previous z.name for this edge
    std::string accepted_z_name =
        environment.nodes[topEdgeElt->z_name_idx].name;

    // Reinit ui.vect, z.name, zi.vect, z.name.idx
    if (environment.verbose) {
      Rcout << "# DO: Add new zi to {ui}: " << topEdgeElt->z_name_idx << endl;
    }
    // move top z_name_idx from zi_vect to ui_vect
    topEdgeElt->ui_list.push_back(topEdgeElt->zi_list[topEdgeElt->z_name_idx]);
    topEdgeElt->zi_list.erase(
        topEdgeElt->zi_list.begin() + topEdgeElt->z_name_idx);
    topEdgeElt->z_name_idx = -1;

    double* v = NULL;
    if (!environment.is_continuous[posX] && !environment.is_continuous[posY] &&
        std::all_of(topEdgeElt->ui_list.cbegin(), topEdgeElt->ui_list.cend(),
            [&environment](int i) { return !environment.is_continuous[i]; })) {
      v = computeEnsInformationNew(environment, posX, posY,
          environment.edges[posX][posY].shared_info->ui_list, vector<int>(),
          environment.cplx);

      topEdgeElt->Ixy_ui = v[1];
      topEdgeElt->Nxy_ui = v[0];
      topEdgeElt->cplx = v[2];
    } else {
      v = computeEnsInformationContinuous(environment, posX, posY,
          environment.edges[posX][posY].shared_info->ui_list, vector<int>(),
          environment.cplx);
      topEdgeElt->Nxy_ui = v[0];
      topEdgeElt->Ixy_ui = v[1];
      topEdgeElt->cplx = v[2];
    }
    delete[] v;
    double topEdgeElt_kxy_ui = topEdgeElt->cplx;

    if (environment.degenerate)
      topEdgeElt_kxy_ui =
          topEdgeElt->cplx + (topEdgeElt->ui_list.size() * log(3.0));

    int nRemainingEdges = environment.numSearchMore + environment.numNoMore;

    if (environment.verbose) {
      Rcout << "# --> nbrEdges L = " << nRemainingEdges << "\n";
      Rcout << "# --> nbrProp P = " << environment.n_nodes << "\n\n";
      Rcout << "topEdgeElt->Ixy_ui " << topEdgeElt->Ixy_ui << "\n";
      Rcout << "topEdgeElt_kxy_ui " << topEdgeElt_kxy_ui << "\n";
      Rcout << "environment.log_eta " << environment.log_eta << "\n";
      Rcout << "IsPhantom? "
           << (topEdgeElt->Ixy_ui - topEdgeElt_kxy_ui - environment.log_eta <= 0)
           << endl;
    }
    if (topEdgeElt->Ixy_ui - topEdgeElt_kxy_ui - environment.log_eta <= 0) {
      // Conditional independence found, remove edge
      if (environment.verbose) {
        Rcout << "# PHANTOM" << environment.nodes[posX].name << ","
             << environment.nodes[posY].name << "\n";
      }
      environment.unsettled_list.erase(
          environment.unsettled_list.begin() + max);
      environment.numSearchMore--;
      // Set the connection to 0 on the adj matrix
      environment.edges[posX][posY].status = 0;
      environment.edges[posY][posX].status = 0;
      // Save the phantom status
      topEdgeElt->connected = 0;
    } else {
      // Reinit Rxyz_ui
      topEdgeElt->Rxyz_ui = environment.thresPc;

      if (environment.verbose) {
        Rcout << "# Do SearchForNewContributingNodeAndItsRank\n";
      }

      if (topEdgeElt->zi_list.size() > 0) {
        SearchForNewContributingNodeAndItsRank(environment, posX, posY);
      }

      if (environment.verbose) {
        if (environment.edges[posX][posY].shared_info->z_name_idx == -1)
          Rcout << "# See topEdgeElt[['z.name']]: NA\n";
        else
          Rcout << "# See topEdgeElt[['z.name']]: "
               << environment.nodes[topEdgeElt->zi_list[topEdgeElt->z_name_idx]]
                      .name
               << "\n";
      }
      //// Update the information about the edge
      if (topEdgeElt->z_name_idx != -1) {
        if (environment.verbose) {
          Rcout << "# Do update myAllEdges$searchMore\n";
        }
        // myGv$allEdges[["searchMore"]][[topEdgeElt[["key"]]]] = topEdgeElt

      } else {
        if (environment.verbose) {
          Rcout << "# Do update myAllEdges$noMore\n";
        }
        // Move this edge from the list searchMore to noMore
        environment.connected_list.push_back(environment.unsettled_list[max]);
        environment.numNoMore++;
        environment.unsettled_list.erase(
            environment.unsettled_list.begin() + max);
        environment.numSearchMore--;
        // Update the status of the edge
        topEdgeElt->connected = 1;
      }
    }

    // Sort all pairs xy with a contributing node z in decreasing order of
    // their ranks, R(xy;z| )
    if (environment.verbose) {
      Rcout << "# Do Sort all pairs by Rxyz_ui\n";
    }

    max = 0;
    for (int i = 0; i < environment.numSearchMore; i++) {
      if (environment
              .edges[environment.unsettled_list[i].i]
                    [environment.unsettled_list[i].j]
              .shared_info->Rxyz_ui >
          environment
              .edges[environment.unsettled_list[max].i]
                    [environment.unsettled_list[max].j]
              .shared_info->Rxyz_ui)
        max = i;
    }
    printProgress(1.0 * (start_numSearchMore - environment.numSearchMore) /
                      (start_numSearchMore),
        loop_start_time, progress_percentile);
  }
  Rcerr << "\n";
  std::sort(
      environment.connected_list.begin(), environment.connected_list.end());
  return (true);
}

bool BCC::isConsistent(int x, int y, const vector<int>& vect_z) const {
  // For each node z in vect_z, check if
  // (1) z lies on a path between node x and node y, and
  // (2) z is a non-child of either x or y.
  if (vect_z.empty()) return true;
  std::set<int> set_z = get_candidate_z(x, y);
  for (auto& z : vect_z) {
    if (set_z.find(z) == set_z.end() ||
        (environment.edges[z][x].status <= 0 &&
            environment.edges[z][y].status <= 0))
      return false;
  }
  return true;
}

void BCC::bcc_aux(int u, int& time, vector<int>& parent, vector<int>& lowest,
    vector<int>& depth, std::stack<std::pair<int, int> >& st) {
  // Auxiliary recurrent method for biconnected component decomposition.
  //
  // @param int u Current node under consideration.
  // @param int& time global time used to set the depth of each vertex.
  // @param vector<int> parent Parent vertex of each vertex in the dfs search.
  // @param vector<int> lowest Lowest point of each vertex.
  // @param vector<int> depth Time when each vertex is visited in the dfs
  // search.
  // @param stack<pair<int> > st Stack for the dfs search.
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

      bcc_aux(v, time, parent, lowest, depth, st);

      lowest[u] = std::min(lowest[u], lowest[v]);
      if ((parent[u] == -1 && children > 1) ||
          (parent[u] != -1 && lowest[v] >= depth[u])) {
        is_cp[u] = 1;
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

void BCC::bcc() {
  // Biconnected components decomposition of the graph contained in the
  // environment, allowing for search of candidate vertices for separation.
  int time = 0;
  int n_nodes = environment.n_nodes;
  vector<int> depth(n_nodes, -1), lowest(n_nodes, -1), parent(n_nodes, -1);
  std::stack<std::pair<int, int> > st;

  for (int u = 0; u < n_nodes; u++) {
    if (depth[u] == -1) bcc_aux(u, time, parent, lowest, depth, st);

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

  int bc_tree_size = count(is_cp.begin(), is_cp.end(), 1) + bcc_list.size();
  bc_tree_adj_list.assign(bc_tree_size, set<int>());
  bc_tree_inverse_index.assign(bc_tree_size, -1);
  bc_tree_node_is_cp.assign(bc_tree_size, 0);

  int bc_tree_index = 0;
  for (size_t index = 0; index < bcc_list.size(); index++) {
    int rep = bc_tree_index++;
    bc_tree_inverse_index[rep] = index;

    for (auto& node : bcc_list[index]) {
      bcc_set_indices[node].insert(index);

      if (is_cp[node]) {
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

vector<int> BCC::bc_tree_bfs(int start, int end) const {
  // Return the shortest path between two nodes in the block-cut tree.
  //
  // @param int start, end Starting, ending nodes.
  // @return Path as a vector of node indices.
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

std::set<int> BCC::get_candidate_z(int x, int y) const {
  // Find and set all candidate Z for a given pair of vertices
  // using biconnected components and block-cut tree.
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
    vector<int> bc_tree_path = bc_tree_bfs(start, end);
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

void BCC::set_candidate_z(int x, int y) {
  auto& vect_z = environment.edges[x][y].shared_info->zi_list;
  auto set_z = get_candidate_z(x, y);
  copy_if(set_z.begin(), set_z.end(), back_inserter(vect_z),
      [this, x, y](int i) {
        return (environment.latent || environment.edges[i][x].status_prev > 0 ||
                environment.edges[i][y].status_prev > 0);
      });
}

}  // namespace reconstruction
}  // namespace miic
