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

using uint = unsigned int;
using std::cout;
using std::endl;
using std::vector;
using namespace miic::computation;
using namespace miic::reconstruction;
using namespace miic::structure;
using namespace miic::utility;

bool SortFunctionNoMore(
    const EdgeID* a, const EdgeID* b, const Environment& environment) {
  return environment.edges[a->i][a->j].shared_info->Ixy_ui >
         environment.edges[b->i][b->j].shared_info->Ixy_ui;
}

class sorterNoMore {
  Environment& environment;

 public:
  sorterNoMore(Environment& env) : environment(env) {}
  bool operator()(EdgeID const* o1, EdgeID const* o2) const {
    return SortFunctionNoMore(o1, o2, environment);
  }
};

bool SortFunction1(
    const EdgeID* a, const EdgeID* b, const Environment& environment) {
  return environment.edges[a->i][a->j].shared_info->Rxyz_ui >
         environment.edges[b->i][b->j].shared_info->Rxyz_ui;
}

class sorter1 {
  Environment& environment;

 public:
  sorter1(Environment& env) : environment(env) {}
  bool operator()(EdgeID const* o1, EdgeID const* o2) const {
    return SortFunction1(o1, o2, environment);
  }
};

template <bool isLatent = false>
void searchAndSetZi(
    Environment& environment, const uint posX, const uint posY) {
  // Search candidate nodes for the separation set of X and Y
  int numZiPos = 0;
  for (uint c = 0; c < environment.numNodes; c++) {
    if (c == posX || c == posY) continue;
    if (!isLatent && !environment.edges[posX][c].status_prev &&
        !environment.edges[posY][c].status_prev)
      continue;
    environment.edges[posX][posY].shared_info->zi_vect_idx.push_back(c);
    numZiPos++;
  }

  if (environment.isVerbose)
    cout << "The number of neighbours is: " << numZiPos << endl;
}

void searchAndSetZi(
    Environment& environment, const uint posX, const uint posY) {
  if (environment.isLatent)
    return searchAndSetZi<true>(environment, posX, posY);
  else
    return searchAndSetZi<false>(environment, posX, posY);
}

namespace miic {
namespace reconstruction {

bool firstStepIteration(Environment& environment, BCC& bcc) {
  // During first step iteration, search for U contributors is not
  // parallelizable see flag "parallelizable" in
  // computeEnsInformationContinuous() l. 1118 in computeEnsInformation.cpp
  environment.firstIterationDone = false;

  for (unsigned i = 0; i < environment.searchMoreAddress.size(); i++)
    delete environment.searchMoreAddress[i];

  for (unsigned i = 0; i < environment.noMoreAddress.size(); i++)
    delete environment.noMoreAddress[i];

  environment.noMoreAddress.clear();
  environment.searchMoreAddress.clear();

  // create and fill the searchMoreAddress struct, that keep track of i and j
  // positions of searchMore Edges
  environment.numSearchMore = 0;
  environment.numNoMore = 0;
  for (uint i = 0; i < environment.numNodes - 1; i++) {
    for (uint j = i + 1; j < environment.numNodes; j++) {
      // Do dot consider edges removed with unconditional independence
      if (!environment.edges[i][j].status) continue;
      environment.edges[i][j].shared_info->reset();
      environment.searchMoreAddress.emplace_back(new EdgeID(i, j));
      environment.numSearchMore++;
    }
  }

  int threadnum = 0;
  bool interrupt = false;
  int prg_numSearchMore = -1;
  cout << "First round of conditional independences :\n";
  if (environment.numSearchMore > 0) {
    if (environment.isVerbose == true)
      cout << "\n# -> searchMore edges, to get zi and noMore...\n";

    for (int i = 0; i < environment.numSearchMore; i++) {
      int posX = environment.searchMoreAddress[i]->i;
      int posY = environment.searchMoreAddress[i]->j;
      if (environment.isVerbose) {
        cout << "\n# --------------------\n# ----> EDGE: "
             << environment.nodes[posX].name << "--"
             << environment.nodes[posY].name << "\n# --------------------";
      }
      if (environment.consistentPhase)
        bcc.set_candidate_z(posX, posY);
      else
        searchAndSetZi(environment, posX, posY);
    }

    if (environment.isVerbose) {
      cout << "SEARCH OF BEST Z: ";
    }

    environment.execTime.startTimeInit = get_wall_time();
#ifdef _OPENMP
#pragma omp parallel for shared(interrupt) firstprivate(threadnum) \
    schedule(dynamic)
#endif
    for (int i = 0; i < environment.numSearchMore; i++) {
      if (interrupt) {
        continue;  // will continue until out of for loop
      }
#ifdef _OPENMP
      threadnum = omp_get_thread_num();
#endif
      if (threadnum == 0) {
        if (checkInterrupt(i / environment.nThreads % 2 == 0)) {
          interrupt = true;
        }
      }
      int posX = environment.searchMoreAddress[i]->i;
      int posY = environment.searchMoreAddress[i]->j;
      if (environment.isVerbose)
        cout << "##  "
             << "XY: " << environment.nodes[posX].name << " "
             << environment.nodes[posY].name << "\n\n";
      if (environment.edges[posX][posY].shared_info->zi_vect_idx.size() > 0) {
        // Search for new contributing node and its rank
        SearchForNewContributingNodeAndItsRank(
            environment, posX, posY, environment.memoryThreads[threadnum]);
      }
      // Dynamic thread allocation makes it so we can't know the end point of
      // thread 0, on average it will be numSearchMore - nThreads/2
      if (threadnum == 0)
        prg_numSearchMore = printProgress(
            1.0 * i / (environment.numSearchMore - environment.nThreads / 2),
            environment.execTime.startTimeInit, prg_numSearchMore);
    }
    // Print finished progress bar
    prg_numSearchMore = printProgress(
        1.0, environment.execTime.startTimeInit, prg_numSearchMore);

    if (interrupt) return false;

    for (int i = 0; i < environment.numSearchMore; i++) {
      int posX = environment.searchMoreAddress[i]->i;
      int posY = environment.searchMoreAddress[i]->j;
      if (environment.edges[posX][posY].shared_info->z_name_idx != -1) {
        if (environment.isVerbose) {
          cout << "## ------!!--> Update the edge element in 'searchMore': "
               << environment
                      .nodes[environment.edges[posX][posY]
                                 .shared_info
                                 ->zi_vect_idx[environment.edges[posX][posY]
                                                   .shared_info->z_name_idx]]
                      .name
               << " is a good zi candidate\n";
        }
      } else {
        if (environment.isVerbose) {
          cout << "## ------!!--> Remove the edge element from searchMore.\n## "
                  "------!!--> Add edge to 'noMore' (no good zi candidate)\n";
        }
        // Put this edge element to "noMore"
        environment.noMoreAddress.push_back(environment.searchMoreAddress[i]);
        environment.numNoMore++;
        // Remove the element from searchMore
        environment.searchMoreAddress.erase(
            environment.searchMoreAddress.begin() + i);
        environment.numSearchMore--;
        i--;
        // Update the status
        environment.edges[posX][posY].shared_info->connected = 1;
      }
      if (environment.isVerbose) cout << "\n";
    }

    if (checkInterrupt()) {
      return false;
    }
    // sort the ranks
    std::sort(environment.searchMoreAddress.begin(),
        environment.searchMoreAddress.end(), sorter1(environment));
  }
  environment.firstIterationDone = true;
  return (true);
}

bool skeletonIteration(Environment& environment) {
  int iIteration_count = 0;
  int max = 0;

  if (environment.isVerbose)
    cout << "Number of numSearchMore: " << environment.numSearchMore << endl;

  cout << "\nSkeleton iteration :\n";
  environment.execTime.startTimeIter = get_wall_time();
  int start_numSearchMore = environment.numSearchMore;

  int prg_numSearchMore = -1;

  while (environment.numSearchMore > 0) {
    if (checkInterrupt()) {
      return (false);
    }
    iIteration_count++;
    if (environment.isVerbose) {
      cout << "\n# Iteration " << iIteration_count << "\n";
    }
    // Get the first edge
    int posX = environment.searchMoreAddress[max]->i;
    int posY = environment.searchMoreAddress[max]->j;

    if (environment.isVerbose)
      cout << "Pos x : " << posX << " , pos y: " << posY << endl;

    auto topEdgeElt = environment.edges[posX][posY].shared_info;

    if (environment.isVerbose) cout << "# Before adding new zi to {ui}: ";

    // Keep the previous z.name for this edge
    std::string accepted_z_name =
        environment.nodes[topEdgeElt->z_name_idx].name;

    // Reinit ui.vect, z.name, zi.vect, z.name.idx
    if (environment.isVerbose) {
      cout << "# DO: Add new zi to {ui}: " << topEdgeElt->z_name_idx << endl;
    }
    // move top z_name_idx from zi_vect to ui_vect
    topEdgeElt->ui_vect_idx.push_back(
        topEdgeElt->zi_vect_idx[topEdgeElt->z_name_idx]);
    topEdgeElt->zi_vect_idx.erase(
        topEdgeElt->zi_vect_idx.begin() + topEdgeElt->z_name_idx);
    topEdgeElt->z_name_idx = -1;

    double* v = NULL;
    if (!environment.is_continuous[posX] && !environment.is_continuous[posY] &&
        std::all_of(topEdgeElt->ui_vect_idx.cbegin(),
            topEdgeElt->ui_vect_idx.cend(),
            [&environment](int i) { return !environment.is_continuous[i]; })) {
      v = computeEnsInformationNew(environment,
          &environment.edges[posX][posY].shared_info->ui_vect_idx[0],
          environment.edges[posX][posY].shared_info->ui_vect_idx.size(), NULL,
          0, -1, posX, posY, environment.cplx, environment.m);

      topEdgeElt->Ixy_ui = v[1];
      topEdgeElt->Nxy_ui = v[0];
      topEdgeElt->cplx = v[2];
      free(v);
    } else {
      v = computeEnsInformationContinuous(environment,
          &environment.edges[posX][posY].shared_info->ui_vect_idx[0],
          environment.edges[posX][posY].shared_info->ui_vect_idx.size(), NULL,
          0, -1, posX, posY, environment.cplx, environment.m);
      topEdgeElt->Nxy_ui = v[0];
      topEdgeElt->Ixy_ui = v[1];
      topEdgeElt->cplx = v[2];

      delete[] v;
    }

    double topEdgeElt_kxy_ui = topEdgeElt->cplx;

    if (environment.isDegeneracy)
      topEdgeElt_kxy_ui =
          topEdgeElt->cplx + (topEdgeElt->ui_vect_idx.size() * log(3));

    int nRemainingEdges = environment.numSearchMore + environment.numNoMore;

    if (environment.isVerbose) {
      cout << "# --> nbrEdges L = " << nRemainingEdges << "\n";
      cout << "# --> nbrProp P = " << environment.numNodes << "\n\n";
      cout << "topEdgeElt->Ixy_ui " << topEdgeElt->Ixy_ui << "\n";
      cout << "topEdgeElt_kxy_ui " << topEdgeElt_kxy_ui << "\n";
      cout << "environment.logEta " << environment.logEta << "\n";
      cout << "IsPhantom? "
           << (topEdgeElt->Ixy_ui - topEdgeElt_kxy_ui - environment.logEta <= 0)
           << endl;
    }
    if (topEdgeElt->Ixy_ui - topEdgeElt_kxy_ui - environment.logEta <= 0) {
      if (environment.isVerbose) {
        cout << "# PHANTOM" << environment.nodes[posX].name << ","
             << environment.nodes[posY].name << "\n";
      }

      // Move this edge from the list searchMore to phantom
      environment.searchMoreAddress.erase(
          environment.searchMoreAddress.begin() + max);
      environment.numSearchMore--;

      // Set the connection to 0 on the adj matrix
      environment.edges[posX][posY].status = 0;
      environment.edges[posY][posX].status = 0;
      // Save the phantom status
      topEdgeElt->connected = 0;
    } else {
      // Reinit Rxyz_ui
      topEdgeElt->Rxyz_ui = environment.thresPc;

      if (environment.isVerbose) {
        cout << "# Do SearchForNewContributingNodeAndItsRank\n";
      }

      if (topEdgeElt->zi_vect_idx.size() > 0) {
        SearchForNewContributingNodeAndItsRank(
            environment, posX, posY, environment.m);
      }

      if (environment.isVerbose) {
        if (environment.edges[posX][posY].shared_info->z_name_idx == -1)
          cout << "# See topEdgeElt[['z.name']]: NA\n";
        else
          cout << "# See topEdgeElt[['z.name']]: "
               << environment
                      .nodes[topEdgeElt->zi_vect_idx[topEdgeElt->z_name_idx]]
                      .name
               << "\n";
      }
      //// Update the information about the edge
      if (topEdgeElt->z_name_idx != -1) {
        if (environment.isVerbose) {
          cout << "# Do update myAllEdges$searchMore\n";
        }
        // myGv$allEdges[["searchMore"]][[topEdgeElt[["key"]]]] = topEdgeElt

      } else {
        if (environment.isVerbose) {
          cout << "# Do update myAllEdges$noMore\n";
        }
        //// Move this edge from the list searchMore to noMore
        environment.noMoreAddress.push_back(environment.searchMoreAddress[max]);
        environment.numNoMore++;
        environment.searchMoreAddress.erase(
            environment.searchMoreAddress.begin() + max);
        environment.numSearchMore--;
        //// Update the status of the edge
        topEdgeElt->connected = 1;
      }
    }

    //// Sort all pairs xy with a contributing node z in decreasing order of
    /// their ranks, R(xy;z| )
    if (environment.isVerbose) {
      cout << "# Do Sort all pairs by Rxyz_ui\n";
    }

    max = 0;
    for (int i = 0; i < environment.numSearchMore; i++) {
      if (environment
              .edges[environment.searchMoreAddress[i]->i]
                    [environment.searchMoreAddress[i]->j]
              .shared_info->Rxyz_ui >
          environment
              .edges[environment.searchMoreAddress[max]->i]
                    [environment.searchMoreAddress[max]->j]
              .shared_info->Rxyz_ui)
        max = i;
    }
    // cout << 1.0*(start_numSearchMore -
    // environment.numSearchMore)/(start_numSearchMore-1) << "\t" <<
    // prg_numSearchMore << "\n" << flush;
    prg_numSearchMore =
        printProgress(1.0 * (start_numSearchMore - environment.numSearchMore) /
                          (start_numSearchMore),
            environment.execTime.startTimeIter, prg_numSearchMore);
  }
  cout << "\n";
  std::sort(environment.noMoreAddress.begin(), environment.noMoreAddress.end(),
      sorterNoMore(environment));
  return (true);
}

bool BCC::is_consistent(int x, int y, const vector<int>& vect_z) const {
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
  int numNodes = environment.numNodes;
  int children = 0;
  depth[u] = lowest[u] = ++time;

  for (int v = 0; v < numNodes; v++) {
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
  int numNodes = environment.numNodes;
  vector<int> depth(numNodes, -1), lowest(numNodes, -1), parent(numNodes, -1);
  std::stack<std::pair<int, int> > st;

  for (int u = 0; u < numNodes; u++) {
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

  for (int i = 0; i < numNodes; i++) {
    for (int j = 0; j < numNodes; j++) {
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
  int numNodes = bc_tree_adj_list.size();
  vector<int> visited(numNodes, 0);

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
  std::insert_iterator<set<int> > insert_it = inserter(set_z, set_z.begin());

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

      const set<int>& bcc_set = bcc_list[bc_tree_inverse_index[node]];
      copy_if(bcc_set.begin(), bcc_set.end(), insert_it,
          [this, x, y](int i) { return i != x && i != y; });
    }
    return set_z;
  } else {
    const set<int>& bcc_set = bcc_list[common_bcc[0]];
    copy_if(bcc_set.begin(), bcc_set.end(), insert_it,
        [this, x, y](int i) { return i != x && i != y; });
    return set_z;
  }
}

void BCC::set_candidate_z(int x, int y) {
  vector<int>& vect_z = environment.edges[x][y].shared_info->zi_vect_idx;
  std::insert_iterator<vector<int> > insert_it =
      inserter(vect_z, vect_z.begin());
  set<int> set_z = get_candidate_z(x, y);
  copy_if(set_z.begin(), set_z.end(), insert_it, [this, x, y](int i) {
    return (environment.isLatent || environment.edges[i][x].status_prev > 0 ||
            environment.edges[i][y].status_prev > 0);
  });
}

}  // namespace reconstruction
}  // namespace miic
