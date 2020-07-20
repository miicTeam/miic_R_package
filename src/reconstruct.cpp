#include "reconstruct.h"

#include <Rcpp.h>
#include <unistd.h>

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "orientation_probability.h"
#include "utilities.h"

using Rcpp::_;
using Rcpp::as;
using Rcpp::DataFrame;
using Rcpp::List;
using Rcpp::Rcout;
using std::string;
using std::vector;
using namespace miic::reconstruction;
using namespace miic::structure;
using namespace miic::utility;

List empty_results() {
  List result;
  result = List::create(_["interrupted"] = true);
  return (result);
}

// [[Rcpp::export]]
List reconstruct(List input_data, List arg_list) {
  Environment environment(input_data, arg_list);

  double startTime = get_wall_time();

  environment.memoryThreads = new MemorySpace[environment.n_threads];
  for (int i = 0; i < environment.n_threads; i++) {
    createMemorySpace(environment, environment.memoryThreads[i]);
  }
  createMemorySpace(environment, environment.m);

  // Initialize skeleton, find unconditional independence
  if (!skeletonInitialization(environment)) {
    List result =
        List::create(
            _["error"]       = "error during skeleton initialization",
            _["interrupted"] = true);
    return result;
  }

  long double spentTime = (get_wall_time() - startTime);
  environment.exec_time.init = spentTime;
  environment.exec_time.init_iter = spentTime;
  environment.exec_time.iter = 0;
  if (environment.verbose) {
    Rcout << "\n# ----> First contributing node elapsed time:" << spentTime
              << "sec\n\n";
  }
  BCC bcc(environment);
  auto cycle_tracker = CycleTracker(environment);
  vector<vector<string>> confVect;
  vector<vector<string>> orientations;
  do {
    if (environment.consistent > 0) bcc.analyse();
    // Save the neighbours in the status_prev structure
    // and revert to the structure at the moment of initialization
    for (int i = 0; i < environment.n_nodes; i++) {
      for (int j = 0; j < environment.n_nodes; j++) {
        environment.edges[i][j].status_prev = environment.edges[i][j].status;
        environment.edges[i][j].status = environment.edges[i][j].status_init;
      }
    }
    // If interrupted
    if (!firstStepIteration(environment, bcc)) return empty_results();

    if (environment.numNoMore == 0 && environment.numSearchMore == 0) {
      if (environment.verbose)
        Rcout << "# ------| Only phantom edges found.\n";
    } else if (environment.numSearchMore > 0) {
      // Search for other Contributing node(s) (possible only for the edges
      /// still in 'searchMore', ie. 2)
      if (environment.verbose) {
        Rcout << "\n# ---- Other Contributing node(s) ----\n\n";
      }
      startTime = get_wall_time();

      // If interrupted
      if (!skeletonIteration(environment)) return (empty_results());

      long double spentTime = (get_wall_time() - startTime);
      environment.exec_time.iter += spentTime;
      environment.exec_time.init_iter += spentTime;
    }

    startTime = get_wall_time();
    if (environment.n_shuffles > 0) {
      Rcout << "Computing confidence cut with permutations..." << std::flush;
      confVect = confidenceCut(environment);
      long double spentTime = (get_wall_time() - startTime);
      environment.exec_time.cut += spentTime;
      Rcout << " done." << std::endl;
    } else {
      environment.exec_time.cut = 0;
    }
    // Oriente edges for non-consistent/orientation consistent algorithm
    if (environment.orientation_phase && environment.numNoMore > 0 &&
        environment.consistent <= 1) {
      orientations = orientationProbability(environment);
    }
    Rcout << "Number of edges: " << environment.numNoMore << std::endl;
  } while (environment.consistent > 0 && !cycle_tracker.hasCycle());

  int union_n_edges = 0;
  for (int i = 1; i < environment.n_nodes; i++) {
    for (int j = 0; j < i; j++) {
      if (environment.edges[i][j].status) {
        union_n_edges++;
      }
    }
  }
  environment.numNoMore = union_n_edges;

  // skeleton consistent algorithm
  if (environment.numNoMore > 0 && environment.consistent == 2) {
    orientations = orientationProbability(environment);
    // Check inconsistency after orientation, add undirected edge to
    // pairs with inconsistent conditional independence.
    bcc.analyse();
    int n_inconsistency = 0;
    std::vector<std::pair<int, int> > inconsistent_edges;
    for (int i = 1; i < environment.n_nodes; i++) {
      for (int j = 0; j < i; j++) {
        const Edge& edge = environment.edges[i][j];
        if (edge.status || bcc.is_consistent(i, j, edge.shared_info->ui_list))
          continue;
        if (environment.verbose) {
          Rcout << environment.nodes[i].name << ",\t"
                    << environment.nodes[j].name << "\t| "
                    << toNameString(environment, edge.shared_info->ui_list)
                    << std::endl;
        }
        inconsistent_edges.emplace_back(i, j);
        ++n_inconsistency;
      }
    }
    for (const auto& k : inconsistent_edges) {
      environment.edges[k.first][k.second].status = 1;
      environment.edges[k.second][k.first].status = 1;
      environment.edges[k.first][k.second].shared_info->setUndirected();
    }
    Rcout << n_inconsistency << " inconsistent conditional independences"
              << " found after orientation." << std::endl;
  }

  vector<double> time;
  time.push_back(environment.exec_time.init);
  time.push_back(environment.exec_time.iter);
  time.push_back(environment.exec_time.cut);
  time.push_back(environment.exec_time.init_iter + environment.exec_time.cut);

  List result;
  result = List::create(
      _["adj_matrix"]        = getAdjMatrix(environment),
      _["edges"]             = getEdgesInfoTable(environment),
      _["orientations.prob"] = orientations,
      _["time"]              = time,
      _["interrupted"]       = false);
  if (environment.n_shuffles > 0) {
    result.push_back(confVect, "confData");
  }
  if (environment.consistent > 0) {
    result.push_back(cycle_tracker.adj_matrices, "adj_matrices");
  }

  for (int i = 0; i < environment.n_threads; i++) {
    deleteMemorySpace(environment, environment.memoryThreads[i]);
  }
  deleteMemorySpace(environment, environment.m);

  return result;
}

namespace miic {
namespace reconstruction {

bool CycleTracker::hasCycle() {
  int n_edge = env_.numNoMore;
  // Before saving the current iteration, search among previous iterations
  // those with the same number of edges
  auto range = edge_index_map_.equal_range(n_edge);
  bool no_cycle_found = range.first == range.second;
  // Indices of iteration that is possibly the end point of a cycle in the
  // backtracking sense. Example: suppose that iteration #1, #3, #6 have the
  // same number of edges as the current iteration #8, and each of them is a
  // possible start point of a cycle, therefore #2, #4, #7 are the corresponding
  // possible end points of the cycle (#8 #7 ... #2), (#8 #7 ... #4), (#8 #7).
  std::deque<int> iter_indices;
  for (auto it = range.first; it != range.second; ++it)
    iter_indices.push_back(it->second + 1);
  saveIteration();
  if (n_saved > env_.max_iteration) {
    Rcout << "Max number of iterations reached: " << env_.max_iteration
              << '\n';
    return true;
  }
  if (no_cycle_found) return false;
  // Backtracking requires starting from the largest index first
  std::sort(iter_indices.begin(), iter_indices.end(), std::greater<int>());
  // Set of edges that are to be marked as connected and undirected
  std::set<int> edges_union;
  // Check if an edge is changed. vector is chosen over map for quicker access
  // and simpler syntax, at the cost of extra memory trace and (possible) extra
  // time complexity (in practice there are very few changes between each pair
  // of iterations).
  vector<int> changed(env_.n_nodes * (env_.n_nodes - 1) / 2, 0);
  // backtracking over iteration to get changed_edges
  int cycle_size = 0;
  for (const auto& iter : iterations_) {
    ++cycle_size;
    for (const auto& k : iter.changed_edges) {
      edges_union.insert(k.first);
      // compare edge status in the previous iteration against the latest edge
      // status
      std::pair<int, int> p = getEdgeIndex2D(k.first);
      changed[k.first] = (k.second != env_.edges[p.first][p.second].status);
    }
    if (iter.index != iter_indices.front()) continue;
    iter_indices.pop_front();
    // if any edge has been changed
    if (std::any_of(
            changed.begin(), changed.end(), [](int j) { return j != 0; })) {
      // no cycle
      if (iter_indices.empty())
        return false;
      else
        continue;
    }
    for (auto& k : edges_union) {
      std::pair<int, int> p = getEdgeIndex2D(k);
      env_.edges[p.first][p.second].status = 1;
      env_.edges[p.second][p.first].status = 1;
      env_.edges[p.first][p.second].shared_info->setUndirected();
    }

    Rcout << "cycle found of size " << cycle_size << std::endl;
    break;
  }
  // fill the adj_matrices
  for (auto i = iterations_.begin(), e = iterations_.begin() + cycle_size;
       i != e; ++i) {
    adj_matrices.push_back(i->adj_matrix_1d);
  }
  return true;
}
}  // namespace reconstruction
}  // namespace miic
