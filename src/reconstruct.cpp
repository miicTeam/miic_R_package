#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>

#include <string>
#include <vector>

#include "biconnected_component.h"
#include "confidence_cut.h"
#include "cycle_tracker.h"
#include "orientation.h"
#include "skeleton.h"
#include "utilities.h"

using Rcpp::_;
using Rcpp::List;
using Rcpp::Rcout;
using std::string;
using std::vector;
using namespace miic::reconstruction;
using namespace miic::structure;
using namespace miic::utility;

List empty_results() { return List::create(_["interrupted"] = true); }

// [[Rcpp::export]]
List reconstruct(List input_data, List arg_list) {
  Environment environment(input_data, arg_list);

  int max_level =
      *std::max_element(begin(environment.levels), end(environment.levels));
  size_t li_alloc_size = 8192;  // in bytes, extra space for fragmentation
  size_t n_integers = 4 * (environment.n_samples + 2) + 6 * (max_level + 1) +
                      3 * ((environment.n_samples + 1) * 7) +
                      4 * environment.n_nodes +
                      2 * (max_level + 1) * (max_level + 1) +
                      2 * max_level * max_level * max_level;  // coarse_counts
  size_t n_doubles = 3 * (max_level + 1) + 2 * environment.n_nodes;
  li_alloc_size += sizeof(int) * n_integers;
  li_alloc_size += sizeof(double) * n_doubles;
#ifdef _OPENMP
#pragma omp parallel  // each thread has its own instance of li_alloc_ptr
#endif
  li_alloc_ptr = std::make_unique<LinearAllocator>(li_alloc_size);

  auto lap_start = getLapStartTime();
  Rcout << "Search all pairs for unconditional independence relations...\n";
  // Initialize skeleton, find unconditional independence
  if (!initializeSkeleton(environment)) return empty_results();
  environment.exec_time.init += getLapInterval(lap_start);

  BiconnectedComponent bcc(environment);
  CycleTracker cycle_tracker(environment);
  vector<vector<string>> orientations;
  int iter_count{0};
  bool is_consistent{false};
  do {
    if (environment.consistent != 0) bcc.analyse();
    // Store current status in status_prev and revert to the structure at the
    // moment of initialization
    for (int i = 0; i < environment.n_nodes; i++) {
      for (int j = 0; j < environment.n_nodes; j++) {
        environment.edges[i][j].status_prev = environment.edges[i][j].status;
        environment.edges[i][j].status = environment.edges[i][j].status_init;
      }
    }
    lap_start = getLapStartTime();
    Rcout << "Search for candidate separating nodes...\n";
    // If interrupted
    if (!setBestContributingNode(environment, bcc)) return empty_results();

    if (!environment.unsettled_list.empty()) {
      Rcout << "Search for conditional independence relations...\n";
      // If interrupted
      if (!searchForConditionalIndependence(environment))
        return (empty_results());
    }
    environment.exec_time.iter += getLapInterval(lap_start);

    if (environment.n_shuffles > 0) {
      lap_start = getLapStartTime();
      Rcout << "Compute confidence cut with permutations..." << std::flush;
      setConfidence(environment);
      size_t n_connected = environment.connected_list.size();
      confidenceCut(environment);
      Rcout << n_connected - environment.connected_list.size()
            << " edges cut.\n";
      environment.exec_time.cut += getLapInterval(lap_start);
    }
    // Oriente edges for non-consistent/orientation consistent algorithm
    if (environment.orientation_phase && !environment.connected_list.empty() &&
        environment.consistent <= 1) {
      lap_start = getLapStartTime();
      Rcout << "Search for edge directions...\n";
      orientations = orientationProbability(environment);
      environment.exec_time.ori += getLapInterval(lap_start);
    }
    if (environment.consistent != 0)
      Rcout << "Iteration " << iter_count << ' ';
    Rcout << "Number of edges: " << environment.numNoMore << '\n';
    is_consistent = cycle_tracker.hasCycle();
    if (is_consistent) {
      Rcout << "cycle found of size " << cycle_tracker.getCycleSize() << '\n';
      break;
    }
    if (++iter_count > environment.max_iteration) {
      Rcout << "Iteration limit " << environment.max_iteration << " reached.\n";
      break;
    }
  } while (environment.consistent != 0);

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
    lap_start = getLapStartTime();
    orientations = orientationProbability(environment);
    environment.exec_time.ori += getLapInterval(lap_start);
    // Check inconsistency after orientation, add undirected edge to
    // pairs with inconsistent conditional independence.
    bcc.analyse();
    int n_inconsistency = 0;
    vector<std::pair<int, int>> inconsistent_edges;
    for (int i = 1; i < environment.n_nodes; i++) {
      for (int j = 0; j < i; j++) {
        const Edge& edge = environment.edges[i][j];
        if (edge.status || bcc.isConsistent(i, j, edge.shared_info->ui_list))
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
          << " found after orientation.\n";
  }

  const auto& time = environment.exec_time;
  List result = List::create(
      _["adj_matrix"]        = getAdjMatrix(environment),
      _["edges"]             = getEdgesInfoTable(environment),
      _["orientations.prob"] = orientations,
      _["time"]              = vector<double>{
        time.init, time.iter, time.cut, time.ori, time.getTotal()},
      _["interrupted"]       = false);
  if (environment.consistent != 0) {
    int size = is_consistent ? cycle_tracker.getCycleSize()
                             : environment.max_iteration;
    result.push_back(cycle_tracker.getAdjMatrices(size), "adj_matrices");
    result.push_back(is_consistent, "is_consistent");
  }
  return result;
}
