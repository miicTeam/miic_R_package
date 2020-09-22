#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>

#include <string>
#include <vector>

#include "biconnected_component.h"
#include "confidence_cut.h"
#include "cycle_tracker.h"
#include "environment.h"
#include "orientation.h"
#include "r_cpp_interface.h"
#include "skeleton.h"
#include "utilities.h"
#include "tmiic.h"

using Rcpp::_;
using Rcpp::as;
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
  // Initialize Environment with mandatory inputs
  Environment environment(as<int>(arg_list["n_samples"]),
      as<int>(arg_list["n_nodes"]), as<vector<int>>(input_data["factor"]),
      as<vector<int>>(input_data["order"]),
      as<vector<int>>(arg_list["is_continuous"]),
      as<vector<int>>(arg_list["levels"]));

  // Set optional parameters
  setEnvironmentFromR(input_data, arg_list, environment);

  size_t li_alloc_size = getLinearAllocatorSize(environment.n_samples,
      environment.n_nodes, environment.maxbins, environment.initbins,
      environment.is_continuous, environment.levels);
#ifdef _OPENMP
#pragma omp parallel  // each thread has its own instance of li_alloc_ptr
#endif
  li_alloc_ptr = std::make_unique<LinearAllocator>(li_alloc_size);

  // Start reconstruction
  auto lap_start = getLapStartTime();
  Rcout << "Search all pairs for unconditional independence relations...\n";
  // Initialize skeleton, find unconditional independence
  if (!initializeSkeleton(environment)) return empty_results();
  environment.exec_time.init += getLapInterval(lap_start);

  BiconnectedComponent bcc(
      environment.edges, environment.consistent, environment.latent);
  CycleTracker cycle_tracker(environment.edges, environment.connected_list);
  vector<vector<string>> orientations;
  int iter_count{0};
  bool is_consistent{false};
  do {
    if (environment.consistent != 0) {
      // In temporal mode, duplicate temporarily edges over history for the
      // consistency assessment
      if (environment.tau_max >= 1)
        tmiic::repeatEdgesOverHistory (environment);
      bcc.analyse();
      if (environment.tau_max >= 1)
        tmiic::dropPastEdges (environment);
    }
    // Store current status in status_prev and revert to the structure at the
    // moment of initialization
    for (int i = 0; i < environment.n_nodes; i++) {
      for (int j = 0; j < environment.n_nodes; j++) {
        auto& edge = environment.edges(i, j);
        edge.status_prev = edge.status;
        edge.status = edge.status_init;
        if (edge.status != 0) edge.proba_head = 0.5;
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
    if (environment.orientation && !environment.connected_list.empty()) {
      lap_start = getLapStartTime();
      Rcout << "Search for edge directions...\n";
      //
      // In temporal mode, when latent variable discovery is activated,
      // we temporarily duplicate edges over history assuming stationarity to
      // increase the number of possible unshielded triples for orientation
      //
      if ( (environment.tau_max >= 1) && (environment.latent_orientation) )
        tmiic::repeatEdgesOverHistory (environment);
      orientations = orientationProbability(environment);
      if ( (environment.tau_max >= 1) && (environment.latent_orientation) )
          tmiic::dropPastEdges (environment);
      environment.exec_time.ori += getLapInterval(lap_start);
    }
    if (environment.consistent != 0)
      Rcout << "Iteration " << iter_count << ' ';
    Rcout << "Number of edges: " << environment.connected_list.size() << '\n';
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

  const auto& time = environment.exec_time;
  List result = List::create(
      _["adj_matrix"]        = getAdjMatrix(environment.edges),
      _["proba_adj_matrix"]  = getProbaAdjMatrix(environment.edges),
      _["edges"]             = getEdgesInfoTable(environment.edges,
                                   environment.nodes),
      _["orientations.prob"] = orientations,
      _["time"]              = vector<double>{time.init, time.iter, time.cut,
                                   time.ori, time.getTotal()},
      _["interrupted"]       = false);
  if (environment.consistent != 0) {
    int size = is_consistent ? cycle_tracker.getCycleSize()
                             : environment.max_iteration;
    result.push_back(cycle_tracker.getAdjMatrices(size), "adj_matrices");
    result.push_back(
        cycle_tracker.getProbaAdjMatrices(size), "proba_adj_matrices");
    result.push_back(is_consistent, "is_consistent");
  }
  return result;
}
