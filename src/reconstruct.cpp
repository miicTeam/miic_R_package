#include "reconstruct.h"

#include <Rcpp.h>
#include <unistd.h>

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "orientation_probability.h"
#include "structure.h"
#include "utilities.h"

using uint = unsigned int;
using Rcpp::_;
using Rcpp::as;
using Rcpp::DataFrame;
using Rcpp::List;
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
List reconstruct(DataFrame input_data, List arg_list) {
  Environment environment;

  environment.data = as<vector<vector<string>>>(input_data);
  auto var_names = as<vector<string>>(arg_list["var_names"]);
  std::transform(var_names.begin(), var_names.end(),
      std::back_inserter(environment.nodes),
      [](string name) { return Node(name); });
  environment.is_continuous = as<vector<int>>(arg_list["is_continuous"]);
  environment.n_nodes = environment.is_continuous.size();
  environment.n_threads = as<int>(arg_list["n_threads"]);
#ifdef _OPENMP
  int threads = environment.n_threads;
  if (threads < 0) threads = omp_get_num_procs();
  omp_set_num_threads(threads);
#endif

  environment.test_mar = as<bool>(arg_list["test_mar"]);

  environment.n_eff = as<int>(arg_list["n_eff"]);
  environment.half_v_structure = as<int>(arg_list["half_v_structure"]);

  string cplx_flag = as<string>(arg_list["cplx"]);
  environment.cplx = 1;  // use nml complexity by default
  if (cplx_flag.compare("mdl") == 0) environment.cplx = 0;

  std::string consistent_flag = as<std::string>(arg_list["consistent"]);
  environment.consistent = 0;
  if (consistent_flag.compare("orientation") == 0)
    environment.consistent = 1;
  if (consistent_flag.compare("skeleton") == 0) environment.consistent = 2;
  environment.max_iteration = as<int>(arg_list["max_iteration"]);

  std::string latent_flag = as<std::string>(arg_list["latent"]);
  environment.latent = false;
  environment.latent_orientation = false;
  if (latent_flag.compare("yes") == 0) environment.latent = true;
  if (latent_flag.compare("orientation") == 0)
    environment.latent_orientation = true;

  environment.is_k23 = as<bool>(arg_list["is_k23"]);
  environment.degenerate = as<bool>(arg_list["degenerate"]);
  bool orientation_phase = as<bool>(arg_list["orientation"]);
  environment.propagation = as<bool>(arg_list["propagation"]);
  environment.no_init_eta = as<bool>(arg_list["no_init_eta"]);

  environment.n_shuffles = as<int>(arg_list["n_shuffles"]);
  environment.conf_threshold = as<double>(arg_list["conf_threshold"]);

  environment.sample_weights = as<vector<double>>(arg_list["sample_weights"]);

  environment.verbose = as<bool>(arg_list["verbose"]);

  double startTime;

  vector<vector<string> > retShuffle;
  vector<vector<string> > retShuffleAvg;
  vector<string> shfVec;

  string filePath;

  // set the environment
  srand(0);
  setEnvironment(environment);
  vector<string> v(as<vector<string>>(arg_list["black_box"]));
  if (v.size() > 1) readBlackbox(v, environment);

  startTime = get_wall_time();

  environment.memoryThreads = new MemorySpace[environment.n_threads];
  for (uint i = 0; i < environment.n_threads; i++) {
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
  if (environment.verbose == true) {
    std::cout << "\n# ----> First contributing node elapsed time:" << spentTime
              << "sec\n\n";
  }
  // Run the skeleton iteration phase if consistency is required.
  BCC bcc(environment);
  auto cycle_tracker = CycleTracker(environment);
  vector<vector<string> > confVect;
  vector<vector<string> > orientations;
  do {
    if (environment.consistent > 0) bcc.analyse();
    // Save the neighbours in the status_prev structure
    // and revert to the structure at the moment of initialization
    for (uint i = 0; i < environment.n_nodes; i++) {
      for (uint j = 0; j < environment.n_nodes; j++) {
        environment.edges[i][j].status_prev = environment.edges[i][j].status;
        environment.edges[i][j].status = environment.edges[i][j].status_init;
      }
    }
    // If interrupted
    if (!firstStepIteration(environment, bcc)) return empty_results();

    if (environment.numNoMore == 0 && environment.numSearchMore == 0) {
      if (environment.verbose)
        std::cout << "# ------| Only phantom edges found.\n";
    } else if (environment.numSearchMore > 0) {
      // Search for other Contributing node(s) (possible only for the edges
      /// still in 'searchMore', ie. 2)
      if (environment.verbose) {
        std::cout << "\n# ---- Other Contributing node(s) ----\n\n";
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
      std::cout << "Computing confidence cut with permutations..." << std::flush;
      confVect = confidenceCut(environment);
      long double spentTime = (get_wall_time() - startTime);
      environment.exec_time.cut += spentTime;
      std::cout << " done." << std::endl;
    } else {
      environment.exec_time.cut = 0;
    }
    // Oriente edges for non-consistent/orientation consistent algorithm
    if (orientation_phase && environment.numNoMore > 0 &&
        environment.consistent <= 1) {
      orientations = orientationProbability(environment);
    }
    std::cout << "Number of edges: " << environment.numNoMore << std::endl;
  } while (environment.consistent > 0 && !cycle_tracker.hasCycle());

  int union_n_edges = 0;
  for (uint i = 1; i < environment.n_nodes; i++) {
    for (uint j = 0; j < i; j++) {
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
    uint n_inconsistency = 0;
    std::vector<std::pair<uint, uint> > inconsistent_edges;
    for (uint i = 1; i < environment.n_nodes; i++) {
      for (uint j = 0; j < i; j++) {
        const Edge& edge = environment.edges[i][j];
        if (edge.status || bcc.is_consistent(i, j, edge.shared_info->ui_list))
          continue;
        if (environment.verbose) {
          std::cout << environment.nodes[i].name << ",\t"
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
    std::cout << n_inconsistency << " inconsistent conditional independences"
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

  for (uint i = 0; i < environment.n_threads; i++) {
    deleteMemorySpace(environment, environment.memoryThreads[i]);
  }
  deleteMemorySpace(environment, environment.m);
  deleteStruct(environment);

  return result;
}

namespace miic {
namespace reconstruction {

bool CycleTracker::hasCycle() {
  uint n_edge = env_.numNoMore;
  // Before saving the current iteration, search among previous iterations
  // those with the same number of edges
  auto range = edge_index_map_.equal_range(n_edge);
  bool no_cycle_found = range.first == range.second;
  // Indices of iteration that is possibly the end point of a cycle in the
  // backtracking sense. Example: suppose that iteration #1, #3, #6 have the
  // same number of edges as the current iteration #8, and each of them is a
  // possible start point of a cycle, therefore #2, #4, #7 are the corresponding
  // possible end points of the cycle (#8 #7 ... #2), (#8 #7 ... #4), (#8 #7).
  std::deque<uint> iter_indices;
  for (auto it = range.first; it != range.second; ++it)
    iter_indices.push_back(it->second + 1);
  saveIteration();
  if (n_saved > env_.max_iteration) {
    std::cout << "Max number of iterations reached: " << env_.max_iteration
              << '\n';
    return true;
  }
  if (no_cycle_found) return false;
  // Backtracking requires starting from the largest index first
  std::sort(iter_indices.begin(), iter_indices.end(), std::greater<uint>());
  // Set of edges that are to be marked as connected and undirected
  std::set<uint> edges_union;
  // Check if an edge is changed. vector is chosen over map for quicker access
  // and simpler syntax, at the cost of extra memory trace and (possible) extra
  // time complexity (in practice there are very few changes between each pair
  // of iterations).
  std::vector<int> changed(env_.n_nodes * (env_.n_nodes - 1) / 2, 0);
  // backtracking over iteration to get changed_edges
  int cycle_size = 0;
  for (const auto& iter : iterations_) {
    ++cycle_size;
    for (const auto& k : iter.changed_edges) {
      edges_union.insert(k.first);
      // compare edge status in the previous iteration against the latest edge
      // status
      std::pair<uint, uint> p = getEdgeIndex2D(k.first);
      changed[k.first] = (k.second != env_.edges[p.first][p.second].status);
    }
    if (iter.index != iter_indices.front()) continue;
    iter_indices.pop_front();
    // if any edge has been changed
    if (std::any_of(
            changed.begin(), changed.end(), [](uint j) { return j != 0; })) {
      // no cycle
      if (iter_indices.empty()) return false;
      continue;
    }
    for (auto& k : edges_union) {
      std::pair<uint, uint> p = getEdgeIndex2D(k);
      env_.edges[p.first][p.second].status = 1;
      env_.edges[p.second][p.first].status = 1;
      env_.edges[p.first][p.second].shared_info->setUndirected();
    }

    std::cout << "cycle found of size " << cycle_size << std::endl;
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
