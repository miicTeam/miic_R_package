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

namespace miic {
namespace reconstruction {

using uint = unsigned int;
using Rcpp::_;
using Rcpp::List;
using std::string;
using std::vector;
using namespace structure;
using namespace utility;

bool skeletonChanged(const Environment& environment) {
  for (uint i = 0; i < environment.numNodes; i++) {
    for (uint j = 0; j < environment.numNodes; j++) {
      if (environment.edges[i][j].status !=
          environment.edges[i][j].status_prev) {
        return true;
      }
    }
  }
  return false;
}

List empty_results() {
  List result;
  result = List::create(_["interrupted"] = true);
  return (result);
}

extern "C" SEXP reconstruct(SEXP inputDataR, SEXP typeOfDataR, SEXP cntVarR,
    SEXP numNodesR, SEXP nThreadsR, SEXP edgefileR, SEXP blackBoxR, SEXP effNR,
    SEXP cplxR, SEXP etaR, SEXP hvsR, SEXP isLatentR, SEXP isTplReuseR,
    SEXP isK23R, SEXP isDegeneracyR, SEXP orientationR, SEXP propagationR,
    SEXP isNoInitEtaR, SEXP confidenceShuffleR, SEXP confidenceThresholdR,
    SEXP sampleWeightsR, SEXP consistentR, SEXP testDistR, SEXP verboseR) {
  vector<vector<double> > outScore;
  vector<vector<string> > edgesMatrix;
  std::stringstream output;

  Environment environment;
  environment.myVersion = "V79";

  environment.myTest = false;

  environment.numNodes = Rcpp::as<int>(numNodesR);
  environment.nThreads = Rcpp::as<int>(nThreadsR);
#ifdef _OPENMP
  int threads = environment.nThreads;
  if (threads < 0) threads = omp_get_num_procs();
  omp_set_num_threads(threads);
#endif

  environment.testDistribution = Rcpp::as<bool>(testDistR);
  environment.vectorData = Rcpp::as<vector<string> >(inputDataR);

  environment.effN = Rcpp::as<int>(effNR);
  environment.typeOfData = Rcpp::as<int>(typeOfDataR);
  environment.halfVStructures = Rcpp::as<int>(hvsR);

  string cplx_flag = Rcpp::as<string>(cplxR);
  environment.cplx = 1;  // use nml complexity by default
  if (cplx_flag.compare("mdl") == 0) environment.cplx = 0;

  std::string consistent_flag = Rcpp::as<std::string>(consistentR);
  environment.consistentPhase = 0;
  if (consistent_flag.compare("orientation") == 0)
    environment.consistentPhase = 1;
  if (consistent_flag.compare("skeleton") == 0) environment.consistentPhase = 2;

  std::string latent_flag = Rcpp::as<std::string>(isLatentR);
  environment.isLatent = false;
  environment.isLatentOnlyOrientation = false;
  if (latent_flag.compare("yes") == 0) environment.isLatent = true;
  if (latent_flag.compare("ort") == 0)
    environment.isLatentOnlyOrientation = true;

  environment.isTplReuse = Rcpp::as<bool>(isTplReuseR);
  environment.isK23 = Rcpp::as<bool>(isK23R);
  environment.isDegeneracy = Rcpp::as<bool>(isDegeneracyR);
  bool orientation_phase = Rcpp::as<bool>(orientationR);
  environment.isPropagation = Rcpp::as<bool>(propagationR);
  environment.isNoInitEta = Rcpp::as<bool>(isNoInitEtaR);

  environment.confidenceThreshold = Rcpp::as<double>(confidenceThresholdR);
  environment.numberShuffles = Rcpp::as<int>(confidenceShuffleR);

  environment.sampleWeightsVec = Rcpp::as<vector<double> >(sampleWeightsR);
  environment.cntVarVec = Rcpp::as<vector<int> >(cntVarR);

  environment.isVerbose = Rcpp::as<bool>(verboseR);

  double startTime;

  vector<vector<string> > retShuffle;
  vector<vector<string> > retShuffleAvg;
  vector<string> shfVec;

  string filePath;

  environment.isAllGaussian = false;

  // set the environment
  srand(0);
  setEnvironment(environment);
  vector<string> v(Rcpp::as<vector<string> >(blackBoxR));
  if (v.size() > 1) readBlackbox(v, environment);

  startTime = get_wall_time();

  environment.memoryThreads = new MemorySpace[environment.nThreads];
  for (uint i = 0; i < environment.nThreads; i++) {
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
  environment.execTime.init = spentTime;
  if (environment.isVerbose == true) {
    std::cout << "\n# ----> First contributing node elapsed time:" << spentTime
              << "sec\n\n";
  }
  // Run the skeleton iteration phase if consistency is required.
  BCC bcc(environment);
  auto cycle_tracker = CycleTracker(environment);
  vector<vector<string> > confVect;
  vector<vector<string> > orientations;
  do {
    if (environment.consistentPhase) bcc.analyse();
    // Save the neighbours in the status_prev structure
    // and revert to the structure at the moment of initialization
    for (uint i = 0; i < environment.numNodes; i++) {
      for (uint j = 0; j < environment.numNodes; j++) {
        environment.edges[i][j].status_prev = environment.edges[i][j].status;
        environment.edges[i][j].status = environment.edges[i][j].status_init;
      }
    }
    // If interrupted
    if (!firstStepIteration(environment, bcc)) return empty_results();

    if (environment.numNoMore == 0 && environment.numSearchMore == 0) {
      if (environment.isVerbose == true)
        std::cout << "# ------| Only phantom edges found.\n";
    } else if (environment.numSearchMore > 0) {
      // Search for other Contributing node(s) (possible only for the edges
      /// still in 'searchMore', ie. 2)
      if (environment.isVerbose == true) {
        std::cout << "\n# ---- Other Contributing node(s) ----\n\n";
      }
      startTime = get_wall_time();

      // If interrupted
      if (!skeletonIteration(environment)) return (empty_results());

      long double spentTime = (get_wall_time() - startTime);
      environment.execTime.iter = spentTime;
      environment.execTime.initIter =
          environment.execTime.init + environment.execTime.iter;
    }

    startTime = get_wall_time();
    if (environment.numberShuffles > 0) {
      std::cout << "Computing confidence cut with permutations..." << std::flush;
      confVect = confidenceCut(environment);
      long double spentTime = (get_wall_time() - startTime);
      environment.execTime.cut = spentTime;
      std::cout << " done." << std::endl;
    } else {
      environment.execTime.cut = 0;
    }
    // Do not test distributions for missing at random for orientation
    environment.testDistribution = false;

    // Oriente edges for non-consistent/orientation consistent algorithm
    if (orientation_phase && environment.numNoMore > 0 &&
        environment.consistentPhase <= 1) {
      orientations = orientationProbability(environment);
    }
    std::cout << "Number of edges: " << environment.numNoMore << std::endl;
  } while (environment.consistentPhase && !cycle_tracker.hasCycle());

  int union_n_edges = 0;
  for (uint i = 1; i < environment.numNodes; i++) {
    for (uint j = 0; j < i; j++) {
      if (environment.edges[i][j].status) {
        union_n_edges++;
      }
    }
  }
  environment.numNoMore = union_n_edges;

  if (environment.numNoMore > 0) {
    // skeleton consistent algorithm
    if (environment.consistentPhase == 2) {
      orientations = orientationProbability(environment);
      // Check inconsistency after orientation, add undirected edge to
      // pairs with inconsistent conditional independence.
      bcc.analyse();
      uint n_inconsistency = 0;
      std::vector<std::pair<uint, uint> > inconsistent_edges;
      for (uint i = 1; i < environment.numNodes; i++) {
        for (uint j = 0; j < i; j++) {
          const Edge& edge = environment.edges[i][j];
          if (edge.status ||
              bcc.is_consistent(i, j, edge.shared_info->ui_vect_idx))
            continue;
          if (environment.isVerbose) {
            std::cout << environment.nodes[i].name << ",\t"
                      << environment.nodes[j].name << "\t| "
                      << vectorToStringNodeName(
                             environment, edge.shared_info->ui_vect_idx)
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
    edgesMatrix = saveEdgesListAsTable(environment);
  }

  vector<vector<string> > adjMatrix;
  adjMatrix = getAdjMatrix(environment);

  vector<double> time;
  time.push_back(environment.execTime.init);
  time.push_back(environment.execTime.iter);
  time.push_back(environment.execTime.cut);
  time.push_back(environment.execTime.initIter + environment.execTime.cut);

  List result;
  if (environment.numberShuffles > 0) {
    result = List::create(
        _["confData"]          = confVect,
        _["adjMatrix"]         = adjMatrix,
        _["edges"]             = edgesMatrix,
        _["orientations.prob"] = orientations,
        _["time"]              = time,
        _["interrupted"]       = false);
  } else {
    result = List::create(
        _["adjMatrix"]         = adjMatrix,
        _["edges"]             = edgesMatrix,
        _["orientations.prob"] = orientations,
        _["time"]              = time,
        _["interrupted"]       = false);
  }
  for (uint i = 0; i < environment.nThreads; i++) {
    deleteMemorySpace(environment, environment.memoryThreads[i]);
  }
  deleteMemorySpace(environment, environment.m);
  deleteStruct(environment);

  return result;
}

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
  if (n_saved > max_iter_size) {
    std::cout << "Max number of iterations reached: " << max_iter_size << '\n';
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
  std::vector<int> changed(env_.numNodes * (env_.numNodes - 1) / 2, 0);
  // backtracking over iteration to get changed_edges
  uint cycle_size = 0;
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
      if (iter_indices.empty()) break;
      continue;
    }
    for (auto& k : edges_union) {
      std::pair<uint, uint> p = getEdgeIndex2D(k);
      env_.edges[p.first][p.second].status = 1;
      env_.edges[p.second][p.first].status = 1;
      env_.edges[p.first][p.second].shared_info->setUndirected();
    }

    std::cout << "cycle found of size " << cycle_size << std::endl;
    return true;
  }
  return false;
}

}  // namespace reconstruction
}  // namespace miic
