#ifndef MIIC_ENVIRONMENT_H_
#define MIIC_ENVIRONMENT_H_

#include <Rcpp.h>

#include <map>
#include <string>
#include <vector>

#include "computation_cache.h"
#include "structure.h"

#define N_COL_NML 1000
namespace miic {
namespace structure {

namespace structure_impl {

using std::string;
using std::vector;
using computation::CompCache;

struct Environment {
  vector<vector<int>> data_numeric;
  int n_samples;
  int n_nodes;
  vector<vector<double>> data_double;
  // data_numeric_idx[i] = index of i'th smallest value in data_double
  vector<vector<int>> data_numeric_idx;
  int* oneLineMatrix;

  vector<int> is_continuous;
  vector<int> levels;
  int n_eff;
  vector<Node> nodes{};
  Edge** edges;
  bool orientation_phase;
  bool propagation;
  // Level of consistency required for the graph
  // 0: no consistency requirement
  // 1: skeleton consistent
  // 2: orientation consistent
  int consistent{0};
  // When consistent > 0, the maximum number of iterations allowed when trying
  // to find a consistent graph.
  int max_iteration;
  // A latent (conditioning) node (w.r.t. node X and Y) is a node that is a
  // neighbor of neither X nor Y.
  // true: consider latent node during the search of conditioning nodes as well
  // as during the orientation.
  bool latent{false};
  // true: consider latent node during the orientation only.
  bool latent_orientation{false};
  // Whether or not do MAR (Missing at random) test using KL-divergence
  bool test_mar;
  // Complexity mode. 0: mdl 1: nml
  int cplx{1};
  // If firstStepIteration is done
  bool first_iter_done{false};
  // List of ids of edge whose status is not yet determined
  vector<EdgeID> unsettled_list;
  // List of ids of edge whose status is sure to be connected
  vector<EdgeID> connected_list;
  int numSearchMore{0};
  int numNoMore{0};

  int n_shuffles;
  double conf_threshold;

  int** iterative_cuts;
  vector<double> sample_weights;
  bool flag_sample_weights;
  double* noise_vec;

  MemorySpace m;
  MemorySpace* memoryThreads;

  double log_eta = 0;
  bool is_k23;
  bool degenerate;
  bool no_init_eta{false};
  int half_v_structure;

  // Set the probability threshold for the rank if the contribution probability
  // is the min value
  int thresPc{0};

  int maxbins;
  int initbins;

  std::map<CacheInfoKey, double> look_scores;
  std::map<CacheInfoKey, CacheScoreValue> look_scores_orientation;

  ExecutionTime exec_time;
  int n_threads;
  CompCache cache;
  bool verbose;

  Environment(const Rcpp::List&, const Rcpp::List&);
  Environment() = default;

  ~Environment() {
    delete[] oneLineMatrix;
    for (int i = 0; i < n_nodes; i++) delete[] edges[i];
    delete[] edges;
    delete[] noise_vec;
  }

  void readBlackbox(const vector<vector<int>>&);
};

}  // namespace structure_impl
using structure_impl::Environment;
}  // namespace structure
}  // namespace miic

#endif  // MIIC_ENVIRONMENT_H_
