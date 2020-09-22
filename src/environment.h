#ifndef MIIC_ENVIRONMENT_H_
#define MIIC_ENVIRONMENT_H_

#include <vector>

#include "computation_cache.h"
#include "structure.h"

namespace miic {
namespace structure {

namespace detail {

using std::vector;
using computation::CompCache;

struct Environment {
  int n_samples;
  int n_nodes;
  // In temporal node, store the number of non lagged nodes 
  int n_nodes_not_lagged=-1;
  Grid2d<int> data_numeric;
  Grid2d<double> data_double;
  // data_numeric_idx(j, i) = index of i'th smallest value in data_double(j, )
  Grid2d<int> data_numeric_idx;

  vector<int> is_continuous;
  vector<int> levels;
  vector<int> has_na;
  int n_eff;
  vector<Node> nodes;
  Grid2d<Edge> edges;
  bool orientation = false;
  double ori_proba_ratio = 1;
  bool propagation = false;
  // Level of consistency required for the graph
  // 0: no consistency requirement
  // 1: skeleton consistent
  // 2: orientation consistent
  int consistent = 0;
  // When consistent > 0, the maximum number of iterations allowed when trying
  // to find a consistent graph.
  int max_iteration = 100;
  // A latent (conditioning) node (w.r.t. node X and Y) is a node that is a
  // neighbor of neither X nor Y.
  // true: consider latent node during the search of conditioning nodes as well
  // as during the orientation.
  bool latent = false;
  // true: consider latent node during the orientation only.
  bool latent_orientation = false;
  // Whether or not do MAR (Missing at random) test using KL-divergence
  bool test_mar = false;
  // Complexity mode. 0: mdl 1: nml
  int cplx = 1;
  // List of ids of edge whose status is not yet determined
  vector<EdgeID> unsettled_list;
  // List of ids of edge whose status is sure to be connected
  vector<EdgeID> connected_list;

  int n_shuffles = 0;
  double conf_threshold = 0;

  vector<double> sample_weights;
  bool flag_sample_weights = false;
  vector<double> noise_vec;

  double log_eta = 0;
  bool is_k23 = true;
  bool degenerate = false;
  bool no_init_eta = false;
  bool half_v_structure = false;

  int maxbins = 50;
  int initbins;

  ExecutionTime exec_time;
  int n_threads = 1;
  CompCache cache;
  bool verbose = false;

  // Maximum lag. Switch miic to temporal mode if >=1
  int tau=-1;
  
  Environment(int n_samples, int n_nodes, vector<int> vec_numeric,
      vector<int> vec_index, vector<int> is_continuous_, vector<int> levels_);
  Environment() = default;

  void readBlackbox(const Grid2d<int>&);
};

}  // namespace detail
using detail::Environment;
}  // namespace structure
}  // namespace miic

#endif  // MIIC_ENVIRONMENT_H_
