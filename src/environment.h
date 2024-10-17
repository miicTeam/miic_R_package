#ifndef MIIC_ENVIRONMENT_H_
#define MIIC_ENVIRONMENT_H_

#include <vector>
#include <string>

#include "computation_cache.h"
#include "structure.h"
#include "layers.h"

namespace miic {
namespace structure {

namespace detail {

using std::vector;
using computation::CompCache;

struct Environment {
  int n_samples;
  int n_nodes;
  // Used only in temporal mode, number of nodes at lag0 or not lagged
  int n_nodes_not_lagged = -1;

  Grid2d<int> data_numeric;
  Grid2d<double> data_double;
  // data_numeric_idx(j, i) = index of i'th smallest value in data_double(j, )
  Grid2d<int> data_numeric_idx;
  // Store the miic's mode: 0=Standard (IID), 1=Temporal (stationary)
  // (use an int to be ready for addition of future modes)
  int mode = 0;
  // As we foresee to have a stationary and a non stationary temporal version,
  // this flag is true if we are in a temporal mode, whatever stationarity
  bool temporal = false;

  // Identify if any node is marked as contextual
  bool any_contextual = false;
  // For each node, contains 0 = not a contextual node or 1 = contextual node
  vector<int> is_contextual;
  // Identify if any node is marked as consequence
  bool any_consequence = false;
  // For each node, contains 0 = not a consequence node or 1 = consequence node
  vector<int> is_consequence;
  vector<int> is_continuous;
  vector<int> levels;
  vector<int> has_na;
  double n_eff;
  vector<Node> nodes;
  Grid2d<Edge> edges;
  bool orientation = false;
  double ort_proba_ratio = 1;
  bool propagation = false;
  // Level of consistency required for the graph
  // 0: no consistency requirement
  // 1: orientation consistent
  // 2: skeleton consistent
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
  // Complexity mode. 0: bic (formerly mdl) 1: nml
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
  bool degenerate = false;
  bool no_init_eta = false;
  bool half_v_structure = false;
  // If true, allow for negative shifted mutual information
  bool negative_info = false;

  int maxbins = 50;
  int initbins;

  ExecutionTime exec_time;
  int n_threads = 1;
  CompCache cache;

  bool verbose = false;
  //
  // Temporal mode
  //
  // Max number of layers
  // Even if not recommended, the number of layers can be different for each
  // variables, layer_max is the maximum number of layers for all variables
  int layer_max = -1;
  // Number of layers for each variable
  vector<int> list_n_layers;
  // Class for each node
  vector<int> nodes_class;
  // Lag for each node
  // Note that we consider contextual variables as very old (= INT_MAX)
  // so from the time point of view, they are never the consequence
  // of another variable
  vector<int> nodes_lags;
  // Store nodes index shift, giving for each node the same lagged node
  // in the next layer
  // (e.g.  variables: x_lag0, ctx_var, y_lag0, x_lag1, y_lag1
  //  => nodes_shifts:   3   ,    0   ,   2   ,   0   ,   0)
  vector<int> nodes_shifts;
  //
  // Multi-layers mode
  //
  bool is_layered = false;
  vector<int> nodes_layers;
  vector<Layer> layers;
  //
  // Constructors
  //
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
