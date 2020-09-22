#include "r_cpp_interface.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include <vector>

namespace miic {
namespace utility {

constexpr double kMagnitudeTies = 0.00005;

using namespace structure;
using Rcpp::as;
using std::vector;

void setEnvironmentFromR(const Rcpp::List& input_data,
    const Rcpp::List& arg_list, Environment& environment) {
  const int n_nodes = environment.n_nodes;
  const int n_samples = environment.n_samples;
  if (input_data.containsElementNamed("double"))
    environment.data_double = Grid2d<double>(
        n_nodes, n_samples, as<vector<double>>(input_data["double"]));

  vector<int> list_delta_taus;
  if ( arg_list.containsElementNamed ("tau") ) {
    environment.list_taus = as<vector<int>> (arg_list["tau"]);
    environment.tau_max = *std::max_element ( environment.list_taus.begin(),
                                              environment.list_taus.end() );
  }
  if (environment.tau_max >= 1) {
    environment.n_nodes_not_lagged = n_nodes;
    for (auto& one_tau : environment.list_taus)
      if (one_tau >= 1)
        environment.n_nodes_not_lagged -= one_tau;

    list_delta_taus = vector<int> (1, environment.n_nodes_not_lagged);
    if ( arg_list.containsElementNamed ("delta_tau") )
      list_delta_taus = as<vector<int>> (arg_list["delta_tau"]);
  }

  if (arg_list.containsElementNamed("n_eff"))
    environment.n_eff = as<double>(arg_list["n_eff"]);
  if (environment.n_eff < 0 || environment.n_eff > n_samples) {
    environment.n_eff = static_cast<double>(n_samples);

    if (environment.tau_max >= 1) {
      //
      // In temporal mode, when delta_tau is > 1, we need to tune the n_eff
      // Note: this implementation is intended for an unique delta_tau
      // (excepting non lagged variables like contextual ones)
      // and needs an update to take into account different delta_taus
      //
      int delta_tau_max = *std::max_element ( list_delta_taus.begin(),
                                              list_delta_taus.end() );
      if (delta_tau_max > 1) {
        environment.n_eff = n_samples / delta_tau_max;
        Rcpp::Rcout << "Note: n_eff has been set to " <<  environment.n_eff
                    << " (nb lagged samples=" << n_samples
                    << " / max delta_tau=" << delta_tau_max << ")\n";
      }
    }
  }
    
  if (arg_list.containsElementNamed("var_names")) {
    auto var_names = as<vector<std::string>>(arg_list["var_names"]);
    std::transform(var_names.begin(), var_names.end(),
        std::back_inserter(environment.nodes),
        [](std::string name) { return Node(name); });
  }

  if (arg_list.containsElementNamed("orientation"))
    environment.orientation = as<bool>(arg_list["orientation"]);

  if (arg_list.containsElementNamed("ori_proba_ratio"))
    environment.ori_proba_ratio = as<double>(arg_list["ori_proba_ratio"]);

  if (arg_list.containsElementNamed("propagation"))
    environment.propagation = as<bool>(arg_list["propagation"]);

  if (arg_list.containsElementNamed("consistent")) {
    auto consistent_flag = as<std::string>(arg_list["consistent"]);
    if (consistent_flag.compare("orientation") == 0)
      environment.consistent = 1;
    else if (consistent_flag.compare("skeleton") == 0)
      environment.consistent = 2;
  }

  if (arg_list.containsElementNamed("max_iteration"))
    environment.max_iteration = as<int>(arg_list["max_iteration"]);

  if (arg_list.containsElementNamed("latent")) {
    auto latent_flag = as<std::string>(arg_list["latent"]);
    if (latent_flag.compare("yes") == 0) {
      environment.latent = true;
      environment.latent_orientation = true;
    } else if (latent_flag.compare("orientation") == 0)
      environment.latent_orientation = true;
  }

  if (arg_list.containsElementNamed("test_mar"))
    environment.test_mar = as<bool>(arg_list["test_mar"]);

  if (arg_list.containsElementNamed("cplx")) {
    if (as<std::string>(arg_list["cplx"]).compare("mdl") == 0)
      environment.cplx = 0;
  }

  if (arg_list.containsElementNamed("n_shuffles"))
    environment.n_shuffles = as<int>(arg_list["n_shuffles"]);

  if (arg_list.containsElementNamed("conf_threshold"))
    environment.conf_threshold = as<double>(arg_list["conf_threshold"]);

  if (arg_list.containsElementNamed("is_contextual")) {
    environment.is_contextual = as<vector<int>>(arg_list["is_contextual"]);
  } else {
    environment.is_contextual.resize(n_nodes, 0);
  }
  //
  // Variables considered as consequence only
  //
  if (arg_list.containsElementNamed("is_consequence")) {
    environment.is_consequence = as<vector<int>>(arg_list["is_consequence"]);
  } else {
    environment.is_consequence.resize(n_nodes, 0);
  }
  environment.any_consequence = std::any_of (environment.is_consequence.begin(),
    environment.is_consequence.end(), [](bool v) { return v; });
  //
  // Remove edges between consequence variables
  //
  if (environment.any_consequence) {
    for (int i = 0; i < n_nodes; i++)
      for (int j = 0; j < n_nodes; j++)
        if (environment.is_consequence[i] && environment.is_consequence[j]) {
          environment.edges(i, j).status = 0;
          environment.edges(i, j).status_init = 0;
          environment.edges(i, j).status_prev = 0;
        }
  }

  if (arg_list.containsElementNamed("sample_weights"))
    environment.sample_weights = as<vector<double>>(arg_list["sample_weights"]);
  if (environment.sample_weights.empty()) {
    double uniform_weight = environment.n_eff / n_samples;
    environment.sample_weights.resize(n_samples, uniform_weight);
  }

  std::generate(begin(environment.noise_vec), end(environment.noise_vec),
      []() { return R::runif(-kMagnitudeTies, kMagnitudeTies); });

  if (arg_list.containsElementNamed("is_k23"))
    environment.is_k23 = as<bool>(arg_list["is_k23"]);

  if (arg_list.containsElementNamed("degenerate"))
    environment.degenerate = as<bool>(arg_list["degenerate"]);

  if (arg_list.containsElementNamed("no_init_eta"))
    environment.no_init_eta = as<bool>(arg_list["no_init_eta"]);

  if (arg_list.containsElementNamed("half_v_structure"))
    environment.half_v_structure = as<bool>(arg_list["half_v_structure"]);

  if (arg_list.containsElementNamed("max_bins"))
    environment.maxbins = as<int>(arg_list["max_bins"]);
  if (environment.maxbins > n_samples)
    environment.maxbins = n_samples;

  if (arg_list.containsElementNamed("n_threads"))
    environment.n_threads = as<int>(arg_list["n_threads"]);
#ifdef _OPENMP
  if (environment.n_threads <= 0) environment.n_threads = omp_get_num_procs();
  omp_set_num_threads(environment.n_threads);
#endif

  if (environment.tau_max >= 1) {
    //
    // Precompute the class of each variable (unique ID for the part before "_lagX")
    //
    for (int node_idx = 0; node_idx < environment.n_nodes_not_lagged; ++node_idx)
      environment.nodes_class.push_back (node_idx);
    for (int tau_idx = 1; tau_idx <= environment.tau_max; ++tau_idx)
      for (int node_idx = 0; node_idx < environment.n_nodes_not_lagged; ++node_idx)
        if (tau_idx <= environment.list_taus[node_idx])
          environment.nodes_class.push_back (node_idx);
    //
    // Precompute the lag of each variable: layer * delta_tau
    //
    environment.nodes_lags.assign (environment.n_nodes_not_lagged, 0);
    for (int tau_idx = 1; tau_idx <= environment.tau_max; ++tau_idx)
      for (int node_idx = 0; node_idx < environment.n_nodes_not_lagged; ++node_idx)
        if (tau_idx <= environment.list_taus[node_idx])
          environment.nodes_lags.push_back (tau_idx * list_delta_taus[node_idx]);
    //
    // For contextual variables, we consider them as very old, 
    // so they are never the consequence of another variable
    //
    for (int ctx_idx = 0; ctx_idx < environment.is_contextual.size();  ++ctx_idx)
      if (environment.is_contextual[ctx_idx])
        environment.nodes_lags[ctx_idx] = INT_MAX;
    //
    // Precompute the index shifts from a var to its next lagged counterpart
    //
    int n_nodes_shifts = environment.n_nodes_not_lagged;
    vector <bool> end_reached (environment.n_nodes_not_lagged, false);
    for (int tau_idx = 1; tau_idx <= environment.tau_max+1; ++tau_idx)
      for (int node_idx = 0; node_idx < environment.n_nodes_not_lagged; ++node_idx)
        if (tau_idx <= environment.list_taus[node_idx])
          environment.nodes_shifts.push_back (n_nodes_shifts);
        else if (!end_reached[node_idx]) {
          end_reached[node_idx] = true;
          environment.nodes_shifts.push_back (0);
          --n_nodes_shifts;
        }
    //
    // In temporal mode, we do not start from a complete graph
    // => Remove all edges not having a node on the layer 0
    //
    for (int i = environment.n_nodes_not_lagged; i < n_nodes; i++) 
      for (int j = environment.n_nodes_not_lagged; j < n_nodes; j++) {
        environment.edges(i, j).status = 0;
        environment.edges(i, j).status_init = 0;
        environment.edges(i, j).status_prev = 0;
        environment.edges(i, j).proba_head = -1;
      }
    //
    // In addition, non lagged variable (i.e.:contextual) can only have edges
    // with nodes of the layer 0 (others can be found by stationarity)
    //
    for (int i = 0; i < environment.n_nodes_not_lagged; i++) 
      if (environment.is_contextual[i])
        for (int j = environment.n_nodes_not_lagged; j < n_nodes; j++) {
          environment.edges(i, j).status = 0;
          environment.edges(i, j).status_init = 0;
          environment.edges(i, j).status_prev = 0;
          environment.edges(i, j).proba_head = -1;
          environment.edges(j, i).status = 0;
          environment.edges(j, i).status_init = 0;
          environment.edges(j, i).status_prev = 0;
          environment.edges(j, i).proba_head = -1;
        }
  }
  
  if (arg_list.containsElementNamed("negative_info"))
    environment.negative_info = as<bool>(arg_list["negative_info"]);

  if (arg_list.containsElementNamed("verbose"))
    environment.verbose = as<bool>(arg_list["verbose"]);

  if (arg_list.containsElementNamed("black_box")) {
    auto black_box_vec = as<vector<int>>(arg_list["black_box"]);
    int n_pairs = black_box_vec.size() / 2;
    environment.readBlackbox(Grid2d<int>(n_pairs, 2, std::move(black_box_vec)));
  }
}
}  // namespace utility
}  // namespace miic
