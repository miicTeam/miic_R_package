#include "environment.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>

#include "compute_info.h"
#include "utilities.h"

using Rcpp::as;
using std::string;
using std::vector;
using namespace miic::utility;
using namespace miic::computation;

#define MAGNITUDE_TIES 0.00005f
#define N_COL_NML 1000

namespace miic {
namespace structure {

Environment::Environment(
    const Rcpp::DataFrame& input_data, const Rcpp::List& arg_list)
    : data(as<vector<vector<string>>>(input_data)),
      n_samples(data.size()),
      n_nodes(data[0].size()),
      data_double(n_nodes),
      data_numeric(n_samples, vector<int>(n_nodes)),
      is_continuous(as<vector<int>>(arg_list["is_continuous"])),
      levels(n_nodes),
      n_eff(as<int>(arg_list["n_eff"])),
      orientation_phase(as<bool>(arg_list["orientation"])),
      propagation(as<bool>(arg_list["propagation"])),
      max_iteration(as<int>(arg_list["max_iteration"])),
      test_mar(as<bool>(arg_list["test_mar"])),
      n_shuffles(as<int>(arg_list["n_shuffles"])),
      conf_threshold(as<double>(arg_list["conf_threshold"])),
      sample_weights(as<vector<double>>(arg_list["sample_weights"])),
      is_k23(as<bool>(arg_list["is_k23"])),
      degenerate(as<bool>(arg_list["degenerate"])),
      no_init_eta(as<bool>(arg_list["no_init_eta"])),
      half_v_structure(as<int>(arg_list["half_v_structure"])),
      initbins(std::min(30, int(0.5 + cbrt(n_samples)))),
      n_threads(as<int>(arg_list["n_threads"])),
      verbose(as<bool>(arg_list["verbose"])) {
  srand(0);
  readFileType();
  if (n_eff == -1 || n_eff > n_samples)
    n_eff = n_samples;
  auto var_names = as<vector<string>>(arg_list["var_names"]);
  std::transform(var_names.begin(), var_names.end(), std::back_inserter(nodes),
      [](string name) { return Node(name); });

  auto consistent_flag = as<std::string>(arg_list["consistent"]);
  if (consistent_flag.compare("orientation") == 0)
    consistent = 1;
  else if (consistent_flag.compare("skeleton") == 0)
    consistent = 2;

  auto latent_flag = as<std::string>(arg_list["latent"]);
  if (latent_flag.compare("yes") == 0)
    latent = true;
  else if (latent_flag.compare("orientation") == 0)
    latent_orientation = true;

  if (as<string>(arg_list["cplx"]).compare("mdl") == 0) cplx = 0;

  if (sample_weights.empty()) {
    double uniform_weight{1};
    if (n_eff != n_samples)
      uniform_weight = static_cast<double>(n_eff) / n_samples;
    sample_weights.resize(n_samples, uniform_weight);
  }

#ifdef _OPENMP
  if (n_threads < 0) n_threads = omp_get_num_procs();
  omp_set_num_threads(n_threads);
#endif
  // for continuous
  auto any_continuous = std::any_of(
      is_continuous.begin(), is_continuous.end(), [](int i) { return i == 1; });
  if (any_continuous) {
    // create the data matrix for factors indexes
    data_numeric_idx = new int*[n_nodes];
    for (int i = 0; i < n_nodes; i++) {
      data_numeric_idx[i] = new int[n_samples];
      for (int j = 0; j < n_samples; j++) data_numeric_idx[i][j] = -1;
    }
  }
  // transform to factors
  for (int i = 0; i < n_nodes; i++) {
    if (!is_continuous[i])
      transformToFactors(i);
    else
      // update data_numeric not taking into account repetition
      transformToFactorsContinuous(i);
  }
  // for continuous
  if (any_continuous) {
    for (int j = 0; j < n_nodes; j++) {
      if (is_continuous[j]) {
        // update data_numeric taking into account repetition
        transformToFactorsContinuousIdx(j);
        transformToFactors(j);
      }
    }
  }
  // Set a variables with all properties name and levels
  setNumberLevels();
  // create the 1000 entries to store c2 values
  c2terms = new double[n_samples + 1];
  for (int i = 0; i < n_samples + 1; i++) {
    c2terms[i] = -1;
  }
  // for mixed, create the log(j) lookup table with j=0..n_samples;
  looklog = new double[n_samples + 2];
  looklog[0] = 0.0;
  for (int i = 1; i < n_samples + 2; i++) {
    looklog[i] = log(1.0 * i);
  }

  lookH = new double[n_samples + 2];
  lookH[0] = 0.0;
  for (int i = 1; i < n_samples + 2; i++) {
    // lookH[i] = i*looklog[i]-(i+1)*looklog[(i+1)];
    lookH[i] = i * looklog[i];
  }
  // Number of levels r for which we want to store the
  // stochastic NML complexity LogC(n,r) for n in [1,N].
  // For r>N_COL_NML LogC() is computed with the normal
  // recurrence and the result is not stored.
  int ncol = N_COL_NML;
  cterms = new double*[ncol];
  for (int K = 0; K < (ncol); K++) {
    cterms[K] = new double[n_samples + 1];
    for (int i = 0; i < (n_samples + 1); i++) {
      if (K == 1)
        cterms[K][i] = 0;
      else if (i == 0)
        cterms[K][i] = 0;
      else
        cterms[K][i] = -1;
    }
  }
  // Initialize the c2 terms
  for (int i = 0; i < (n_samples + 1); i++) {
    computeLogC(i, 2, looklog, cterms);
  }

  lookchoose = new double*[ncol];
  for (int K = 0; K < (ncol); K++) {
    lookchoose[K] = new double[n_samples + 1];
    for (int i = 0; i < (n_samples + 1); i++) {
      lookchoose[K][i] = -1;
    }
  }

  // create the edge structure and keep track of how many searchMore we have
  edges = new Edge*[n_nodes];
  for (int i = 0; i < n_nodes; i++)
    edges[i] = new Edge[n_nodes];

  for (int i = 0; i < n_nodes; i++) {
    for (int j = 0; j < n_nodes; j++) {
      if ((!is_continuous[i] && levels[i] == n_samples) ||
          (!is_continuous[j] && levels[j] == n_samples)) {
        // If a node is discrete with as many levels as there are samples, its
        // information with other nodes is null.
        edges[i][j].status = 0;
        edges[i][j].status_prev = 0;
      } else {
        // Initialise all other edges.
        edges[i][j].status = 1;
        edges[i][j].status_prev = 1;
      }
    }
  }

  for (int i = 0; i < n_nodes; i++) {
    edges[i][i].status = 0;
    edges[i][i].status_prev = 0;
  }

  noise_vec = new double[2 * n_samples];
  for (int i = 0; i < 2 * n_samples; i++) {
    noise_vec[i] =
        std::rand() / ((RAND_MAX + 1u) / MAGNITUDE_TIES) - MAGNITUDE_TIES / 2;
  }

  readBlackbox(as<vector<vector<int>>>(arg_list["black_box"]));
}

void Environment::transformToFactors(int i) {
  // create a dictionary to store the factors of the strings
  std::map<string, int> myMap;
  myMap["NA"] = -1;
  myMap[""] = -1;
  int factor = 0;

  for (int j = 0; j < n_samples; j++) {
    auto it = myMap.find(data[j][i]);
    if (it != myMap.end()) {
      data_numeric[j][i] = it->second;
    } else {
      myMap[data[j][i]] = factor;
      data_numeric[j][i] = factor;
      factor++;
    }
  }
}

void Environment::transformToFactorsContinuous(int i) {
  std::multimap<double, int> myMap;

  vector<double> clmn;
  for (int j = 0; j < n_samples; j++) {
    string entry = data[j][i];
    if (entry.compare("NA") != 0 && entry.compare("") != 0)
      clmn.push_back(atof(entry.c_str()));
  }

  sort(clmn.begin(), clmn.end());

  for (size_t j = 0; j < clmn.size(); j++) {
    myMap.insert(std::pair<double, int>(clmn[j], j));
  }

  for (int j = 0; j < n_samples; ++j) {
    string entry = data[j][i];
    if (entry.compare("NA") != 0 && entry.compare("") != 0) {
      data_numeric[j][i] = myMap.find(atof(entry.c_str()))->second;

      auto iterpair = myMap.equal_range(atof(entry.c_str()));
      for (auto it = iterpair.first; it != iterpair.second; ++it) {
        if (it->second == data_numeric[j][i]) {
          myMap.erase(it);
          break;
        }
      }
    } else {
      data_numeric[j][i] = -1;
    }
  }
}

void Environment::transformToFactorsContinuousIdx(int i) {
  std::map<int, int> myMap;

  int entry;
  for (int j = 0; j < n_samples; j++) {
    entry = data_numeric[j][i];
    if (entry != -1) myMap[entry] = j;
  }

  int j = 0;
  for (auto it = myMap.begin(); it != myMap.end(); ++it) {
    data_numeric_idx[i][j] = it->second;
    ++j;
  }
}

// Set the number of levels for each node (the maximum level of each column)
void Environment::setNumberLevels() {
  for (int i = 0; i < n_nodes; i++) {
    int max = 0;
    for (int j = 0; j < n_samples; j++) {
      if (data_numeric[j][i] > max) max = data_numeric[j][i];
    }
    levels[i] = max + 1;
  }
}

void Environment::readFileType() {
  for (int j = 0; j < n_nodes; j++) {
    if (is_continuous[j]) {
      data_double[j].resize(n_samples);
      for (int i = 0; i < n_samples; i++) {
        if (data[i][j].compare("NA") == 0 || data[i][j].compare("") == 0) {
          data_double[j][i] = std::numeric_limits<double>::quiet_NaN();
        } else {
          data_double[j][i] = atof(data[i][j].c_str());
        }
      }
    }
  }
}

void Environment::readBlackbox(const vector<vector<int>>& node_list) {
  for (const auto& pair : node_list) {
    auto x = pair[0], y = pair[1];
      edges[x][y].status = 0;
      edges[y][x].status = 0;
  }
}

}  // namespace structure
}  // namespace miic
