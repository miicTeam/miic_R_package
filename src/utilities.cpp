#include "utilities.h"

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iterator>  // std::ostream_iterator
#include <unordered_map>
#include <sstream>  //std::stringstream
#include <string>
#include <vector>

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "nanoflann.h"

#define KNN_K 5

namespace miic {
namespace utility {

using std::string;
using std::vector;
using namespace miic::computation;
using namespace miic::structure;

namespace {

// Check if the i'th sample of node X and Y and ui in ui_list has no NA
// (-1) in data_numeric
bool SampleHasNoNA(int X, int Y, const vector<int>& ui_list,
    const vector<vector<int>>& data_numeric, int i) {
  return data_numeric[i][X] != -1 && data_numeric[i][Y] != -1 &&
         std::all_of(begin(ui_list), end(ui_list),
             [&data_numeric, i](int u) { return data_numeric[i][u] != -1; });
}

double kl(const TempGrid2d<double>& freqs1, const TempGrid2d<double>& freqs2) {
  double kl_div = 0;
  for (size_t i = 0; i < freqs1.n_rows(); i++) {
    for (size_t j = 0; j < freqs1.n_cols(); j++) {
      double freq1 = freqs1(i, j);
      if (freq1 != 0) kl_div += freq1 * log(freq1 / freqs2(i, j));
    }
  }
  return kl_div;
}

// count1 and count2 are of the same size (= nlevels of the variable).
// The space of count1 is a subspace of that of count2.
// n_samples1&2 are counts of non NA samples.
double kl(const TempVector<int>& count1, const TempVector<int>& count2,
    int n_samples1, int n_samples2) {
  double kl_div = 0;
  for (size_t i = 0; i < count1.size(); ++i) {
    double freq1 = static_cast<double>(count1[i]) / n_samples1;
    double freq2 = static_cast<double>(count2[i]) / n_samples2;
    if (freq1 != 0) {
      kl_div += freq1 * log(freq1 / freq2);
    }
  }
  return kl_div;
}

void getJointMixed(const Environment& environment, int i, int j,
    TempVector<int>& mixedDiscrete, TempVector<double>& mixedContinuous,
    TempVector<int>& curr_sample_is_not_NA) {
  int discrete_pos = environment.is_continuous[i] ? j : i;
  int continuous_pos = environment.is_continuous[i] ? i : j;
  // Fill marginal distributions
  int n_samples_non_na = 0;
  const auto& ui_list = environment.edges(i, j).shared_info->ui_list;
  for (int k = 0; k < environment.n_samples; k++) {
    curr_sample_is_not_NA[k] = 0;
    if (SampleHasNoNA(i, j, ui_list, environment.data_numeric, k)) {
      curr_sample_is_not_NA[k] = 1;
      mixedContinuous[n_samples_non_na] =
          environment.data_double[continuous_pos][k];
      mixedDiscrete[n_samples_non_na] =
          environment.data_numeric[k][discrete_pos];
      n_samples_non_na++;
    }
  }
}

using my_kd_tree_t =
    KDTreeVectorOfVectorsAdaptor<vector<vector<double>>, double>;
// Compute the distance to the kth nearest neigbhour of the given point in the
// given space. The point object is the coordinates of the point with has as
// many dimensions as the space.
//
// Args
// 	-point		: a ndims-dimensional vector containing the coordinates
// of the given point. 	-kdTree		:
// 	-k			: the rank of the nearest neigbhour's distance
// to return
double compute_k_nearest_distance(
    vector<double> point, my_kd_tree_t::index_t* index, int k) {
  TempAllocatorScope scope;

  TempVector<size_t> ret_indexes(k);
  TempVector<double> out_dists_sqr(k);
  nanoflann::KNNResultSet<double> resultSet(k);
  resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

  index->findNeighbors(resultSet, &point[0], nanoflann::SearchParams(10));

  return (sqrt(out_dists_sqr[k - 1]));
}

// Computes the Kullback-Leibler divergence between two joint (2D) distributions
// of real values based on the KNN estimation (F. Perez-Cruz 2004).
// <space 1> is the subsampling of <space 2> after removing NAs.
double compute_kl_divergence_continuous(vector<vector<double>>& space1,
    vector<vector<double>>& space2, int n1, int n2, int ndims, int k,
    const TempVector<bool>& flag_break_ties, const TempVector<int>& map_samples,
    const vector<double>& noise_vec) {
  double sumlog = 0;
  double noise;
  int i_map;
  for (int j = 0; j < ndims; j++) {
    i_map = 0;
    for (int i = 0; i < n2; i++) {
      if (flag_break_ties[j]) {
        noise = noise_vec[(j * n2) + i];
        space2[i][j] += noise;
        if (map_samples[i] == 1) {
          space1[i_map][j] += noise;
          i_map++;
        }
      }
    }
  }
  // construct a kd-tree index:
  // Dimensionality set at run-time (default: L2)
  my_kd_tree_t mat_index1(ndims, space1, 10);
  mat_index1.index->buildIndex();
  my_kd_tree_t mat_index2(ndims, space2, 10);
  mat_index2.index->buildIndex();

  for (int i = 0; i < n1; i++) {
    vector<double> point(ndims);
    for (int j = 0; j < ndims; j++) {
      point[j] = space1[i][j];
    }

    sumlog += log(compute_k_nearest_distance(point, mat_index2.index, k) /
                  compute_k_nearest_distance(point, mat_index1.index, k));
  }

  return ndims * (sumlog / n1) + log(1.0 * (n2 - 1) / (n1 - 1));
}

void getJointSpace(const Environment& environment, int i, int j,
    vector<vector<double>>& jointSpace,
    TempVector<int>& curr_sample_is_not_NA) {
  int n_samples_non_na = 0;
  const auto& ui_list = environment.edges(i, j).shared_info->ui_list;
  for (int k = 0; k < environment.n_samples; k++) {
    curr_sample_is_not_NA[k] = 0;
    if (SampleHasNoNA(i, j, ui_list, environment.data_numeric, k)) {
      curr_sample_is_not_NA[k] = 1;
      jointSpace[n_samples_non_na][0] = environment.data_double[i][k];
      jointSpace[n_samples_non_na][1] = environment.data_double[j][k];
      n_samples_non_na++;
    }
  }
}

TempGrid2d<double> getJointFreqs(const Environment& environment, int i, int j,
    const TempVector<int>& sample_is_not_NA = TempVector<int>()) {
  TempGrid2d<double> joint_freqs(
      environment.levels[i], environment.levels[j], 0);

  int n_samples_non_na = 0;
  const auto& ui_list = environment.edges(i, j).shared_info->ui_list;
  for (int k = 0; k < environment.n_samples; k++) {
    if ((!sample_is_not_NA.empty() && sample_is_not_NA[k]) ||
        (sample_is_not_NA.empty() &&
            SampleHasNoNA(i, j, ui_list, environment.data_numeric, k))) {
      ++joint_freqs(
          environment.data_numeric[k][i], environment.data_numeric[k][j]);
      ++n_samples_non_na;
    }
  }
  for (auto& freq : joint_freqs)
    freq /= n_samples_non_na;

  return joint_freqs;
}

}  // anonymous namespace

vector<vector<int>> getAdjMatrix(const Environment& env) {
  vector<vector<int>> adj_matrix(env.n_nodes, vector<int>(env.n_nodes, 0));
  for (int i = 1; i < env.n_nodes; i++) {
    for (int j = 0; j < i; j++) {
      adj_matrix[i][j] = env.edges(i, j).status;
      adj_matrix[j][i] = env.edges(j, i).status;
    }
  }
  return adj_matrix;
}

string toNameString(const Environment& env, const vector<int>& vec) {
  if (vec.empty()) {
    return "NA";
  } else {
    std::stringstream ss;
    std::transform(vec.begin(), vec.end() - 1,
        std::ostream_iterator<string>(ss, ","),
        [&env](int i) { return env.nodes[i].name; });
    ss << env.nodes[vec.back()].name;
    return ss.str();
  }
}

vector<vector<string>> getEdgesInfoTable(const Environment& env) {
  vector<EdgeID> edge_list;
  for (int i = 0; i < env.n_nodes - 1; i++) {
    for (int j = i + 1; j < env.n_nodes; j++) {
      edge_list.emplace_back(i, j, env.edges(i, j));
    }
  }
  std::sort(edge_list.begin(), edge_list.end());

  vector<vector<string>> table;
  table.emplace_back(std::initializer_list<string>{"x", "y", "z.name",
      "ai.vect", "zi.vect", "Ixy", "Ixy_ai", "cplx", "Rxyz_ai", "category",
      "Nxy_ai", "confidence"});
  for (const auto& edge : edge_list) {
    auto info = edge.getEdge().shared_info;
    double confidence = -1;
    if (info->exp_shuffle != -1)
      confidence = exp(info->cplx - info->Ixy_ui) / info->exp_shuffle;

    using std::to_string;
    table.emplace_back(std::initializer_list<string>{
        env.nodes[edge.X].name,
        env.nodes[edge.Y].name,
        info->top_z == -1 ? "NA" : env.nodes[info->top_z].name,
        toNameString(env, info->ui_list),
        toNameString(env, info->zi_list),
        to_string(info->Ixy),
        to_string(info->Ixy_ui),
        to_string(info->cplx),
        to_string(info->Rxyz_ui),
        to_string(info->connected),
        to_string(info->Nxy_ui),
        to_string(confidence)
    });
  }

  return table;
}

static void chkIntFn(void* dummy) { R_CheckUserInterrupt(); }
bool checkInterrupt() { return (R_ToplevelExec(chkIntFn, NULL) == FALSE); }

void printProgress(double percent, TimePoint start_time, int& percentile_prev) {
  constexpr int bar_length_total{40};
  if (std::isnan(percent) || std::isinf(percent) || percent < 0 || percent > 1)
    return;
  int percentile = static_cast<int>(percent * 100);
  // Only update progress bar if the percentile has changed
  if (percentile == percentile_prev) return;
  percentile_prev = percentile;
  int lpad = static_cast<int>(percent * bar_length_total);
  int rpad = bar_length_total - lpad;
  double sec_elapsed = getLapInterval(start_time);
  double sec_remaining = sec_elapsed / percent * (1 - percent);
  std::stringstream eta;
  if (std::isinf(sec_remaining)) {
    eta << "--";
  } else {
    if (sec_remaining > 60) {
      int minutes = sec_remaining / 60;
      if (minutes > 60) {
        int hours = minutes / 60;
        eta << hours << "h";
      }
      eta << minutes % 60 << "m";
    }
    eta << static_cast<int>(sec_remaining) % 60 << "s";
  }
  string lpad_str = string(bar_length_total, '=');
  string rpad_str = ">" + string(bar_length_total - 1, '-');
  // To stderr
  REprintf("\r[%.*s%.*s] %3d%% eta: %-10s", lpad, lpad_str.c_str(), rpad,
      rpad_str.c_str(), percentile, eta.str().c_str());
  R_FlushConsole();
}

double compute_kl_divergence(int X, int Y, Environment& environment,
    int samplesNotNA, const TempVector<int>& AllLevels_red,
    const TempVector<int>& sample_is_not_NA) {
  TempAllocatorScope scope;

  if (!environment.is_continuous[X] && !environment.is_continuous[Y]) {
    // 1 - XY discrete
    // Joint freqs X,Y after adding the new contributor (Z)
    auto freqs1 = getJointFreqs(environment, X, Y, sample_is_not_NA);
    // Joint freqs X,Y before adding the new contributor (with the current
    // conditioning ui)
    auto freqs2 = getJointFreqs(environment, X, Y);

    return samplesNotNA * kl(freqs1, freqs2);
  } else if (environment.is_continuous[X] && environment.is_continuous[Y]) {
    int current_samplesNotNA = getNumSamplesNonNA(environment, X, Y);
    // 2 - XY continuous
    // Retrieve marginal distibutions with the current conditioning Us
    vector<vector<double>> joint_base(current_samplesNotNA, vector<double>(2));
    TempVector<int> curr_sample_is_not_NA(environment.n_samples);
    getJointSpace(environment, X, Y, joint_base, curr_sample_is_not_NA);

    TempVector<int> map_samples(current_samplesNotNA);
    int i_map = 0;
    for (int i = 0; i < environment.n_samples; i++) {
      if (curr_sample_is_not_NA[i] == 1) {  // sample i is present in X;Y|U
        map_samples[i_map] = 0;
        if (sample_is_not_NA[i] == 1) {  // sample i is also present in X;Y|U,Z
          map_samples[i_map] = 1;
        }
        i_map++;
      }
    }
    vector<vector<double>> joint_non_na(samplesNotNA, vector<double>(2));
    int i_non_na = 0;
    for (int i = 0; i < environment.n_samples; i++) {
      if (sample_is_not_NA[i] == 1) {
        joint_non_na[i_non_na][0] = environment.data_double[X][i];
        joint_non_na[i_non_na][1] = environment.data_double[Y][i];
        i_non_na++;
      }
    }
    TempVector<bool> flag_break_ties(2);
    for (int k = 0; k < 2; k++) {
      flag_break_ties[k] = false || (AllLevels_red[k] != samplesNotNA) ||
                           (AllLevels_red[k] != current_samplesNotNA);
    }

    return samplesNotNA * compute_kl_divergence_continuous(joint_non_na,
                               joint_base, samplesNotNA, current_samplesNotNA,
                               2, KNN_K, flag_break_ties, map_samples,
                               environment.noise_vec);
  } else {
    // 3 - One discrete and one continuous
    int discrete_pos, continuous_pos, continuous_pos_binary;
    if (!environment.is_continuous[X]) {
      discrete_pos = X;
      continuous_pos = Y;
      continuous_pos_binary = 1;
    } else {
      discrete_pos = Y;
      continuous_pos = X;
      continuous_pos_binary = 0;
    }
    int n_discrete_levels = environment.levels[discrete_pos];
    // Full and reduced data may not have the same number of unique levels

    // Retrieve marginal distibutions with the current conditioning Us
    int current_samplesNotNA = getNumSamplesNonNA(environment, X, Y);
    TempVector<int> mixedDiscrete(current_samplesNotNA);
    TempVector<double> mixedContinuous(current_samplesNotNA);
    TempVector<int> curr_sample_is_not_NA(environment.n_samples);
    getJointMixed(environment, X, Y, mixedDiscrete, mixedContinuous,
        curr_sample_is_not_NA);

    // Create count vectors for the discrete variable
    TempVector<int> count_non_na(n_discrete_levels, 0);
    TempVector<int> count_base(n_discrete_levels, 0);
    for (int i = 0; i < current_samplesNotNA; i++) {
      count_base[mixedDiscrete[i]]++;
    }
    for (int i = 0; i < environment.n_samples; i++) {
      // Make sure to use environment data so that the levels match (may be
      // recoded in reduced data)
      if(sample_is_not_NA[i])
        count_non_na[environment.data_numeric[i][discrete_pos]]++;
    }
    double kldiv = 0;
    // Compute the sum count(y) * KL(X_nonNA|y || X|y) over all values of Y y
    for (int level = 0; level < n_discrete_levels; level++) {
      TempAllocatorScope scope;

      TempVector<int> map_level(count_base[level]);
      int i_level = 0;
      for (int i = 0; i < environment.n_samples; i++) {
        if (environment.data_numeric[i][discrete_pos] == level) {
          if (curr_sample_is_not_NA[i]) {
            map_level[i_level] = 0;
            if (sample_is_not_NA[i] == 1) {
              map_level[i_level] = 1;
            }
            i_level++;
          }
        }
      }

      vector<vector<double>> continuous_base(
          count_base[level], vector<double>(1));
      i_level = 0;
      for (int i = 0; i < current_samplesNotNA; i++) {
        if (mixedDiscrete[i] == level) {
          continuous_base[i_level][0] = mixedContinuous[i];
          i_level++;
        }
      }

      vector<vector<double>> continuous_non_na(
          count_non_na[level], vector<double>(1));
      int i_level_non_na = 0;
      for (int i = 0; i < environment.n_samples; i++) {
        if (sample_is_not_NA[i] == 1 &&
            environment.data_numeric[i][discrete_pos] == level) {
          continuous_non_na[i_level_non_na][0] =
              environment.data_double[continuous_pos][i];
          i_level_non_na++;
        }
      }
      TempVector<bool> flag_break_ties{
          (AllLevels_red[continuous_pos_binary] != samplesNotNA) ||
          (AllLevels_red[continuous_pos_binary] != current_samplesNotNA)};

      if (count_non_na[level] > KNN_K) {
        double temp =
            count_non_na[level] *
            compute_kl_divergence_continuous(continuous_non_na, continuous_base,
                count_non_na[level], count_base[level], 1, KNN_K,
                flag_break_ties, map_level, environment.noise_vec);
        if (temp > 0) kldiv += temp;
      }
    }  // level loop
    // add KL(Y(!NA) || Y)
    kldiv += samplesNotNA *
             kl(count_non_na, count_base, samplesNotNA, current_samplesNotNA);
    return kldiv;
  }
}

int getNumSamplesNonNA(const Environment& environment, int i, int j) {
  int n_samples_non_na = 0;
  const auto& ui_list = environment.edges(i, j).shared_info->ui_list;
  for (int k = 0; k < environment.n_samples; k++) {
    if (SampleHasNoNA(i, j, ui_list, environment.data_numeric, k))
      ++n_samples_non_na;
  }
  return n_samples_non_na;
}

// Counts and marks the rows that contain "NA"s for columns X, Y, Z, ui
// \return The number of non NA samples and modifies the vectors
// sample_is_not_na and NAs_count
int countNonNA(int X, int Y, int Z, const vector<int>& ui_list,
    const vector<vector<int>>& data_numeric, TempVector<int>& sample_is_not_na,
    TempVector<int>& NAs_count) {
  int n_samples = data_numeric.size();
  int na_count = 0;
  for (int i = 0; i < n_samples; i++) {
    bool has_na = (Z != -1 && data_numeric[i][Z] == -1) ||
                  !SampleHasNoNA(X, Y, ui_list, data_numeric, i);
    if (has_na) ++na_count;

    sample_is_not_na[i] = !has_na;
    NAs_count[i] = na_count;
  }

  return n_samples - na_count;
}

/**
 * Reduces the full data to its "reduced" form without NA rows given the edge
 * X-Y | U,(Z). In the case of discrete variables, levels are remapped so that
 * the reduced data still has levels that start at 0. For continuous variables,
 * the ranks are updated (\a data_numeric_idx_red).
 *
 * \return A boolean telling whether or not sample weights are needed on the
 * reduced data. More importantly, modifies the vectors passed as argument :
 * \a all_levels_red contains the number of levels of the reduced variables
 * \a is_continuous_red tells whether the reduced variables are continuous
 * \a posArray_red contains the positions of variables in the reduced data
 * \a data_numeric_red contains the discrete levels of the reduced data
 * \a data_numeric_idx_red contains the continuous ranks of the reduced data
 * \a sample_weights_red contains the sample weights of the reduced data
 */
bool filterNA(int X, int Y, int Z, const vector<int>& ui_list,
    const vector<vector<int>>& data_numeric,
    const vector<vector<int>>& data_numeric_idx, const vector<int>& levels,
    const vector<int>& is_continuous, const vector<double>& sample_weights,
    const TempVector<int>& sample_is_not_NA, const TempVector<int>& NAs_count,
    TempGrid2d<int>& data_numeric_red, TempGrid2d<int>& data_numeric_idx_red,
    TempVector<int>& levels_red, TempVector<int>& is_continuous_red,
    TempVector<int>& posArray_red, TempVector<double>& sample_weights_red,
    bool any_na) {
  TempAllocatorScope scope;

  int n_samples = data_numeric.size();
  int n_ui = ui_list.size();
  TempVector<int> posArray(n_ui + 3, -1);
  posArray[0] = X;
  posArray[1] = Y;
  std::copy(begin(ui_list), end(ui_list), begin(posArray) + 2);
  posArray.back() = Z;

  // Map to make sure that the levels of the reduced data start at zero
  std::unordered_map<int, int> new_levels;

  bool flag_sample_weights{false};
  for (int j = 0; j < n_ui + 3; j++) {
    int index = posArray[j];
    if (index == -1) continue;

    posArray_red[j] = j;
    is_continuous_red[j] = is_continuous[index];
    int k1 = 0;         // position in the full data
    int k2 = 0;         // position in the reduced data
    int nnr = 0;        // number of non repeated values
    int prev_val = -1;  // previous value
    int updated_discrete_level = 0;
    new_levels.clear();

    for (int i = 0; i < n_samples; i++) {
      if (!any_na) {
        data_numeric_red(j, i) = data_numeric[i][index];
        if (is_continuous_red[j])
          data_numeric_idx_red(j, i) = data_numeric_idx[index][i];
        if (j == 0)
          sample_weights_red[i] = sample_weights[i];
      } else {
        if (sample_is_not_NA[i] == 1) {
          // Row at index i does not contain any NA
          int old_val = data_numeric[i][index];
          auto it = new_levels.find(old_val);
          if (it == end(new_levels)) {
            // If level has not already been seen add it to the map
            // and increment the number of unique levels in reduced data
            new_levels.insert({old_val, updated_discrete_level});
            data_numeric_red(j, k1) = updated_discrete_level++;
          } else {
            data_numeric_red(j, k1) = it->second;
          }
          if (j == 0) {
            sample_weights_red[k1] = sample_weights[i];
            if (sample_weights_red[k1] != 1.0) flag_sample_weights = true;
          }
          ++k1;
        }
        if (is_continuous_red[j] == 0) continue;
        // Variable j is continuous
        int si = data_numeric_idx[index][i];  // position of ith sample (order)
        if (si == -1 || sample_is_not_NA[si] == 0) continue;
        // Row at position si does not contain any NA, rank is updated
        // taking into account the number of NAs up to si.
        data_numeric_idx_red(j, k2) = si - NAs_count[si];
        ++k2;
        // check whether is a different values or repeated
        if (data_numeric[si][index] != prev_val) {
          ++nnr;
          prev_val = data_numeric[si][index];
        }
      }
    }
    if (!any_na) {
      levels_red[j] = levels[index];
    } else {
      // Update with the effective number of levels
      if (is_continuous_red[j] == 1) {
        if (nnr < 3) is_continuous_red[j] = 0;
        levels_red[j] = nnr;
      } else
        levels_red[j] = updated_discrete_level;
    }
  }
  return flag_sample_weights;
}

// Calculate maximal required memory in bytes for linear allocator
// DO NOT MODIFY unless you know what you are doing
size_t getLinearAllocatorSize(int n_samples, int n_nodes, int maxbins,
    int initbins, const vector<int>& is_continuous, const vector<int>& levels) {
  using std::max;
  constexpr size_t s_int = sizeof(int), s_double = sizeof(double);
  bool all_discrete = std::all_of(
      begin(is_continuous), end(is_continuous), [](int i) { return i == 0; });
  int max_level{0};  // max discrete number of levels
  for (int i = 0; i < n_nodes; ++i) {
    if (!is_continuous[i] && levels[i] > max_level) max_level = levels[i];
  }
  size_t m_discrete =
      s_int * (4 * n_samples + 6 * max_level + max_level * max_level) +
      s_double * max_level;
  size_t m_discretize =
      s_int * (3 * maxbins + (maxbins + 2) * max(initbins + 1, 2 * n_samples) +
                  n_samples) +
      s_double * (3 * maxbins + 2 * n_nodes);
  size_t m_mutual_info = s_int * (36 * n_samples);
  size_t m_continuous = s_int * (3 * n_nodes + 2 * maxbins * n_nodes);
  m_continuous += s_int * (n_samples * n_nodes + 7 + 4 * n_samples + n_nodes) +
                  s_double * 2 * /*STEPMAX1*/ 50;
  m_continuous += max(m_discretize, m_mutual_info);
  size_t m_computation = m_discrete;
  if (!all_discrete) m_computation = max(m_discrete, m_continuous);
  size_t m_get_info =
      s_int * (2 * n_samples + 3 * n_nodes + 2 * n_samples * n_nodes) +
      s_double * n_samples;
  size_t m_utility = s_int * max(max_level * max_level, n_nodes);
  m_utility = max(m_utility,
      s_int * (4 * n_samples + 2 * max_level + 2) + s_double * n_samples);
  size_t m_max = m_get_info + max(m_utility, m_computation);
  return 4096 + m_max;  // Some extra space for fragmentation
}

}  // namespace utility
}  // namespace miic
