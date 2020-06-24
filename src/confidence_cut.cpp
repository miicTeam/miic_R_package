#include <math.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "compute_ens_information.h"
#include "reconstruct.h"
#include "structure.h"
#include "utilities.h"

using std::string;
using std::vector;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;

void shuffle_lookup(int* array, int* array2, size_t n) {
  if (n <= 1) return;
  for (size_t i = 0; i < n - 1; i++) {
    size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
    int t = array[j];
    array[j] = array[i];
    array[i] = t;
    array2[t] = i;
  }
}

bool SortFunctionNoMore2(
    const EdgeID* a, const EdgeID* b, const Environment& environment) {
  return (environment.edges[a->i][a->j].shared_info->Ixy_ui >
          environment.edges[b->i][b->j].shared_info->Ixy_ui);
}

class sorterNoMore2 {
  Environment& environment;

 public:
  sorterNoMore2(Environment& env) : environment(env) {}
  bool operator()(EdgeID const* o1, EdgeID const* o2) const {
    return SortFunctionNoMore2(o1, o2, environment);
  }
};

vector<vector<string> > miic::reconstruction::confidenceCut(
    Environment& environment) {
  int** safe_state;
  int** safe_stateIdx;

  int* lookup = new int[environment.n_samples];
  int* lookup2 = new int[environment.n_samples];
  for (int i = 0; i < environment.n_samples; i++) lookup[i] = i;

  double noMore = environment.numNoMore;
  // Allocate the true edges table
  int** inferredEdges_tab;
  inferredEdges_tab = new int*[environment.numNoMore];
  for (int i = 0; i < environment.numNoMore; i++)
    inferredEdges_tab[i] = new int[2];

  int pos = 0;

  for (int i = 0; i < environment.n_nodes - 1; i++) {
    for (int j = i + 1; j < environment.n_nodes; j++) {
      if (environment.edges[i][j].status) {
        inferredEdges_tab[pos][0] = i;
        inferredEdges_tab[pos][1] = j;
        pos++;
      }
    }
  }

  // Create a back up of the data, for later randomization
  safe_state = new int*[environment.n_samples];
  for (int i = 0; i < environment.n_samples; i++)
    safe_state[i] = new int[environment.n_nodes];

  auto any_continuous = std::any_of(environment.is_continuous.begin(),
      environment.is_continuous.end(), [](int i) { return i == 1; });
  if (any_continuous) {
    safe_stateIdx = new int*[environment.n_nodes];
    for (int i = 0; i < environment.n_nodes; i++)
      safe_stateIdx[i] = new int[environment.n_samples];
  }
  // copy to safe state
  for (int i = 0; i < environment.n_samples; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      safe_state[i][j] = environment.data_numeric[i][j];
      if (environment.is_continuous[j])
        safe_stateIdx[j][i] = environment.data_numeric_idx[j][i];
    }
  }

  int* nodes_toShf = new int[environment.n_nodes];
  for (int i = 0; i < environment.n_nodes; i++) nodes_toShf[i] = 0;
  // indexes of nodes to shuffle
  for (int i = 0; i < environment.numNoMore; i++) {
    nodes_toShf[inferredEdges_tab[i][0]] = 1;
  }

  double* confVect = new double[environment.numNoMore];
  for (int nb = 0; nb < environment.numNoMore; nb++) confVect[nb] = 0;

  int* ptrVarIdx = new int[2];

  // loop on the number of shuffling
  for (int nb = 1; nb <= environment.n_shuffles; nb++) {
    // Shuffle the dataset only for the variables present in nodes_toShf
    for (int col = 0; col < environment.n_nodes; col++) {
      if (nodes_toShf[col] == 1) {
        int row2 = 0;
        shuffle_lookup(lookup, lookup2, environment.n_samples);
        if (environment.is_continuous[col]) {
          for (int i = 0; i < environment.n_samples; i++) {
            lookup2[i] = i;
          }

          sort2arraysConfidence(environment.n_samples, lookup, lookup2);
        }
        for (int row = 0; row < environment.n_samples; row++) {
          environment.data_numeric[row][col] = safe_state[lookup[row]][col];
        }

        for (int row = 0; row < environment.n_samples; row++) {
          if (environment.is_continuous[col]) {
            if (environment.data_numeric[lookup2[row]][col] != -1) {
              environment.data_numeric_idx[col][row2] = lookup2[row];
              row2++;
            }
          }
        }
        if (environment.is_continuous[col]) {
          while (row2 < environment.n_samples) {
            environment.data_numeric_idx[col][row2] = -1;
            row2++;
          }
        }
      }
    }

    for (int i = 0; i < environment.n_samples; i++) {
      for (int j = 0; j < environment.n_nodes; j++) {
        environment.oneLineMatrix[j * environment.n_samples + i] =
            environment.data_numeric[i][j];
      }
    }

    int X, Y;
    double NIxy_ui, k_xy_ui;
    // evaluate the mutual information for every edge
    for (int i = 0; i < environment.numNoMore; i++) {
      X = inferredEdges_tab[i][0];
      Y = inferredEdges_tab[i][1];

      double* res;
      if (!environment.is_continuous[X] &&
          !environment.is_continuous[Y]) {
        // discrete case
        res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, X, Y,
            environment.cplx, environment.m);
        NIxy_ui = res[1];
        k_xy_ui = res[2];
        free(res);
      } else {
        // mixed case
        res = computeEnsInformationContinuous(environment, NULL, 0, NULL, 0, -1,
            X, Y, environment.cplx, environment.m);
        NIxy_ui = res[1];
        k_xy_ui = res[2];
        free(res);
      }

      double ni = NIxy_ui - k_xy_ui;
      if (ni <= 0) {
        ni = 0;
      }
      confVect[i] += exp(-ni);
    }
  }
  // evaluate the average confidence
  for (int nb = 0; nb < environment.numNoMore; nb++) {
    confVect[nb] /= environment.n_shuffles;
  }
  // put values > 1 to 1
  for (int nb = 0; nb < environment.numNoMore; nb++)
    if (confVect[nb] > 1) confVect[nb] = 1;
  // remove edges based on confidence cut
  double confidence;
  vector<int> toDelete;

  for (int i = 0; i < environment.numNoMore; i++) {
    int X = inferredEdges_tab[i][0];
    int Y = inferredEdges_tab[i][1];
    confidence = exp(-(environment.edges[X][Y].shared_info->Ixy_ui -
                       environment.edges[X][Y].shared_info->cplx));
    confVect[i] = confidence / confVect[i];
    if (confVect[i] > environment.conf_threshold) {
      environment.edges[X][Y].shared_info->connected = 0;
      environment.edges[X][Y].status = 0;
      environment.edges[Y][X].status = 0;
      toDelete.push_back(i);
    }
  }
  std::cout << "# -- number of edges cut: " << toDelete.size() << "\n";
  // Delete from vector
  environment.connected_list.clear();
  for (int i = 0; i < environment.numNoMore; i++) {
    if (!(std::find(toDelete.begin(), toDelete.end(), i) != toDelete.end())) {
      int X = inferredEdges_tab[i][0];
      int Y = inferredEdges_tab[i][1];
      environment.connected_list.emplace_back(new EdgeID(X, Y));
    }
  }

  for (int X = 0; X < environment.n_nodes - 1; X++) {
    for (int Y = X + 1; Y < environment.n_nodes; Y++) {
      if (environment.edges[X][Y].status == -2 ||
          environment.edges[X][Y].status == 2 ||
          environment.edges[X][Y].status == 6) {
        environment.edges[X][Y].status = 1;
        environment.edges[Y][X].status = 1;
      }
    }
  }
  // Copy data back
  for (int i = 0; i < environment.n_samples; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      environment.data_numeric[i][j] = safe_state[i][j];
    }
  }

  for (int i = 0; i < environment.n_samples; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      environment.oneLineMatrix[j * environment.n_samples + i] =
          environment.data_numeric[i][j];
    }
  }

  if (any_continuous) {
    // Create the data matrix for factors indexes
    for (int i = 0; i < environment.n_nodes; i++) {
      for (int j = 0; j < environment.n_samples; j++)
        environment.data_numeric_idx[i][j] = -1;
    }

    for (int j = 0; j < environment.n_nodes; j++) {
      if (environment.is_continuous[j]) {
        transformToFactorsContinuous(environment, j);
        transformToFactorsContinuousIdx(environment, j);
        transformToFactors(environment, j);
      }
    }
  }  // End copy data back

  std::sort(environment.connected_list.begin(),
      environment.connected_list.end(), sorterNoMore2(environment));
  environment.numNoMore = environment.connected_list.size();

  delete[] ptrVarIdx;

  for (int i = 0; i < environment.n_samples; i++) delete safe_state[i];
  delete[] safe_state;

  delete[] lookup;
  delete[] lookup2;
  delete[] nodes_toShf;

  vector<vector<string> > confVect1;
  confVect1.emplace_back(
      std::initializer_list<string>{"x", "y", "confidence_ratio"});
  for (int i = 0; i < noMore; i++) {
    confVect1.emplace_back(std::initializer_list<string>{
        environment.nodes[inferredEdges_tab[i][0]].name,
        environment.nodes[inferredEdges_tab[i][1]].name,
        std::to_string(confVect[i])});
  }

  for (int i = 0; i < noMore; i++) delete inferredEdges_tab[i];
  delete[] inferredEdges_tab;
  delete[] confVect;

  return confVect1;
}
