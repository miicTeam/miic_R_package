#include "proba_orientation.h"

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>  // std::iota
#include <tuple>  // std::tie
#include <regex>

#include "debug.h"
#define _DEBUG  1

namespace miic {
namespace reconstruction {

using std::vector;
using std::fabs;

namespace {

constexpr double kEps = 1.0e-12;
constexpr double kEpsDiff = 1.0e-12;
constexpr int kRemoveTripleMark = -1;

struct TripleComparator {
  const vector<double>& scores;

  TripleComparator(const vector<double>& val_vec) : scores(val_vec) {}

  bool operator()(int i1, int i2) { return scores[i1] > scores[i2]; }
};

bool isLikely(double proba) { return proba > (0.5 + kEps); }
bool isUnlikely(double proba) { return proba < (0.5 - kEps); }
bool isOriented(double proba) { return isLikely(proba) || isUnlikely(proba); }

// Given a Triple t0 (X -- Z -- Y) with the highest log_score, and another
// Triple t who shares one edge with t0, and whose probabilities can be updated
// by those of t0. Suppose the shared edge is W -- Z (W can be X or Y), in the
// triple t, Z is not necessarily the middle node.
//
// Input w2z, z2w the probabilities of t0 (W *-> Z and W <-* Z respectively),
// where * can be head (< or >) or tail (-)
// Input/Output probas, probas2, the ProbaArrays of t to be updated
//
// Denote by M (can be W or Z) the middle node of t, by U (can be Z or W) the
// other node in the shared edge, and by V the third node, so t is either (U
// -- M -- V) or (V -- M -- U)
//
// Input w2z_index, z2w_index value in [0, 1, 2, 3]. Index of the probabilites
// in the ProbaArray of t that corresponds to w2z, z2w, If t is (U -- M -- V),
// then (w2z_index, z2w_index) is (0, 1) or (1, 0), if t is (V -- M -- U), then
// (w2z_index, z2w_index) is (2, 3) or (3, 2)
void updateProba(int w2z_index, int z2w_index, double w2z, double z2w,
    bool latent, bool propagation, ProbaArray& probas, ProbaArray& probas2,
    bool& need_propagation, bool& remove_triple) {
  probas[w2z_index] = w2z;
  if ((!latent && isLikely(probas[w2z_index])) ||
      (propagation && isUnlikely(probas[w2z_index])))
    probas[z2w_index] = z2w;
  double v2m = std::max(w2z_index, z2w_index) == 1 ? probas[2] : probas[1];
  double u2m = (w2z_index == 1 || w2z_index == 2) ? w2z : z2w;
  if (isOriented(u2m) && isOriented(v2m)) {
    remove_triple = true;
  } else {
    probas2[w2z_index] = w2z;
    if ((!latent && isLikely(probas[w2z_index])) ||
        (propagation && isUnlikely(probas[w2z_index])))
      probas2[z2w_index] = z2w;
    need_propagation = true;
  }
}

// For an unshielded Triple X -- Z -- Y, given three point mutual infomation
// I3 = I(X;Y;Z|ui) and probability w2z > 0.5 of orientation from one side node
// W (X or Y) to Z being a head, return the induced probability v2z (and log
// proba) of orientation from the other node V (Y or X) to Z
std::pair<double, double> getInducedProbability(double I3, double w2z) {
  // Denote p = 1.0 / (1 + exp(-abs(I3)))
  // If I3 > 0, then p is the probability of non-v-structure W *-> Z --* V, and
  // v2z = w2z * p is the conditional probability Pr(v2z is tail | w2z is head)
  // if I3 < 0, then p is the probability of v-structure W *-> Z <-* V, and
  // v2z = w2z * p is the conditional probability Pr(v2z is head | w2z is head)
  // Calculate log version first to retain precision in case of large abs(I3)
  double log_v2z = log1p(w2z - 1) - log1p(exp(-fabs(I3)));
  return std::make_pair(expm1(log_v2z) + 1, log_v2z);
}

// Update score performing putative propagation
// See Proposition 1.ii and 2 of Verny et al., 2017 (Supplementary Text)
void propagate(bool latent, bool propagation, double I3,
    const ProbaArray& probas, ProbaArray& probas2, double& score,
    double& log_score) {
  double x2z = probas[1], y2z = probas[2];
  if (I3 > 0) {  // Proposition 2
    // Define score in case of no true propagation below
    score = fmax(
        fmin(probas2[1], 1 - probas2[2]), fmin(probas2[2], 1 - probas2[1]));
    log_score = log1p(score - 1);
    double p_tail{0.5}, logp_tail;
    if (isLikely(x2z)) {  // Try from X *-> Z to Z --* Y
      std::tie(p_tail, logp_tail) = getInducedProbability(I3, x2z);
      // Propagate to Z [-]-* Y if no previously higher putative propagation
      if (isLikely(p_tail) && (probas2[2] > (1 - p_tail + kEps))) {
        log_score = logp_tail;
        score = p_tail;
        probas2[2] = 1 - p_tail;
        // Propagate to Z --[>] Y if no previously higher putative propagation
        if (propagation && probas2[3] < (1 - probas2[2] - kEps))
          probas2[3] = 1 - probas2[2];
      }
    } else if (isLikely(y2z)) {  // Try from Z <-* Y to X *-- Z
      std::tie(p_tail, logp_tail) = getInducedProbability(I3, y2z);
      // Propagate to X *-[-] Z if no previously higher putative propagation
      if (isLikely(p_tail) && (probas2[1] > (1 - p_tail + kEps))) {
        log_score = logp_tail;
        score = p_tail;
        probas2[1] = 1 - p_tail;
        // Propagate to X [<]-- Z if no previously higher putative propagation
        if (propagation && probas2[0] < (1 - probas2[1] - kEps))
          probas2[0] = 1 - probas2[1];
      }
    }
  } else if (I3 < 0) {  // Proposition 1.ii
    // define score in case of no true propagation below
    if (fabs(probas2[1] - probas2[2]) > kEpsDiff) {
      score = fmin(probas2[1], probas2[2]);
      log_score = log1p(score - 1);
    }
    double p_head{0.5}, logp_head;
    if (isLikely(x2z)) {  // Try from X *-> Z to Z <-* Y
      std::tie(p_head, logp_head) = getInducedProbability(I3, x2z);
      // Propagate to Z [<]-* Y if no previously higher putative propagation
      if (isLikely(p_head) && probas2[2] < (p_head - kEps)) {
        log_score = logp_head;
        score = p_head;
        probas2[2] = p_head;
        // Propagate to Z <-[-] Y if no previously higher putative propagation
        if (!latent && (probas2[3] > (1 - probas2[2] + kEps)))
          probas2[3] = 1 - probas2[2];
      }
    } else if (isLikely(y2z)) {  // Try from Z <-* Y to X *-> Z
      std::tie(p_head, logp_head) = getInducedProbability(I3, y2z);
      // Propagate to X *-[>] Z if no previously higher putative propagation
      if (isLikely(p_head) && (probas2[1] < (p_head - kEps))) {
        log_score = logp_head;
        score = p_head;
        probas2[1] = p_head;
        // Propagate to X [-]-> Z if no previously higher putative propagation
        if (!latent && (probas2[0] > (1 - probas2[1] + kEps)))
          probas2[0] = 1 - probas2[1];
      }
    }
  }
}

}  // anonymous namespace

// Iteratively converge towards partially oriented graphs including possible
// latent variables and Propagation/Non-Propagation rules.
// param triples list of unshielded Triple (X -- Z -- Y)
// param I3_list the 3-point mutual info (N * I'(X;Y;Z|{ui})) of each Triple
// return vector<ProbaArray> Each ProbaArray is bound to a unshielded Triple
vector<ProbaArray> getOriProbasList(const vector<Triple>& triples,
    const vector<double>& I3_list, int tau, 
    const vector<miic::structure::Node>& nodes, 
    bool latent, bool latent_orientation, bool degenerate,
    bool propagation, bool half_v_structure) {
  int n_triples = triples.size();
  vector<ProbaArray> probas_list(n_triples);  // to be returned
  for (auto& p_array : probas_list) p_array.fill(0.5);
  //
  // In temporal mode, use time for orientation
  //
  if (tau >= 1)
    {
    if (_DEBUG)
      {
      Rcpp::Rcout << "\ngetOriProbasList, triples received:\n";
      for (int i = 0; i < n_triples; i++) 
        {
        Rcpp::Rcout << nodes[ triples[i][0] ].name << "-" 
                    << nodes[ triples[i][1] ].name << "-" 
                    << nodes[ triples[i][2] ].name << " ";
        Rcpp::Rcout << probas_list[i][0] << " " 
                    << probas_list[i][1] << " " 
                    << probas_list[i][2] << " " 
                    << probas_list[i][3] << "\n";
        }
      Rcpp::Rcout << "\nModified triples:\n";
      }
    // When orientation can be deduced from time, The probability for the tail end 
    // depends on latent variable discovery: 0 if not activated, 0.5 when activated
    //
    double proba_tail = 0;
    if (latent || latent_orientation)
      proba_tail = 0.5;
    //
    // Orient each triple using time when possible
    //
    std::regex lag_expr(".*lag");
    for (int i = 0; i < n_triples; i++) 
      {
      int nodeX_lag = stoi (std::regex_replace( nodes[ triples[i][0] ].name, lag_expr, "" ) );
      int nodeZ_lag = stoi (std::regex_replace( nodes[ triples[i][1] ].name, lag_expr, "" ) );
      int nodeY_lag = stoi (std::regex_replace( nodes[ triples[i][2] ].name, lag_expr, "" ) );
                              
      if (nodeX_lag < nodeZ_lag)
        {
        probas_list[i][0] = 1;
        probas_list[i][1] = proba_tail;
        }
      else if (nodeX_lag > nodeZ_lag)
        {
        probas_list[i][0] = proba_tail;
        probas_list[i][1] = 1;
        }
      
      if (nodeZ_lag < nodeY_lag)
        {
        probas_list[i][2] = 1;
        probas_list[i][3] = proba_tail;
        }
      else if (nodeZ_lag > nodeY_lag)
        {
        probas_list[i][2] = proba_tail;
        probas_list[i][3] = 1;
        }
      if (_DEBUG)
        {
        Rcpp::Rcout << nodes[ triples[i][0] ].name << "-" 
                    << nodes[ triples[i][1] ].name << "-" 
                    << nodes[ triples[i][2] ].name << " ";
        Rcpp::Rcout << probas_list[i][0] << " " 
                    << probas_list[i][1] << " " 
                    << probas_list[i][2] << " " 
                    << probas_list[i][3] << "\n";
        }
      }
    }

  auto probas_list2 = probas_list;  // copy
  // For a Triple (X, Z, Y), score is initialized as the probability of the
  // arrowhead Pr(X *-> Z)
  vector<double> score(n_triples, 0.5);
  // For dataset with large n_samples, use log score for numerical presicion
  vector<double> log_score(n_triples);
  vector<int> orderTpl(n_triples);
  std::iota(begin(orderTpl), end(orderTpl), 0);
  // Initialize score and log_score
  for (int i = 0; i < n_triples; i++) {
    if (I3_list[i] >= 0) {
      log_score[i] = log(score[i]);  // log(0.5)
    } else {
      if (!degenerate) {
        // if I3_list < 0 (likely a v-structure),
        // Pr(X *-> Z) = Pr(Y *-> Z) = (1 + exp(I3)) / (1 + 3 * exp(I3))
        // See Proposition 1.i of Verny et al., 2017 (Supplementary Text)
        // use log1p and expm1 to accommodate large N(n_samples) case
        log_score[i] = log1p(exp(I3_list[i])) - log1p(3 * exp(I3_list[i]));
        score[i] = expm1(log_score[i]) + 1;
      } else {
        // larger than p without degenerate
        score[i] = (3 - 2 * exp(I3_list[i])) / (3 - exp(I3_list[i]));
      }
      probas_list2[i][1] = score[i];  // Pr(X *-> Z)
      probas_list2[i][2] = score[i];  // Pr(Y *-> Z)
      if (!latent_orientation) {
        probas_list2[i][0] = 1 - probas_list2[i][1];
        probas_list2[i][3] = 1 - probas_list2[i][2];
      }
    }
  }
  while (!orderTpl.empty()) {
    // Order triples in decreasing log score
    std::sort(begin(orderTpl), end(orderTpl), TripleComparator(log_score));
    int max_idx = orderTpl[0];
    double maxscoreTpl = score[max_idx];
    if (!(maxscoreTpl > 0.5 + kEps)) break;
    orderTpl.erase(begin(orderTpl));  // Remove the triple with max score

    auto& max_triple = triples[max_idx];
    auto& max_probas = probas_list[max_idx];
    auto& max_probas2 = probas_list2[max_idx];

    int X{-1}, Z{-1}, Y{-1};
    // Correspond to ProbaArray[0-3]: *2+, proba of arrowhead from * to +
    double z2x{0.5}, x2z{0.5}, y2z{0.5}, z2y{0.5};
    // if arrowhead/tail on Z (x 0-*1 z 2-3 y) is not already established
    // through an earlier propagation
    if ((fabs(max_probas[1] - 0.5) < (fabs(max_probas2[1] - 0.5) - kEps)) &&
        (half_v_structure || I3_list[max_idx] > 0 ||
            max_probas[2] > (0.5 - kEps))) {
      // establish arrowhead/tail final proba on 1 (x 0-*1 z 2-3 y)
      X = max_triple[0];
      Z = max_triple[1];
      z2x = max_probas2[0];
      x2z = max_probas2[1];
      max_probas[1] = max_probas2[1];
      if ((!latent_orientation && isLikely(max_probas[1])) ||
          (propagation && isUnlikely(max_probas[1]))) {
        // establish arrowhead/tail if no latent_orientation or
        // arrowhead final proba on 0 (x 0<-1 z 2-3 y)
        max_probas[0] = max_probas2[0];
      }
    }
    if ((fabs(max_probas[2] - 0.5) < (fabs(max_probas2[2] - 0.5) - kEps)) &&
        (half_v_structure || I3_list[max_idx] > 0 ||
            max_probas[1] > (0.5 - kEps))) {
      // establish arrowhead/tail final proba on 2 (x 0-1 z 2*-3 y)
      Z = max_triple[1];
      Y = max_triple[2];
      y2z = max_probas2[2];
      z2y = max_probas2[3];
      max_probas[2] = max_probas2[2];
      if ((!latent_orientation && isLikely(max_probas[2])) ||
          (propagation && isUnlikely(max_probas[2]))) {
        // establish arrowhead/tail if no latent_orientation or
        // arrowhead final proba on 3 (x 0-1 z 2->3 y)
        max_probas[3] = max_probas2[3];
      }
    }

    for (auto i : orderTpl) {
      bool need_propagation = false;
      bool remove_triple = false;
      if (X != -1) {
        if (triples[i][0] == X && triples[i][1] == Z) {
          updateProba(1, 0, x2z, z2x, latent_orientation, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][0] == Z && triples[i][1] == X) {
          updateProba(0, 1, x2z, z2x, latent_orientation, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][2] == X && triples[i][1] == Z) {
          updateProba(2, 3, x2z, z2x, latent_orientation, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Z && triples[i][1] == X) {
          updateProba(3, 2, x2z, z2x, latent_orientation, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        }
        if (need_propagation)
          propagate(latent_orientation, propagation, I3_list[i], probas_list[i],
              probas_list2[i], score[i], log_score[i]);
      }    // if (X != -1)
      need_propagation = false;
      if (Y != -1) {
        if (triples[i][0] == Y && triples[i][1] == Z) {
          updateProba(1, 0, y2z, z2y, latent_orientation, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][0] == Z && triples[i][1] == Y) {
          updateProba(0, 1, y2z, z2y, latent_orientation, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Y && triples[i][1] == Z) {
          updateProba(2, 3, y2z, z2y, latent_orientation, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Z && triples[i][1] == Y) {
          updateProba(3, 2, y2z, z2y, latent_orientation, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        }
        if (need_propagation)
          propagate(latent_orientation, propagation, I3_list[i], probas_list[i],
              probas_list2[i], score[i], log_score[i]);
      }
      if (remove_triple) i = kRemoveTripleMark;
    }  // for (auto i : orderTpl)
    orderTpl.erase(remove(begin(orderTpl), end(orderTpl), kRemoveTripleMark),
        end(orderTpl));
  }

  if (_DEBUG)
    {
    Rcpp::Rcout << "\nEnd of fuction getOriProbasList:\n";
    for (int i = 0; i < n_triples; i++) 
      {
      Rcpp::Rcout << nodes[ triples[i][0] ].name << "-" 
                  << nodes[ triples[i][1] ].name << "-" 
                  << nodes[ triples[i][2] ].name << " ";
      Rcpp::Rcout << probas_list[i][0] << " " 
                  << probas_list[i][1] << " " 
                  << probas_list[i][2] << " " 
                  << probas_list[i][3] << "\n";
      }
    }
  return probas_list;
}

}  // namespace reconstruction
}  // namespace miic
