#include "proba_orientation.h"

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>  // std::iota
#include <tuple>  // std::tie

namespace miic {
namespace reconstruction {

using std::vector;
using std::fabs;

namespace {

constexpr double kEps = 1.0e-12;
constexpr double kEpsDiff = 1.0e-12;
constexpr int kRemoveTripleMark = -1;

bool isHead(double proba) { return proba > (0.5 + kEps); }
bool isTail(double proba) { return proba < (0.5 - kEps); }
bool isOriented(double proba) { return isHead(proba) || isTail(proba); }
// Given an edge (X px -- py Y) with probabilities px and py, and the value of
// px(py) is determined, then the value of py(px) can be inferred if
// 1. px(py) is likely an arrowhead but latent variable (bi-directed edge) is
//    not allowed.
// 2. px(py) is likely an arrowtail and propagation of orientation is allowed.
bool isInferable(double proba, bool latent, bool propagation) {
  return (isHead(proba) && !latent) || (isTail(proba) && propagation);
}

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
    bool latent, bool propagation, ProbaArray& final, ProbaArray& current,
    bool& need_propagation, bool& remove_triple) {
  const bool inferable = isInferable(w2z, latent, propagation);

  final[w2z_index] = w2z;
  if (inferable) final[z2w_index] = z2w;
  // If both probabilites of the mid node are already oriented, the triple will
  // not be able to update other triples of lower score in the list.
  if (isOriented(final[1]) && isOriented(final[2])) {
    remove_triple = true;
  } else {
    need_propagation = true;
    // Keep probabilites as the current best, to be later compared with the
    // result of putative propagation.
    current[w2z_index] = w2z;
    if (inferable) current[z2w_index] = z2w;
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
    const ProbaArray& final, ProbaArray& current, double& score,
    double& log_score) {
  double x2z = final[1], y2z = final[2];
  if (I3 > 0) {  // Proposition 2
    // Define score in case of no true propagation below
    score = fmax(
        fmin(current[1], 1 - current[2]), fmin(current[2], 1 - current[1]));
    log_score = log1p(score - 1);
    double p_tail{0.5}, logp_tail;
    if (isHead(x2z)) {  // Try from X *-> Z to Z --* Y
      std::tie(p_tail, logp_tail) = getInducedProbability(I3, x2z);
      // Propagate to Z [-]-* Y if no previously higher putative propagation
      if ((p_tail > 0.5 + kEps) && (current[2] > (1 - p_tail + kEps))) {
        log_score = logp_tail;
        score = p_tail;
        current[2] = 1 - p_tail;
        // Propagate to Z --[>] Y if no previously higher putative propagation
        if (propagation && current[3] < (1 - current[2] - kEps))
          current[3] = 1 - current[2];
      }
    } else if (isHead(y2z)) {  // Try from Z <-* Y to X *-- Z
      std::tie(p_tail, logp_tail) = getInducedProbability(I3, y2z);
      // Propagate to X *-[-] Z if no previously higher putative propagation
      if ((p_tail > 0.5 + kEps) && (current[1] > (1 - p_tail + kEps))) {
        log_score = logp_tail;
        score = p_tail;
        current[1] = 1 - p_tail;
        // Propagate to X [<]-- Z if no previously higher putative propagation
        if (propagation && current[0] < (1 - current[1] - kEps))
          current[0] = 1 - current[1];
      }
    }
  } else if (I3 < 0) {  // Proposition 1.ii
    // define score in case of no true propagation below
    if (fabs(current[1] - current[2]) > kEpsDiff) {
      score = fmin(current[1], current[2]);
      log_score = log1p(score - 1);
    }
    double p_head{0.5}, logp_head;
    if (isHead(x2z)) {  // Try from X *-> Z to Z <-* Y
      std::tie(p_head, logp_head) = getInducedProbability(I3, x2z);
      // Propagate to Z [<]-* Y if no previously higher putative propagation
      if ((p_head > 0.5 + kEps) && current[2] < (p_head - kEps)) {
        log_score = logp_head;
        score = p_head;
        current[2] = p_head;
        // Propagate to Z <-[-] Y if no previously higher putative propagation
        if (!latent && (current[3] > (1 - current[2] + kEps)))
          current[3] = 1 - current[2];
      }
    } else if (isHead(y2z)) {  // Try from Z <-* Y to X *-> Z
      std::tie(p_head, logp_head) = getInducedProbability(I3, y2z);
      // Propagate to X *-[>] Z if no previously higher putative propagation
      if ((p_head > 0.5 + kEps) && (current[1] < (p_head - kEps))) {
        log_score = logp_head;
        score = p_head;
        current[1] = p_head;
        // Propagate to X [-]-> Z if no previously higher putative propagation
        if (!latent && (current[0] > (1 - current[1] + kEps)))
          current[0] = 1 - current[1];
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
    const vector<double>& I3_list, const vector<int>& is_contextual,
    bool latent, bool degenerate, bool propagation, bool half_v_structure) {
  int n_triples = triples.size();
  vector<ProbaArray> probas_final(n_triples);  // to be returned
  for (auto& p_array : probas_final) p_array.fill(0.5);
  // Set initial probability for triples involving contextual variables
  for (int i = 0; i < n_triples; ++i) {
    int X = triples[i][0], Z = triples[i][1], Y = triples[i][2];
    if (is_contextual[X]) {
      // X --> Z, forced orientation
      probas_final[i][1] = 1.0;
      probas_final[i][0] = 0.0;
    }
    if (is_contextual[Z]) {
      // X <-- Z, forced orientation
      probas_final[i][0] = 1.0;
      probas_final[i][1] = 0.0;
      // Z --> Y, forced orientation
      probas_final[i][3] = 1.0;
      probas_final[i][2] = 0.0;
    }
    if (is_contextual[Y]) {
      // Z <-- Y, forced orientation
      probas_final[i][2] = 1.0;
      probas_final[i][3] = 0.0;
    }
  }
  // Keep the current best probas of each triple during the orientation process
  // Best means most extreme w.r.t. (away from) 0.5
  auto probas_current = probas_final;  // copy
  // For a Triple (X, Z, Y), score is initialized as the probability of the
  // arrowhead Pr(X *-> Z)
  vector<double> score(n_triples, 0.5);
  // For dataset with large n_samples, use log score for numerical presicion
  vector<double> log_score(n_triples);
  vector<int> orderTpl(n_triples);
  std::iota(begin(orderTpl), end(orderTpl), 0);
  // Initialize score and log_score
  for (int i = 0; i < n_triples; i++) {
    // Try to propagate any already initialized arrowhead
    if (isHead(probas_final[i][1]) || isHead(probas_final[i][2])) {
      propagate(latent, propagation, I3_list[i], probas_final[i],
          probas_current[i], score[i], log_score[i]);
    } else {
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
        probas_current[i][1] = score[i];  // Pr(X *-> Z)
        probas_current[i][2] = score[i];  // Pr(Y *-> Z)
        if (!latent) {
          probas_current[i][0] = 1 - probas_current[i][1];
          probas_current[i][3] = 1 - probas_current[i][2];
        }
      }
    }
  }

  auto compareTriples = [&log_score, &I3_list](int a, int b) {
    // log scores are non-positive, when the score (proba) is close to 1, i.e.,
    // when the abs(NI3) is large enough, log score can be subnormal.
    if (log_score[a] != log_score[b]) {
      return log_score[a] > log_score[b];
    } else {
      return fabs(I3_list[a]) > fabs(I3_list[b]);
    }
  };

  while (!orderTpl.empty()) {
    // Order triples in decreasing log score
    std::sort(begin(orderTpl), end(orderTpl), compareTriples);
    int max_idx = orderTpl[0];
    if (!(score[max_idx] > 0.5 + kEps)) break;
    orderTpl.erase(begin(orderTpl));  // Remove the triple with max score

    const auto& max_triple = triples[max_idx];
    const auto& max_current = probas_current[max_idx];
    auto& max_final = probas_final[max_idx];

    int X{-1}, Z{-1}, Y{-1};
    // Correspond to ProbaArray[0-3]: *2+, proba of arrowhead from * to +
    double z2x{0.5}, x2z{0.5}, y2z{0.5}, z2y{0.5};
    // if arrowhead/tail on Z (x 0-*1 z 2-3 y) is not already established
    // through an earlier propagation
    if ((fabs(max_final[1] - 0.5) < (fabs(max_current[1] - 0.5) - kEps)) &&
        (half_v_structure || I3_list[max_idx] > 0 ||
            max_final[2] > (0.5 - kEps))) {
      // establish arrowhead/tail final proba on 1 (x 0-*1 z 2-3 y)
      max_final[1] = max_current[1];
      if (isInferable(max_final[1], latent, propagation))
        max_final[0] = max_current[0];

      X = max_triple[0];
      Z = max_triple[1];
      z2x = max_current[0];
      x2z = max_current[1];
    }
    if ((fabs(max_final[2] - 0.5) < (fabs(max_current[2] - 0.5) - kEps)) &&
        (half_v_structure || I3_list[max_idx] > 0 ||
            max_final[1] > (0.5 - kEps))) {
      // establish arrowhead/tail final proba on 2 (x 0-1 z 2*-3 y)
      max_final[2] = max_current[2];
      if (isInferable(max_final[2], latent, propagation))
        max_final[3] = max_current[3];

      Z = max_triple[1];
      Y = max_triple[2];
      y2z = max_current[2];
      z2y = max_current[3];
    }
    // No new final orientation found, goto next triple
    if (X == -1 && Y == -1) continue;
    // Update probas of all remaining triples with the newly found orientation
    for (auto i : orderTpl) {
      bool need_propagation = false;
      bool remove_triple = false;
      if (X != -1) {
        if (triples[i][0] == X && triples[i][1] == Z) {
          updateProba(1, 0, x2z, z2x, latent, propagation, probas_final[i],
              probas_current[i], need_propagation, remove_triple);
        } else if (triples[i][0] == Z && triples[i][1] == X) {
          updateProba(0, 1, x2z, z2x, latent, propagation, probas_final[i],
              probas_current[i], need_propagation, remove_triple);
        } else if (triples[i][2] == X && triples[i][1] == Z) {
          updateProba(2, 3, x2z, z2x, latent, propagation, probas_final[i],
              probas_current[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Z && triples[i][1] == X) {
          updateProba(3, 2, x2z, z2x, latent, propagation, probas_final[i],
              probas_current[i], need_propagation, remove_triple);
        }
        if (need_propagation)
          propagate(latent, propagation, I3_list[i], probas_final[i],
              probas_current[i], score[i], log_score[i]);
      }  // if (X != -1)
      need_propagation = false;
      if (Y != -1) {
        if (triples[i][0] == Y && triples[i][1] == Z) {
          updateProba(1, 0, y2z, z2y, latent, propagation, probas_final[i],
              probas_current[i], need_propagation, remove_triple);
        } else if (triples[i][0] == Z && triples[i][1] == Y) {
          updateProba(0, 1, y2z, z2y, latent, propagation, probas_final[i],
              probas_current[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Y && triples[i][1] == Z) {
          updateProba(2, 3, y2z, z2y, latent, propagation, probas_final[i],
              probas_current[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Z && triples[i][1] == Y) {
          updateProba(3, 2, y2z, z2y, latent, propagation, probas_final[i],
              probas_current[i], need_propagation, remove_triple);
        }
        if (need_propagation)
          propagate(latent, propagation, I3_list[i], probas_final[i],
              probas_current[i], score[i], log_score[i]);
      }  // if (Y != -1)
      if (remove_triple) i = kRemoveTripleMark;
    }  // for (auto i : orderTpl)
    orderTpl.erase(remove(begin(orderTpl), end(orderTpl), kRemoveTripleMark),
        end(orderTpl));
  }

  return probas_final;
}

}  // namespace reconstruction
}  // namespace miic
