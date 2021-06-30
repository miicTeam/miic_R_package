#include "proba_orientation.h"

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>  // std::numeric_limits
#include <numeric>  // std::iota
#include <tuple>  // std::tie

namespace miic {
namespace reconstruction {

using std::vector;
using std::fabs;

namespace {

constexpr int kRemoveTripleMark = -1;
constexpr double kScoreLowest = std::numeric_limits<double>::lowest();

bool isHead(double score) { return score > 0; }  // proba > 0.5
bool isTail(double score) { return score < 0; }  // proba < 0.5
bool isOriented(double score) { return isHead(score) || isTail(score); }
// Given an edge (X px -- py Y) with probabilities px and py, and the value of
// px(py) is determined, then the value of py(px) can be inferred if
// 1. px(py) is likely an arrowhead but latent variable (bi-directed edge) is
//    not allowed.
// 2. px(py) is likely an arrowtail and propagation of orientation is allowed.
bool isInducible(double score, bool latent, bool propagation) {
  return (isHead(score) && !latent) || (isTail(score) && propagation);
}

// Given a Triple t0 (X -- Z -- Y) with the highest score, and another
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
void updateScore(int w2z_index, int z2w_index, double w2z, double z2w,
    bool inducible, ProbaArray& best, ProbaArray& current,
    bool& need_induction, bool& remove_triple) {
  best[w2z_index] = w2z;
  if (inducible) best[z2w_index] = z2w;
  // If both probabilites of the mid node are already oriented, the triple will
  // not be able to update other triples of lower score in the list.
  if (isOriented(best[1]) && isOriented(best[2])) {
    remove_triple = true;
  } else {
    need_induction = true;
    // Keep probabilites as the current best, to be later compared with the
    // result of putative propagation.
    current[w2z_index] = w2z;
    if (inducible) current[z2w_index] = z2w;
  }
}

// For an unshielded Triple X -- Z -- Y, given three point mutual infomation
// I3 = I(X;Y;Z|ui) and probability w2z > 0.5 of orientation from one side node
// W (X or Y) to Z being a head, return the induced probability score v2z
// of orientation from the other node V (Y or X) to Z
double getInducedScore(double I3, double w2z) {
  // Denote p = 1.0 / (1 + exp(-abs(I3)))
  // If I3 > 0, then p is the probability of non-v-structure W *-> Z --* V, and
  // v2z = w2z * p is the probability Pr(v2z is tail), as score.
  // if I3 < 0, then p is the probability of v-structure W *-> Z <-* V, and
  // v2z = w2z * p is the probability Pr(v2z is head), as score.
  double s_min, s_max;  // structured binding is only available after C++17
  std::tie(s_min, s_max) = std::minmax(fabs(I3), w2z);
  return s_min - std::log1p(std::exp(s_min - s_max) + std::exp(-s_max));
}

// Given X *- (x2z) Z (y2z) -* Y, try to induce x2z from y2z, or vice versa.
// See Proposition 1.ii and 2 of Verny et al., 2017 (Supplementary Text).
void induceScore(bool latent, bool propagation, double I3,
    const ProbaArray& best, ProbaArray& current, double& max_score) {
  double x2z = best[1], y2z = best[2];
  if (I3 > 0) {  // Proposition 2
    // Define score in case of no update below
    max_score = std::fmax(
        std::fmin(current[1], -current[2]), std::fmin(current[2], -current[1]));
    if (isHead(x2z)) {  // Try from X *-> Z to Z --* Y
      double s_tail = getInducedScore(I3, x2z);
      // Update Z [-]-* Y if current proba is less likely a tail
      if (s_tail > 0 && current[2] > -s_tail) {
        max_score = s_tail;
        current[2] = -s_tail;
        // Update Z --[>] Y if current proba is less likely a head
        if (propagation && current[3] < -current[2])
          current[3] = -current[2];
      }
    } else if (isHead(y2z)) {  // Try from Z <-* Y to X *-- Z
      double s_tail = getInducedScore(I3, y2z);
      // Update X *-[-] Z if current proba is less likely a tail
      if (s_tail > 0 && current[1] > -s_tail) {
        max_score = s_tail;
        current[1] = -s_tail;
        // Update X [<]-- Z if current proba is less likely a head
        if (propagation && current[0] < -current[1])
          current[0] = -current[1];
      }
    }
  } else if (I3 < 0) {  // Proposition 1.ii
    // Define score in case of no update below
    if (fabs(current[1] - current[2]) > 0) {
      max_score = std::fmin(current[1], current[2]);
    }
    if (isHead(x2z)) {  // Try from X *-> Z to Z <-* Y
      double s_head = getInducedScore(I3, x2z);
      // Update Z [<]-* Y if current proba is less likely a head
      if (s_head > 0 && current[2] < s_head) {
        max_score = s_head;
        current[2] = s_head;
        // Update Z <-[-] Y if current proba is less likely a tail
        if (!latent && current[3] > -current[2])
          current[3] = -current[2];
      }
    } else if (isHead(y2z)) {  // Try from Z <-* Y to X *-> Z
      double s_head = getInducedScore(I3, y2z);
      // Update X *-[>] Z if current proba is less likely a head
      if (s_head > 0 && current[1] < s_head) {
        max_score = s_head;
        current[1] = s_head;
        // Update X [-]-> Z if current proba is less likely a tail
        if (!latent && current[0] > -current[1])
          current[0] = -current[1];
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
  // A score is a quantity almost proportional to abs(I3). All probabilities
  // can be expressed in the form of 1 / (1 + exp(-score)), they suffer from
  // loss of numerical precision for high score due to the exponential term.
  // Since the mapping between probability and score is bijective, throughout
  // the orientation process, we will work with scores only.
  //     score   |--> probability
  // (-inf, inf) |-->   [0, 1]
  vector<ProbaArray> scores_best(n_triples);
  for (auto& s_array : scores_best) s_array.fill(0);  // 0 score <-> 0.5 proba
  // Set initial probability for triples involving contextual variables
  for (int i = 0; i < n_triples; ++i) {
    int X = triples[i][0], Z = triples[i][1], Y = triples[i][2];
    if (is_contextual[X]) {  // X --* Z, X cannot be the child of Z
      scores_best[i][0] = kScoreLowest;
    }
    if (is_contextual[Z]) {  // X *-- Z --* Y, Z cannot be the child of X or Y
      scores_best[i][1] = kScoreLowest;
      scores_best[i][2] = kScoreLowest;
    }
    if (is_contextual[Y]) {  // Z *-- Y, Y cannot be the child of Z
      scores_best[i][3] = kScoreLowest;
    }
  }
  // Keep the current scores of each triple during the orientation process
  auto scores_current = scores_best;  // copy
  // Keep the highest absolute value of score of each triple for ordering.
  vector<double> max_score(n_triples, 0);
  for (int i = 0; i < n_triples; i++) {
    // Try to induce score with already established arrowhead
    if (isHead(scores_best[i][1]) || isHead(scores_best[i][2])) {
      induceScore(latent, propagation, I3_list[i], scores_best[i],
          scores_current[i], max_score[i]);
    } else if (I3_list[i] < 0) {
      // Initialize scores of v-structure
      if (!degenerate) {
        // if I3_list < 0 (likely a v-structure),
        // Pr(X *-> Z) = Pr(Y *-> Z) = (1 + exp(I3)) / (1 + 3 * exp(I3))
        // See Proposition 1.i of Verny et al., 2017 (Supplementary Text)
        max_score[i] = -I3_list[i] - log(2) + std::log1p(std::exp(I3_list[i]));
      } else {
        // larger than p without degenerate
        max_score[i] = -I3_list[i] + log(3 - 2 * std::exp(I3_list[i]));
      }
      scores_current[i][1] = max_score[i];  // Pr(X *-> Z)
      scores_current[i][2] = max_score[i];  // Pr(Y *-> Z)
      if (!latent) {
        // (p, 1-p) <-> (score, -score)
        scores_current[i][0] = -scores_current[i][1];
        scores_current[i][3] = -scores_current[i][2];
      }
    }
  }

  auto compareTriples = [&max_score](int a, int b) {
    return max_score[a] > max_score[b];
  };

  vector<int> orderTpl(n_triples);
  std::iota(begin(orderTpl), end(orderTpl), 0);
  while (!orderTpl.empty()) {
    // Order triples in decreasing log score
    std::sort(begin(orderTpl), end(orderTpl), compareTriples);
    int top_idx = orderTpl[0];
    if (max_score[top_idx] <= 0) break;  // top proba <= 0.5
    orderTpl.erase(begin(orderTpl));  // Remove the top triple

    const auto& top_triple = triples[top_idx];
    const auto& top_current = scores_current[top_idx];
    auto& top_best = scores_best[top_idx];

    int X{-1}, Z{-1}, Y{-1};
    // Correspond to ProbaArray[0-3]: *2+, proba of arrowhead from * to +
    double z2x{0}, x2z{0}, y2z{0}, z2y{0};
    // Try updating X 0 -- 1 Z
    if ((fabs(top_best[1]) < fabs(top_current[1])) &&
        (half_v_structure || I3_list[top_idx] > 0 || top_best[2] >= 0)) {
      // establish arrowhead/tail final proba on 1 (x 0-*1 z 2-3 y)
      top_best[1] = top_current[1];
      if (isInducible(top_best[1], latent, propagation))
        top_best[0] = top_current[0];

      X = top_triple[0];
      Z = top_triple[1];
      z2x = top_current[0];
      x2z = top_current[1];
    }
    // Try updating Z 2 -- 3 Y
    if ((fabs(top_best[2]) < fabs(top_current[2])) &&
        (half_v_structure || I3_list[top_idx] > 0 || top_best[1] >= 0)) {
      // establish arrowhead/tail final proba on 2 (x 0-1 z 2*-3 y)
      top_best[2] = top_current[2];
      if (isInducible(top_best[2], latent, propagation))
        top_best[3] = top_current[3];

      Z = top_triple[1];
      Y = top_triple[2];
      y2z = top_current[2];
      z2y = top_current[3];
    }
    // No score updated, goto next triple
    if (X == -1 && Y == -1) continue;
    // Update scores of all remaining triples with scores of the top triple
    for (auto i : orderTpl) {
      bool remove_triple = false;
      if (X != -1) {
        int x2z_index{-1}, z2x_index{-1};  // shared edge (scores) markers
        if (triples[i][0] == X && triples[i][1] == Z) {
          x2z_index = 1;
          z2x_index = 0;
        } else if (triples[i][0] == Z && triples[i][1] == X) {
          x2z_index = 0;
          z2x_index = 1;
        } else if (triples[i][2] == X && triples[i][1] == Z) {
          x2z_index = 2;
          z2x_index = 3;
        } else if (triples[i][2] == Z && triples[i][1] == X) {
          x2z_index = 3;
          z2x_index = 2;
        }
        if (x2z_index != -1 && z2x_index != -1) {  // found shared edge
          const bool inducible = isInducible(x2z, latent, propagation);
          bool need_induction{false};
          updateScore(x2z_index, z2x_index, x2z, z2x, inducible, scores_best[i],
              scores_current[i], need_induction, remove_triple);
          if (need_induction)
            induceScore(latent, propagation, I3_list[i], scores_best[i],
                scores_current[i], max_score[i]);
        }
      }  // if (X != -1)
      if (Y != -1) {
        int y2z_index{-1}, z2y_index{-1};  // shared edge (scores) markers
        if (triples[i][0] == Y && triples[i][1] == Z) {
          y2z_index = 1;
          z2y_index = 0;
        } else if (triples[i][0] == Z && triples[i][1] == Y) {
          y2z_index = 0;
          z2y_index = 1;
        } else if (triples[i][2] == Y && triples[i][1] == Z) {
          y2z_index = 2;
          z2y_index = 3;
        } else if (triples[i][2] == Z && triples[i][1] == Y) {
          y2z_index = 3;
          z2y_index = 2;
        }
        if (y2z_index != -1 && z2y_index != -1) {  // found shared edge
          const bool inducible = isInducible(y2z, latent, propagation);
          bool need_induction{false};
          updateScore(y2z_index, z2y_index, y2z, z2y, inducible, scores_best[i],
              scores_current[i], need_induction, remove_triple);
          if (need_induction)
            induceScore(latent, propagation, I3_list[i], scores_best[i],
                scores_current[i], max_score[i]);
        }
      }  // if (Y != -1)
      if (remove_triple) i = kRemoveTripleMark;
    }  // for (auto i : orderTpl)
    orderTpl.erase(remove(begin(orderTpl), end(orderTpl), kRemoveTripleMark),
        end(orderTpl));
  }

  vector<ProbaArray> probas_list(n_triples);
  for (int i=0; i < n_triples; ++i) {
    std::transform(begin(scores_best[i]), end(scores_best[i]),
        begin(probas_list[i]),
        [](double score) { return 1.0 / (1.0 + std::exp(-score)); });
  }

  return probas_list;
}

}  // namespace reconstruction
}  // namespace miic
