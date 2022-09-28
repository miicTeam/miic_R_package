#include "proba_orientation.h"

#include <algorithm>  // std::sort
#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>   // std::numeric_limits
#include <numeric>  // std::iota
#include <tuple>    // std::tie

namespace miic {
namespace reconstruction {

using std::vector;

namespace {

constexpr int kRemovalMark = -1;
constexpr double kScoreLowest = std::numeric_limits<double>::lowest();
constexpr double kScoreMax = std::numeric_limits<double>::max();

bool isHead(double score) { return score > 0; }  // proba > 0.5
bool isTail(double score) { return score < 0; }  // proba < 0.5
// Given an edge (X px -- py Y) with probabilities px and py, and the value of
// px(py) is determined, then the value of py(px) can be induced if
// 1. px(py) is likely an arrowhead but latent variable (bi-directed edge) is
//    not allowed.
// 2. px(py) is likely an arrowtail and propagation of orientation is allowed.
bool isInducible(double score, bool latent, bool propagation) {
  return (isHead(score) && !latent) || (isTail(score) && propagation);
}

// Given a Triple t0 (W -- Z -- V) and another Triple t who shares the edge
// (W -- Z) with t0, update the scores of t by those of t0.
// t can be one of the following:
//   (W -- Z -- U), (Z -- W -- U), (U -- W -- Z), (U -- Z -- W)
//
// Input w2z, z2w the scores of t0 (W *-> Z and W <-* Z respectively),
//   where * can be head (< or >) or tail (-)
// Input w2z_index, z2w_index value in [0, 1, 2, 3]. Index of the scores of t
//   that corresponds to w2z, z2w
// Input/Output score ProbaScore of t to be updated
bool updateScore(int w2z_index, int z2w_index, double w2z, double z2w,
    bool inducible, ScoreArray& score) {
  score[w2z_index] = ProbaScore{w2z, true};
  if (inducible) {
    score[z2w_index] = ProbaScore{z2w, true};
  } else if (!score[z2w_index].settled && (z2w_index == 0 || z2w_index == 3)) {
    // When (propagation == true || latent == false), score[0] and score[3]
    // might be set in a previous induction. It's not inducible now, reset them.
    score[z2w_index].value = 0;
  }
  // If both scores of the mid node are settled, the triple t will not be able
  // to update other triples in the list, and will be removed from the list.
  return score[1].settled && score[2].settled;
}

// Input I3 = N * I'(W;V;Z|{ui}) 3-point infomation of a Triple W -- Z -- V
// Input w2z > 0 arrowhead probability score from W to Z
// Return v2z induced head/tail probability score from V to Z such that
//   1 / (1 + exp(-v2z)) = {1 / (1 + exp(-abs(I3)))} * {1 / (1 + exp(-w2z))}
//
// Denote p = 1 / (1 + exp(-abs(I3))).
// If I3 > 0, then p is the probability of non-v-structure W *-> Z --* V, and
// P(v2z) = P(w2z) * p is the probability Pr(v2z is tail).
// If I3 < 0, then p is the probability of v-structure W *-> Z <-* V, and
// P(v2z) = P(w2z) * p is the probability Pr(v2z is head).
double getInducedScore(double I3, double w2z) {
  double s_min, s_max;  // structured binding is only available after C++17
  std::tie(s_min, s_max) = std::minmax(std::fabs(I3), w2z);
  return s_min - std::log1p(std::exp(s_min - s_max) + std::exp(-s_max));
}

// Given X *- (x2z) Z (y2z) -* Y, try to induce y2z from x2z.
// See Proposition 1.ii and 2 of Verny et al., 2017 (Supplementary Text).
void induceScore(
    bool latent, bool propagation, double I3, ScoreArray& score, double& rank) {
  int x2z_index{-1};                 // Index of settled score
  int y2z_index{-1}, z2y_index{-1};  // Indices of scores to be induced
  if (isHead(score[1].value) && score[1].settled) {  // X *-> Z [</-]-* Y
    x2z_index = 1;
    y2z_index = 2;
    z2y_index = 3;
  } else if (isHead(score[2].value) && score[2].settled) {  // Y *-[-/>] Z <-* X
    x2z_index = 2;
    y2z_index = 1;
    z2y_index = 0;
  }
  if (I3 == 0 || x2z_index == -1) return;

  double induced_score = getInducedScore(I3, score[x2z_index].value);
  auto& y2z = score[y2z_index];
  auto& z2y = score[z2y_index];
  if (induced_score <= 0 || std::fabs(y2z.value) >= induced_score) return;

  rank = induced_score;  // successful induction
  if (I3 > 0) {  // non-v-structure, induced_score is the tail score
    y2z.value = -induced_score;  // tail Z [-]-* Y
    if (propagation && !z2y.settled && z2y.value < -y2z.value)
      z2y.value = -y2z.value;  // head Z --[>] Y
  } else {  // I3 < 0, v-structure, induced_score is the head score
    y2z.value = induced_score;  // head Z [<]-* Y
    if (!latent && !z2y.settled && z2y.value > -y2z.value)
      z2y.value = -y2z.value;  // tail Z <-[-] Y
  }
}

}  // anonymous namespace

// Iteratively converge towards partially directed graphs including possible
// latent variables and Propagation/Non-Propagation rules.
// Input triples list of unshielded Triple (X -- Z -- Y)
// Input I3_list the 3-point information N * I'(X;Y;Z|{ui}) of each Triple
// return vector<ProbaArray> Each ProbaArray is bound to an unshielded Triple
vector<ProbaArray> getOriProbasList(const vector<Triple>& triples,
    const vector<double>& I3_list, const vector<int>& is_contextual,
    const vector<int>& is_consequence, bool latent, bool degenerate,
    bool propagation, bool half_v_structure) {
  // A score is a quantity almost proportional to abs(I3). All probabilities can
  // be expressed in the form of 1 / (1 + exp(-score)), and they suffer from
  // loss of numerical precision for high score due to the exponential term.
  // Since the mapping between probability and score is bijective, throughout
  // the orientation process, we will work with scores only.
  //     score   |--> probability
  // (-inf, inf) |-->   [0, 1]
  int n_triples = triples.size();
  // Each ScoreArray contains 4 ProbaScores, and is associated with an
  // unshielded triple: ScoreArray[0, 1, 2, 3] <> X 0 -- 1 Z 2 -- 3 Y.
  // Each ProbaScore{score, settled} is initialized to {0, false}.
  vector<ScoreArray> scores(n_triples);
  // Rank is the highest absolute value of unsettled score for each triple.
  vector<double> rank(n_triples, 0);
  for (int i = 0; i < n_triples; ++i) {
    // Initialize scores of v-structure
    if (I3_list[i] < 0) {
      if (!degenerate) {
        // Pr(X *-> Z) = Pr(Y *-> Z) = (1 + exp(I3)) / (1 + 3 * exp(I3))
        // See Proposition 1.i of Verny et al., 2017 (Supplementary Text)
        rank[i] = -I3_list[i] - log(2) + std::log1p(std::exp(I3_list[i]));
      } else {
        // larger than p without degenerate
        rank[i] = -I3_list[i] + log(3 - 2 * std::exp(I3_list[i]));
      }
      scores[i][1].value = rank[i];  // Pr(X *-> Z)
      scores[i][2].value = rank[i];  // Pr(Y *-> Z)
      if (!latent) {
        // (p, 1-p) <> (score, -score)
        scores[i][0].value = -scores[i][1].value;
        scores[i][3].value = -scores[i][2].value;
      }
    }
    // Initialize scores of triples involving contextual variables
    int X = triples[i][0], Z = triples[i][1], Y = triples[i][2];
    if (is_contextual[X]) {  // X --* Z, X cannot be the child of Z
      scores[i][0] = ProbaScore{kScoreLowest, true};
    }
    if (is_contextual[Z]) {  // X *-- Z --* Y, Z cannot be the child of X or Y
      scores[i][1] = ProbaScore{kScoreLowest, true};
      scores[i][2] = ProbaScore{kScoreLowest, true};
    }
    if (is_contextual[Y]) {  // Z *-- Y, Y cannot be the child of Z
      scores[i][3] = ProbaScore{kScoreLowest, true};
    }
    //
    // Initialize scores of triples involving consequence variables
    //
    if (is_consequence[X]) {  // X <-* Z, X cannot be the parent of Z
      scores[i][0] = ProbaScore{kScoreMax, true};
      if (!latent)
        scores[i][1] = ProbaScore{kScoreLowest, true};
    }
    if (is_consequence[Y]) {  // Z *-> Y, Y cannot be the parent of Z
      scores[i][3] = ProbaScore{kScoreMax, true};
      if (!latent)
        scores[i][2] = ProbaScore{kScoreLowest, true};
    }
    if (is_consequence[Z]) {  // X *-> Z <-* Y, Z cannot be the parent of X or Y
      scores[i][1] = ProbaScore{kScoreMax, true};
      scores[i][2] = ProbaScore{kScoreMax, true};
      if (!latent) {
        scores[i][0] = ProbaScore{kScoreLowest, true};
        scores[i][3] = ProbaScore{kScoreLowest, true};
      }
    }

    // Try to induce score with already (manually) settled arrowhead
    if (isHead(scores[i][1].value) || isHead(scores[i][2].value))
      induceScore(latent, propagation, I3_list[i], scores[i], rank[i]);
  }

  auto compareTriples = [&rank](int a, int b) { return rank[a] > rank[b]; };
  // indices of triples
  vector<int> order(n_triples);
  std::iota(begin(order), end(order), 0);
  while (!order.empty()) {
    // Order triples in decreasing rank
    std::sort(begin(order), end(order), compareTriples);
    int top = order[0];
    if (rank[top] <= 0) break;  // top proba <= 0.5
    order.erase(begin(order));  // Remove the top triple from the list

    int X{-1}, Z{-1}, Y{-1};
    ProbaScore& z2x = scores[top][0];
    ProbaScore& x2z = scores[top][1];
    ProbaScore& y2z = scores[top][2];
    ProbaScore& z2y = scores[top][3];
    // Try updating final score on X -- Z
    if (!x2z.settled) {
      // Conflict between I3 (<= 0, v-structure) and y2z (< 0, non-v-structure),
      // reset x2z (and z2x if unsettled) if half_v_structure is not allowed.
      if (y2z.value < 0 && I3_list[top] <= 0 && !half_v_structure) {
        x2z.value = 0;
        if (!z2x.settled) z2x.value = 0;
      } else {
        X = triples[top][0];
        Z = triples[top][1];
      }
    }
    // Try updating final score on Z -- Y
    if (!y2z.settled) {
      if (x2z.value < 0 && I3_list[top] <= 0 && !half_v_structure) {
        y2z.value = 0;
        if (!z2y.settled) z2y.value = 0;
      } else {
        Z = triples[top][1];
        Y = triples[top][2];
      }
    }
    if (X == -1 && Y == -1) continue;  // No score updated, goto next triple
    // Update scores of all triples in the list sharing edge with the top triple
    for (auto i : order) {
      int w2z_index{-1}, z2w_index{-1};   // shared edge markers of triples[i]
      double w2z_value{0}, z2w_value{0};  // final scores of the shared edge
      if (X != -1) {
        if (triples[i][0] == X && triples[i][1] == Z) {
          w2z_index = 1;
          z2w_index = 0;
        } else if (triples[i][0] == Z && triples[i][1] == X) {
          w2z_index = 0;
          z2w_index = 1;
        } else if (triples[i][2] == X && triples[i][1] == Z) {
          w2z_index = 2;
          z2w_index = 3;
        } else if (triples[i][2] == Z && triples[i][1] == X) {
          w2z_index = 3;
          z2w_index = 2;
        }
      }  // if (X != -1)
      if (w2z_index != -1 && z2w_index != -1) {  // shared edge found
        w2z_value = x2z.value;
        z2w_value = z2x.value;
      } else if (Y != -1) {
        if (triples[i][0] == Y && triples[i][1] == Z) {
          w2z_index = 1;
          z2w_index = 0;
        } else if (triples[i][0] == Z && triples[i][1] == Y) {
          w2z_index = 0;
          z2w_index = 1;
        } else if (triples[i][2] == Y && triples[i][1] == Z) {
          w2z_index = 2;
          z2w_index = 3;
        } else if (triples[i][2] == Z && triples[i][1] == Y) {
          w2z_index = 3;
          z2w_index = 2;
        }
        if (w2z_index != -1 && z2w_index != -1) {  // shared edge found
          w2z_value = y2z.value;
          z2w_value = z2y.value;
        }
      }  // if (Y != -1)
      if (w2z_index != -1 && z2w_index != -1) {  // shared edge found
        const bool inducible = isInducible(w2z_value, latent, propagation);
        const bool remove_triple = updateScore(
            w2z_index, z2w_index, w2z_value, z2w_value, inducible, scores[i]);
        if (remove_triple)
          i = kRemovalMark;
        else
          induceScore(latent, propagation, I3_list[i], scores[i], rank[i]);
      }
    }  // for (auto i : order)
    order.erase(remove(begin(order), end(order), kRemovalMark), end(order));
  }
  // Get proba = 1 / (1 + exp(-score)) from scores
  vector<ProbaArray> probas_list(n_triples);
  for (int i = 0; i < n_triples; ++i) {
    std::transform(begin(scores[i]), end(scores[i]), begin(probas_list[i]),
        [](ProbaScore& score) { return 1.0 / (1.0 + std::exp(-score.value)); });
  }
  return probas_list;
}

}  // namespace reconstruction
}  // namespace miic
