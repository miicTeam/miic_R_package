#include "proba_orientation.h"

#include <Rcpp.h>

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>

namespace miic {
namespace reconstruction {

using std::vector;

namespace {

constexpr double kEps = 1.0e-12;
constexpr double kEpsDiff = 1.0e-12;
constexpr int kRemoveTripleMark = -1;
#define _MY_PRINT_ 0

struct TripleComparator {
  const vector<double>& scores;

  TripleComparator(const vector<double>& val_vec) : scores(val_vec) {}

  bool operator()(int i1, int i2) { return scores[i1] > scores[i2]; }
};

bool isHead(double proba) { return proba > (0.5 + kEps); }
bool isTail(double proba) { return proba < (0.5 - kEps); }
bool isOriented(double proba) { return isHead(proba) || isTail(proba); }

// Given a Triple t0 (X -- Z -- Y) with the highest log_score, and another
// Triple t who shares one edge with t0, and whose probabilities can be updated
// by those of t0. Suppose the shared edge is W -- Z (W can be X or Y), in the
// triple t, Z is not necessarily the middle node.
//
// Input w2z, z2w the probabilities of t0 (W -* Z and W *- Z respectively),
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
  if ((!latent && isHead(probas[w2z_index])) ||
      (propagation && isTail(probas[w2z_index])))
    probas[z2w_index] = z2w;
  double v2m = std::max(w2z_index, z2w_index) == 1 ? probas[2] : probas[1];
  double u2m = (w2z_index == 1 || w2z_index == 2) ? w2z : z2w;
  if (isOriented(u2m) && isOriented(v2m)) {
    remove_triple = true;
  } else {
    probas2[w2z_index] = w2z;
    if ((!latent && isHead(probas[w2z_index])) ||
        (propagation && isTail(probas[w2z_index])))
      probas2[z2w_index] = z2w;
    need_propagation = true;
  }
}

// See Proposition 7.ii and 8 of Affeldt & Isambert, 2015
// For an unshielded Triple X -- Z -- Y, given three point mutual infomation
// I3 = I(X;Y;Z|ui) and probability w2z of orientation from one side node W (X
// or Y) to Z, return the induced proba (and log proba) of orientation.
// Suppose the other node is V (Y or X), if I3 > 0, then returned probability is
// Z -> V (Proposition 8), otherwise it is V -> Z (Proposition 7.ii)
std::pair<double, double> getInducedProbability(double I3, double w2z) {
  // Calculate log version to retain precision in case of large abs(I3)
  // FIXME: experimental
  // Only keep w2z * (1.0 / (1 + exp(-abs(I3)))) term, drop 0.5 * (1 - w2z) to
  // reduce the accumulation of orientation error if there is any
  double logp_induced = log1p(w2z - 1) - log1p(exp(-abs(I3)));
  return std::make_pair(expm1(logp_induced) + 1, logp_induced);
}

// See Proposition 7.ii and 8 of Affeldt & Isambert, 2015
// Update score performing putative propagation
void propagate(bool latent, bool propagation, double I3, double& score,
    double& log_score, ProbaArray& probas, ProbaArray& probas2) {
  double p_induced, logp_induced;
  if (I3 > 0) {  // Proposition 8
    // Define score in case of no true propagation below
    score = fmax(
        fmin(probas2[1], 1 - probas2[2]), fmin(probas2[2], 1 - probas2[1]));
    log_score = log1p(score - 1);
    // Try from X -> Z to Z -> Y
    double x2z = probas[1];
    std::tie(p_induced, logp_induced) = getInducedProbability(I3, x2z);
    // Propagate to Z -> Y if (x2z > p_induced > 0.5) and no previously higher
    // putative propagation
    if (isHead(x2z) && isHead(p_induced) &&
        (probas2[2] > (1 - p_induced + kEps))) {
      log_score = logp_induced;
      score = p_induced;
      probas2[2] = 1 - p_induced;
      if (propagation && probas2[3] < (p_induced - kEps))
        probas2[3] = p_induced;
    } else {  // Try from Y -> Z to Z -> X
      double y2z = probas[2];
      std::tie(p_induced, logp_induced) = getInducedProbability(I3, y2z);
      if (isHead(y2z) && isHead(p_induced) &&
          (probas2[1] > (1 - p_induced + kEps))) {
        // Propagate to Z -> X if (y2z > p_induced > 0.5) and no previously
        // higher putative propagation
        log_score = logp_induced;
        score = p_induced;
        probas2[1] = 1 - p_induced;
        if (propagation && probas2[0] < (p_induced - kEps))
          probas2[0] = p_induced;
      }
    }
  } else if (I3 < 0) {  // Proposition 7.ii
    // define score in case of no true propagation below
    if (fabs(probas2[1] - probas2[2]) > kEpsDiff) {
      score = fmin(probas2[1], probas2[2]);
      log_score = log1p(score - 1);
    }
    // Try from X -> Z to Y -> Z
    double x2z = probas[1];
    std::tie(p_induced, logp_induced) = getInducedProbability(I3, x2z);
    // Propagate to Y -> Z if (x2z > p_induced > 0.5) and no previously higher
    // putative propagation, update score that decreased due to < 0 propagation!
    if (isHead(p_induced) && probas2[2] < (p_induced - kEps)) {
      log_score = logp_induced;
      score = p_induced;
      probas2[2] = p_induced;
      if (!latent && (probas2[3] > (1 - p_induced + kEps)))
        probas2[3] = 1 - probas2[2];
    } else {  // Try from Y -> Z to X -> Z
      double y2z = probas[2];
      std::tie(p_induced, logp_induced) = getInducedProbability(I3, y2z);
      if (isHead(p_induced) && (probas2[1] < (p_induced - kEps))) {
        // Propagate to X -> Z if (y2z > p_induced > 0.5) and no previously
        // higher putative propagation, update score that decreased due to < 0
        // propagation!
        log_score = logp_induced;
        score = p_induced;
        probas2[1] = p_induced;
        if (!latent && (probas2[0] > (1 - p_induced + kEps)))
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
// return vector<ProbaArray> Each ProbaArray is bound to a unshielded Triple:
// ProbaArray[i] > 0.5 means an arrowhead (< or >), otherwise it is a tail (-).
// ProbaArray[0]: probability of the orientation X *- Z
// ProbaArray[1]: probability of the orientation X -* Z
// ProbaArray[2]: probability of the orientation Z *- Y
// ProbaArray[3]: probability of the orientation Z -* Y
// where * can be head (< or >) or tail (-)
// X [0]--[1] Z [2]--[3] Y
vector<ProbaArray> getOriProbasList(const vector<Triple>& triples,
    const vector<double>& I3_list, bool latent, bool degenerate,
    bool propagation, bool half_v_structure) {
  int n_triples = triples.size();
  vector<ProbaArray> probas_list(n_triples);  // to be returned
  for (auto& p_array : probas_list) p_array.fill(0.5);

  auto probas_list2 = probas_list;  // copy
  // For a Triple X -- Z -- Y, if I3_list < 0, gives the probability of the
  // orientation Pr(X -> Z) = Pr(Y -> Z) = (1 + exp(I3)) / (1 + 3 * exp(I3))
  // See Proposition 7 of S. Affeldt & H. Isambert, 2015
  vector<double> score(n_triples, 0.5);
  // For dataset with large n_samples, use log score for numerical presicion
  vector<double> log_score(n_triples);
  vector<int> orderTpl(n_triples);
  std::iota(begin(orderTpl), end(orderTpl), 0);
  // Initialize ScoreTpl
  for (int i = 0; i < n_triples; i++) {
    if (I3_list[i] >= 0) {
      log_score[i] = log(score[i]);  // log(0.5)
    } else {
      if (!degenerate) {
        // use log1p and expm1 to accommodate large N(n_samples) case
        log_score[i] = log1p(exp(I3_list[i])) - log1p(3 * exp(I3_list[i]));
        score[i] = expm1(log_score[i]) + 1;
      } else {
        // larger than p without degenerate
        score[i] = (3 - 2 * exp(I3_list[i])) / (3 - exp(I3_list[i]));
      }
      probas_list2[i][1] = score[i];  // Pr(X -> Z)
      probas_list2[i][2] = score[i];  // Pr(Y -> Z)
      if (!latent) {
        probas_list2[i][0] = 1 - probas_list2[i][1];
        probas_list2[i][3] = 1 - probas_list2[i][2];
      }
    }
  }

  // FIXME: Hard cap to avoid infinite loop, make no "real sense"
  int count{0};
  while (++count < 200000 && !orderTpl.empty()) {
    // Order triples in decreasing log score
    std::sort(begin(orderTpl), end(orderTpl), TripleComparator(log_score));
    int max_idx = orderTpl[0];
    double maxscoreTpl = score[max_idx];
    if (!(maxscoreTpl > 0.5 + kEps)) break;
    orderTpl.erase(begin(orderTpl));  // Remove the triple with max score

    auto& max_triple = triples[max_idx];
    auto& max_probas = probas_list[max_idx];
    auto& max_probas2 = probas_list2[max_idx];

#if _MY_PRINT_
    for (size_t i = 0; i < orderTpl.size(); i++) {
      Rprintf("count=%i scoreTpl[orderTpl[%i]=%i]=%g", count, i, orderTpl[i],
          score[orderTpl[i]]);
    }
    Rprintf(" maxTpl=%i  maxscoreTpl=%g\n", max_idx, maxscoreTpl);

    Rprintf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
    for (int i = 0; i < n_triples; i++) {
      Rprintf(
          "!!! %i triples %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n",
          i, triples[i][0], probas_list[i][0], probas_list[i][1], triples[i][1],
          probas_list[i][2], probas_list[i][3], triples[i][2], I3_list[i],
          score[i]);
    }
    Rprintf("\n");
    for (int i = 0; i < n_triples; i++) {
      Rprintf(
          "!!! %i Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n", i,
          triples[i][0], probas_list2[i][0], probas_list2[i][1], triples[i][1],
          probas_list2[i][2], probas_list2[i][3], triples[i][2], I3_list[i],
          score[i]);
    }

    Rprintf(
        "!INTER latent=%i count=%i degenerate=%i maxTpl=%i maxscoreTpl=%g \n",
        latent, count, degenerate, max_idx, maxscoreTpl);

    Rprintf("maxTpl=%i X=%i  %g(%g)--%g(%g)  n2=%i  %g(%g)--%g(%g)  Y=%i \n",
        max_idx, max_triple[0], max_probas[0], max_probas2[0], max_probas[1],
        max_probas2[1], max_triple[1], max_probas[2], max_probas2[2],
        max_probas[3], max_probas2[3], max_triple[2]);
#endif  // _MY_PRINT_

    int X{-1}, Z{-1}, Y{-1};
    // Correspond to ProbaArray[0-3]: *2+, proba of orientation from * to +
    double z2x{0.5}, x2z{0.5}, y2z{0.5}, z2y{0.5};
    // if arrowhead/tail on Z (x 0-*1 z 2-3 y) is not already established
    // through an earlier propagation
    if ((abs(max_probas[1] - 0.5) < (abs(max_probas2[1] - 0.5) - kEps)) &&
        (half_v_structure || I3_list[max_idx] > 0 ||
            max_probas[2] > (0.5 - kEps))) {
      // establish arrowhead/tail final proba on 1 (x 0-*1 z 2-3 y)
      X = max_triple[0];
      Z = max_triple[1];
      z2x = max_probas2[0];
      x2z = max_probas2[1];
      max_probas[1] = max_probas2[1];
      if ((!latent && isHead(max_probas[1])) ||
          (propagation && isTail(max_probas[1]))) {
        // establish arrowhead/tail if no latent or
        // arrowhead final proba on 0 (x 0<-1 z 2-3 y)
        max_probas[0] = max_probas2[0];
      }
    }
    if ((abs(max_probas[2] - 0.5) < (abs(max_probas2[2] - 0.5) - kEps)) &&
        (half_v_structure || I3_list[max_idx] > 0 ||
            max_probas[1] > (0.5 - kEps))) {
      // establish arrowhead/tail final proba on 2 (x 0-1 z 2*-3 y)
      Z = max_triple[1];
      Y = max_triple[2];
      y2z = max_probas2[2];
      z2y = max_probas2[3];
      max_probas[2] = max_probas2[2];
      if ((!latent && isHead(max_probas[2])) ||
          (propagation && isTail(max_probas[2]))) {
        // establish arrowhead/tail if no latent or
        // arrowhead final proba on 3 (x 0-1 z 2->3 y)
        max_probas[3] = max_probas2[3];
      }
    }

    for (auto i : orderTpl) {
      bool need_propagation = false;
      bool remove_triple = false;
      if (X != -1) {
        if (triples[i][0] == X && triples[i][1] == Z) {
          updateProba(1, 0, x2z, z2x, latent, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][0] == Z && triples[i][1] == X) {
          updateProba(0, 1, x2z, z2x, latent, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][2] == X && triples[i][1] == Z) {
          updateProba(2, 3, x2z, z2x, latent, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Z && triples[i][1] == X) {
          updateProba(3, 2, x2z, z2x, latent, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        }
        if (need_propagation)
          propagate(latent, propagation, I3_list[i], score[i], log_score[i],
              probas_list[i], probas_list2[i]);
      }    // if (X != -1)
      need_propagation = false;
      if (Y != -1) {
        if (triples[i][0] == Y && triples[i][1] == Z) {
          updateProba(1, 0, y2z, z2y, latent, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][0] == Z && triples[i][1] == Y) {
          updateProba(0, 1, y2z, z2y, latent, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Y && triples[i][1] == Z) {
          updateProba(2, 3, y2z, z2y, latent, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        } else if (triples[i][2] == Z && triples[i][1] == Y) {
          updateProba(3, 2, y2z, z2y, latent, propagation, probas_list[i],
              probas_list2[i], need_propagation, remove_triple);
        }
        if (need_propagation)
          propagate(latent, propagation, I3_list[i], score[i], log_score[i],
              probas_list[i], probas_list2[i]);
      }
      if (remove_triple) i = kRemoveTripleMark;
    }  // for (auto i : orderTpl)
    orderTpl.erase(remove(begin(orderTpl), end(orderTpl), kRemoveTripleMark),
        end(orderTpl));
  }

#if _MY_PRINT_
  Rprintf("\n\n consist-neg-pos\n");
  Rprintf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
  for (int i = 0; i < n_triples; i++) {
    Rprintf("!!! triples %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n",
        triples[i][0], probas_list[i][0], probas_list[i][1], triples[i][1],
        probas_list[i][2], probas_list[i][3], triples[i][2], I3_list[i],
        score[i]);
  }
  for (int i = 0; i < n_triples; i++) {
    Rprintf("!!! Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n",
        triples[i][0], probas_list2[i][0], probas_list2[i][1], triples[i][1],
        probas_list2[i][2], probas_list2[i][3], triples[i][2], I3_list[i],
        score[i]);
  }

  if (!orderTpl.empty()) {
    int max_idx = orderTpl[0];
    double maxscoreTpl = score[max_idx];
    Rprintf("!END latent=%i count=%i degenerate=%i maxTpl=%i maxscoreTpl=%g \n",
        latent, count, degenerate, max_idx, maxscoreTpl);
  }
#endif  // _MY_PRINT_

  return probas_list;
}

}  // namespace reconstruction
}  // namespace miic
