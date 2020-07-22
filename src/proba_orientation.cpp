#include "proba_orientation.h"

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

namespace miic {
namespace reconstruction {

using std::vector;
using namespace structure;

constexpr double kEps = 1.0e-12;
constexpr double kEpsDiff = 1.0e-12;
#define _MY_PRINT_ 0

bool isHead(double proba) { return proba > (0.5 + kEps); }
bool isTail(double proba) { return proba < (0.5 - kEps); }
bool isOriented(double proba) { return isHead(proba) || isTail(proba); }

struct TripleComparator {
  const vector<double>& scores;

  TripleComparator(const vector<double>& val_vec) : scores(val_vec) {}

  bool operator()(int i1, int i2) { return scores[i1] > scores[i2]; }
};

double logF2(double scoreTpl, double I3) {
  double scoreN;
  if (scoreTpl < (1 - kEps))
    scoreN = log(-2 / (scoreTpl * (1 / (1 + exp(-fabs(I3))) - 0.5) - 0.5) - 3);
  else
    scoreN = fabs(I3) + log1p(-2 * exp(-fabs(I3)));
  return scoreN;
}

bool updateProba(int z, int w, double w2z, double z2w, bool latent,
    bool propagation, ProbaArray& probas, ProbaArray& probas2) {
  probas[z] = w2z;
  if ((!latent && isHead(probas[z])) || (propagation && isTail(probas[z])))
    probas[w] = z2w;
  // if (z, w) is (0, 1) or (1, 0), otherside2z is 2,
  // if (z, w) is (2, 3) or (3, 2), otherside2z is 1,
  int proba_otherside2z = std::max(z, w) == 1 ? probas[2] : probas[1];
  double proba_2z = (z == 1 || z == 2) ? w2z : z2w;
  if (isOriented(proba_2z) && isOriented(proba_otherside2z)) {
    return false;
  } else {
    probas2[z] = w2z;
    if ((!latent && isHead(probas[z])) || (propagation && isTail(probas[z])))
      probas2[w] = z2w;
    return true;
  }
}

// Iteratively converge towards partially oriented graphs including possible
// latent variables and Propagation/Non-Propagation rules.
// param triples list of unshielded Triple (X -- Z -- Y)
// param I3_list the 3-point mutual info (N * I'(X;Y;Z|{ui})) of each Triple
// return vector<ProbaArray> Each ProbaArray is bound to a unshielded Triple:
// If ProbaArray[i] > 0.5, then it is an arrowhead, otherwise it is a tail.
// ProbaArray[0]: probability that there is an arrowhead X <- Z
// ProbaArray[1]: probability that there is an arrowhead X -> Z
// ProbaArray[2]: probability that there is an arrowhead Z <- Y
// ProbaArray[3]: probability that there is an arrowhead Z -> Y
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
  vector<double> scoreTpl(n_triples, 0.5);
  // For dataset with large n_samples, use log score for numerical presicion
  vector<double> logscoreTpl(n_triples);
  // If I3_list < 0, gives abs(I3_list)
  vector<double> scoresN(n_triples, 0);
  vector<int> orderTpl(n_triples);
  std::iota(begin(orderTpl), end(orderTpl), 0);
  // Initialize ScoreTpl
  for (int i = 0; i < n_triples; i++) {
    if (I3_list[i] >= 0) {
      logscoreTpl[i] = log(scoreTpl[i]);  // log(0.5)
    } else {
      if (!degenerate) {
        // use log1p and expm1 to accommodate large N(n_samples) case
        logscoreTpl[i] = log1p(exp(I3_list[i])) - log1p(3 * exp(I3_list[i]));
        scoreTpl[i] = expm1(logscoreTpl[i]) + 1;
        scoresN[i] = -I3_list[i];
      } else {
        // larger than p without degenerate
        scoreTpl[i] = (3 - 2 * exp(I3_list[i])) / (3 - exp(I3_list[i]));
      }
      probas_list2[i][1] = scoreTpl[i];  // Pr(X -> Z)
      probas_list2[i][2] = scoreTpl[i];  // Pr(Y -> Z)
      if (!latent) {
        probas_list2[i][0] = 1 - probas_list2[i][1];
        probas_list2[i][3] = 1 - probas_list2[i][2];
      }
    }
  }
  // Order triples in decreasing log score
  std::sort(begin(orderTpl), end(orderTpl), TripleComparator(logscoreTpl));

  int count{0};
  int max_idx = orderTpl[0];
  double maxscoreTpl = scoreTpl[max_idx];

#if _MY_PRINT_
  Rprintf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
  for (int i = 0; i < n_triples; i++) {
    Rprintf("!!! triples %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n",
        triples[i][0], probas_list[i][0], probas_list[i][1],
        triples[i][1], probas_list[i][2], probas_list[i][3],
        triples[i][2], I3_list[i], scoreTpl[i]);
  }
  for (int i = 0; i < n_triples; i++) {
    Rprintf("!!! Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n",
        triples[i][0], probas_list2[i][0], probas_list2[i][1],
        triples[i][1], probas_list2[i][2], probas_list2[i][3],
        triples[i][2], I3_list[i], scoreTpl[i]);
  }
  Rprintf("!START latent=%i count=%i degenerate=%i maxTpl=%i maxscoreTpl=%g \n",
      latent, count, degenerate, max_idx, maxscoreTpl);

  for (int i = 0; i < n_triples; i++) {
    Rprintf(" scoreTpl[orderTpl[%i]=%i]=%g  \n", i, orderTpl[i],
        scoreTpl[orderTpl[i]]);
  }
  Rprintf(" maxTpl=%i  maxscoreTpl=%g\n", max_idx, maxscoreTpl);
#endif  // _MY_PRINT_

  do {
    count++;
#if _MY_PRINT_
    Rprintf("! count=%i degenerate=%i maxTpl=%i maxscoreTpl=%g \n", count,
        degenerate, max_idx, maxscoreTpl);
#endif  // _MY_PRINT_

    auto& max_triple = triples[max_idx];
    auto& max_probas = probas_list[max_idx];
    auto& max_probas2 = probas_list2[max_idx];
    logscoreTpl[max_idx] = std::numeric_limits<double>::lowest();
    scoreTpl[max_idx] = -1;
    scoresN[max_idx] = -1;

#if _MY_PRINT_
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
#if _MY_PRINT_
    Rprintf(
        "!X=%i z2x=%g x2z=%g count=%i degenerate=%i maxTpl=%i maxscoreTpl=%g "
        "\n",
        X, z2x, x2z, count, degenerate, max_idx, maxscoreTpl);
#endif  // _MY_PRINT_
    if (X != -1) {
      for (const auto i : orderTpl) {
        if (i == max_idx) continue;
        bool need_update = true;
        if (triples[i][0] == X && triples[i][1] == Z) {
          //need_update = updateProba(1, 0, x2z, z2x, latent, propagation, probas_list[i], probas_list2[i]);
          // sharing fist edge, same symmetry
          probas_list[i][1] = x2z;
          if ((!latent && isHead(probas_list[i][1])) ||
              (propagation && isTail(probas_list[i][1])))
            probas_list[i][0] = z2x;
          if (isOriented(x2z) && isOriented(probas_list[i][2])) {
            logscoreTpl[i] = std::numeric_limits<double>::lowest();
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused triples
            scoresN[i] = -1;
            need_update = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][1] = x2z;
            if ((!latent && isHead(probas_list[i][1])) ||
                (propagation && isTail(probas_list[i][1])))
              probas_list2[i][0] = z2x;
          }
        } else if (triples[i][0] == Z && triples[i][1] == X) {
          //need_update = updateProba(0, 1, x2z, z2x, latent, propagation, probas_list[i], probas_list2[i]);
          // sharing fist edge, antisymmetry
          probas_list[i][0] = x2z;
          if ((!latent && isHead(probas_list[i][0])) ||
              (propagation && isTail(probas_list[i][0])))
            probas_list[i][1] = z2x;
          if (isOriented(z2x) && isOriented(probas_list[i][2])) {
            logscoreTpl[i] = std::numeric_limits<double>::lowest();
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused triples
            scoresN[i] = -1;
            need_update = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][0] = x2z;
            if ((!latent && isHead(probas_list[i][0])) ||
                (propagation && isTail(probas_list[i][0])))
              probas_list2[i][1] = z2x;
          }
        } else if (triples[i][2] == X && triples[i][1] == Z) {
          //need_update = updateProba(2, 3, x2z, z2x, latent, propagation, probas_list[i], probas_list2[i]);
          // sharing second edge, antisymmetry
          probas_list[i][2] = x2z;
          if ((!latent && isHead(probas_list[i][2])) ||
              (propagation && isTail(probas_list[i][2])))
            probas_list[i][3] = z2x;
          if (isOriented(x2z) && isOriented(probas_list[i][1])) {
            logscoreTpl[i] = std::numeric_limits<double>::lowest();
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused triples
            scoresN[i] = -1;
            need_update = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][2] = x2z;
            if ((!latent && isHead(probas_list[i][2])) ||
                (propagation && isTail(probas_list[i][2])))
              probas_list2[i][3] = z2x;
          }
        } else if (triples[i][2] == Z && triples[i][1] == X) {
          //need_update = updateProba(3, 2, x2z, z2x, latent, propagation, probas_list[i], probas_list2[i]);
          // sharing second edge, same symmetry
          probas_list[i][3] = x2z;
          if ((!latent && isHead(probas_list[i][3])) ||
              (propagation && isTail(probas_list[i][3])))
            probas_list[i][2] = z2x;
          if (isOriented(z2x) && isOriented(probas_list[i][1])) {
            logscoreTpl[i] = std::numeric_limits<double>::lowest();
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused triples
            scoresN[i] = -1;
            need_update = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][3] = x2z;
            if ((!latent && isHead(probas_list[i][3])) ||
                (propagation && isTail(probas_list[i][3])))
              probas_list2[i][2] = z2x;
          }
        } else {
          need_update = false;  // non edited tpl
        }

        if (!need_update) continue;
        // update score performing putative propagation
        double p, pp, logpp;
        if (I3_list[i] > 0) {
          // define score in case of no true propagation below
          scoreTpl[i] = fmin(probas_list2[i][1], 1 - probas_list2[i][2]);
          scoreTpl[i] = fmax(
              scoreTpl[i], fmin(probas_list2[i][2], 1 - probas_list2[i][1]));
          logscoreTpl[i] = log1p(scoreTpl[i] - 1);
          scoresN[i] = logF2(scoreTpl[i], I3_list[i]);

          p = probas_list[i][1];
          // only p*(1.0/(1+exp(-I3_list[i])) term
          logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
          pp = expm1(logpp) + 1;
          // 2->3  condition of propagation (p > pp > 0.5) and no previously
          // higher putative propagation
          if (p > (0.5 + kEps) && pp > (0.5 + kEps) &&
              (probas_list2[i][2] > (1 - pp + kEps))) {
            logscoreTpl[i] = logpp;
            scoreTpl[i] = pp;
            scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
            probas_list2[i][2] = 1 - pp;
            if (propagation && probas_list2[i][3] < (pp - kEps))
              probas_list2[i][3] = pp;
          } else {  // other direction? 2->1
            p = probas_list[i][2];
            // only p*(1.0/(1+exp(-I3_list[i])) term
            logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
            pp = expm1(logpp) + 1;
            if (p > (0.5 + kEps) && pp > (0.5 + kEps) &&
                (probas_list2[i][1] > (1 - pp + kEps))) {
              // 2->1   condition of propagation (p > pp > 0.5) and no
              // previously higher putative propagation
              logscoreTpl[i] = logpp;
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
              probas_list2[i][1] = 1 - pp;
              if (propagation && probas_list2[i][0] < (pp - kEps))
                probas_list2[i][0] = pp;
            }
          }
        } else if (I3_list[i] < 0) {
          // define score in case of no true propagation below
          if (fabs(probas_list2[i][1] - probas_list2[i][2]) > kEpsDiff) {
            scoreTpl[i] = fmin(probas_list2[i][1], probas_list2[i][2]);
            logscoreTpl[i] = log1p(scoreTpl[i] - 1);
            scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
          }
          p = probas_list[i][1];
          // only p*(1.0/(1+exp(I3_list[i])) term
          logpp = log1p(p - 1) - log1p(exp(I3_list[i]));
          pp = expm1(logpp) + 1;
          if (pp > (0.5 + kEps) && (probas_list2[i][2] < (pp - kEps))) {
            // 2->3 condition of propagation (p > pp > 0.5) and no previously
            // higher putative propagation
            // update score which has decreased due to < 0 propagation!
            logscoreTpl[i] = logpp;
            scoreTpl[i] = pp;
            scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
            probas_list2[i][2] = pp;
            if (!latent && (probas_list2[i][3] > (1 - pp + kEps)))
              probas_list2[i][3] = 1 - probas_list2[i][2];
          } else {
            p = probas_list[i][2];
            // only p*(1.0/(1+exp(I3_list[i])) term
            logpp = log1p(p - 1) - log1p(exp(I3_list[i]));
            pp = expm1(logpp) + 1;
            if (pp > (0.5 + kEps) && (probas_list2[i][1] < (pp - kEps))) {
              // 2->1 update score which has decreased due to < 0 propagation!
              logscoreTpl[i] = logpp;
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
              probas_list2[i][1] = pp;
              if (!latent && (probas_list2[i][0] > (1 - pp + kEps)))
                probas_list2[i][0] = 1 - probas_list2[i][1];
            }
          }
        }  // if (I3_list[i])
      }  // for
    }  // if (X != -1)

    if (Y != -1) {
      for (const auto i : orderTpl) {
        if (i == max_idx) continue;
        bool need_update = true;
        // sharing fist edge, same symmetry
        if (triples[i][0] == Y && triples[i][1] == Z) {
          //need_update = updateProba(1, 0, y2z, z2y, latent, propagation, probas_list[i], probas_list2[i]);
          probas_list[i][1] = y2z;
          if ((!latent && isHead(probas_list[i][1])) ||
              (propagation && isTail(probas_list[i][1])))
            probas_list[i][0] = z2y;
          if (isOriented(y2z) && isOriented(probas_list[i][2])) {
            logscoreTpl[i] = std::numeric_limits<double>::lowest();
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused triples
            scoresN[i] = -1;
            need_update = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][1] = y2z;
            if ((!latent && isHead(probas_list[i][1])) ||
                (propagation && isTail(probas_list[i][1])))
              probas_list2[i][0] = z2y;
          }
        } else if (triples[i][0] == Z && triples[i][1] == Y) {
          //need_update = updateProba(0, 1, y2z, z2y, latent, propagation, probas_list[i], probas_list2[i]);
          // sharing fist edge, antisymmetry
          probas_list[i][0] = y2z;
          if ((!latent && isHead(probas_list[i][0])) ||
              (propagation && isTail(probas_list[i][0])))
            probas_list[i][1] = z2y;
          if (isOriented(z2y) && isOriented(probas_list[i][2])) {
            logscoreTpl[i] = std::numeric_limits<double>::lowest();
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused triples
            scoresN[i] = -1;
            need_update = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][0] = y2z;
            if ((!latent && isHead(probas_list[i][0])) ||
                (propagation && isTail(probas_list[i][0])))
              probas_list2[i][1] = z2y;
          }
        } else if (triples[i][2] == Y && triples[i][1] == Z) {
          //need_update = updateProba(2, 3, y2z, z2y, latent, propagation, probas_list[i], probas_list2[i]);
          // sharing second edge, antisymmetry
          probas_list[i][2] = y2z;
          if ((!latent && isHead(probas_list[i][2])) ||
              (propagation && isTail(probas_list[i][2])))
            probas_list[i][3] = z2y;
          if (isOriented(y2z) && isOriented(probas_list[i][1])) {
            logscoreTpl[i] = std::numeric_limits<double>::lowest();
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused triples
            scoresN[i] = -1;
            need_update = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][2] = y2z;
            if ((!latent && isHead(probas_list[i][2])) ||
                (propagation && isTail(probas_list[i][2])))
              probas_list2[i][3] = z2y;
          }
        } else if (triples[i][2] == Z && triples[i][1] == Y) {
          //need_update = updateProba(3, 2, y2z, z2y, latent, propagation, probas_list[i], probas_list2[i]);
          // sharing second edge, same symmetry
          probas_list[i][3] = y2z;
          if ((!latent && isHead(probas_list[i][3])) ||
              (propagation && isTail(probas_list[i][3])))
            probas_list[i][2] = z2y;
          if (isOriented(z2y) && isOriented(probas_list[i][1])) {
            logscoreTpl[i] = std::numeric_limits<double>::lowest();
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused triples
            scoresN[i] = -1;
            need_update = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][3] = y2z;
            if ((!latent && isHead(probas_list[i][3])) ||
                (propagation && isTail(probas_list[i][3])))
              probas_list2[i][2] = z2y;
          }
        } else {
          need_update = false;  // non edited tpl
        }

        if (!need_update) continue;
        // update score performing putative propagation
        double p, pp, logpp;
        if (I3_list[i] > 0) {
          // define score in case of no true propagation below
          scoreTpl[i] = fmin(probas_list2[i][1], 1 - probas_list2[i][2]);
          scoreTpl[i] = fmax(
              scoreTpl[i], fmin(probas_list2[i][2], 1 - probas_list2[i][1]));
          logscoreTpl[i] = log1p(scoreTpl[i] - 1);
          scoresN[i] = logF2(scoreTpl[i], I3_list[i]);

          p = probas_list[i][1];
          // only p*(1.0/(1+exp(-I3_list[i])) term
          logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
          pp = expm1(logpp) + 1;
          // 2->3  condition of propagation (p > pp > 0.5) and no previously
          // higher putative propagation
          if (p > (0.5 + kEps) && pp > (0.5 + kEps) &&
              (probas_list2[i][2] > (1 - pp + kEps))) {
            logscoreTpl[i] = logpp;
            scoreTpl[i] = pp;
            scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
            probas_list2[i][2] = 1 - pp;
            if (propagation && (probas_list2[i][3] < (pp - kEps)))
              probas_list2[i][3] = pp;
          } else {  // other direction? 2->1
            p = probas_list[i][2];
            // only p*(1.0/(1+exp(-I3_list[i])) term
            logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
            pp = expm1(logpp) + 1;
            if (p > (0.5 + kEps) && pp > (0.5 + kEps) &&
                (probas_list2[i][1] > (1 - pp + kEps))) {
              // 2->1   condition of propagation (p > pp > 0.5) and no
              // previously higher putative propagation
              logscoreTpl[i] = logpp;
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
              probas_list2[i][1] = 1 - pp;
              if (propagation && (probas_list2[i][0] < (pp - kEps)))
                probas_list2[i][0] = pp;
            }
          }
        } else if (I3_list[i] < 0) {
          // define score in case of no true propagation below
          if (fabs(probas_list2[i][1] - probas_list2[i][2]) > kEpsDiff) {
            scoreTpl[i] = fmin(probas_list2[i][1], probas_list2[i][2]);
            logscoreTpl[i] = log1p(scoreTpl[i] - 1);
            scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
          }
          p = probas_list[i][1];
          // only p*(1.0/(1+exp(I3_list[i])) term
          logpp = log1p(p - 1) - log1p(exp(I3_list[i]));
          pp = expm1(logpp) + 1;
          if (pp > (0.5 + kEps) && (probas_list2[i][2] < (pp - kEps))) {
            // 2->3 condition of propagation (p > pp > 0.5) and no previously
            // higher putative propagation
            // update score which has decreased due to < 0 propagation!
            logscoreTpl[i] = logpp;
            scoreTpl[i] = pp;
            scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
            probas_list2[i][2] = pp;
            if (!latent && (probas_list2[i][3] > (1 - pp + kEps)))
              probas_list2[i][3] = 1 - probas_list2[i][2];
          } else {
            p = probas_list[i][2];
            // only p*(1.0/(1+exp(I3_list[i])) term
            logpp = log1p(p - 1) - log1p(exp(I3_list[i]));
            pp = expm1(logpp) + 1;
            if (pp > (0.5 + kEps) && (probas_list2[i][1] < (pp - kEps))) {
              // 2->1 update score which has decreased due to < 0 propagation!
              logscoreTpl[i] = logpp;
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
              probas_list2[i][1] = pp;
              if (!latent && (probas_list2[i][0] > (1 - pp + kEps)))
                probas_list2[i][0] = 1 - probas_list2[i][1];
            }
          }
        }
      }
    }
    //  Order triples in increasing ScoreTpl of scoresN
    std::sort(begin(orderTpl), end(orderTpl), TripleComparator(logscoreTpl));
    max_idx = orderTpl[0];
    maxscoreTpl = scoreTpl[max_idx];
#if _MY_PRINT_
    for (int i = 0; i < n_triples; i++) {
      Rprintf("count=%i scoreTpl[orderTpl[%i]=%i]=%g", count, i, orderTpl[i],
          scoreTpl[orderTpl[i]]);
    }
    Rprintf(" maxTpl=%i  maxscoreTpl=%g\n", max_idx, maxscoreTpl);

    Rprintf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
    for (int i = 0; i < n_triples; i++) {
      Rprintf("!!! %i triples %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n", i,
          triples[i][0], probas_list[i][0],
          probas_list[i][1], triples[i][1],
          probas_list[i][2], probas_list[i][3],
          triples[i][2], I3_list[i], scoreTpl[i]);
    }
    Rprintf("\n");
    for (int i = 0; i < n_triples; i++) {
      Rprintf("!!! %i Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n", i,
          triples[i][0], probas_list2[i][0],
          probas_list2[i][1], triples[i][1],
          probas_list2[i][2], probas_list2[i][3],
          triples[i][2], I3_list[i], scoreTpl[i]);
    }

    Rprintf(
        "!INTER latent=%i count=%i degenerate=%i maxTpl=%i maxscoreTpl=%g \n",
        latent, count, degenerate, max_idx, maxscoreTpl);
#endif  // _MY_PRINT_
  } while ((maxscoreTpl > (0.5 + kEps)) && count < 200000);

#if _MY_PRINT_
  Rprintf("\n\n consist-neg-pos\n");
  Rprintf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
  for (int i = 0; i < n_triples; i++) {
    Rprintf("!!! triples %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n",
        triples[i][0], probas_list[i][0],
        probas_list[i][1], triples[i][1],
        probas_list[i][2], probas_list[i][3],
        triples[i][2], I3_list[i], scoreTpl[i]);
  }
  for (int i = 0; i < n_triples; i++) {
    Rprintf("!!! Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3_list=%g scoreTpl=%g\n",
        triples[i][0], probas_list2[i][0],
        probas_list2[i][1], triples[i][1],
        probas_list2[i][2], probas_list2[i][3],
        triples[i][2], I3_list[i], scoreTpl[i]);
  }

  Rprintf("!END latent=%i count=%i degenerate=%i maxTpl=%i maxscoreTpl=%g \n",
      latent, count, degenerate, max_idx, maxscoreTpl);
#endif  // _MY_PRINT_

  return probas_list;
}

}  // namespace reconstruction
}  // namespace miic
