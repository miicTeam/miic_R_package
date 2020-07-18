#include "proba_orientation.h"

#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

namespace miic {
namespace reconstruction {

using std::vector;
using namespace structure;

constexpr double eps = 1.0e-12;
constexpr double eps_diff = 1.0e-12;
#define _MY_PRINT_ 0

struct TripletComparator {
  const vector<double>& scores;

  TripletComparator(const vector<double>& val_vec) : scores(val_vec) {}

  bool operator()(int i1, int i2) { return scores[i1] < scores[i2]; }
};

double logF2(double scoreTpl, double I3) {
  double scoreN;
  if (scoreTpl < (1 - eps))
    scoreN = log(-2 / (scoreTpl * (1 / (1 + exp(-fabs(I3))) - 0.5) - 0.5) - 3);
  else
    scoreN = fabs(I3) + log1p(-2 * exp(-fabs(I3)));
  return scoreN;
}

// Iteratively converge towards partially oriented graphs including possible
// latent variables and Propagation/Non-Propagation rules.
// param triples list of unshielded Triple (X -- Z -- Y)
// param I3_list the I3 (I'(X;Y;Z|{ui})) of each Triple
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

  int i, j, maxTpl, count = 0;
  bool ok;
  double maxscoreTpl;
  double p, pp;
  double logpp;

  vector<double> scoreTpl(n_triples);
  vector<double> logscoreTpl(n_triples);
  vector<double> orderTpl(n_triples);
  vector<double> scoresN(n_triples);

  // Initialization ScoreTpl
  for (i = 0; i < n_triples; i++) {
    if (I3_list[i] < 0) {
      if (!degenerate) {
        logscoreTpl[i] = log1p(exp(I3_list[i])) - log1p(3 * exp(I3_list[i]));
        scoreTpl[i] = expm1(logscoreTpl[i]) + 1;
        scoresN[i] = -(I3_list[i]);
      } else {
        // larger than p without deg
        scoreTpl[i] = (3 - 2 * exp(I3_list[i])) / (3 - exp(I3_list[i]));
      }
      probas_list2[i][1] = scoreTpl[i];
      probas_list2[i][2] = scoreTpl[i];
      if (!latent) {
        probas_list2[i][0] = 1 - probas_list2[i][1];
        probas_list2[i][3] = 1 - probas_list2[i][2];
      }
    } else {  // I3[i] non-negative
      logscoreTpl[i] = log(0.5);
      scoreTpl[i] = 0.5;
      scoresN[i] = 0;
    }
  }

#if _MY_PRINT_
  Rprintf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
  for (i = 0; i < NbTpl; i++) {
    Rprintf("!!! Tpl %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n",
        Tpl[i][0], ProbArrowhead[i][0], ProbArrowhead[i][1],
        Tpl[i][1], ProbArrowhead[i][2], ProbArrowhead[i][3],
        Tpl[i][2], I3[i], scoreTpl[i]);
  }
  for (i = 0; i < NbTpl; i++) {
    Rprintf("!!! Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n",
        Tpl[i][0], ProbArrowhead2[i][0], ProbArrowhead2[i][1],
        Tpl[i][1], ProbArrowhead2[i][2], ProbArrowhead2[i][3],
        Tpl[i][2], I3[i], scoreTpl[i]);
  }

  Rprintf("!START latent=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", latent,
      count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

  // Order Tpl in increasing ScoreTpl
  std::iota(begin(orderTpl), end(orderTpl), 0);
  std::sort(begin(orderTpl), end(orderTpl), TripletComparator(logscoreTpl));

  maxTpl = orderTpl[n_triples - 1];
  maxscoreTpl = scoreTpl[maxTpl];

#if _MY_PRINT_
  for (i = 0; i < NbTpl; i++) {
    Rprintf(" scoreTpl[orderTpl[%i]=%i]=%g  \n", i, orderTpl[i],
        scoreTpl[orderTpl[i]]);
  }
  Rprintf(" maxTpl=%i  maxscoreTpl=%g\n", maxTpl, maxscoreTpl);
#endif

  do {
    count++;
#if _MY_PRINT_
    Rprintf("! count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", count, deg, maxTpl,
        maxscoreTpl);
#endif  // _MY_PRINT_

    // assign topTpl
    int X{-1}, Z{-1}, Y{-1};
    double p1{0.5}, p2{0.5}, p3{0.5}, p4{0.5};

    i = maxTpl;
    logscoreTpl[i]=-__DBL_MAX__;
    scoreTpl[i] = -1;
    scoresN[i] = -1;

#if _MY_PRINT_
    Rprintf("maxTpl=%i n1=%i  %g(%g)--%g(%g)  n2=%i  %g(%g)--%g(%g)  n3=%i \n",
        maxTpl, Tpl[i][0], ProbArrowhead[i][0],
        ProbArrowhead2[i][0], ProbArrowhead[i][1],
        ProbArrowhead2[i][1], Tpl[i][1],
        ProbArrowhead[i][2], ProbArrowhead2[i][2],
        ProbArrowhead[i][3], ProbArrowhead2[i][3],
        Tpl[i][2]);
#endif  // _MY_PRINT_

    // arrowhead proba
    p = fmax(probas_list2[i][1], 1 - probas_list2[i][1]);
    // if arrowhead/tail on 1 (x 0-*1 z 2-3 y) is not already established
    // through an earlier propagation
    if ((probas_list[i][1] < (p - eps)) &&
        (probas_list[i][1] > (1 - p + eps)) &&
        (half_v_structure || I3_list[i] > 0 ||
            probas_list[i][2] > (0.5 - eps))) {
      // establish arrowhead/tail final proba on 1 (x 0-*1 z 2-3 y)
      probas_list[i][1] = probas_list2[i][1];

      X = triples[i][0];
      Z = triples[i][1];
      // not sure this is useful;) //change 20151023
      p1 = probas_list2[i][0];
      // if no latent or tail on 1  (x 0<-1 z 2-3 y) // change 20151022 herve
      if ((!latent && probas_list[i][1] > (0.5 + eps)) ||
          (propagation && probas_list[i][1] < (0.5 - eps))) {
        // establish arrowhead/tail if no latent or
        // arrowhead final proba on 0 (x 0<-1 z 2-3 y)
        probas_list[i][0] = probas_list2[i][0];
      }
      p2 = probas_list[i][1];
    }

    p = fmax(probas_list2[i][2], 1 - probas_list2[i][2]);

    if ((probas_list[i][2] < (p - eps)) &&
        (probas_list[i][2] > (1 - p + eps)) &&
        (half_v_structure || I3_list[i] > 0 ||
            probas_list[i][1] > (0.5 - eps))) {
      // establish arrowhead/tail final proba on 2 (x 0-1 z 2*-3 y)
      probas_list[i][2] = probas_list2[i][2];

      Z = triples[i][1];
      Y = triples[i][2];
      p3 = probas_list[i][2];
      // not sure this is useful;) //change 20151023
      p4 = probas_list2[i][3];
      // if no latent or tail on 2  (x 0-1 z 2->3 y) // change 20151022 herve
      if ((!latent && probas_list[i][2] > (0.5 + eps)) ||
          (propagation && probas_list[i][2] < (0.5 - eps))) {
        // establish arrowhead/tail if no latent or
        // arrowhead final proba on 3 (x 0-1 z 2->3 y)
        probas_list[i][3] = probas_list2[i][3];
      }
    }

#if _MY_PRINT_
    Rprintf("!n1=%i p1=%g p2=%g count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", n1,
        p1, p2, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

    maxscoreTpl = 0.0;

    if (X != -1) {
      // only loop on relevant Tpl i != maxTpl, thanks to the sorting
      // of their scoreTpl
      for (j = n_triples - 2; j >= 0; j--) {
        i = orderTpl[j];
        ok = true;
        // sharing fist edge, same symmetry
        if (triples[i][0] == X && triples[i][1] == Z) {
          probas_list[i][1] = p2;
          // if(!latent) ProbArrowhead[0*NbTpl+i] = p1;
          if ((!latent && probas_list[i][1] > (0.5 + eps)) ||
              (propagation && probas_list[i][1] < (0.5 - eps)))
            probas_list[i][0] = p1;
          // if p2!=0.5 && other edge already oriented // change 20151101
          if (((p2 > (0.5 + eps)) || (p2 < (0.5 - eps))) &&
              ((probas_list[i][2] > (0.5 + eps)) ||
                  (probas_list[i][2] < (0.5 - eps)))) {
            logscoreTpl[i] = -__DBL_MAX__;
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = false;
          } else {  // copy final proba on putative proba too
            if ((!latent && probas_list[i][1] > (0.5 + eps)) ||
                (propagation && probas_list[i][1] < (0.5 - eps))) {
              probas_list2[i][0] = p1;
            }
            probas_list2[i][1] = p2;
          }
        } else if (triples[i][0] == Z && triples[i][1] == X) {
          // sharing fist edge, antisymmetry
          probas_list[i][0] = p2;
          if ((!latent && probas_list[i][0] > (0.5 + eps)) ||
              (propagation && probas_list[i][0] < (0.5 - eps)))
            probas_list[i][1] = p1;
          // if other edge already oriented
          // if p1!=0.5 && other edge already oriented
          if (((p1 > (0.5 + eps)) || (p1 < (0.5 - eps))) &&
              ((probas_list[i][2] > (0.5 + eps)) ||
                  (probas_list[i][2] < (0.5 - eps)))) {
            logscoreTpl[i]=-__DBL_MAX__;
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][0] = p2;
            if ((!latent && probas_list[i][0] > (0.5 + eps)) ||
                (propagation && probas_list[i][0] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              probas_list2[i][1] = p1;
          }
        } else if (triples[i][2] == X && triples[i][1] == Z) {
          // sharing second edge, antisymmetry
          probas_list[i][2] = p2;

          if ((!latent && probas_list[i][2] > (0.5 + eps)) ||
              (propagation && probas_list[i][2] < (0.5 - eps)))
            // edit HI 20150222  // change 20151022 herve
            probas_list[i][3] = p1;
          // if p2!=0.5 && other edge already oriented // change 20151101
          if (((p2 > (0.5 + eps)) || (p2 < (0.5 - eps))) &&
              ((probas_list[i][1] > (0.5 + eps)) ||
                  (probas_list[i][1] < (0.5 - eps)))) {
            logscoreTpl[i]=-__DBL_MAX__;
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = false;
          } else {  // copy final proba on putative proba too
            if ((!latent && probas_list[i][2] > (0.5 + eps)) ||
                (propagation && probas_list[i][2] < (0.5 - eps)))
              // edit HI 20150222  // change 20151022 herve
              probas_list2[i][3] = p1;

            probas_list2[i][2] = p2;
          }
        } else if (triples[i][2] == Z && triples[i][1] == X) {
          // sharing second edge, same symmetry
          probas_list[i][3] = p2;
          if ((!latent && probas_list[i][3] > (0.5 + eps)) ||
              (propagation && probas_list[i][3] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            probas_list[i][2] = p1;
          // if p1!=0.5 && other edge already oriented // change 20151101
          if (((p1 > (0.5 + eps)) || (p1 < (0.5 - eps))) &&
              ((probas_list[i][1] > (0.5 + eps)) ||
                  (probas_list[i][1] < (0.5 - eps)))) {
            logscoreTpl[i]=-__DBL_MAX__;
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][3] = p2;
            if ((!latent && probas_list[i][3] > (0.5 + eps)) ||
                (propagation && probas_list[i][3] < (0.5 - eps)))
              probas_list2[i][2] = p1;
          }
        } else {
          ok = false;  // non edited tpl
        }

#if _MY_PRINT_
        Rprintf(
            "!n1=%i i=%i j=%i ok=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g "
            "\n",
            n1, i, j, ok, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

        // update score performing putative propagation
        if (ok) {
          if (I3_list[i] > 0) {
            // define score in case of no true propagation below
            scoreTpl[i] = fmin(probas_list2[i][1], 1 - probas_list2[i][2]);
            scoreTpl[i] = fmax(
                scoreTpl[i], fmin(probas_list2[i][2], 1 - probas_list2[i][1]));
            logscoreTpl[i] = log1p(scoreTpl[i] - 1);
            scoresN[i] = logF2(scoreTpl[i], I3_list[i]);

            p = probas_list[i][1];
            // only p*(1.0/(1+exp(-I3[i])) term
            logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
            pp = expm1(logpp) + 1;
            // 2->3  condition of propagation (p > pp > 0.5) and no
            // previously higher putative propagation
            if (p > (0.5 + eps) && pp > (0.5 + eps) &&
                (probas_list2[i][2] > (1 - pp + eps))) {
              logscoreTpl[i]=logpp;
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
              probas_list2[i][2] = 1 - pp;
              if (propagation && probas_list2[i][3] < (pp - eps))
                probas_list2[i][3] = pp;
            } else {  // other direction? 2->1
              p = probas_list[i][2];
              // only p*(1.0/(1+exp(-I3[i])) term
              logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
              pp = expm1(logpp) + 1;
              if (p > (0.5 + eps) && pp > (0.5 + eps) &&
                  (probas_list2[i][1] > (1 - pp + eps))) {
                logscoreTpl[i] = logpp;
                scoreTpl[i] = pp;
                scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
                probas_list2[i][1] = 1 - pp;
                if (propagation && probas_list2[i][0] < (pp - eps))
                  probas_list2[i][0] = pp;
              }
            }
          } else if (I3_list[i] < 0) {
            if (fabs(probas_list2[i][1] - probas_list2[i][2]) > eps_diff) {
              scoreTpl[i] = fmin(probas_list2[i][1], probas_list2[i][2]);
              logscoreTpl[i] = log1p(scoreTpl[i] - 1);
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
            }
            p = probas_list[i][1];
            // only p*(1.0/(1+exp(I3[i])) term
            logpp = log1p(p - 1) - log1p(exp(I3_list[i]));
            pp = expm1(logpp) + 1;

            if (pp > (0.5 + eps) && (probas_list2[i][2] < (pp - eps))) {
              // 2->3 condition of propagation (p > pp > 0.5) and no previously
              // higher putative propagation

              // update score which has decreased due to < 0 propagation!
              logscoreTpl[i] = logpp;
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);

              probas_list2[i][2] = pp;
              if (!latent && (probas_list2[i][3] > (1 - pp + eps)))
                probas_list2[i][3] = 1 - probas_list2[i][2];
            } else {
              p = probas_list[i][2];
              // only p*(1.0/(1+exp(I3[i])) term
              logpp = log1p(p - 1) - log1p(exp(I3_list[i]));
              pp = expm1(logpp) + 1;
              if (pp > (0.5 + eps) && (probas_list2[i][1] < (pp - eps))) {
                // update score which has decreased due to < 0 propagation!
                logscoreTpl[i] = logpp;
                scoreTpl[i] = pp;
                scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
                probas_list2[i][1] = pp;
                if (!latent && (probas_list2[i][0] > (1 - pp + eps)))
                  probas_list2[i][0] = 1 - probas_list2[i][1];
              }
            }
          }
        }
      }
    }

#if _MY_PRINT_
    Rprintf("!n3=%i p3=%g p4=%g count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", n3,
        p3, p4, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

    if (Y != -1) {
      // only loop on relevant Tpl i != maxTpl, thanks to the sorting of their
      // scoreTpl
      for (j = n_triples - 2; j >= 0; j--) {
        i = orderTpl[j];
        ok = true;
        // sharing fist edge, same symmetry
        if (triples[i][0] == Y && triples[i][1] == Z) {
          probas_list[i][1] = p3;
          if ((!latent && probas_list[i][1] > (0.5 + eps)) ||
              (propagation && probas_list[i][1] < (0.5 - eps)))
            // edit HI 20150222  // change 20151022 herve
            probas_list[i][0] = p4;
          // if p3!=0.5 && other edge already oriented // change 20151101
          if (((p3 > (0.5 + eps)) || (p3 < (0.5 - eps))) &&
              ((probas_list[i][2] > (0.5 + eps)) ||
                  (probas_list[i][2] < (0.5 - eps)))) {
            logscoreTpl[i] = -__DBL_MAX__;
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = false;
          } else {  // copy final proba on putative proba too
            if ((!latent && probas_list[i][1] > (0.5 + eps)) ||
                (propagation && probas_list[i][1] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              probas_list2[i][0] = p4;

            probas_list2[i][1] = p3;
          }
        } else if (triples[i][0] == Z && triples[i][1] == Y) {
          // sharing fist edge, antisymmetry
          probas_list[i][0] = p3;
          if ((!latent && probas_list[i][0] > (0.5 + eps)) ||
              (propagation && probas_list[i][0] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            probas_list[i][1] = p4;
          // if p4!=0.5 && other edge already oriented // change 20151101
          if (((p4 > (0.5 + eps)) || (p4 < (0.5 - eps))) &&
              ((probas_list[i][2] > (0.5 + eps)) ||
                  (probas_list[i][2] < (0.5 - eps)))) {
            logscoreTpl[i] = -__DBL_MAX__;
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][0] = p3;
            if ((!latent && probas_list[i][0] > (0.5 + eps)) ||
                (propagation && probas_list[i][0] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              probas_list2[i][1] = p4;
          }
        } else if (triples[i][2] == Y && triples[i][1] == Z) {
          // sharing second edge, antisymmetry
          probas_list[i][2] = p3;
          if ((!latent && probas_list[i][2] > (0.5 + eps)) ||
              (propagation && probas_list[i][2] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            probas_list[i][3] = p4;
          // if p3!=0.5 && other edge already oriented // change 20151101
          if (((p3 > (0.5 + eps)) || (p3 < (0.5 - eps))) &&
              ((probas_list[i][1] > (0.5 + eps)) ||
                  (probas_list[i][1] < (0.5 - eps)))) {
            logscoreTpl[i] = -__DBL_MAX__;
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = false;
          } else {  // copy final proba on putative proba too
            if ((!latent && probas_list[i][2] > (0.5 + eps)) ||
                (propagation && probas_list[i][2] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              probas_list2[i][3] = p4;

            probas_list2[i][2] = p3;
          }
        } else if (triples[i][2] == Z && triples[i][1] == Y) {
          // sharing second edge, same symmetry
          probas_list[i][3] = p3;
          if ((!latent && probas_list[i][3] > (0.5 + eps)) ||
              (propagation && probas_list[i][3] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            probas_list[i][2] = p4;
          // if p4!=0.5 && other edge already oriented // change 20151101
          if (((p4 > (0.5 + eps)) || (p4 < (0.5 - eps))) &&
              ((probas_list[i][1] > (0.5 + eps)) ||
                  (probas_list[i][1] < (0.5 - eps)))) {
            logscoreTpl[i] = -__DBL_MAX__;
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = false;
          } else {  // copy final proba on putative proba too
            probas_list2[i][3] = p3;
            if ((!latent && probas_list[i][3] > (0.5 + eps)) ||
                (propagation && probas_list[i][3] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              probas_list2[i][2] = p4;
          }
        } else {
          ok = false;  // non edited tpl
        }

#if _MY_PRINT_
        Rprintf(
            "!n3=%i i=%i j=%i ok=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g "
            "\n",
            n3, i, j, ok, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_
        // update score performing putative propagation
        if (ok) {
          if (I3_list[i] > 0) {
            // define score in case of no true propagation below
            scoreTpl[i] = fmin(probas_list2[i][1], 1 - probas_list2[i][2]);
            scoreTpl[i] = fmax(
                scoreTpl[i], fmin(probas_list2[i][2], 1 - probas_list2[i][1]));
            logscoreTpl[i] = log1p(scoreTpl[i] - 1);
            scoresN[i] = logF2(scoreTpl[i], I3_list[i]);

            p = probas_list[i][1];
            // only p*(1.0/(1+exp(-I3[i])) term
            logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
            pp = expm1(logpp) + 1;
            // change 20151031
            if (p > (0.5 + eps) && pp > (0.5 + eps) &&
                (probas_list2[i][2] > (1 - pp + eps))) {
              // 2->3  condition of propagation (p > pp > 0.5) and no previously
              // higher putative propagation
              logscoreTpl[i] = logpp;
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
              probas_list2[i][2] = 1 - pp;
              if (propagation && (probas_list2[i][3] < (pp - eps)))
                probas_list2[i][3] = pp;  // change 20151031
            } else {                      // other direction? 2->1
              p = probas_list[i][2];
              // only p*(1.0/(1+exp(-I3[i])) term
              logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
              pp = expm1(logpp) + 1;
              // change 20151031
              if (p > (0.5 + eps) && pp > (0.5 + eps) &&
                  (probas_list2[i][1] > (1 - pp + eps))) {
                // 2->1   condition of propagation (p > pp > 0.5) and no
                // previously higher putative propagation
                logscoreTpl[i] = logpp;
                scoreTpl[i] = pp;
                scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
                probas_list2[i][1] = 1 - pp;
                if (propagation && (probas_list2[i][0] < (pp - eps)))
                  probas_list2[i][0] = pp;
              }
            }
          } else if (I3_list[i] < 0) {
            // define score in case of no true propagation below
            if (fabs(probas_list2[i][1] - probas_list2[i][2]) > eps_diff) {
              scoreTpl[i] = fmin(probas_list2[i][1], probas_list2[i][2]);
              logscoreTpl[i] = log1p(scoreTpl[i] - 1);
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
            }
            p = probas_list[i][1];
            // only p*(1.0/(1+exp(-I3[i])) term
            logpp = log1p(p - 1) - log1p(exp(-I3_list[i]));
            pp = expm1(logpp) + 1;
            if (pp > (0.5 + eps) && (probas_list2[i][2] < (pp - eps))) {
              // 2->3 condition of propagation (p > pp > 0.5) and no previously
              // higher putative propagation
              logscoreTpl[i] = logpp;
              scoreTpl[i] = pp;  // update score which has decreased due to < 0
                                 // propagation!
              scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
              probas_list2[i][2] = pp;
              if (!latent && (probas_list2[i][3] > (1 - pp + eps)))
                probas_list2[i][3] = 1 - probas_list2[i][2];  // change 20151031
            } else {
              p = probas_list[i][2];
              // only p*(1.0/(1+exp(I3[i])) term
              logpp = log1p(p - 1) - log1p(exp(I3_list[i]));
              pp = expm1(logpp) + 1;
              if (pp > (0.5 + eps) && (probas_list2[i][1] < (pp - eps))) {
                // 2->1
                logscoreTpl[i] = logpp;
                // update score which has decreased due to < 0 propagation!
                scoreTpl[i] = pp;
                scoresN[i] = logF2(scoreTpl[i], I3_list[i]);
                probas_list2[i][1] = pp;
                if (!latent && (probas_list2[i][0] > (1 - pp + eps)))
                  probas_list2[i][0] = 1 - probas_list2[i][1];
              }
            }
          }
        }
      }
    }
#if _MY_PRINT_
    Rprintf("!n3=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", n3, count, deg,
        maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

    //  Order Tpl in increasing ScoreTpl of scoresN
    std::sort(begin(orderTpl), end(orderTpl), TripletComparator(scoreTpl));

    maxTpl = orderTpl[n_triples - 1];
    maxscoreTpl = scoreTpl[maxTpl];
#if _MY_PRINT_
    for (i = 0; i < NbTpl; i++) {
      Rprintf("count=%i scoreTpl[orderTpl[%i]=%i]=%g", count, i, orderTpl[i],
          scoreTpl[orderTpl[i]]);
    }
    Rprintf(" maxTpl=%i  maxscoreTpl=%g\n", maxTpl, maxscoreTpl);

    Rprintf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
    for (i = 0; i < NbTpl; i++) {
      Rprintf("!!! %i Tpl %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n", i,
          Tpl[i][0], ProbArrowhead[i][0],
          ProbArrowhead[i][1], Tpl[i][1],
          ProbArrowhead[i][2], ProbArrowhead[i][3],
          Tpl[i][2], I3[i], scoreTpl[i]);
    }
    Rprintf("\n");
    for (i = 0; i < NbTpl; i++) {
      Rprintf("!!! %i Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n", i,
          Tpl[i][0], ProbArrowhead2[i][0],
          ProbArrowhead2[i][1], Tpl[i][1],
          ProbArrowhead2[i][2], ProbArrowhead2[i][3],
          Tpl[i][2], I3[i], scoreTpl[i]);
    }

    Rprintf("!INTER latent=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n",
        latent, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

    i = maxTpl;
  } while ((maxscoreTpl > (0.5 + eps)) && count < 200000);

#if _MY_PRINT_
  Rprintf("\n\n consist-neg-pos\n");
  Rprintf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
  for (i = 0; i < NbTpl; i++) {
    Rprintf("!!! Tpl %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n",
        Tpl[i][0], ProbArrowhead[i][0],
        ProbArrowhead[i][1], Tpl[i][1],
        ProbArrowhead[i][2], ProbArrowhead[i][3],
        Tpl[i][2], I3[i], scoreTpl[i]);
  }
  for (i = 0; i < NbTpl; i++) {
    Rprintf("!!! Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n",
        Tpl[i][0], ProbArrowhead2[i][0],
        ProbArrowhead2[i][1], Tpl[i][1],
        ProbArrowhead2[i][2], ProbArrowhead2[i][3],
        Tpl[i][2], I3[i], scoreTpl[i]);
  }

  Rprintf("!END latent=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", latent,
      count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

  return probas_list;
}

}  // namespace reconstruction
}  // namespace miic
