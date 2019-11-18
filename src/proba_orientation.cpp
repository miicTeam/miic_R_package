#include "proba_orientation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>

#include "utilities.h"

using namespace std;

//#define M_PI		3.14159265358979323846	/* pi */
#define eps 0.000000000001
#define _MY_PRINT_ 0

static double dmaxarg1, dmaxarg2;
#define DMAX(a, b)                 \
  (dmaxarg1 = (a), dmaxarg2 = (b), \
      (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

static double dminarg1, dminarg2;
#define DMIN(a, b)                 \
  (dminarg1 = (a), dminarg2 = (b), \
      (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))

class sort_indices_3points {
  double* mparr;

 public:
  sort_indices_3points(double* parr) : mparr(parr) {}
  bool operator()(int i, int j) const { return mparr[i] < mparr[j]; }
};

double logF2(double scoreTpl, double I3) {
  double scoreN;
  if (scoreTpl < (1 - eps))
    scoreN = log(-2 / (scoreTpl * (1 / (1 + exp(-fabs(I3))) - 0.5) - 0.5) - 3);
  else
    scoreN = fabs(I3) + log1p(-2 * exp(-fabs(I3)));
  return scoreN;
}

int miic::reconstruction::OrientTpl_LV_Deg_Propag(int NbTpl, int* Tpl,
    double* I3, double* ProbArrowhead, int LV, int deg, int Propag,
    int HALFVSTRUCT) {
  int i, j, n1, n2, n3, TRUE = 1, FALSE = 0, *orderTpl, maxTpl, ok, count = 0;
  double maxscoreTpl;
  double p, pp, p1, p2, p3, p4, *ProbArrowhead2, *scoreTpl, *scoresN;

  ProbArrowhead2 = new double[4 * NbTpl];
  scoreTpl = new double[NbTpl];
  orderTpl = new int[NbTpl];
  scoresN = new double[NbTpl];

  // Initialization ProbArrowhead2
  for (i = 0; i < 4 * NbTpl; i++) {
    ProbArrowhead2[i] = ProbArrowhead[i];
  }
  // Initialization ScoreTpl for I3<0
  for (i = 0; i < NbTpl; i++) {
    if (I3[i] < 0) {       // Tpl negative
      if (deg == FALSE) {  // no degenerescence
        scoreTpl[i] =
            (1 + exp(I3[i])) / (1 + 3 * exp(I3[i]));  // update 20150228
        scoresN[i] = -(I3[i]);
      } else {  // degenerescence
        scoreTpl[i] = (3 - 2 * exp(I3[i])) /
                      (3 - exp(I3[i]));  // larger than p without deg
      }
      // modif 20150228
      ProbArrowhead2[1 * NbTpl + i] = scoreTpl[i];
      ProbArrowhead2[2 * NbTpl + i] = scoreTpl[i];

      if (LV == FALSE) {  // no change 20151023
        ProbArrowhead2[0 * NbTpl + i] = 1 - ProbArrowhead2[1 * NbTpl + i];
        ProbArrowhead2[3 * NbTpl + i] = 1 - ProbArrowhead2[2 * NbTpl + i];
      }

    } else {  // Tpl positive
      scoreTpl[i] = 0.5;
      scoresN[i] = 0;
    }
  }

#if _MY_PRINT_
  printf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
  for (i = 0; i < NbTpl; i++) {
    printf("!!! Tpl %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n",
        Tpl[0 * NbTpl + i], ProbArrowhead[0 * NbTpl + i],
        ProbArrowhead[1 * NbTpl + i], Tpl[1 * NbTpl + i],
        ProbArrowhead[2 * NbTpl + i], ProbArrowhead[3 * NbTpl + i],
        Tpl[2 * NbTpl + i], I3[i], scoreTpl[i]);
  }
  for (i = 0; i < NbTpl; i++) {
    printf("!!! Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n",
        Tpl[0 * NbTpl + i], ProbArrowhead2[0 * NbTpl + i],
        ProbArrowhead2[1 * NbTpl + i], Tpl[1 * NbTpl + i],
        ProbArrowhead2[2 * NbTpl + i], ProbArrowhead2[3 * NbTpl + i],
        Tpl[2 * NbTpl + i], I3[i], scoreTpl[i]);
  }

  printf("!START LV=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", LV, count,
      deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

  //  Order Tpl in increasing ScoreTpl /////////////////////
  for (i = 0; i < NbTpl; i++) {
    orderTpl[i] = i;
  }
  std::sort(orderTpl, orderTpl + NbTpl, sort_indices_3points(scoreTpl));
  double max = -10000;
  int maxpos = 0;
  for (i = 0; i < NbTpl; i++) {
    if (scoreTpl[i] > max) {
      maxpos = i;
      max = scoreTpl[i];
    }
  }

  maxTpl = maxpos;
  maxscoreTpl = scoreTpl[maxTpl];

#if _MY_PRINT_
  for (i = 0; i < NbTpl; i++) {
    printf(" scoreTpl[orderTpl[%i]=%i]=%g  \n", i, orderTpl[i],
        scoreTpl[orderTpl[i]]);
  }
  printf(" maxTpl=%i  maxscoreTpl=%g\n", maxTpl, maxscoreTpl);
#endif

  do {
    count++;
#if _MY_PRINT_
    printf("! count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", count, deg, maxTpl,
        maxscoreTpl);
#endif  // _MY_PRINT_

    // assign topTpl ///////////////////////////////////////////:
    n1 = 0;
    n3 = 0;
    p1 = 0.5;
    p2 = 0.5;
    p3 = 0.5;
    p4 = 0.5;

    i = maxTpl;
    scoreTpl[i] = -1;
    scoresN[i] = -1;

#if _MY_PRINT_
    printf("maxTpl=%i n1=%i  %g(%g)--%g(%g)  n2=%i  %g(%g)--%g(%g)  n3=%i \n",
        maxTpl, Tpl[0 * NbTpl + i], ProbArrowhead[0 * NbTpl + i],
        ProbArrowhead2[0 * NbTpl + i], ProbArrowhead[1 * NbTpl + i],
        ProbArrowhead2[1 * NbTpl + i], Tpl[1 * NbTpl + i],
        ProbArrowhead[2 * NbTpl + i], ProbArrowhead2[2 * NbTpl + i],
        ProbArrowhead[3 * NbTpl + i], ProbArrowhead2[3 * NbTpl + i],
        Tpl[2 * NbTpl + i]);
#endif  // _MY_PRINT_

    // arrowhead proba
    p = DMAX(ProbArrowhead2[1 * NbTpl + i],
        1 - ProbArrowhead2[1 * NbTpl + i]);  // change 20151031
    // if arrowhead/tail on 1 (x 0-*1 z 2-3 y) is not already established
    // (through an earlier propagation)
    if ((ProbArrowhead[1 * NbTpl + i] < (p - eps)) &&
        (ProbArrowhead[1 * NbTpl + i] > (1 - p + eps)) &&
        (HALFVSTRUCT == TRUE || I3[i] > 0 ||
            ProbArrowhead[2 * NbTpl + i] > (0.5 - eps))) {  // change 20151024
      // establish arrowhead/tail final proba on 1 (x 0-*1 z 2-3 y)
      ProbArrowhead[1 * NbTpl + i] = ProbArrowhead2[1 * NbTpl + i];

      n1 = Tpl[0 * NbTpl + i];
      n2 = Tpl[1 * NbTpl + i];
      // not sure this is useful;) //change 20151023
      p1 = ProbArrowhead2[0 * NbTpl + i];
      // if no LV or tail on 1  (x 0<-1 z 2-3 y) // change 20151022 herve
      if ((LV == FALSE && ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
          (Propag == TRUE && ProbArrowhead[1 * NbTpl + i] < (0.5 - eps))) {
        // establish arrowhead/tail if no LV or
        // arrowhead final proba on 0 (x 0<-1 z 2-3 y)
        ProbArrowhead[0 * NbTpl + i] = ProbArrowhead2[0 * NbTpl + i];
      }
      p2 = ProbArrowhead[1 * NbTpl + i];
    }

    p = DMAX(ProbArrowhead2[2 * NbTpl + i],
        1 - ProbArrowhead2[2 * NbTpl + i]);  // change 20151031

    // change 20151024
    if ((ProbArrowhead[2 * NbTpl + i] < (p - eps)) &&
        (ProbArrowhead[2 * NbTpl + i] > (1 - p + eps)) &&
        (HALFVSTRUCT == TRUE || I3[i] > 0 ||
            ProbArrowhead[1 * NbTpl + i] > (0.5 - eps))) {
      // establish arrowhead/tail final proba on 2 (x 0-1 z 2*-3 y)
      ProbArrowhead[2 * NbTpl + i] = ProbArrowhead2[2 * NbTpl + i];

      n2 = Tpl[1 * NbTpl + i];
      n3 = Tpl[2 * NbTpl + i];
      p3 = ProbArrowhead[2 * NbTpl + i];
      // not sure this is useful;) //change 20151023
      p4 = ProbArrowhead2[3 * NbTpl + i];
      // if no LV or tail on 2  (x 0-1 z 2->3 y) // change 20151022 herve
      if ((LV == FALSE && ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
          (Propag == TRUE && ProbArrowhead[2 * NbTpl + i] < (0.5 - eps))) {
        // establish arrowhead/tail if no LV or
        // arrowhead final proba on 3 (x 0-1 z 2->3 y)
        ProbArrowhead[3 * NbTpl + i] = ProbArrowhead2[3 * NbTpl + i];
      }
    }

#if _MY_PRINT_
    printf("!n1=%i p1=%g p2=%g count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", n1,
        p1, p2, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

    maxscoreTpl = 0.0;

    if (n1 != 0) {
      // only loop on relevant Tpl i != maxTpl, thanks to the sorting
      // of their scoreTpl
      for (j = NbTpl - 2; j >= 0; j--) {
        i = orderTpl[j];
        ok = TRUE;
        // sharing fist edge, same symmetry
        if (Tpl[0 * NbTpl + i] == n1 && Tpl[1 * NbTpl + i] == n2) {
          ProbArrowhead[1 * NbTpl + i] = p2;
          // if(LV==FALSE) ProbArrowhead[0*NbTpl+i] = p1;
          if ((LV == FALSE && ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
              (Propag == TRUE && ProbArrowhead[1 * NbTpl + i] < (0.5 - eps)))
            // edit HI 20150222 change 20151022 herve
            ProbArrowhead[0 * NbTpl + i] = p1;
          // if p2!=0.5 && other edge already oriented // change 20151101
          if (((p2 > (0.5 + eps)) || (p2 < (0.5 - eps))) &&
              ((ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
                  (ProbArrowhead[2 * NbTpl + i] < (0.5 - eps)))) {
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = FALSE;
          } else {  // copy final proba on putative proba too
            // edit HI 20150222 // change 20151022 herve
            if ((LV == FALSE && ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
                (Propag == TRUE && ProbArrowhead[1 * NbTpl + i] < (0.5 - eps)))
              ProbArrowhead2[0 * NbTpl + i] = p1;

            ProbArrowhead2[1 * NbTpl + i] = p2;
          }
        } else if (Tpl[0 * NbTpl + i] == n2 && Tpl[1 * NbTpl + i] == n1) {
          // sharing fist edge, antisymmetry
          ProbArrowhead[0 * NbTpl + i] = p2;
          if ((LV == FALSE && ProbArrowhead[0 * NbTpl + i] > (0.5 + eps)) ||
              (Propag == TRUE && ProbArrowhead[0 * NbTpl + i] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            ProbArrowhead[1 * NbTpl + i] = p1;
          // if other edge already oriented
          // if p1!=0.5 && other edge already oriented // change 20151101
          if (((p1 > (0.5 + eps)) || (p1 < (0.5 - eps))) &&
              ((ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
                  (ProbArrowhead[2 * NbTpl + i] < (0.5 - eps)))) {
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = FALSE;
          } else {  // copy final proba on putative proba too
            ProbArrowhead2[0 * NbTpl + i] = p2;
            if ((LV == FALSE && ProbArrowhead[0 * NbTpl + i] > (0.5 + eps)) ||
                (Propag == TRUE && ProbArrowhead[0 * NbTpl + i] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              ProbArrowhead2[1 * NbTpl + i] = p1;
          }
        } else if (Tpl[2 * NbTpl + i] == n1 && Tpl[1 * NbTpl + i] == n2) {
          // sharing second edge, antisymmetry
          ProbArrowhead[2 * NbTpl + i] = p2;

          if ((LV == FALSE && ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
              (Propag == TRUE && ProbArrowhead[2 * NbTpl + i] < (0.5 - eps)))
            // edit HI 20150222  // change 20151022 herve
            ProbArrowhead[3 * NbTpl + i] = p1;
          // if p2!=0.5 && other edge already oriented // change 20151101
          if (((p2 > (0.5 + eps)) || (p2 < (0.5 - eps))) &&
              ((ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
                  (ProbArrowhead[1 * NbTpl + i] < (0.5 - eps)))) {
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = FALSE;
          } else {  // copy final proba on putative proba too
            if ((LV == FALSE && ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
                (Propag == TRUE && ProbArrowhead[2 * NbTpl + i] < (0.5 - eps)))
              // edit HI 20150222  // change 20151022 herve
              ProbArrowhead2[3 * NbTpl + i] = p1;

            ProbArrowhead2[2 * NbTpl + i] = p2;
          }
        } else if (Tpl[2 * NbTpl + i] == n2 && Tpl[1 * NbTpl + i] == n1) {
          // sharing second edge, same symmetry
          ProbArrowhead[3 * NbTpl + i] = p2;
          if ((LV == FALSE && ProbArrowhead[3 * NbTpl + i] > (0.5 + eps)) ||
              (Propag == TRUE && ProbArrowhead[3 * NbTpl + i] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            ProbArrowhead[2 * NbTpl + i] = p1;
          // if p1!=0.5 && other edge already oriented // change 20151101
          if (((p1 > (0.5 + eps)) || (p1 < (0.5 - eps))) &&
              ((ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
                  (ProbArrowhead[1 * NbTpl + i] < (0.5 - eps)))) {
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = FALSE;
          } else {  // copy final proba on putative proba too
            ProbArrowhead2[3 * NbTpl + i] = p2;
            // edit HI 20150222 // change 20151022 herve
            if ((LV == FALSE && ProbArrowhead[3 * NbTpl + i] > (0.5 + eps)) ||
                (Propag == TRUE && ProbArrowhead[3 * NbTpl + i] < (0.5 - eps)))
              ProbArrowhead2[2 * NbTpl + i] = p1;
          }
        } else {
          ok = FALSE;  // non edited tpl
        }

#if _MY_PRINT_
        printf(
            "!n1=%i i=%i j=%i ok=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g "
            "\n",
            n1, i, j, ok, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

        // update score performing putative propagation
        if (ok) {
          if (I3[i] > 0) {
            // define score in case of no true propagation below
            scoreTpl[i] = DMIN(ProbArrowhead2[1 * NbTpl + i],
                1 - ProbArrowhead2[2 * NbTpl + i]);
            scoreTpl[i] =
                DMAX(scoreTpl[i], DMIN(ProbArrowhead2[2 * NbTpl + i],
                                      1 - ProbArrowhead2[1 * NbTpl + i]));
            scoresN[i] = logF2(scoreTpl[i], I3[i]);

            p = ProbArrowhead[1 * NbTpl + i];
            pp = 0.5 + p * (1.0 / (1 + exp(-I3[i])) - 0.5);  // update 20150228
            // change 20151031
            // 2->3  condition of propagation (p > pp > 0.5) and no
            // previously higher putative propagation
            if (p > (0.5 + eps) && pp > (0.5 + eps) &&
                (ProbArrowhead2[2 * NbTpl + i] > (1 - pp + eps))) {
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3[i]);

              ProbArrowhead2[2 * NbTpl + i] = 1 - pp;
              if (Propag == TRUE && ProbArrowhead2[3 * NbTpl + i] < (pp - eps))
                ProbArrowhead2[3 * NbTpl + i] = pp;  // change 20151031
            } else {  // other direction? 2->1
              p = ProbArrowhead[2 * NbTpl + i];
              // update 20150228
              pp = 0.5 + p * (1.0 / (1 + exp(-I3[i])) - 0.5);
              // change 20151031
              if (p > (0.5 + eps) && pp > (0.5 + eps) &&
                  (ProbArrowhead2[1 * NbTpl + i] > (1 - pp + eps))) {
                scoreTpl[i] = pp;
                scoresN[i] = logF2(scoreTpl[i], I3[i]);
                ProbArrowhead2[1 * NbTpl + i] = 1 - pp;
                // change 20151031
                if (Propag == TRUE &&
                    ProbArrowhead2[0 * NbTpl + i] < (pp - eps))
                  ProbArrowhead2[0 * NbTpl + i] = pp;
              }
            }
          } else if (I3[i] < 0) {
            if (fabs(ProbArrowhead2[1 * NbTpl + i] -
                     ProbArrowhead2[2 * NbTpl + i]) > eps) {
              scoreTpl[i] = DMIN(
                  ProbArrowhead2[1 * NbTpl + i], ProbArrowhead2[2 * NbTpl + i]);
              scoresN[i] = logF2(scoreTpl[i], I3[i]);
            }
            p = ProbArrowhead[1 * NbTpl + i];
            pp = 0.5 + p * (1.0 / (1 + exp(I3[i])) - 0.5);

            if (pp > (0.5 + eps) &&
                (ProbArrowhead2[2 * NbTpl + i] < (pp - eps))) {
              // 2->3 condition of propagation (p > pp > 0.5) and no previously
              // higher putative propagation

              // update score which has decreased due to < 0 propagation!
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3[i]);

              ProbArrowhead2[2 * NbTpl + i] = pp;
              // change 20151031
              if (LV == FALSE &&
                  (ProbArrowhead2[3 * NbTpl + i] > (1 - pp + eps)))
                ProbArrowhead2[3 * NbTpl + i] =
                    1 - ProbArrowhead2[2 * NbTpl + i];
            } else {
              p = ProbArrowhead[2 * NbTpl + i];
              pp = 0.5 + p * (1.0 / (1 + exp(I3[i])) - 0.5);
              // change 20151031
              if (pp > (0.5 + eps) &&
                  (ProbArrowhead2[1 * NbTpl + i] < (pp - eps))) {
                // update score which has decreased due to < 0 propagation!
                scoreTpl[i] = pp;
                scoresN[i] = logF2(scoreTpl[i], I3[i]);
                ProbArrowhead2[1 * NbTpl + i] = pp;
                // change 20151031
                if (LV == FALSE &&
                    (ProbArrowhead2[0 * NbTpl + i] > (1 - pp + eps)))
                  ProbArrowhead2[0 * NbTpl + i] =
                      1 - ProbArrowhead2[1 * NbTpl + i];
              }
            }
          }
        }
      }
    }

#if _MY_PRINT_
    printf("!n3=%i p3=%g p4=%g count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", n3,
        p3, p4, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

    if (n3 != 0) {
      // only loop on relevant Tpl i != maxTpl, thanks to the sorting of their
      // scoreTpl
      for (j = NbTpl - 2; j >= 0; j--) {
        i = orderTpl[j];
        ok = TRUE;
        // sharing fist edge, same symmetry
        if (Tpl[0 * NbTpl + i] == n3 && Tpl[1 * NbTpl + i] == n2) {
          ProbArrowhead[1 * NbTpl + i] = p3;
          if ((LV == FALSE && ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
              (Propag == TRUE && ProbArrowhead[1 * NbTpl + i] < (0.5 - eps)))
            // edit HI 20150222  // change 20151022 herve
            ProbArrowhead[0 * NbTpl + i] = p4;
          // if p3!=0.5 && other edge already oriented // change 20151101
          if (((p3 > (0.5 + eps)) || (p3 < (0.5 - eps))) &&
              ((ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
                  (ProbArrowhead[2 * NbTpl + i] < (0.5 - eps)))) {
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = FALSE;
          } else {  // copy final proba on putative proba too
            if ((LV == FALSE && ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
                (Propag == TRUE && ProbArrowhead[1 * NbTpl + i] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              ProbArrowhead2[0 * NbTpl + i] = p4;

            ProbArrowhead2[1 * NbTpl + i] = p3;
          }
        } else if (Tpl[0 * NbTpl + i] == n2 && Tpl[1 * NbTpl + i] == n3) {
          // sharing fist edge, antisymmetry
          ProbArrowhead[0 * NbTpl + i] = p3;
          if ((LV == FALSE && ProbArrowhead[0 * NbTpl + i] > (0.5 + eps)) ||
              (Propag == TRUE && ProbArrowhead[0 * NbTpl + i] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            ProbArrowhead[1 * NbTpl + i] = p4;
          // if p4!=0.5 && other edge already oriented // change 20151101
          if (((p4 > (0.5 + eps)) || (p4 < (0.5 - eps))) &&
              ((ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
                  (ProbArrowhead[2 * NbTpl + i] < (0.5 - eps)))) {
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = FALSE;
          } else {  // copy final proba on putative proba too
            ProbArrowhead2[0 * NbTpl + i] = p3;
            if ((LV == FALSE && ProbArrowhead[0 * NbTpl + i] > (0.5 + eps)) ||
                (Propag == TRUE && ProbArrowhead[0 * NbTpl + i] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              ProbArrowhead2[1 * NbTpl + i] = p4;
          }
        } else if (Tpl[2 * NbTpl + i] == n3 && Tpl[1 * NbTpl + i] == n2) {
          // sharing second edge, antisymmetry
          ProbArrowhead[2 * NbTpl + i] = p3;
          if ((LV == FALSE && ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
              (Propag == TRUE && ProbArrowhead[2 * NbTpl + i] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            ProbArrowhead[3 * NbTpl + i] = p4;
          // if p3!=0.5 && other edge already oriented // change 20151101
          if (((p3 > (0.5 + eps)) || (p3 < (0.5 - eps))) &&
              ((ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
                  (ProbArrowhead[1 * NbTpl + i] < (0.5 - eps)))) {
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = FALSE;
          } else {  // copy final proba on putative proba too
            if ((LV == FALSE && ProbArrowhead[2 * NbTpl + i] > (0.5 + eps)) ||
                (Propag == TRUE && ProbArrowhead[2 * NbTpl + i] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              ProbArrowhead2[3 * NbTpl + i] = p4;

            ProbArrowhead2[2 * NbTpl + i] = p3;
          }
        } else if (Tpl[2 * NbTpl + i] == n2 && Tpl[1 * NbTpl + i] == n3) {
          // sharing second edge, same symmetry
          ProbArrowhead[3 * NbTpl + i] = p3;
          if ((LV == FALSE && ProbArrowhead[3 * NbTpl + i] > (0.5 + eps)) ||
              (Propag == TRUE && ProbArrowhead[3 * NbTpl + i] < (0.5 - eps)))
            // edit HI 20150222 // change 20151022 herve
            ProbArrowhead[2 * NbTpl + i] = p4;
          // if p4!=0.5 && other edge already oriented // change 20151101
          if (((p4 > (0.5 + eps)) || (p4 < (0.5 - eps))) &&
              ((ProbArrowhead[1 * NbTpl + i] > (0.5 + eps)) ||
                  (ProbArrowhead[1 * NbTpl + i] < (0.5 - eps)))) {
            scoreTpl[i] = -1;  // remove tpl from the list of yet-unused Tpl
            scoresN[i] = -1;
            ok = FALSE;
          } else {  // copy final proba on putative proba too
            ProbArrowhead2[3 * NbTpl + i] = p3;
            if ((LV == FALSE && ProbArrowhead[3 * NbTpl + i] > (0.5 + eps)) ||
                (Propag == TRUE && ProbArrowhead[3 * NbTpl + i] < (0.5 - eps)))
              // edit HI 20150222 // change 20151022 herve
              ProbArrowhead2[2 * NbTpl + i] = p4;
          }
        } else {
          ok = FALSE;  // non edited tpl
        }

#if _MY_PRINT_
        printf(
            "!n3=%i i=%i j=%i ok=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g "
            "\n",
            n3, i, j, ok, count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_
        // update score performing putative propagation
        if (ok) {
          if (I3[i] > 0) {
            // define score in case of no true propagation below
            scoreTpl[i] = DMIN(ProbArrowhead2[1 * NbTpl + i],
                1 - ProbArrowhead2[2 * NbTpl + i]);
            scoreTpl[i] =
                DMAX(scoreTpl[i], DMIN(ProbArrowhead2[2 * NbTpl + i],
                                      1 - ProbArrowhead2[1 * NbTpl + i]));
            scoresN[i] = logF2(scoreTpl[i], I3[i]);

            p = ProbArrowhead[1 * NbTpl + i];
            pp = 0.5 + p * (1.0 / (1 + exp(-I3[i])) - 0.5);  // update 20150228
            // change 20151031
            if (p > (0.5 + eps) && pp > (0.5 + eps) &&
                (ProbArrowhead2[2 * NbTpl + i] > (1 - pp + eps))) {
              // 2->3  condition of propagation (p > pp > 0.5) and no previously
              // higher putative propagation
              scoreTpl[i] = pp;
              scoresN[i] = logF2(scoreTpl[i], I3[i]);
              ProbArrowhead2[2 * NbTpl + i] = 1 - pp;
              if (Propag == TRUE &&
                  (ProbArrowhead2[3 * NbTpl + i] < (pp - eps)))
                ProbArrowhead2[3 * NbTpl + i] = pp;  // change 20151031
            } else {  // other direction? 2->1
              p = ProbArrowhead[2 * NbTpl + i];
              pp = 0.5 + p * (1.0 / (1 + exp(-I3[i])) - 0.5);
              // change 20151031
              if (p > (0.5 + eps) && pp > (0.5 + eps) &&
                  (ProbArrowhead2[1 * NbTpl + i] > (1 - pp + eps))) {
                // 2->1   condition of propagation (p > pp > 0.5) and no
                // previously higher putative propagation
                scoreTpl[i] = pp;
                scoresN[i] = logF2(scoreTpl[i], I3[i]);
                ProbArrowhead2[1 * NbTpl + i] = 1 - pp;
                if (Propag == TRUE &&
                    (ProbArrowhead2[0 * NbTpl + i] < (pp - eps)))
                  ProbArrowhead2[0 * NbTpl + i] = pp;  // change 20151031
              }
            }
          } else if (I3[i] < 0) {
            // define score in case of no true propagation below
            if (fabs(ProbArrowhead2[1 * NbTpl + i] -
                     ProbArrowhead2[2 * NbTpl + i]) > eps) {  // modif 20150228
              scoreTpl[i] = DMIN(
                  ProbArrowhead2[1 * NbTpl + i], ProbArrowhead2[2 * NbTpl + i]);
              scoresN[i] = logF2(scoreTpl[i], I3[i]);
            }
            p = ProbArrowhead[1 * NbTpl + i];
            pp = 0.5 + p * (1.0 / (1 + exp(I3[i])) - 0.5);
            // change 20151031
            if (pp > (0.5 + eps) &&
                (ProbArrowhead2[2 * NbTpl + i] < (pp - eps))) {
              // 2->3 condition of propagation (p > pp > 0.5) and no previously
              // higher putative propagation
              scoreTpl[i] = pp;  // update score which has decreased due to < 0
                                 // propagation!
              scoresN[i] = logF2(scoreTpl[i], I3[i]);
              ProbArrowhead2[2 * NbTpl + i] = pp;
              if (LV == FALSE &&
                  (ProbArrowhead2[3 * NbTpl + i] > (1 - pp + eps)))
                ProbArrowhead2[3 * NbTpl + i] =
                    1 - ProbArrowhead2[2 * NbTpl + i];  // change 20151031
            } else {
              p = ProbArrowhead[2 * NbTpl + i];
              pp = 0.5 + p * (1.0 / (1 + exp(I3[i])) - 0.5);
              // change 20151031
              if (pp > (0.5 + eps) &&
                  (ProbArrowhead2[1 * NbTpl + i] < (pp - eps))) {
                // 2->1
                scoreTpl[i] = pp;  // update score which has decreased due to <0
                                   // propagation!
                scoresN[i] = logF2(scoreTpl[i], I3[i]);
                ProbArrowhead2[1 * NbTpl + i] = pp;
                if (LV == FALSE &&
                    (ProbArrowhead2[0 * NbTpl + i] > (1 - pp + eps)))
                  ProbArrowhead2[0 * NbTpl + i] =
                      1 - ProbArrowhead2[1 * NbTpl + i];  // change 20151031
              }
            }
          }
        }
      }
    }
#if _MY_PRINT_
    printf("!n3=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", n3, count, deg,
        maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

    //  Order Tpl in increasing ScoreTpl of scoresN
    std::sort(orderTpl, orderTpl + NbTpl, sort_indices_3points(scoreTpl));

    for (i = 0; i < NbTpl; i++) {
      orderTpl[i] = i;
    }

    double max = -10000;
    int maxpos = 0;
    for (i = 0; i < NbTpl; i++) {
      if (scoreTpl[i] > max) {
        maxpos = i;
        max = scoreTpl[i];
      }
    }
    maxTpl = maxpos;
    maxscoreTpl = scoreTpl[maxTpl];
#if _MY_PRINT_
    for (i = 0; i < NbTpl; i++) {
      printf("count=%i scoreTpl[orderTpl[%i]=%i]=%g", count, i, orderTpl[i],
          scoreTpl[orderTpl[i]]);
    }
    printf(" maxTpl=%i  maxscoreTpl=%g\n", maxTpl, maxscoreTpl);

    printf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
    for (i = 0; i < NbTpl; i++) {
      printf("!!! %i Tpl %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n", i,
          Tpl[0 * NbTpl + i], ProbArrowhead[0 * NbTpl + i],
          ProbArrowhead[1 * NbTpl + i], Tpl[1 * NbTpl + i],
          ProbArrowhead[2 * NbTpl + i], ProbArrowhead[3 * NbTpl + i],
          Tpl[2 * NbTpl + i], I3[i], scoreTpl[i]);
    }
    printf("\n");
    for (i = 0; i < NbTpl; i++) {
      printf("!!! %i Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n", i,
          Tpl[0 * NbTpl + i], ProbArrowhead2[0 * NbTpl + i],
          ProbArrowhead2[1 * NbTpl + i], Tpl[1 * NbTpl + i],
          ProbArrowhead2[2 * NbTpl + i], ProbArrowhead2[3 * NbTpl + i],
          Tpl[2 * NbTpl + i], I3[i], scoreTpl[i]);
    }

    printf("!INTER LV=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", LV,
        count, deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

    i = maxTpl;
  } while ((maxscoreTpl > (0.5 + eps)) && count < 200000);

#if _MY_PRINT_
  printf("\n\n consist-neg-pos\n");
  printf("\n\n Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
  for (i = 0; i < NbTpl; i++) {
    printf("!!! Tpl %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n",
        Tpl[0 * NbTpl + i], ProbArrowhead[0 * NbTpl + i],
        ProbArrowhead[1 * NbTpl + i], Tpl[1 * NbTpl + i],
        ProbArrowhead[2 * NbTpl + i], ProbArrowhead[3 * NbTpl + i],
        Tpl[2 * NbTpl + i], I3[i], scoreTpl[i]);
  }
  for (i = 0; i < NbTpl; i++) {
    printf("!!! Tpl2 %i (%g-%g) %i (%g-%g) %i -- I3=%g scoreTpl=%g\n",
        Tpl[0 * NbTpl + i], ProbArrowhead2[0 * NbTpl + i],
        ProbArrowhead2[1 * NbTpl + i], Tpl[1 * NbTpl + i],
        ProbArrowhead2[2 * NbTpl + i], ProbArrowhead2[3 * NbTpl + i],
        Tpl[2 * NbTpl + i], I3[i], scoreTpl[i]);
  }

  printf("!END LV=%i count=%i deg=%i maxTpl=%i maxscoreTpl=%g \n", LV, count,
      deg, maxTpl, maxscoreTpl);
#endif  // _MY_PRINT_

  delete (ProbArrowhead2);
  delete (scoreTpl);
  delete (orderTpl);

  return 0;
}

double* miic::reconstruction::getOrientTplLVDegPropag(int nbrTpl,
    int* ptrAllTpl, double* ptrAllI3, int LV, int isDeg, int isPropag,
    int halfVStructures) {
  int nbrRetProbaValues = -1;  // Nbr proba to return
  double*
      ptrRetProbValues;  // To return ProbArrowhead
                         // >> ProbArrowhead[1][i] 1<-->2 ProbArrowhead[2][i]
                         // >> ProbArrowhead[3][i] 2<-->3 ProbArrowhead[4][i]
                         //
                         // ProbArrowhead > 0.5: an arrow head  (>)
                         // ProbArrowhead < 0.5: an arrow tail  (-)
  nbrRetProbaValues = (4 * nbrTpl);
  ptrRetProbValues = new double[nbrRetProbaValues];
  // Initialise the arrowhead probabilities to 0.5
  for (int i = 0; i < nbrTpl; i++) {
    for (int j = 0; j < 4; j++) {
      ptrRetProbValues[i + j * nbrTpl] = 0.5;
    }
  }
  // Iteratively converge towards partially oriented graphs including possible
  // latent variables and Propagation/Non-Propagation rules.
  OrientTpl_LV_Deg_Propag(nbrTpl, ptrAllTpl, ptrAllI3, ptrRetProbValues, LV,
      isDeg, isPropag, halfVStructures);

  return ptrRetProbValues;
}
