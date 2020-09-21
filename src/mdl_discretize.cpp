#include <Rcpp.h>

#include <algorithm>
#include <cfloat>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

double ramanujan(int n) {
  // Returns log(fac(n)) with Ramanujan's approximation.
  if (n == 0) {
    return (0);
  }
  double N = n * log(1.0 * n) - n +
             log(1.0 * n * (1 + 4 * n * (1 + 2 * n))) / 6 + log(M_PI) / 2L;
  return N;
}

double compute_parametric_complexity(int n, int K, double** sc_look) {
  if (sc_look[n - 1][K - 1] != 0) {
    return (sc_look[n - 1][K - 1]);
  }

  double res;
  if (K == 1) {
    res = 1;
  } else if (K == 2) {
    if (n < 1000) {
      res = 0;
      for (int i = 0; i <= n; i++) {
        int h1 = i;
        int h2 = n - h1;
        res = res + exp(ramanujan(n) - ramanujan(h1) - ramanujan(h2)) *
                        pow((1.0 * h1 / n), h1) * pow((1.0 * h2 / n), h2);
      }
    } else {
      res = sqrt((n * M_PI) / 2) * exp(sqrt(8. / (9. * n * M_PI)) +
                                       (3 * M_PI - 16) / (36. * n * M_PI));
    }
  } else {
    res = compute_parametric_complexity(n, K - 1, sc_look) +
          1.0 * n / (K - 2) * compute_parametric_complexity(n, K - 2, sc_look);
  }

  sc_look[n - 1][K - 1] = res;
  return (res);
}

/* Copyright (c) 2011, Paul Mineiro
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without
 *    modification, are permitted provided that the following conditions
 *    are met:
 *
 *    * Redistributions of source code must retain the above copyright notice,
 *        this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *    * Neither the name of Paul Mineiro nor the names of its contributors
 *        may be used to endorse or promote products derived from this software
 *        without specific prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 *    THE POSSIBILITY OF SUCH DAMAGE.
 *
 * https://github.com/etheory/fastapprox/blob/master/fastapprox */
static inline float fastlog2(float x) {
  union {
    float f;
    uint32_t i;
  } vx = {x};
  union {
    uint32_t i;
    float f;
  } mx = {(vx.i & 0x007FFFFF) | 0x3f000000};
  float y = vx.i;
  y *= 1.1920928955078125e-7f;

  return y - 124.22551499f - 1.498030302f * mx.f -
         1.72587999f / (0.3520887068f + mx.f);
}
static inline float fastlog(float x) { return 0.69314718f * fastlog2(x); }
/* End fastapprox */

double* compute_B(int K, int e, double* candidate_cut_points, double* myDist,
    int n, double epsilon, double xmin, double* B_previous_K, double** sc_look,
    int* count_for_e, int E) {
  double* res = new double[2];

  int ne = count_for_e[e - 1];
  double B;
  if (K == 1) {
    B = -ne * (log(epsilon * ne /
                   ((candidate_cut_points[e - 1] - (xmin - epsilon / 2)) * n)));
    res[0] = 0;
    res[1] = B;
    return (res);
  } else {
    double R_ne_K = compute_parametric_complexity(ne, K, sc_look);
    double min = DBL_MAX;
    int minidx(0);
    int nep(0);
    // Loop on ep
    for (int ep = K - 1; ep < e - 1; ep++) {
      nep = count_for_e[ep - 1];
      double R_nep_prevK = compute_parametric_complexity(nep, K - 1, sc_look);
      double B = B_previous_K[ep - 1] -
                 (ne - nep) * fastlog((epsilon * (ne - nep)) /
                                      ((candidate_cut_points[e - 1] -
                                           candidate_cut_points[ep - 1]) *
                                          n)) +
                 fastlog((R_ne_K / R_nep_prevK) * ((E - K + 2) / (K - 1)));

      if (B < min) {
        min = B;
        minidx = ep;
      }
    }
    int best_ep = minidx;
    double B_current_K = min;
    res[0] = best_ep;
    res[1] = B_current_K;
    return (res);
  }
}

}  // anonymous namespace

// [[Rcpp::export]]
List mydiscretizeMDL(SEXP RmyDist, SEXP RmaxBins) {
  std::vector<double> myDistVec = Rcpp::as<std::vector<double> >(RmyDist);
  int maxBins = Rcpp::as<int>(RmaxBins);
  double epsilon = 0.001;
  int n = myDistVec.size();
  sort(myDistVec.begin(), myDistVec.end());

  double* myDist = &myDistVec[0];

  double* middle_points = new double[n - 1];
  int n_unique_points(0);

  for (int i = 1; i < n - 1; i++) {
    if (myDist[i] != myDist[i + 1]) {
      middle_points[n_unique_points] = (myDist[i] + myDist[i + 1]) / 2;
      n_unique_points++;
    }
  }
  int n_candidate_cut_points = n_unique_points + 1;
  double* candidate_cut_points = new double[n_candidate_cut_points];

  for (int i = 0; i < n_unique_points; i++) {
    candidate_cut_points[i] = middle_points[i];
  }

  int* count_for_e = new int[n_candidate_cut_points];
  for (int e = 0; e < n_unique_points; e++) {
    int count = 0;
    for (int i = 0; i < n; i++) {
      if (myDist[i] < candidate_cut_points[e]) {
        count++;
      } else
        break;
    }
    count_for_e[e] = count;
  }

  double xmin = myDist[0] - epsilon;
  double xmax = myDist[n - 1] + epsilon;

  double** sc_look = new double*[n];
  for (int i = 0; i < n; i++) {
    sc_look[i] = new double[maxBins];
    for (int k = 0; k < maxBins; k++) {
      sc_look[i][k] = 0.;
    }
  }

  double** dyntable = new double*[n_candidate_cut_points];
  int** dyntable_trace = new int*[n_candidate_cut_points];
  for (int i = 0; i < n_candidate_cut_points; i++) {
    dyntable[i] = new double[maxBins];
    dyntable_trace[i] = new int[maxBins];
    for (int k = 0; k < maxBins; k++) {
      dyntable[i][k] = 0.;
      dyntable_trace[i][k] = 0;
    }
  }

  for (int K = 1; K <= maxBins; K++) {
    double* previous_column = new double[n_candidate_cut_points];
    for (int i = 0; i < n_candidate_cut_points;
         i++) {  // Retrieve previous column
      if (K == 1)
        previous_column[i] = 0;
      else
        previous_column[i] = dyntable[i][K - 2];
    }

    for (int e = K; e < n_candidate_cut_points; e++) {
      double* res = compute_B(K, e, candidate_cut_points, myDist, n, epsilon,
          xmin, previous_column, sc_look, count_for_e, n_candidate_cut_points);
      dyntable[e - 1][K - 1] = res[1];
      dyntable_trace[e - 1][K - 1] = res[0];
    }
  }
  int best_K = 0;
  double best_score = DBL_MAX;
  for (int i = 0; i < maxBins; i++) {
    if (dyntable[n_candidate_cut_points - 2][i] < best_score) {
      best_score = dyntable[n_candidate_cut_points - 2][i];
      best_K = i + 1;
    }
  }
  double* cut_points = new double[best_K + 2];
  int K = best_K;
  int curr_idx = dyntable_trace[n_candidate_cut_points - 2][K - 1];
  cut_points[K] = candidate_cut_points[curr_idx - 1];
  int prev_idx = curr_idx;
  K--;
  while (K > 0) {
    curr_idx = dyntable_trace[prev_idx][K];
    cut_points[K] = candidate_cut_points[curr_idx - 1];
    prev_idx = curr_idx;
    K--;
  }
  cut_points[0] = xmin;
  cut_points[best_K + 1] = xmax;

  std::vector<double> values(cut_points, cut_points + best_K + 2);

  NumericMatrix Rdyntable(n_candidate_cut_points, maxBins);
  NumericMatrix Rdyntable_trace(n_candidate_cut_points, maxBins);
  for (int i = 0; i < n_candidate_cut_points; i++) {
    for (int j = 0; j < maxBins; j++) {
      Rdyntable[i + j * n_candidate_cut_points] = dyntable[i][j];
      Rdyntable_trace[i + j * n_candidate_cut_points] = dyntable_trace[i][j];
    }
  }

  // structure the output
  List result = List::create(
      _["cutpoints"] = values
  );

  return result;
}
