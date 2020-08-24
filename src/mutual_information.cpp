#include "mutual_information.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <iostream>

#include "compute_info.h"
#include "utilities.h"

//#define _MY_DEBUG_NEW 1
//#define _MY_DEBUG_NEW_UI_2 1
//#define _MY_DEBUG_NEW_UI 1
//#define _MY_DEBUG_NEW_OPTFUN 1
//#define DEBUG_JOINT 1
//#define _MY_DEBUG_NEW_UI_NML 1
//#define _MY_DEBUG_NEW_3p 1
//#define _MY_DEBUG_NEW_3p_f 1
//#define DEBUG 1

namespace miic {
namespace computation {

using std::vector;
using namespace miic::structure;
using namespace miic::utility;

// INPUT:
// memory_cuts: vector length n with recorded recursvely best cuts,
// possible values 0...n :
// 0->stops (one bins); -k -> stops (two bin ([0 k-1][k ..];
// k->continute [.. k-1][k ...])
// OUTPUT:
// r : number cuts
// cut: vector with cuts point-> [0 cut[0]][cut[0]+1 cut[1]]...[cut[r-2]
// cut[r-1]]
int reconstruction_cut_coarse(const TempVector<int>& memory_cuts,
    const TempVector<int>& memory_cuts2, int np, int n, int *cut) {
  int ncuts = 0;
  int l, s;
  if (memory_cuts[np - 1] == 0) {
    cut[0] = n - 1;
    return 1;
  }

  l = memory_cuts[np - 1];
  while (l > 0) {
    ncuts++;
    l = memory_cuts[l - 1];
  }
  if (l < 0) ncuts++;

  cut[ncuts] = n - 1;  // conventional last cut
  l = ncuts - 1;
  s = memory_cuts[np - 1];
  cut[l] = memory_cuts2[np - 1];  // truly last cut
  l--;
  while (s > 0 && l >= 0) {
    cut[l] = memory_cuts2[s - 1];
    s = memory_cuts[s - 1];
    l--;
  }
  return ncuts + 1;  // number of levels (r)
}

// update datafactors of a variable from the cut positions vector <cut>
// INPUT:
// d: index of variable in datafactors
// varidx: index of variable in sortidx
void update_datafactors(const vector<vector<int>>& sortidx, int varidx,
    int** datafactors, int d, int n, int** cut) {
  int uu = 0;
  for (int j = 0; j < n; ++j) {
    int sjj = sortidx[varidx][j];
    if (j > cut[d][uu]) uu++;
    datafactors[d][sjj] = uu;
  }
  return;
}

// INPUT
// datafactors: [0, ]: x, [1, ]: y, [2 ... <n_ui> + 1, ]: {ui}
// exclude: ignore index <exclude> in the {ui}
// n: Number of samples
// n_ui: number of {ui}
// r: number of levels of each variable [x, y, {ui}]
// OUTPUT
// uiyxfactors: Joint datafactors. [ui, uiy, uix, uixy]
// ruiyx[0,1,2,3]: Number of joint levels. 0: u, 1: uy, 2: ux, 3: uyx
void jointfactors_uiyx(int** datafactors, int exclude, int n, int n_ui, int* r,
    int** uiyxfactors, int* ruiyx) {
  TempAllocatorScope scope;
  // Compute unique hash values for each sample in each of the joint spaces
  TempVector<int> hash_ui(n);
  TempVector<int> hash_uiy(n);
  TempVector<int> hash_uix(n);
  TempVector<int> hash_uiyx(n);
  for (int i = 0; i < n; ++i) {
    hash_ui[i] = 0;
    hash_uiy[i] = datafactors[1][i];
    hash_uix[i] = datafactors[0][i];
    hash_uiyx[i] = hash_uix[i] + hash_uiy[i] * r[0];

    int Pbin_ui = 1;
    for (int l = n_ui + 1; l >= 2; --l) {
      if (l == exclude) continue;

      int df = datafactors[l][i] * Pbin_ui;
      hash_uiyx[i] += df * r[1] * r[0];
      hash_uix[i] += df * r[0];
      hash_uiy[i] += df * r[1];
      hash_ui[i] += df;
      Pbin_ui *= r[l];
    }
  }

  bool too_many_levels = false;
  int n_joint_levels = 1;
  for (int i = 0; i < n_ui + 2 && !too_many_levels; ++i) {
    n_joint_levels *= r[i];
    if (n_joint_levels > 8 * n) {
      too_many_levels = true;  // Too large for the sparse vectors
    }
  }

  ruiyx[0] = 0;  // ui
  ruiyx[1] = 0;  // uiy
  ruiyx[2] = 0;  // uix
  ruiyx[3] = 0;  // uiyx
  if (!too_many_levels) {
    // Use large sparse vectors to have O(n) time complexity (no sort)
    TempVector<int> levels_ui(n_joint_levels);
    TempVector<int> levels_uiy(n_joint_levels);
    TempVector<int> levels_uix(n_joint_levels);
    TempVector<int> levels_uiyx(n_joint_levels);
    for (int i = 0; i < n; ++i) {
      levels_ui[hash_ui[i]] = 1;
      levels_uiy[hash_uiy[i]] = 1;
      levels_uix[hash_uix[i]] = 1;
      levels_uiyx[hash_uiyx[i]] = 1;
    }
    // Use ruiyx[0-3] as level indices, whose final values are the total numbers
    // of joint levels. Order of the levels follow the order of the hash values,
    // which are sorted automatically (as indices) with sparse vectors.
    for (int i = 0; i < n_joint_levels; ++i) {
      if (levels_ui[i] == 1) levels_ui[i] = ruiyx[0]++;
      if (levels_uiy[i] == 1) levels_uiy[i] = ruiyx[1]++;
      if (levels_uix[i] == 1) levels_uix[i] = ruiyx[2]++;
      if (levels_uiyx[i] == 1) levels_uiyx[i] = ruiyx[3]++;
    }

    for (int i = 0; i < n; ++i) {
      uiyxfactors[0][i] = levels_ui[hash_ui[i]];      // ui
      uiyxfactors[1][i] = levels_uiy[hash_uiy[i]];    // uiy
      uiyxfactors[2][i] = levels_uix[hash_uix[i]];    // uix
      uiyxfactors[3][i] = levels_uiyx[hash_uiyx[i]];  // uiyx
    }
  } else {
    // Fall back to O(nlog(n)) time complexity (sort)
    TempVector<int> orderSample_uix(n);
    std::iota(begin(orderSample_uix), end(orderSample_uix), 0);  // [0 to n - 1]
    TempVector<int> orderSample_uiyx(orderSample_uix);           // copy

    std::sort(begin(orderSample_uix), end(orderSample_uix),
        [&hash_uix](int a, int b) { return hash_uix[a] < hash_uix[b]; });
    std::sort(begin(orderSample_uiyx), end(orderSample_uiyx),
        [&hash_uiyx](int a, int b) { return hash_uiyx[a] < hash_uiyx[b]; });

    // hash_uix[a] < hash_uix[b] -> hash_ui[a] <= hash_ui[b]
    int hash_ui_prev = hash_ui[orderSample_uix[0]];
    int hash_uix_prev = hash_uix[orderSample_uix[0]];
    for (const auto index : orderSample_uix) {
      auto hash_ui_current = hash_ui[index];
      auto hash_uix_current = hash_uix[index];
      if (hash_ui_current > hash_ui_prev) ++ruiyx[0];
      if (hash_uix_current > hash_uix_prev) ++ruiyx[2];

      uiyxfactors[0][index] = ruiyx[0];  // ui
      uiyxfactors[2][index] = ruiyx[2];  // uix
      hash_ui_prev = hash_ui_current;
      hash_uix_prev = hash_uix_current;
    }
    // hash_uixy[a] < hash_uixy[b] -> hash_uiy[a] <= hash_uiy[b]
    int hash_uiy_prev = hash_uiy[orderSample_uiyx[0]];
    int hash_uiyx_prev = hash_uiyx[orderSample_uiyx[0]];
    for (const auto index : orderSample_uiyx) {
      auto hash_uiy_current = hash_uiy[index];
      auto hash_uiyx_current = hash_uiyx[index];
      if (hash_uiy_current > hash_uiy_prev) ++ruiyx[1];
      if (hash_uiyx_current > hash_uiyx_prev) ++ruiyx[3];

      uiyxfactors[1][index] = ruiyx[1];  // uiy
      uiyxfactors[3][index] = ruiyx[3];  // uiyx
      hash_uiy_prev = hash_uiy_current;
      hash_uiyx_prev = hash_uiyx_current;
    }
    // number of joint levels
    ++ruiyx[0];  // ui
    ++ruiyx[1];  // uiy
    ++ruiyx[2];  // uix
    ++ruiyx[3];  // uiyx
  }
  return;
}

// INPUT:
// datarank, datafactors, cut
// OUTPUT
// return joint datafactors ui , with number of levels rui
// entropy term Hui
void jointfactors_u(int** datafactors, int* ptrIdx, int n, int n_ui, int* r,
    int* ufactors, int* ru) {
  if (n_ui == 1) {
    for (int i = 0; i < n; ++i) {
      ufactors[i] = datafactors[ptrIdx[0]][i];
    }
    *ru = r[ptrIdx[0]];
    return;
  }
  TempAllocatorScope scope;
  // Compute unique hash value for each sample in the joint space
  TempVector<int> hash_u(n, 0);
  for (int i = 0; i < n; ++i) {
    int Pbin_ui = 1;
    for (int l = n_ui - 1; l >= 0; --l) {
      hash_u[i] += datafactors[ptrIdx[l]][i] * Pbin_ui;
      Pbin_ui *= r[ptrIdx[l]];
    }
  }

  bool too_many_levels = false;
  int n_joint_levels = 1;
  for (int i = 0; i < n_ui && !too_many_levels; ++i) {
    n_joint_levels *= r[ptrIdx[i]];
    if (n_joint_levels > 8 * n) {
      too_many_levels = true;  // Too large for the sparse vectors
    }
  }

  *ru = 0;
  if (!too_many_levels) {
    // Use large sparse vectors to have O(n) time complexity (no sort)
    TempVector<int> levels_ui(n_joint_levels);
    for (const auto h : hash_u)
      levels_ui[h] = 1;
    // Use *ru as level indices, whose final value is the total numbers
    // of joint levels. Order of the levels follow the order of the hash values,
    // which are sorted automatically (as indices) with the sparse vector.
    for (auto& l : levels_ui)
      if (l == 1) l = (*ru)++;

    for (int i = 0; i < n; ++i)
      ufactors[i] = levels_ui[hash_u[i]];
  } else {
    // Fall back to O(nlog(n)) time complexity (sort)
    TempVector<int> orderSample_u(n);
    std::iota(begin(orderSample_u), end(orderSample_u), 0);
    std::sort(begin(orderSample_u), end(orderSample_u),
        [&hash_u](int a, int b) { return hash_u[a] < hash_u[b]; });

    int hash_u_prev = hash_u[orderSample_u[0]];
    for (const auto index : orderSample_u) {
      auto hash_u_current = hash_u[index];
      if (hash_u_current > hash_u_prev) ++*ru;

      ufactors[index] = *ru;
      hash_u_prev = hash_u_current;
    }
    ++*ru;
  }
  return;
}

// rux -> 0:x,1;u,2:ux
vector<double> computeMI_knml(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, int n_eff, std::vector<double> sample_weights,
    std::shared_ptr<CtermCache> cache, int flag) {
  vector<double> I(2);

  int j, x, u, ux;

  double Hux = 0, Hu = 0, Hx = 0, SC = 0;

  double *nx = (double *)calloc(rux[0], sizeof(double));
  double *nu = (double *)calloc(rux[1], sizeof(double));
  double *nux = (double *)calloc(rux[2], sizeof(double));

  for (j = 0; j < n; j++) {
    nx[xfactors[j]] += sample_weights[j];
    nu[ufactors[j]] += sample_weights[j];
    nux[uxfactors[j]] += sample_weights[j];
  }

  for (x = 0; x < rux[0]; x++) {
    if (nx[x] > 0) {
      Hx -= nx[x] * log(nx[x]);
      if (flag == 0 || flag == 2)
        SC += cache->getLogC(fmax(1,int(nx[x]+0.5)), rux[1]);
    }
  }
  for (u = 0; u < rux[1]; u++) {
    if (nu[u] > 0){
      Hu -= nu[u] * log(nu[u]);
      if (flag == 0 || flag == 1)
        SC += cache->getLogC(fmax(1,int(nu[u]+0.5)), rux[0]);
    }
  }

  for (ux = 0; ux < rux[2]; ux++) {
    if (nux[ux] > 0) Hux -= nux[ux] * log(nux[ux]);
  }

  if (flag == 0) SC -= cache->getLogC(n_eff, rux[0]);
  if (flag == 0) SC -= cache->getLogC(n_eff, rux[1]);

  I[0] = cache->getLog(n_eff) + (Hu + Hx - Hux) / n_eff;

  if (flag == 0)
    I[1] = I[0] - 0.5 * SC / n_eff;
  else
    I[1] = I[0] - SC / n_eff;

  free(nx);
  free(nu);
  free(nux);

  return I;
}

vector<double> computeMI_knml(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, std::shared_ptr<CtermCache> cache, int flag) {
  vector<double> I(2);

  int j, x, u, ux;

  double Hux = 0, Hu = 0, Hx = 0, SC = 0;

  int *nx = (int *)calloc(rux[0], sizeof(int));
  int *nu = (int *)calloc(rux[1], sizeof(int));
  int *nux = (int *)calloc(rux[2], sizeof(int));

  for (j = 0; j < n; j++) {
    nx[xfactors[j]]++;
    nu[ufactors[j]]++;
    nux[uxfactors[j]]++;
  }

  for (x = 0; x < rux[0]; x++) {
    if (nx[x] > 0) Hx -= nx[x] * cache->getLog(nx[x]);
    if (flag == 0 || flag == 2)
      SC += cache->getLogC(nx[x], rux[1]);
  }
  for (u = 0; u < rux[1]; u++) {
    if (nu[u] > 0) Hu -= nu[u] * cache->getLog(nu[u]);
    if (flag == 0 || flag == 1)
      SC += cache->getLogC(nu[u], rux[0]);
  }

  for (ux = 0; ux < rux[2]; ux++) {
    if (nux[ux] > 0) Hux -= nux[ux] * cache->getLog(nux[ux]);
  }

  if (flag == 0) SC -= cache->getLogC(n, rux[0]);
  if (flag == 0) SC -= cache->getLogC(n, rux[1]);

  I[0] = cache->getLog(n) + (Hu + Hx - Hux) / n;

  if (flag == 0)
    I[1] = I[0] - 0.5 * SC / n;
  else
    I[1] = I[0] - SC / n;

  free(nx);
  free(nu);
  free(nux);

  return I;
}

// rux -> 0:x,1;u,2:ux
vector<double> computeMI_kmdl(int *xfactors, int *ufactors, int *uxfactors, int *rux,
    int n, std::shared_ptr<CtermCache> cache, int flag) {
  vector<double> I(2);

  int j, x, u, ux;

  double Hux = 0, Hu = 0, Hx = 0, SC = 0;

  int *nx = (int *)calloc(rux[0], sizeof(int));
  int *nu = (int *)calloc(rux[1], sizeof(int));
  int *nux = (int *)calloc(rux[2], sizeof(int));

  for (j = 0; j < n; j++) {
    nx[xfactors[j]]++;
    nu[ufactors[j]]++;
    nux[uxfactors[j]]++;
  }

  for (x = 0; x < rux[0]; x++) {
    if (nx[x] > 0) Hx -= nx[x] * cache->getLog(nx[x]);
  }
  for (u = 0; u < rux[1]; u++) {
    if (nu[u] > 0) Hu -= nu[u] * cache->getLog(nu[u]);
  }

  for (ux = 0; ux < rux[2]; ux++) {
    if (nux[ux] > 0) Hux -= nux[ux] * cache->getLog(nux[ux]);
  }

  SC = 0.5 * cache->getLog(n);
  if (flag == 0 || flag == 1) SC *= (rux[0] - 1);
  if (flag == 0 || flag == 2) SC *= (rux[1] - 1);

  I[0] = cache->getLog(n) + (Hu + Hx - Hux) / n;

  I[1] = I[0] - SC / n;

  free(nx);
  free(nu);
  free(nux);

  return I;
}

}  // namespace computation
}  // namespace miic
