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
int reconstruction_cut_coarse(TempVector<int>& memory_cuts,
    TempVector<int>& memory_cuts2, int np, int n, int *cut) {
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
void update_datafactors(
    vector<vector<int> > &sortidx, int varidx, int** datafactors,
    int d, int n, int **cut) {

  int j, sjj, uu = 0;
  for (j = 0; j <= n - 1; j++) {
    sjj = sortidx[varidx][j];
    if (j > cut[d][uu]) uu++;
    datafactors[d][sjj] = uu;
  }
  return;
}

// INPUT:
// datarank, datafactors
// take 0 1 as x y and 2..<Mui>+1 for {u} indexes
// ignore index <dui> in the u
// OUTPUT
// return joint datafactors ui and uiy uix uixy, with number of levels
// ruiyx[0,1,2,3] ruiyx -> 0:u,1:uy,2:ux,3:uyx
void jointfactors_uiyx(int **datafactors, int dui, int n, int Mui, int *r,
    int **uiyxfactors, int *ruiyx) {
  TempAllocatorScope scope;

  bool tooManyLevels = false;
  // From single datafactors to joint datafactors ui; (ui,y);(ui,x);(ui,y,x)
  int n_joint_levels = 1;
  for (int i = 0; i < Mui + 2 && !tooManyLevels; ++i) {
    n_joint_levels *= r[i];
    if (n_joint_levels > 8 * n) {
      tooManyLevels = true;
    }
  }

  TempVector<int> vecZeroOnesui;
  TempVector<int> vecZeroOnesuiy;
  TempVector<int> vecZeroOnesuix;
  TempVector<int> vecZeroOnesuiyx;
  // declaration
  if (!tooManyLevels) {
    vecZeroOnesui   = TempVector<int>(n_joint_levels);
    vecZeroOnesuiy  = TempVector<int>(n_joint_levels);
    vecZeroOnesuix  = TempVector<int>(n_joint_levels);
    vecZeroOnesuiyx = TempVector<int>(n_joint_levels);
  }

  TempVector<int> datauix(n);
  TempVector<int> datauiy(n);
  TempVector<int> datauiyx(n);
  TempVector<int> dataui(n, 0);
  TempVector<int> orderSample_ux(n);
  TempVector<int> orderSample_uyx(n);
  for (int i = 0; i < n; ++i) {
    datauiyx[i] = datafactors[0][i];  // yx
    datauiyx[i] += datafactors[1][i] * r[0];
    datauix[i] = datafactors[0][i];  // x
    datauiy[i] = datafactors[1][i];  // y

    int Pbin_ui = 1;
    for (int l = Mui + 1; l >= 2; --l) {
      if (l == dui) continue;

      int df = datafactors[l][i] * Pbin_ui;
      datauiyx[i] += df * r[1] * r[0];
      datauix[i] += df * r[0];
      datauiy[i] += df * r[1];
      dataui[i] += df;
      Pbin_ui *= r[l];
    }
    // init
    if (!tooManyLevels) {
      vecZeroOnesui[dataui[i]] = 1;
      vecZeroOnesuiy[datauiy[i]] = 1;
      vecZeroOnesuix[datauix[i]] = 1;
      vecZeroOnesuiyx[datauiyx[i]] = 1;
    } else {
      orderSample_ux[i] = i;
      orderSample_uyx[i] = i;
    }
  }

  if (!tooManyLevels) {
    // create the vector for storing the position
    int pos = 0;
    int pos1 = 0;
    int pos2 = 0;
    int pos3 = 0;
    for (int i = 0; i < n_joint_levels; ++i) {
      if (vecZeroOnesui[i] == 1)
        vecZeroOnesui[i] = pos++;

      if (vecZeroOnesuiy[i] == 1)
        vecZeroOnesuiy[i] = pos1++;

      if (vecZeroOnesuix[i] == 1)
        vecZeroOnesuix[i] = pos2++;

      if (vecZeroOnesuiyx[i] == 1)
        vecZeroOnesuiyx[i] = pos3++;
    }

    ruiyx[0] = pos;   // ui
    ruiyx[1] = pos1;  // uiy
    ruiyx[2] = pos2;  // uix
    ruiyx[3] = pos3;  // uiyx

    for (int i = 0; i < n; ++i) {
      uiyxfactors[0][i] = vecZeroOnesui[dataui[i]];  // ui
      uiyxfactors[1][i] = vecZeroOnesuiy[datauiy[i]];  // uiy
      uiyxfactors[2][i] = vecZeroOnesuix[datauix[i]];  // uix
      uiyxfactors[3][i] = vecZeroOnesuiyx[datauiyx[i]];  // uiyx
    }
  } else {
    std::sort(begin(orderSample_ux), end(orderSample_ux),
        [&datauix](int a, int b) { return datauix[a] < datauix[b]; });
    std::sort(begin(orderSample_uyx), end(orderSample_uyx),
        [&datauiyx](int a, int b) { return datauiyx[a] < datauiyx[b]; });

    // joint datafactors without gaps
    // compute term constant H(y,Ui) and H(Ui)
    int ix, iyx, iix, iiyx;

    ruiyx[0] = 0;  // ui
    ruiyx[1] = 0;  // uiy
    ruiyx[2] = 0;  // uix
    ruiyx[3] = 0;  // uiyx

    ix = orderSample_ux[0];
    iyx = orderSample_uyx[0];
    uiyxfactors[0][ix] = 0;   // ui
    uiyxfactors[1][iyx] = 0;  // uiy
    uiyxfactors[2][ix] = 0;   // uix
    uiyxfactors[3][iyx] = 0;  // uiyx

    for (int i = 1; i < n; ++i) {
      iix = ix;
      iiyx = iyx;
      ix = orderSample_ux[i];
      iyx = orderSample_uyx[i];

      if (dataui[ix] > dataui[iix]) {
        ruiyx[0]++;
      }
      uiyxfactors[0][ix] = ruiyx[0];  // ui

      if (datauiy[iyx] > datauiy[iiyx]) {
        ruiyx[1]++;
      }
      uiyxfactors[1][iyx] = ruiyx[1];  // uiy

      if (datauix[ix] > datauix[iix]) {
        ruiyx[2]++;
      }
      uiyxfactors[2][ix] = ruiyx[2];  // uix

      if (datauiyx[iyx] > datauiyx[iiyx]) {
        ruiyx[3]++;
      }
      uiyxfactors[3][iyx] = ruiyx[3];  // uiyx
    }
    // number joint levels
    ruiyx[0]++;  // ui
    ruiyx[1]++;  // uiy
    ruiyx[2]++;  // uix
    ruiyx[3]++;  // uiyx
  }
  return;
}

// INPUT:
// datarank, datafactors, cut
// OUTPUT
// return joint datafactors ui , with number of levels rui
// entropy term Hui
void jointfactors_u(int **datafactors, int *ptrIdx, int n, int Mui, int *r,
    int *ufactors, int *ru) {
  TempAllocatorScope scope;
  // from cuts to joint datafactors
  int jj;
  int j, l;
  if (Mui == 1) {
    for (jj = 0; jj <= n - 1; jj++) {
      ufactors[jj] = datafactors[ptrIdx[0]][jj];
    }
    *ru = r[ptrIdx[0]];
    return;
  }
  // update joint datafactors (with gaps) ui and (ui,y)
  int df, Pbin_ui;
  TempVector<int> datau(n, 0);

  bool tooManyLevels = false;

  int nbrLevelsJoint = 1;
  for (l = 0; l < Mui && !tooManyLevels; l++) {
    nbrLevelsJoint *= r[ptrIdx[l]];
    if (nbrLevelsJoint > 8 * n) tooManyLevels = true;
  }
  // decl
  TempVector<int> vecZeroOnesui;
  TempVector<int> orderSample_u(n);
  if (!tooManyLevels) {
    vecZeroOnesui = TempVector<int>(nbrLevelsJoint + 1);
  }

  for (jj = 0; jj <= n - 1; jj++) {
    Pbin_ui = 1;
    for (l = Mui - 1; l >= 0; l--) {
      df = datafactors[ptrIdx[l]][jj] * Pbin_ui;
      datau[jj] += df;
      Pbin_ui *= r[ptrIdx[l]];
    }
    // init
    if (!tooManyLevels) {
      vecZeroOnesui[datau[jj]] = 1;
    } else {
      orderSample_u[jj] = jj;
    }
  }

  if (!tooManyLevels) {
    // create the vector for storing the position
    int pos = 0;

    for (jj = 0; jj <= nbrLevelsJoint; jj++) {
      if (vecZeroOnesui[jj] == 1) {
        vecZeroOnesui[jj] = pos;
        pos += 1;
      }
    }

    *ru = pos;  // ui

    for (j = 0; j <= n - 1; j++) {
      ufactors[j] = vecZeroOnesui[datau[j]];  // ui
    }

  } else {
    std::sort(begin(orderSample_u), end(orderSample_u),
        [&datau](int a, int b) { return datau[a] < datau[b]; });

    int ix, iix;
    ix = orderSample_u[0];
    *ru = 0;
    ufactors[ix] = 0;  // ui

    for (j = 0; j < n - 1; j++) {
      iix = ix;
      ix = orderSample_u[j + 1];

      if (datau[ix] > datau[iix]) {
        *ru = *ru + 1;
      }
      ufactors[ix] = *ru;  // ui
    }
    *ru = *ru + 1;
  }

#if DEBUG_JOINT
  Rprintf("\nj -> datau[j]\n");
  for (j = 0; j <= n - 1; j++) {
    Rprintf("%d -> %d \n", j, datau[j]);
  }
  R_FlushConsole();
#endif

#if DEBUG_JOINT
  Rprintf("j,orderSample[j+1],sampleKey[j+1],datau[orderSample[j+1]]\n");
  for (j = 0; j <= n - 1; j++) {
    Rprintf("%d: %d %d, %d \n", j, orderSample[j + 1], sampleKey[j + 1],
        datau[orderSample[j + 1]]);
  }
  R_FlushConsole();
#endif

  // joint datafactors without gaps
  // compute term constant H(y,Ui) and H(Ui)

#if DEBUG
  Rprintf("Hu=%lf\n", *Hu);
#endif

#if DEBUG_JOINT
  Rprintf("j,uiyfactors[0][i] (ru=%d)\n", *ru);
  for (j = 0; j <= n - 1; j++) {
    Rprintf("%d-> %d\n", j, ufactors[j]);
  }
  R_FlushConsole();
#endif

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
