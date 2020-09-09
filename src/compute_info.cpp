#include "compute_info.h"

#include <algorithm>  // std::sort, std::minmax
#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>  // std::numeric_limits
#include <tuple>   // std::tie
#include <unordered_set>

#include "structure.h"
#include "utilities.h"

constexpr int MDL = 0;

// cplx = 0 --> MDL, cplx = 1 --> NML
namespace miic {
namespace computation {

using namespace miic::structure;
using namespace miic::utility;
using std::vector;

double* computeCondMutualInfoDiscrete(const TempGrid2d<int>& data,
    const TempVector<int>& levels, const TempVector<int>& var_idx,
    const TempVector<double>& weights, int cplx,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_samples = data.n_cols();
  int n_nodes = data.n_rows();
  int id_x = var_idx[0], id_y = var_idx[1];
  int rx = levels[id_x], ry = levels[id_y];

  TempVector<int> hash_ui(n_samples + 1);
  TempVector<int> hash_uiy(n_samples + 1);
  TempVector<int> hash_uiyx(n_samples + 1);
  // Compute unique hash values for each sample in each of the joint spaces
  for (int k = 0; k < n_samples; k++) {
    hash_uiyx[k] = data(id_y, k) * rx + data(id_x, k);
    hash_uiy[k] = data(id_y, k) * rx;
    hash_ui[k] = 0;

    int PBin = rx;        // rx
    for (int j = 2; j < n_nodes; j++) {
      PBin *= levels[var_idx[j - 1]];
      int increment = data(var_idx[j], k) * PBin;
      hash_uiyx[k] += increment;     // Lxyui
      hash_uiy[k] += increment;      // Lyui
      hash_ui[k] += increment;       // Lui
    }
  }
  int n_joint_levels = 1;
  for (int i = 0; i < n_nodes; ++i) {
    n_joint_levels *= levels[var_idx[i]];
  }
  // For termination of counts
  hash_uiyx[n_samples] = n_joint_levels;
  hash_uiy[n_samples] = n_joint_levels;
  hash_ui[n_samples] = n_joint_levels;

  TempVector<int> order(n_samples + 1);
  // [0 to n_samples]
  std::iota(begin(order), end(order), 0);
  std::sort(begin(order), end(order),
      [&hash_uiyx](int a, int b) { return hash_uiyx[a] < hash_uiyx[b]; });

  int Nyui{0}, Nui{0}, Ntot{0};
  double info_xui_y{0}, info_yui_x{0}, info_ui_y{0}, info_ui_x{0};
  double logC_xui_y{0}, logC_yui_x{0}, logC_ui_y{0}, logC_ui_x{0};

  TempVector<int> Nxui(rx, 0);
  TempVector<int> Nx(rx, 0);
  TempVector<int> Ny(ry, 0);
  // initialization of counts and mutual infos & logCs
  int Lxyui = hash_uiyx[order[0]];  // previous hash_xyui
  int Lyui = hash_uiy[order[0]];    // previous hash_yui
  int Lui = hash_ui[order[0]];      // previous hash_ui
  int X = data(id_x, order[0]);
  int Y = data(id_y, order[0]);
  double Pxyui = 0;
  // make the counts and compute mutual infos & logCs
  for (int k = 0; k <= n_samples; k++) {
    int i = order[k];
    if (hash_uiyx[i] <= Lxyui) {
      Pxyui += weights[i];
      continue;
    }
    Lxyui = hash_uiyx[i];
    // X has changed, conclude on previous X
    int Nxyui = (int)Pxyui;
    if (Nxyui > 0) {
      Nui += Nxyui;
      Nyui += Nxyui;
      Nxui[X] += Nxyui;
      Nx[X] += Nxyui;
      Ny[Y] += Nxyui;
      Ntot += Nxyui;

      double NlogN = Nxyui * cache->getLog(Nxyui);
      info_xui_y += NlogN;
      info_yui_x += NlogN;
    }
    if (k < n_samples) {
      Pxyui = weights[i];  // reset cumulative weight to current sample
      X = data(id_x, i);
      Y = data(id_y, i);
    }

    if (hash_uiy[i] <= Lyui) continue;
    Lyui = hash_uiy[i];
    // Y has changed, conclude on previous Y
    if (Nyui > 0) {
      double NlogN = Nyui * cache->getLog(Nyui);
      info_yui_x -= NlogN;
      info_ui_y += NlogN;
      if (cplx != MDL) {
        logC_yui_x += cache->getLogC(Nyui, rx);
      }
      Nyui = 0;
    }

    if (hash_ui[i] <= Lui) continue;
    Lui = hash_ui[i];
    // ui has changed, conclude on previous ui
    for (int j = 0; j < rx; j++) {
      int Nxuij = Nxui[j];
      if (Nxuij > 0) {
        double NlogN = Nxuij * cache->getLog(Nxuij);
        info_xui_y -= NlogN;
        info_ui_x += NlogN;
        if (cplx != MDL) {
          logC_xui_y += cache->getLogC(Nxuij, ry);
        }
        Nxui[j] = 0;
      }
    }

    if (Nui > 0) {
      double NlogN = Nui * cache->getLog(Nui);
      info_ui_y -= NlogN;
      info_ui_x -= NlogN;
      if (cplx != MDL) {
        logC_ui_x += cache->getLogC(Nui, rx);
        logC_ui_y += cache->getLogC(Nui, ry);
      }
      Nui = 0;
    }
  }
  // increment for info for Nx[X] and Ny[Y] contributions
  for (int j = 0; j < rx; j++) {
    int Nxj = Nx[j];
    if (Nxj > 0) {
      double NlogN = Nxj * log(Nxj / (1.0 * Ntot));
      info_yui_x -= NlogN;
      info_ui_x -= NlogN;
    }
  }
  for (int j = 0; j < ry; j++) {
    int Nyj = Ny[j];
    if (Nyj > 0) {
      double NlogN = Nyj * log(Nyj / (1.0 * Ntot));
      info_xui_y -= NlogN;
      info_ui_y -= NlogN;
    }
  }

  double info3xy_ui = 0.5 * (info_xui_y + info_yui_x);
  double info2xy_ui = 0.5 * (info_ui_y + info_ui_x);

  if (cplx == MDL) {
    int Prui = 1;
    double logN = cache->getLog(Ntot);
    for (int j = 2; j < n_nodes; j++)
      Prui *= levels[var_idx[j]];
    logC_xui_y = 0.5 * (ry - 1) * (rx * Prui - 1) * logN;
    logC_yui_x = 0.5 * (rx - 1) * (ry * Prui - 1) * logN;
    logC_ui_y = 0.5 * (ry - 1) * (Prui - 1) * logN;
    logC_ui_x = 0.5 * (rx - 1) * (Prui - 1) * logN;
  }

  double logC3xy_ui = 0.5 * (logC_xui_y + logC_yui_x);
  double logC2xy_ui = 0.5 * (logC_ui_y + logC_ui_x);

  double* res = new double[3];
  res[0] = Ntot;
  res[1] = info3xy_ui - info2xy_ui;
  res[2] = logC3xy_ui - logC2xy_ui;

  return res;
}

double* computeInfo3PointAndScoreDiscrete(const TempGrid2d<int>& data,
    const TempVector<int>& levels, const TempVector<int>& var_idx,
    const TempVector<double>& weights, int cplx,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_samples = data.n_cols();
  int n_nodes_xyui = data.n_rows() - 1;
  int id_x = var_idx[0], id_y = var_idx[1], id_z = var_idx.back();
  int rx = levels[id_x], ry = levels[id_y], rz = levels[id_z];

  // [N+1][4] [sample id (in the list with na), Lxyui, Lyui, Lui]
  TempGrid2d<int> sample(n_samples + 1, 4);
  for (int i = 0; i < n_samples; i++) {
    sample(i, 0) = i;  // sample number
  }

  // compute Lxyui, Lyui, Lui indices for later counting purpose
  for (int k = 0; k < n_samples; k++) {
    int PBin = rx;        // rx
    sample(k, 1) = data(id_y, k) * rx + data(id_x, k); // Lxyui
    sample(k, 2) = data(id_y, k) * rx;   // Lyui initialisation
    sample(k, 3) = 0;           // Lui initialisation

    for (int j = 2; j < n_nodes_xyui; j++) {
      PBin *= levels[var_idx[j - 1]];
      int increment = data(var_idx[j], k) * PBin;
      sample(k, 1) += increment;  // Lxyui
      sample(k, 2) += increment;  // Lyui
      sample(k, 3) += increment;  // Lui
    }
  }
  int n_joint_levels = 1;
  for (int i = 0; i < n_nodes_xyui; ++i) {
    n_joint_levels *= levels[var_idx[i]];
  }

  // extra sample id (useful for termination of counts)
  sample(n_samples, 0) = n_samples;
  // max Lxyui (useful for termination of counts)
  sample(n_samples, 1) = n_joint_levels;
  // max Lyui  (useful for termination of counts)
  sample(n_samples, 2) = n_joint_levels;
  // max Lui   (useful for termination of counts)
  sample(n_samples, 3) = n_joint_levels;

  TempVector<int> order(n_samples + 1);
  // [0 to n_samples_non_na]
  std::iota(begin(order), end(order), 0);
  std::sort(begin(order), end(order),
      [&sample](int a, int b) { return sample(a, 1) < sample(b, 1); });

  // Initialze variables
  int Lxyui = sample(order[0], 1);  // min Lxyui
  int Lyui = sample(order[0], 2);   // min Lyui
  int Lui = sample(order[0], 3);    // min Lui

  int X = data(id_x, order[0]);
  int Y = data(id_y, order[0]);
  int Z = data(id_z, order[0]);
  // to terminate loop properly below
  sample(order[n_samples], 0) = order[0];

  double info_xui_y{0}, info_yui_x{0}, info_ui_y{0}, info_ui_x{0};
  double logC_xui_y{0}, logC_yui_x{0}, logC_ui_y{0}, logC_ui_x{0};

  double info_xuiz_y{0}, info_yuiz_x{0}, info_uiz_y{0}, info_uiz_x{0};
  double logC_xuiz_y{0}, logC_yuiz_x{0}, logC_uiz_y{0}, logC_uiz_x{0};

  double info_xui_z{0}, info_yui_z{0}, info_ui_z{0};
  double logC_xui_z{0}, logC_yui_z{0}, logC_ui_z{0};

  double info_xy_ui{0}, info_xy_uiz{0}, info_xz_ui{0}, info_yz_ui{0};
  double logC_xy_ui{0}, logC_xy_uiz{0}, logC_xz_ui{0}, logC_yz_ui{0};

  // 11 test variables for counts, should all sum to nSample0
  int Nxyuis{0}, Nyuis{0}, Nxuis{0}, Nuis{0}, Nxs{0}, Nys{0};
  int Nzs{0}, Nuizs{0}, Nyuizs{0}, Nxyuizs{0}, Nxuizs{0};

  int Nyui{0}, Nui{0};
  int n_samples_zi{0};
  TempVector<int> Nxui(rx, 0);
  TempVector<int> Nx(rx, 0);
  TempVector<int> Ny(ry, 0);
  TempVector<double> Pxyuiz(rz, 0);
  TempVector<int> Nyuiz(rz, 0);
  TempVector<int> Nuiz(rz, 0);
  TempVector<int> Nz(rz, 0);
  Pxyuiz[Z] = weights[order[0]];
  TempGrid2d<int> Nxuiz(rx, rz, 0);  // [X][Z]
  // make the counts and compute mutual infos & logCs
  for (int k = 1; k <= n_samples; k++) {
    // check whether zi contains NA
    int i = order[k];
    if (k < n_samples) {
      Z = data(id_z, i);
      if (Z == -1) continue;
    }

    if (sample(i, 1) <= Lxyui) {
      Pxyuiz[Z] += weights[i];
      continue;
    }
    Lxyui = sample(i, 1);

    int Nxyuiz = 0;
    for (int l = 0; l < rz; l++) {
      double Pxyuizl = Pxyuiz[l];
      if (Pxyuizl <= 0) continue;

      int Nxyuizl = (int)Pxyuizl;

      if (Nxyuizl > 0) {
        double NlogN = Nxyuizl * cache->getLog(Nxyuizl);
        info_xuiz_y += NlogN;
        info_yuiz_x += NlogN;
        Nxyuiz += Nxyuizl;

        Nz[l] += Nxyuizl;
        Nuiz[l] += Nxyuizl;
        Nyuiz[l] += Nxyuizl;
        Nxuiz(X, l) += Nxyuizl;
      }
      Pxyuiz[l] = 0;
    }
    Pxyuiz[Z] = weights[i];

    if (Nxyuiz > 0) {
      double NlogN = Nxyuiz * cache->getLog(Nxyuiz);
      info_xui_y += NlogN;
      info_yui_x += NlogN;

      n_samples_zi += Nxyuiz;

      Nxyuizs += Nxyuiz;

      Nui += Nxyuiz;
      Nyui += Nxyuiz;
      Nxui[X] += Nxyuiz;

      Nx[X] += Nxyuiz;
      Ny[Y] += Nxyuiz;

      Nxyuis += Nxyuiz;
    }

    if (k < n_samples) X = data(id_x, i);
    if (k < n_samples) Y = data(id_y, i);

    if (sample(i, 2) == Lyui) continue;
    Lyui = sample(i, 2);

    if (Nyui > 0) {
      double NlogN = Nyui * cache->getLog(Nyui);
      info_yui_x -= NlogN;
      info_yui_z -= NlogN;
      info_ui_y += NlogN;
      if (cplx != MDL) {
        logC_yui_x += cache->getLogC(Nyui, rx);
        logC_yui_z += cache->getLogC(Nyui, rz);
      }

      for (int l = 0; l < rz; l++) {
        int Nyuizl = Nyuiz[l];
        if (Nyuizl > 0) {
          double NlogN = Nyuizl * cache->getLog(Nyuizl);
          info_yui_z += NlogN;
          info_uiz_y += NlogN;
          info_yuiz_x -= NlogN;
          if (cplx != MDL) {
            logC_yuiz_x += cache->getLogC(Nyuizl, rx);
          }
          Nyuizs += Nyuizl;
          Nyuiz[l] = 0;
        }
      }
      Nyuis += Nyui;
      Nyui = 0;
    }

    if (sample(i, 3) == Lui) continue;
    Lui = sample(i, 3);

    if (Nui <= 0) continue;

    double NlogN = Nui * cache->getLog(Nui);
    info_ui_x -= NlogN;
    info_ui_y -= NlogN;
    info_ui_z -= NlogN;
    if (cplx != MDL) {
      logC_ui_x += cache->getLogC(Nui, rx);
      logC_ui_y += cache->getLogC(Nui, ry);
      logC_ui_z += cache->getLogC(Nui, rz);
    }
    Nuis += Nui;
    Nui = 0;

    for (int l = 0; l < rz; l++) {
      int Nuizl = Nuiz[l];
      if (Nuizl > 0) {
        double NlogN = Nuizl * cache->getLog(Nuizl);
        info_ui_z += NlogN;
        info_uiz_x -= NlogN;
        info_uiz_y -= NlogN;
        if (cplx != MDL) {
          logC_uiz_x += cache->getLogC(Nuizl, rx);
          logC_uiz_y += cache->getLogC(Nuizl, ry);
        }
        Nuizs += Nuizl;
        Nuiz[l] = 0;
      }
    }

    for (int j = 0; j < rx; j++) {
      int Nxuij = Nxui[j];
      if (Nxuij == 0) continue;

      double NlogN = Nxuij * cache->getLog(Nxuij);
      info_xui_y -= NlogN;
      info_xui_z -= NlogN;
      info_ui_x += NlogN;
      if (cplx != MDL) {
        logC_xui_y += cache->getLogC(Nxuij, ry);
        logC_xui_z += cache->getLogC(Nxuij, rz);
      }
      Nxuis += Nxuij;
      Nxui[j] = 0;

      for (int l = 0; l < rz; l++) {
        int Nxuizjl = Nxuiz(j, l);
        if (Nxuizjl == 0) continue;

        double NlogN = Nxuizjl * cache->getLog(Nxuizjl);
        info_xui_z += NlogN;
        info_uiz_x += NlogN;
        info_xuiz_y -= NlogN;
        if (cplx != MDL) {
          logC_xuiz_y += cache->getLogC(Nxuizjl, ry);
        }
        Nxuizs += Nxuizjl;
        Nxuiz(j, l) = 0;
      }
    }
  }
  // increment info with Nx[X], Ny[Y] and Nz[Z] contributions
  for (int j = 0; j < rx; j++) {
    int Nxj = Nx[j];
    if (Nxj > 0) {
      double NlogN = Nxj * log(Nxj / (1.0 * n_samples_zi));
      info_yui_x -= NlogN;
      info_ui_x -= NlogN;
      info_uiz_x -= NlogN;
      info_yuiz_x -= NlogN;
      Nxs += Nxj;
      Nx[j] = 0;
    }
  }
  for (int j = 0; j < ry; j++) {
    int Nyj = Ny[j];
    if (Nyj > 0) {
      double NlogN = Nyj * log(Nyj / (1.0 * n_samples_zi));
      info_xui_y -= NlogN;
      info_ui_y -= NlogN;
      info_uiz_y -= NlogN;
      info_xuiz_y -= NlogN;
      Nys += Nyj;
      Ny[j] = 0;
    }
  }
  for (int l = 0; l < rz; l++) {
    int Nzl = Nz[l];
    if (Nzl > 0) {
      double NlogN = Nzl * log(Nzl / (1.0 * n_samples_zi));
      info_xui_z -= NlogN;
      info_yui_z -= NlogN;
      info_ui_z -= NlogN;
      Nzs += Nzl;
      Nz[l] = 0;
    }
  }
  // NI(xy|ui)
  double info3xy_ui = 0.5 * (info_xui_y + info_yui_x);
  double info2xy_ui = 0.5 * (info_ui_y + info_ui_x);

  // check maximum mutual infos - cplx terms
  int Prui{1};
  double logN;
  if (cplx == MDL) {
    Prui = 1;
    logN = cache->getLog(n_samples_zi);
    for (int j = 2; j < n_nodes_xyui; j++)
      Prui *= levels[var_idx[j]];
    logC_xui_y = 0.5 * (ry - 1) * (rx * Prui - 1) * logN;
    logC_yui_x = 0.5 * (rx - 1) * (ry * Prui - 1) * logN;
    logC_ui_y = 0.5 * (ry - 1) * (Prui - 1) * logN;
    logC_ui_x = 0.5 * (rx - 1) * (Prui - 1) * logN;
    logC_uiz_y = 0.5 * (ry - 1) * (rz * Prui - 1) * logN;
    logC_yui_z = 0.5 * (rz - 1) * (ry * Prui - 1) * logN;
    logC_ui_z = 0.5 * (rz - 1) * (Prui - 1) * logN;
  }
  double logC3xy_ui = 0.5 * (logC_xui_y + logC_yui_x);
  double logC2xy_ui = 0.5 * (logC_ui_y + logC_ui_x);
  // NI(yz|ui)
  double info3yz_ui = 0.5 * (info_uiz_y + info_yui_z);
  double info2yz_ui = 0.5 * (info_ui_y + info_ui_z);

  double logC3yz_ui = 0.5 * (logC_uiz_y + logC_yui_z);
  double logC2yz_ui = 0.5 * (logC_ui_y + logC_ui_z);
  // NI(xz|ui)
  double info3xz_ui = 0.5 * (info_xui_z + info_uiz_x);
  double info2xz_ui = 0.5 * (info_ui_z + info_ui_x);

  if (cplx == MDL) {
    logC_uiz_x = 0.5 * (rx - 1) * (rz * Prui - 1) * logN;
    logC_xui_z = 0.5 * (rz - 1) * (rx * Prui - 1) * logN;
  }
  double logC3xz_ui = 0.5 * (logC_uiz_x + logC_xui_z);
  double logC2xz_ui = 0.5 * (logC_ui_x + logC_ui_z);
  // NI(xy|uiz)
  double info3xy_uiz = 0.5 * (info_xuiz_y + info_yuiz_x);
  double info2xy_uiz = 0.5 * (info_uiz_y + info_uiz_x);

  if (cplx == MDL) {
    Prui *= rz;
    logC_xuiz_y = 0.5 * (ry - 1) * (rx * Prui - 1) * logN;
    logC_yuiz_x = 0.5 * (rx - 1) * (ry * Prui - 1) * logN;
    logC_uiz_y = 0.5 * (ry - 1) * (Prui - 1) * logN;
    logC_uiz_x = 0.5 * (rx - 1) * (Prui - 1) * logN;
  }
  double logC3xy_uiz = 0.5 * (logC_xuiz_y + logC_yuiz_x);
  double logC2xy_uiz = 0.5 * (logC_uiz_y + logC_uiz_x);
  info_xy_ui = info3xy_ui - info2xy_ui;     // info NI(xy|ui)
  logC_xy_ui = logC3xy_ui - logC2xy_ui;     // cplx k_(xy|ui)
  info_yz_ui = info3yz_ui - info2yz_ui;     // info NI(yz|ui)
  logC_yz_ui = logC3yz_ui - logC2yz_ui;     // cplx k_(yz|ui)
  info_xz_ui = info3xz_ui - info2xz_ui;     // info NI(xz|ui)
  logC_xz_ui = logC3xz_ui - logC2xz_ui;     // cplx k_(xz|ui)
  info_xy_uiz = info3xy_uiz - info2xy_uiz;  // info NI(xy|uiz)
  logC_xy_uiz = logC3xy_uiz - logC2xy_uiz;  // cplx k_(xy|uiz)
  // compute score and store z with max score
  double xz = (info_xz_ui - info_xy_ui) - (logC_xz_ui - logC_xy_ui);
  double yz  = (info_yz_ui - info_xy_ui) - (logC_yz_ui - logC_xy_ui);
  double xyz = (info_xy_ui - info_xy_uiz) - (logC_xy_ui - logC_xy_uiz);

  double lower{0}, higher{0};
  std::tie(lower, higher) = std::minmax(xz, yz);
  double dpi = lower - log1p(exp(lower - higher));

  double Rzi = xyz < dpi ? xyz : dpi;

  double* ptrRetValues = new double[3];
  ptrRetValues[0] = Rzi;
  ptrRetValues[1] = info_xy_ui - info_xy_uiz;
  ptrRetValues[2] = logC_xy_ui - logC_xy_uiz;
  return ptrRetValues;
}

}  // namespace computation
}  // namespace miic
