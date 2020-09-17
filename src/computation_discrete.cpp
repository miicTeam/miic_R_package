#include "computation_discrete.h"

#include <algorithm>  // std::min
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>  // std::iota

#include "linear_allocator.h"
#include "mutual_information.h"
#include "structure.h"

constexpr int MDL = 0;
namespace miic {
namespace computation {

using namespace miic::structure;
using miic::utility::TempAllocatorScope;
using std::log;
using std::vector;

InfoBlock computeCondMutualInfoDiscrete(const TempGrid2d<int>& data,
    const TempVector<int>& r_list, const TempVector<int>& var_idx,
    const TempVector<double>& weights, int cplx,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_samples = data.n_cols();
  int id_x = var_idx[0], id_y = var_idx[1];
  int rx = r_list[id_x], ry = r_list[id_y];
  TempVector<int> ui_list(begin(var_idx) + 2, end(var_idx));

  TempVector<int> order = getDataOrder(data, r_list, var_idx);
  TempVector<int> hash_u(n_samples, 0);
  int rui = fillHashList(data, r_list, ui_list, hash_u);

  int Nyui{0}, Nui{0}, Ntot{0};
  double info_xui_y{0}, info_yui_x{0}, info_ui_y{0}, info_ui_x{0};
  double logC_xui_y{0}, logC_yui_x{0}, logC_ui_y{0}, logC_ui_x{0};
  // initialization of counts and mutual infos & logCs
  TempVector<int> Nxui(rx, 0);
  TempVector<int> Nx(rx, 0);
  TempVector<int> Ny(ry, 0);

  // Sentinels whose change of value (compared to the next sample) indicates
  // that the related counts should be added to mutual information (as NlogN)
  int X{data(id_x, order[0])};
  int X_next{-1};
  int Y{data(id_y, order[0])};
  int Y_next{-1};
  int Lui = hash_u[order[0]];
  int Lui_next{-1};

  double Pxyui = 0;
  // make the counts and compute mutual infos & logCs
  for (int k = 0; k < n_samples; ++k) {
    Pxyui += weights[order[k]];
    int i_next = k + 1 < n_samples ? order[k + 1] : -1;
    if (i_next != -1) {
      X_next = data(id_x, i_next);
      Y_next = data(id_y, i_next);
      Lui_next = hash_u[i_next];
      if (X_next == X && Y_next == Y && Lui_next == Lui) continue;
    }
    // Conclude on current count
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
    Pxyui = 0;  // reset cumulative weight
    if (i_next != -1) {
      X = X_next;
      if (Y_next == Y && Lui_next == Lui) continue;
    }
    // Conclude on current count
    if (Nyui > 0) {
      double NlogN = Nyui * cache->getLog(Nyui);
      info_yui_x -= NlogN;
      info_ui_y += NlogN;
      if (cplx != MDL) {
        logC_yui_x += cache->getLogC(Nyui, rx);
      }
      Nyui = 0;
    }
    if (i_next != -1) {
      Y = Y_next;
      if (Lui_next == Lui) continue;
      Lui = Lui_next;
    }
    // Conclude on current count
    for (auto& Nxuij : Nxui) {
      if (Nxuij > 0) {
        double NlogN = Nxuij * cache->getLog(Nxuij);
        info_xui_y -= NlogN;
        info_ui_x += NlogN;
        if (cplx != MDL) {
          logC_xui_y += cache->getLogC(Nxuij, ry);
        }
        Nxuij = 0;  // reset counter
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
    double logN = cache->getLog(Ntot);
    logC_xui_y = 0.5 * (ry - 1) * (rx * rui - 1) * logN;
    logC_yui_x = 0.5 * (rx - 1) * (ry * rui - 1) * logN;
    logC_ui_y = 0.5 * (ry - 1) * (rui - 1) * logN;
    logC_ui_x = 0.5 * (rx - 1) * (rui - 1) * logN;
  }

  double logC3xy_ui = 0.5 * (logC_xui_y + logC_yui_x);
  double logC2xy_ui = 0.5 * (logC_ui_y + logC_ui_x);

  return InfoBlock{Ntot, info3xy_ui - info2xy_ui, logC3xy_ui - logC2xy_ui};
}

Info3PointBlock computeInfo3PointAndScoreDiscrete(const TempGrid2d<int>& data,
    const TempVector<int>& r_list, const TempVector<int>& var_idx,
    const TempVector<double>& weights, int cplx,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_samples = data.n_cols();
  int n_nodes = data.n_rows() - 1;  // excluding Z
  int id_x = var_idx[0], id_y = var_idx[1], id_z = var_idx.back();
  int rx = r_list[id_x], ry = r_list[id_y], rz = r_list[id_z];
  TempVector<int> ui_list(begin(var_idx) + 2, begin(var_idx) + n_nodes);
  TempVector<int> xyui_list(begin(var_idx), begin(var_idx) + n_nodes);

  TempVector<int> order = getDataOrder(data, r_list, xyui_list);
  TempVector<int> hash_u(n_samples, 0);
  int rui = fillHashList(data, r_list, ui_list, hash_u);

  double info_xui_y{0}, info_yui_x{0}, info_ui_y{0}, info_ui_x{0};
  double logC_xui_y{0}, logC_yui_x{0}, logC_ui_y{0}, logC_ui_x{0};

  double info_xuiz_y{0}, info_yuiz_x{0}, info_uiz_y{0}, info_uiz_x{0};
  double logC_xuiz_y{0}, logC_yuiz_x{0}, logC_uiz_y{0}, logC_uiz_x{0};

  double info_xui_z{0}, info_yui_z{0}, info_ui_z{0};
  double logC_xui_z{0}, logC_yui_z{0}, logC_ui_z{0};

  int Nyui{0}, Nui{0}, Ntot{0};
  TempVector<int> Nxui(rx, 0);
  TempVector<int> Nx(rx, 0);
  TempVector<int> Ny(ry, 0);
  TempVector<int> Nyuiz(rz, 0);
  TempVector<int> Nuiz(rz, 0);
  TempVector<int> Nz(rz, 0);
  TempGrid2d<int> Nxuiz(rx, rz, 0);  // [X][Z]
  // Sentinels whose change of value (compared to the next sample) indicates
  // that the related counts should be added to mutual information (as NlogN)
  int X{data(id_x, order[0])};
  int X_next{-1};
  int Y{data(id_y, order[0])};
  int Y_next{-1};
  int Lui = hash_u[order[0]];
  int Lui_next{-1};

  TempVector<double> Pxyuiz(rz, 0);
  // make the counts and compute mutual infos & logCs
  for (int k = 0; k < n_samples; ++k) {
    int Z = data(id_z, order[k]);
    Pxyuiz[Z] += weights[order[k]];
    int i_next = k + 1 < n_samples ? order[k + 1] : -1;
    if (i_next != -1) {
      X_next = data(id_x, i_next);
      Y_next = data(id_y, i_next);
      Lui_next = hash_u[i_next];
      if (X_next == X && Y_next == Y && Lui_next == Lui) continue;
    }
    // Conclude on current count
    int Nxyui = 0;
    for (int l = 0; l < rz; l++) {
      double Pxyuizl = Pxyuiz[l];
      if (Pxyuizl == 0) continue;

      int Nxyuizl = (int)Pxyuizl;
      if (Nxyuizl > 0) {
        double NlogN = Nxyuizl * cache->getLog(Nxyuizl);
        info_xuiz_y += NlogN;
        info_yuiz_x += NlogN;

        Nxyui += Nxyuizl;
        Nz[l] += Nxyuizl;
        Nuiz[l] += Nxyuizl;
        Nyuiz[l] += Nxyuizl;
        Nxuiz(X, l) += Nxyuizl;
      }
      Pxyuiz[l] = 0;
    }

    if (Nxyui > 0) {
      double NlogN = Nxyui * cache->getLog(Nxyui);
      info_xui_y += NlogN;
      info_yui_x += NlogN;

      Ntot += Nxyui;
      Nui += Nxyui;
      Nyui += Nxyui;
      Nxui[X] += Nxyui;
      Nx[X] += Nxyui;
      Ny[Y] += Nxyui;
    }
    if (i_next != -1) {
      X = X_next;
      if (Y_next == Y && Lui_next == Lui) continue;
    }
    // Conclude on current count
    if (Nyui > 0) {
      double NlogN = Nyui * cache->getLog(Nyui);
      info_yui_x -= NlogN;
      info_yui_z -= NlogN;
      info_ui_y += NlogN;
      if (cplx != MDL) {
        logC_yui_x += cache->getLogC(Nyui, rx);
        logC_yui_z += cache->getLogC(Nyui, rz);
      }
      for (auto& Nyuizl : Nyuiz) {
        if (Nyuizl > 0) {
          double NlogN = Nyuizl * cache->getLog(Nyuizl);
          info_yui_z += NlogN;
          info_uiz_y += NlogN;
          info_yuiz_x -= NlogN;
          if (cplx != MDL) {
            logC_yuiz_x += cache->getLogC(Nyuizl, rx);
          }
          Nyuizl = 0;
        }
      }
      Nyui = 0;
    }
    if (i_next != -1) {
      Y = Y_next;
      if (Lui_next == Lui) continue;
      Lui = Lui_next;
    }
    if (Nui == 0) continue;

    double NlogN = Nui * cache->getLog(Nui);
    info_ui_x -= NlogN;
    info_ui_y -= NlogN;
    info_ui_z -= NlogN;
    if (cplx != MDL) {
      logC_ui_x += cache->getLogC(Nui, rx);
      logC_ui_y += cache->getLogC(Nui, ry);
      logC_ui_z += cache->getLogC(Nui, rz);
    }
    Nui = 0;

    for (auto& Nuizl : Nuiz) {
      if (Nuizl > 0) {
        double NlogN = Nuizl * cache->getLog(Nuizl);
        info_ui_z += NlogN;
        info_uiz_x -= NlogN;
        info_uiz_y -= NlogN;
        if (cplx != MDL) {
          logC_uiz_x += cache->getLogC(Nuizl, rx);
          logC_uiz_y += cache->getLogC(Nuizl, ry);
        }
        Nuizl = 0;
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
        Nxuiz(j, l) = 0;
      }
    }
  }
  // increment info with Nx[X], Ny[Y] and Nz[Z] contributions
  for (int j = 0; j < rx; j++) {
    int Nxj = Nx[j];
    if (Nxj > 0) {
      double NlogN = Nxj * log(Nxj / (1.0 * Ntot));
      info_yui_x -= NlogN;
      info_ui_x -= NlogN;
      info_uiz_x -= NlogN;
      info_yuiz_x -= NlogN;
      Nx[j] = 0;
    }
  }
  for (int j = 0; j < ry; j++) {
    int Nyj = Ny[j];
    if (Nyj > 0) {
      double NlogN = Nyj * log(Nyj / (1.0 * Ntot));
      info_xui_y -= NlogN;
      info_ui_y -= NlogN;
      info_uiz_y -= NlogN;
      info_xuiz_y -= NlogN;
      Ny[j] = 0;
    }
  }
  for (int l = 0; l < rz; l++) {
    int Nzl = Nz[l];
    if (Nzl > 0) {
      double NlogN = Nzl * log(Nzl / (1.0 * Ntot));
      info_xui_z -= NlogN;
      info_yui_z -= NlogN;
      info_ui_z -= NlogN;
      Nz[l] = 0;
    }
  }
  // NI(xy|ui)
  double info3xy_ui = 0.5 * (info_xui_y + info_yui_x);
  double info2xy_ui = 0.5 * (info_ui_y + info_ui_x);

  double logN{0};
  // check maximum mutual infos - cplx terms
  if (cplx == MDL) {
    logN = cache->getLog(Ntot);
    logC_xui_y = 0.5 * (ry - 1) * (rx * rui - 1) * logN;
    logC_yui_x = 0.5 * (rx - 1) * (ry * rui - 1) * logN;
    logC_ui_y = 0.5 * (ry - 1) * (rui - 1) * logN;
    logC_ui_x = 0.5 * (rx - 1) * (rui - 1) * logN;
    logC_uiz_y = 0.5 * (ry - 1) * (rz * rui - 1) * logN;
    logC_yui_z = 0.5 * (rz - 1) * (ry * rui - 1) * logN;
    logC_ui_z = 0.5 * (rz - 1) * (rui - 1) * logN;
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
    logC_uiz_x = 0.5 * (rx - 1) * (rz * rui - 1) * logN;
    logC_xui_z = 0.5 * (rz - 1) * (rx * rui - 1) * logN;
  }
  double logC3xz_ui = 0.5 * (logC_uiz_x + logC_xui_z);
  double logC2xz_ui = 0.5 * (logC_ui_x + logC_ui_z);
  // NI(xy|uiz)
  double info3xy_uiz = 0.5 * (info_xuiz_y + info_yuiz_x);
  double info2xy_uiz = 0.5 * (info_uiz_y + info_uiz_x);

  if (cplx == MDL) {
    rui *= rz;
    logC_xuiz_y = 0.5 * (ry - 1) * (rx * rui - 1) * logN;
    logC_yuiz_x = 0.5 * (rx - 1) * (ry * rui - 1) * logN;
    logC_uiz_y = 0.5 * (ry - 1) * (rui - 1) * logN;
    logC_uiz_x = 0.5 * (rx - 1) * (rui - 1) * logN;
  }
  double logC3xy_uiz = 0.5 * (logC_xuiz_y + logC_yuiz_x);
  double logC2xy_uiz = 0.5 * (logC_uiz_y + logC_uiz_x);
  double info_xy_ui = info3xy_ui - info2xy_ui;
  double logC_xy_ui = logC3xy_ui - logC2xy_ui;
  double info_yz_ui = info3yz_ui - info2yz_ui;
  double logC_yz_ui = logC3yz_ui - logC2yz_ui;
  double info_xz_ui = info3xz_ui - info2xz_ui;
  double logC_xz_ui = logC3xz_ui - logC2xz_ui;
  double info_xy_uiz = info3xy_uiz - info2xy_uiz;
  double logC_xy_uiz = logC3xy_uiz - logC2xy_uiz;

  double xz = (info_xz_ui - info_xy_ui) - (logC_xz_ui - logC_xy_ui);
  double yz  = (info_yz_ui - info_xy_ui) - (logC_yz_ui - logC_xy_ui);
  // Data processing inequality
  double dpi = std::fmin(xz, yz) - log1p(exp(-std::fabs(xz - yz)));

  double Ixyz_ui = info_xy_ui - info_xy_uiz;
  double kxyz_ui = logC_xy_ui - logC_xy_uiz;
  double Ikxyz_ui = Ixyz_ui - kxyz_ui;

  double Rscore = std::fmin(Ikxyz_ui, dpi);

  return Info3PointBlock{Rscore, Ixyz_ui, kxyz_ui};
}

}  // namespace computation
}  // namespace miic
