#include "computation_discrete.h"

#include <algorithm>  // std::sort, std::minmax
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>  // std::iota
#include <tuple>    // std::tie

#include "linear_allocator.h"
#include "structure.h"

constexpr int MDL = 0;

// cplx = 0 --> MDL, cplx = 1 --> NML
namespace miic {
namespace computation {

using namespace miic::structure;
using miic::utility::TempAllocatorScope;
using std::vector;

TempGrid2d<int> getHashTable(const TempGrid2d<int>& data,
    const TempVector<int>& r_list, const TempVector<int>& var_idx, int n_vars) {
  int n_samples = data.n_cols();
  int id_x = var_idx[0], id_y = var_idx[1];
  int rx = r_list[id_x];

  TempGrid2d<int> hash_table(n_samples, 3);
  // Compute unique hash values for each sample in each of the joint spaces
  for (int k = 0; k < n_samples; k++) {
    int r_joint = rx;
    hash_table(k, 0) = data(id_y, k) * rx + data(id_x, k);  // Lxyui
    hash_table(k, 1) = data(id_y, k) * rx;                  // Lyui
    hash_table(k, 2) = 0;                                   // Lui

    for (int j = 2; j < n_vars; j++) {
      r_joint *= r_list[var_idx[j - 1]];
      int increment = data(var_idx[j], k) * r_joint;
      hash_table(k, 0) += increment;
      hash_table(k, 1) += increment;
      hash_table(k, 2) += increment;
    }
  }
  return hash_table;
}

InfoBlock computeCondMutualInfoDiscrete(const TempGrid2d<int>& data,
    const TempVector<int>& r_list, const TempVector<int>& var_idx,
    const TempVector<double>& weights, int cplx,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_samples = data.n_cols();
  int n_nodes = data.n_rows();
  int id_x = var_idx[0], id_y = var_idx[1];
  int rx = r_list[id_x], ry = r_list[id_y];

  // n_rows = n_samples, n_cols = 3 [0: Lxyui, 1: Lyui, 2: Lui]
  TempGrid2d<int> hash_table = getHashTable(data, r_list, var_idx, n_nodes);
  TempVector<int> order(n_samples);
  // [0 to n_samples]
  std::iota(begin(order), end(order), 0);
  std::sort(begin(order), end(order), [&hash_table](int a, int b) {
    return hash_table(a, 0) < hash_table(b, 0);
  });

  int Nyui{0}, Nui{0}, Ntot{0};
  double info_xui_y{0}, info_yui_x{0}, info_ui_y{0}, info_ui_x{0};
  double logC_xui_y{0}, logC_yui_x{0}, logC_ui_y{0}, logC_ui_x{0};
  // initialization of counts and mutual infos & logCs
  TempVector<int> Nxui(rx, 0);
  TempVector<int> Nx(rx, 0);
  TempVector<int> Ny(ry, 0);
  // Sentinels whose change of value (compared to the next sample) indicates
  // that the related counts should be added to mutual information (as NlogN)
  int Lxyui = hash_table(order[0], 0);
  int Lyui = hash_table(order[0], 1);
  int Lui = hash_table(order[0], 2);
  int X = data(id_x, order[0]);
  int Y = data(id_y, order[0]);

  double Pxyui = 0;
  // make the counts and compute mutual infos & logCs
  for (int k = 0; k < n_samples; ++k) {
    Pxyui += weights[order[k]];
    int i_next = k + 1 < n_samples ? order[k + 1] : -1;
    if (i_next != -1) {
      int next = hash_table(i_next, 0);
      if (next == Lxyui) continue;
      Lxyui = next;
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
      X = data(id_x, i_next);
      int next = hash_table(i_next, 1);
      if (next == Lyui) continue;
      Lyui = next;
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
      Y = data(id_y, i_next);
      int next = hash_table(i_next, 2);
      if (next == Lui) continue;
      Lui = next;
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
    int joint_rui = 1;
    double logN = cache->getLog(Ntot);
    for (int j = 2; j < n_nodes; j++)
      joint_rui *= r_list[var_idx[j]];
    logC_xui_y = 0.5 * (ry - 1) * (rx * joint_rui - 1) * logN;
    logC_yui_x = 0.5 * (rx - 1) * (ry * joint_rui - 1) * logN;
    logC_ui_y = 0.5 * (ry - 1) * (joint_rui - 1) * logN;
    logC_ui_x = 0.5 * (rx - 1) * (joint_rui - 1) * logN;
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

  // n_rows = n_samples, n_cols = 3 [0: Lxyui, 1: Lyui, 2: Lui]
  TempGrid2d<int> hash_table = getHashTable(data, r_list, var_idx, n_nodes);
  TempVector<int> order(n_samples);
  // [0 to n_samples_non_na]
  std::iota(begin(order), end(order), 0);
  std::sort(begin(order), end(order), [&hash_table](int a, int b) {
    return hash_table(a, 0) < hash_table(b, 0);
  });

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
  int Lxyui = hash_table(order[0], 0);
  int Lyui = hash_table(order[0], 1);
  int Lui = hash_table(order[0], 2);
  int X = data(id_x, order[0]);
  int Y = data(id_y, order[0]);

  TempVector<double> Pxyuiz(rz, 0);
  // make the counts and compute mutual infos & logCs
  for (int k = 0; k < n_samples; ++k) {
    int Z = data(id_z, order[k]);
    Pxyuiz[Z] += weights[order[k]];
    int i_next = k + 1 < n_samples ? order[k + 1] : -1;
    if (i_next != -1) {
      int next = hash_table(i_next, 0);
      if (next == Lxyui) continue;
      Lxyui = next;
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
      X = data(id_x, i_next);
      int next = hash_table(i_next, 1);
      if (next == Lyui) continue;
      Lyui = next;
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
      Y = data(id_y, i_next);
      int next = hash_table(i_next, 2);
      if (next == Lui) continue;
      Lui = next;
    }
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

  int joint_rui{1};
  double logN{0};
  // check maximum mutual infos - cplx terms
  if (cplx == MDL) {
    for (int j = 2; j < n_nodes; j++)
      joint_rui *= r_list[var_idx[j]];
    logN = cache->getLog(Ntot);
    logC_xui_y = 0.5 * (ry - 1) * (rx * joint_rui - 1) * logN;
    logC_yui_x = 0.5 * (rx - 1) * (ry * joint_rui - 1) * logN;
    logC_ui_y = 0.5 * (ry - 1) * (joint_rui - 1) * logN;
    logC_ui_x = 0.5 * (rx - 1) * (joint_rui - 1) * logN;
    logC_uiz_y = 0.5 * (ry - 1) * (rz * joint_rui - 1) * logN;
    logC_yui_z = 0.5 * (rz - 1) * (ry * joint_rui - 1) * logN;
    logC_ui_z = 0.5 * (rz - 1) * (joint_rui - 1) * logN;
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
    logC_uiz_x = 0.5 * (rx - 1) * (rz * joint_rui - 1) * logN;
    logC_xui_z = 0.5 * (rz - 1) * (rx * joint_rui - 1) * logN;
  }
  double logC3xz_ui = 0.5 * (logC_uiz_x + logC_xui_z);
  double logC2xz_ui = 0.5 * (logC_ui_x + logC_ui_z);
  // NI(xy|uiz)
  double info3xy_uiz = 0.5 * (info_xuiz_y + info_yuiz_x);
  double info2xy_uiz = 0.5 * (info_uiz_y + info_uiz_x);

  if (cplx == MDL) {
    joint_rui *= rz;
    logC_xuiz_y = 0.5 * (ry - 1) * (rx * joint_rui - 1) * logN;
    logC_yuiz_x = 0.5 * (rx - 1) * (ry * joint_rui - 1) * logN;
    logC_uiz_y = 0.5 * (ry - 1) * (joint_rui - 1) * logN;
    logC_uiz_x = 0.5 * (rx - 1) * (joint_rui - 1) * logN;
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
  // compute score
  double xz = (info_xz_ui - info_xy_ui) - (logC_xz_ui - logC_xy_ui);
  double yz  = (info_yz_ui - info_xy_ui) - (logC_yz_ui - logC_xy_ui);
  double xyz = (info_xy_ui - info_xy_uiz) - (logC_xy_ui - logC_xy_uiz);

  double lower{0}, higher{0};
  std::tie(lower, higher) = std::minmax(xz, yz);
  double dpi = lower - log1p(exp(lower - higher));

  double Rscore = xyz < dpi ? xyz : dpi;
  double Ixyz_ui = info_xy_ui - info_xy_uiz;
  double kxyz_ui = logC_xy_ui - logC_xy_uiz;

  return Info3PointBlock{Rscore, Ixyz_ui, kxyz_ui};
}

}  // namespace computation
}  // namespace miic
