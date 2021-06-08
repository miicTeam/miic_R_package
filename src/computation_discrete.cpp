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

using std::log;
using std::lround;
using namespace miic::structure;
using miic::utility::TempAllocatorScope;

InfoBlock computeCondMutualInfoDiscrete(const TempGrid2d<int>& data,
    const TempVector<int>& r_list, const TempVector<int>& var_idx,
    const TempVector<double>& weights, int cplx,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_samples = data.n_cols();
  int id_x = var_idx[0], id_y = var_idx[1];
  int rx = r_list[id_x], ry = r_list[id_y];
  double n_eff = accumulate(begin(weights), end(weights), 0.0);

  if (var_idx.size() == 2) {
    TempVector<int> xy_factors(n_samples);
    int rxy = setJointFactors(data, r_list, var_idx, xy_factors);

    TempVector<int> r_temp{rx, ry, rxy};
    return computeMI(data.getConstRow(id_x), data.getConstRow(id_y), xy_factors,
        r_temp, n_eff, weights, cache, cplx, 0);
  }

  TempVector<int> ui_list(begin(var_idx) + 2, end(var_idx));
  TempVector<int> order = getDataOrder(data, r_list, var_idx);
  TempVector<int> hash_u(n_samples, 0);
  int ru = fillHashList(data, r_list, ui_list, hash_u);
  // Entropy terms
  double Hu{0}, Huy{0}, Hux{0}, Huyx{0};
  // Complexity terms
  double logC_ux_y{0}, logC_uy_x{0}, logC_u_y{0}, logC_u_x{0};
  // Counting variables
  double Nu{0}, Nuy{0}, Nuyx{0}, N_total{0};
  TempVector<double> Nux_list(rx, 0);
  // Sentinels whose change of value (compared to the next sample) indicates
  // that the related counts should be added to mutual information (as NlogN)
  int X{data(id_x, order[0])};
  int X_next{-1};
  int Y{data(id_y, order[0])};
  int Y_next{-1};
  int Lui = hash_u[order[0]];
  int Lui_next{-1};

  // make the counts and compute mutual infos & logCs
  for (int k = 0; k < n_samples; ++k) {
    Nuyx += weights[order[k]];
    int i_next = k + 1 < n_samples ? order[k + 1] : -1;
    if (i_next != -1) {
      X_next = data(id_x, i_next);
      Y_next = data(id_y, i_next);
      Lui_next = hash_u[i_next];
      if (X_next == X && Y_next == Y && Lui_next == Lui) continue;
    }
    // Conclude on current count
    if (Nuyx > 0) {
      Nu += Nuyx;
      Nuy += Nuyx;
      Nux_list[X] += Nuyx;
      N_total += Nuyx;

      Huyx -= Nuyx * log(Nuyx);
    }
    Nuyx = 0;  // reset cumulative weight
    if (i_next != -1) {
      X = X_next;
      if (Y_next == Y && Lui_next == Lui) continue;
    }
    // Conclude on current count
    if (Nuy > 0) {
      Huy -= Nuy * log(Nuy);
      if (cplx != MDL) {
        logC_uy_x += cache->getLogC(lround(Nuy), rx);
      }
      Nuy = 0;
    }
    if (i_next != -1) {
      Y = Y_next;
      if (Lui_next == Lui) continue;
      Lui = Lui_next;
    }
    // Conclude on current count
    for (auto& Nxu : Nux_list) {
      if (Nxu > 0) {
        Hux -= Nxu * log(Nxu);
        if (cplx != MDL) {
          logC_ux_y += cache->getLogC(lround(Nxu), ry);
        }
        Nxu = 0;  // reset counter
      }
    }
    if (Nu > 0) {
      Hu -= Nu * log(Nu);
      if (cplx != MDL) {
        auto Nu_long = lround(Nu);
        logC_u_x += cache->getLogC(Nu_long, rx);
        logC_u_y += cache->getLogC(Nu_long, ry);
      }
      Nu = 0;
    }
  }

  if (cplx == MDL) {
    double logN = log(N_total);
    logC_ux_y = 0.5 * (ry - 1) * (rx * ru - 1) * logN;
    logC_uy_x = 0.5 * (rx - 1) * (ry * ru - 1) * logN;
    logC_u_y = 0.5 * (ry - 1) * (ru - 1) * logN;
    logC_u_x = 0.5 * (rx - 1) * (ru - 1) * logN;
  }

  double Ixy_ui = Hux + Huy - Hu - Huyx;
  double kxy_ui = 0.5 * (logC_ux_y - logC_u_y + logC_uy_x - logC_u_x);

  return InfoBlock{N_total, Ixy_ui, kxy_ui};
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
  int ru = fillHashList(data, r_list, ui_list, hash_u);
  // Entropy terms
  double Hu{0}, Huy{0}, Hux{0}, Huyx{0};
  double Hzu{0}, Hzuy{0}, Hzux{0}, Hzuyx{0};
  // Complexity terms
  double logC_ux_y{0}, logC_uy_x{0}, logC_u_y{0}, logC_u_x{0};
  double logC_zux_y{0}, logC_zuy_x{0}, logC_zu_y{0}, logC_zu_x{0};
  double logC_ux_z{0}, logC_uy_z{0}, logC_u_z{0};
  // Counting variables
  double Nuy{0}, Nu{0}, N_total{0};
  TempVector<double> Nux_list(rx, 0);
  TempVector<double> Nzu_list(rz, 0);
  TempVector<double> Nzuy_list(rz, 0);
  TempGrid2d<double> Nzux_list(rx, rz, 0);  // [X][Z]
  TempVector<double> Nzuyx_list(rz, 0);
  // Sentinels whose change of value (compared to the next sample) indicates
  // that the related counts should be added to mutual information (as NlogN)
  int X{data(id_x, order[0])};
  int X_next{-1};
  int Y{data(id_y, order[0])};
  int Y_next{-1};
  int Lui = hash_u[order[0]];
  int Lui_next{-1};
  // make the counts and compute entropy and logCs
  for (int k = 0; k < n_samples; ++k) {
    int Z = data(id_z, order[k]);
    Nzuyx_list[Z] += weights[order[k]];
    int i_next = k + 1 < n_samples ? order[k + 1] : -1;
    if (i_next != -1) {
      X_next = data(id_x, i_next);
      Y_next = data(id_y, i_next);
      Lui_next = hash_u[i_next];
      if (X_next == X && Y_next == Y && Lui_next == Lui) continue;
    }
    // Conclude on current count
    double Nuyx = 0;
    for (int l = 0; l < rz; l++) {
      double Nzuyx = Nzuyx_list[l];
      if (Nzuyx == 0) continue;

      Hzuyx -= Nzuyx * log(Nzuyx);

      Nuyx += Nzuyx;
      Nzu_list[l] += Nzuyx;
      Nzuy_list[l] += Nzuyx;
      Nzux_list(X, l) += Nzuyx;

      Nzuyx_list[l] = 0;
    }

    if (Nuyx > 0) {
      Huyx -= Nuyx * log(Nuyx);

      N_total += Nuyx;
      Nu += Nuyx;
      Nuy += Nuyx;
      Nux_list[X] += Nuyx;
    }
    if (i_next != -1) {
      X = X_next;
      if (Y_next == Y && Lui_next == Lui) continue;
    }
    // Conclude on current count
    if (Nuy > 0) {
      Huy -= Nuy * log(Nuy);
      if (cplx != MDL) {
        auto Nuy_long = lround(Nuy);
        logC_uy_x += cache->getLogC(Nuy_long, rx);
        logC_uy_z += cache->getLogC(Nuy_long, rz);
      }
      for (auto& Nzuy : Nzuy_list) {
        if (Nzuy > 0) {
          Hzuy -= Nzuy * log(Nzuy);
          if (cplx != MDL) {
            logC_zuy_x += cache->getLogC(lround(Nzuy), rx);
          }
          Nzuy = 0;
        }
      }
      Nuy = 0;
    }
    if (i_next != -1) {
      Y = Y_next;
      if (Lui_next == Lui) continue;
      Lui = Lui_next;
    }
    if (Nu == 0) continue;

    Hu -= Nu * log(Nu);
    if (cplx != MDL) {
      auto Nu_long = lround(Nu);
      logC_u_x += cache->getLogC(Nu_long, rx);
      logC_u_y += cache->getLogC(Nu_long, ry);
      logC_u_z += cache->getLogC(Nu_long, rz);
    }
    Nu = 0;

    for (auto& Nzu : Nzu_list) {
      if (Nzu > 0) {
        Hzu -= Nzu * log(Nzu);
        if (cplx != MDL) {
          auto Nzu_long = lround(Nzu);
          logC_zu_x += cache->getLogC(Nzu_long, rx);
          logC_zu_y += cache->getLogC(Nzu_long, ry);
        }
        Nzu = 0;
      }
    }

    for (int j = 0; j < rx; j++) {
      double Nux = Nux_list[j];
      if (Nux == 0) continue;

      Hux -= Nux * log(Nux);
      if (cplx != MDL) {
        auto Nux_long = lround(Nux);
        logC_ux_y += cache->getLogC(Nux_long, ry);
        logC_ux_z += cache->getLogC(Nux_long, rz);
      }
      Nux_list[j] = 0;

      for (int l = 0; l < rz; l++) {
        double Nzux = Nzux_list(j, l);
        if (Nzux == 0) continue;

        Hzux -= Nzux * log(Nzux);
        if (cplx != MDL) {
          logC_zux_y += cache->getLogC(lround(Nzux), ry);
        }
        Nzux_list(j, l) = 0;
      }
    }
  }

  // check maximum mutual infos - cplx terms
  if (cplx == MDL) {
    double logN = log(N_total);
    logC_ux_y = 0.5 * (ry - 1) * (rx * ru - 1) * logN;
    logC_uy_x = 0.5 * (rx - 1) * (ry * ru - 1) * logN;
    logC_u_y = 0.5 * (ry - 1) * (ru - 1) * logN;
    logC_u_x = 0.5 * (rx - 1) * (ru - 1) * logN;
    logC_zu_y = 0.5 * (ry - 1) * (rz * ru - 1) * logN;
    logC_uy_z = 0.5 * (rz - 1) * (ry * ru - 1) * logN;
    logC_u_z = 0.5 * (rz - 1) * (ru - 1) * logN;

    logC_zu_x = 0.5 * (rx - 1) * (rz * ru - 1) * logN;
    logC_ux_z = 0.5 * (rz - 1) * (rx * ru - 1) * logN;

    int rzu = ru * rz;
    logC_zux_y = 0.5 * (ry - 1) * (rx * rzu - 1) * logN;
    logC_zuy_x = 0.5 * (rx - 1) * (ry * rzu - 1) * logN;
    logC_zu_y = 0.5 * (ry - 1) * (rzu - 1) * logN;
    logC_zu_x = 0.5 * (rx - 1) * (rzu - 1) * logN;
  }

  double info_xy_ui = Hux + Huy - Hu - Huyx;
  double logC_xy_ui = 0.5 * (logC_ux_y - logC_u_y + logC_uy_x - logC_u_x);
  double info_yz_ui = Huy + Hzu - Hu - Hzuy;
  double logC_yz_ui = 0.5 * (logC_zu_y - logC_u_y + logC_uy_z - logC_u_z);
  double info_xz_ui = Hux + Hzu - Hu - Hzux;
  double logC_xz_ui = 0.5 * (logC_zu_x - logC_u_x + logC_ux_z - logC_u_z);
  double info_xy_uiz = Hzux + Hzuy - Hzu - Hzuyx;
  double logC_xy_uiz = 0.5 * (logC_zux_y - logC_zu_y + logC_zuy_x - logC_zu_x);

  double xz = (info_xz_ui - logC_xz_ui) - (info_xy_ui - logC_xy_ui);
  double yz  = (info_yz_ui - logC_yz_ui) - (info_xy_ui - logC_xy_ui);
  // Data processing inequality
  double dpi = std::fmin(xz, yz) - std::log1p(exp(-std::fabs(xz - yz)));

  double Ixyz_ui = info_xy_ui - info_xy_uiz;
  double kxyz_ui = logC_xy_ui - logC_xy_uiz;
  double Ikxyz_ui = Ixyz_ui - kxyz_ui;

  double Rscore = std::fmin(Ikxyz_ui, dpi);

  return Info3PointBlock{Rscore, Ixyz_ui, kxyz_ui};
}

}  // namespace computation
}  // namespace miic
