#include "compute_info.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>  // std::numeric_limits
#include <unordered_set>

#include "structure.h"
#include "utilities.h"

constexpr int MDL = 0;
constexpr int randomrescaling = 1;

// cplx = 0 --> MDL, cplx = 1 --> NML
// If nbrZi== 0, return nSample[0]     & nSample[0]*I(xy|{ui})      & k_xy_ui
// If nbrZi > 0, return nSample[z_top] & nSample[z_top]*I(xy|{ui})  & k_xy_ui
//                      z_top          & nSample[z_top]*I(xy|{ui}z) & k_xy_uiz
//                      R_top          & nSample[z_top]*I(xyz|{ui}) & k_xyz_ui

namespace miic {
namespace computation {

using namespace miic::structure;
using namespace miic::utility;
using std::vector;

double* getAllInfoNEW(const vector<vector<int>>& data,
    const vector<int>& all_levels, int id_x, int id_y,
    const vector<int>& ui_list, const TempVector<int>& zi_list,
    int n_samples_eff, int cplx, int k23, const vector<double>& weights,
    const TempGrid2d<double>& freqs1, bool test_mar,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_ui = ui_list.size();
  int n_zi = zi_list.size();
  int n_samples = data.size();

  // Output pointer
  int nbrRetValues = 3;
  // If no zi, return nSample0 & nSample0*I(xy|{ui}) & k_xy_ui
  // If nbrZi > 0, return
  //   nSample[z_top] & nSample[z_top]*I(xy|{ui})  & k_xy_ui
  //   z_top          & nSample[z_top]*I(xy|{ui}z) & k_xy_uiz
  //   R_top          & nSample[z_top]*I(xyz|{ui}) & k_xyz_ui
  if (n_zi > 0) nbrRetValues = 9;

  double* ptrRetValues = new double[nbrRetValues];

  ptrRetValues[0] = -1;
  ptrRetValues[1] = -1;
  ptrRetValues[2] = -1;
  if (n_zi > 0) {
    ptrRetValues[3] = -1;
    ptrRetValues[4] = -1;
    ptrRetValues[5] = -1;
    ptrRetValues[6] = -1;
    ptrRetValues[7] = -1;
    ptrRetValues[8] = -1;
  }

  // Define the total number of variables (xi, xj, {ui}, {zi})
  int n_nodes_xyui = n_ui + 2;
  // [N+1][7] [sample id (in the list with na), Lxyui, Lyui, Lui, X, Y, Z]
  TempGrid2d<int> sample(n_samples + 1, 7);
  // find samples without NA in x,y,ui and store their id in sample(k, 0)
  TempVector<int> var_idx(n_nodes_xyui);
  var_idx[0] = id_x;
  var_idx[1] = id_y;
  std::copy(begin(ui_list), end(ui_list), begin(var_idx) + 2);
  int n_samples_non_na = 0;
  for (int i = 0; i < n_samples; i++) {
    bool ok = true;
    for (const auto node : var_idx) {
      // check that X, Y, Uj do not contain NA
      if (data[i][node] == -1) {
        ok = false;
        break;
      }
    }
    if (ok) {
      sample(n_samples_non_na++, 0) = i;  // sample number
    }
  }
  int rx_reduced = all_levels[id_x];
  int ry_reduced = all_levels[id_y];
  // update dbin
  if (n_samples_non_na < n_samples) {
    std::unordered_set<int> sx, sy;
    for (int k = 0; k < n_samples_non_na; k++) {
      int i = sample(k, 0);
      sx.insert(data[i][id_x]);
      sy.insert(data[i][id_y]);
    }
    rx_reduced = sx.size();
    ry_reduced = sy.size();
  }

  if (n_samples_non_na == 0) return ptrRetValues;

  TempVector<int> r_list(n_nodes_xyui + 1);  // plus z
  // initialisation of n_levels for continuous variables in x,y, ui
  for (int i = 0; i < n_nodes_xyui; i++) {
    r_list[i] = all_levels[var_idx[i]];
  }
  int PBin;
  // compute Lxyui, Lyui, Lui indices for later counting purpose
  for (int k = 0; k < n_samples_non_na; k++) {
    int i_sample = sample(k, 0);  // sample number

    PBin = r_list[0];        // rx
    sample(k, 4) = data[i_sample][id_x];  // X
    sample(k, 5) = data[i_sample][id_y];  // Y
    sample(k, 1) = sample(k, 5) * r_list[0] + sample(k, 4); // Lxyui
    sample(k, 2) = sample(k, 5) * r_list[0];   // Lyui initialisation
    sample(k, 3) = 0;           // Lui initialisation

    for (int j = 2; j < n_nodes_xyui; j++) {
      PBin *= r_list[j - 1];
      int increment = data[i_sample][var_idx[j]] * PBin;
      sample(k, 1) += increment;  // Lxyui
      sample(k, 2) += increment;  // Lyui
      sample(k, 3) += increment;  // Lui
    }
  }
  int n_joint_levels = PBin * r_list[n_nodes_xyui - 1];

  // extra sample id (useful for termination of counts)
  sample(n_samples_non_na, 0) = n_samples_non_na;
  // max Lxyui (useful for termination of counts)
  sample(n_samples_non_na, 1) = n_joint_levels;
  // max Lyui  (useful for termination of counts)
  sample(n_samples_non_na, 2) = n_joint_levels;
  // max Lui   (useful for termination of counts)
  sample(n_samples_non_na, 3) = n_joint_levels;

  TempVector<int> orderSample(n_samples_non_na + 1);
  // [0 to n_samples_non_na]
  std::iota(begin(orderSample), end(orderSample), 0);
  std::sort(begin(orderSample), end(orderSample),
      [&sample](int a, int b) { return sample(a, 1) < sample(b, 1); });

  TempGrid2d<int> sortedSample(sample.n_rows(), sample.n_cols());
  for (int k = 0; k < n_samples_non_na + 1; k++) {
    int i = orderSample[k];
    for (int j = 0; j < 6; j++) {
      sortedSample(k, j) = sample(i, j);
    }
  }

  int Nyui = 0;
  int Nui = 0;

  int X = sortedSample(0, 4);
  int Y = sortedSample(0, 5);

  double info_xui_y{0}, info_yui_x{0}, info_ui_y {0}, info_ui_x {0};
  double logC_xui_y{0}, logC_yui_x{0}, logC_ui_y {0}, logC_ui_x {0};

  // 6 test variables for counts, should all sum to nSample0
  int Nxyuis{0}, Nyuis{0}, Nxuis{0}, Nuis{0}, Nxs{0}, Nys{0}, Ntot{0};

  TempVector<int> Nxui(r_list[0], 0);
  TempVector<int> Nx(r_list[0], 0);
  TempVector<int> Ny(r_list[1], 0);
  // initialization of counts and mutual infos & logCs
  int Lxyui = sortedSample(0, 1);  // min Lxyui
  int Lyui = sortedSample(0, 2);   // min Lyui
  int Lui = sortedSample(0, 3);    // min Lui
  double Pxyui = weights[sortedSample(0, 0)];
  // make the counts and compute mutual infos & logCs
  for (int k = 1; k <= n_samples_non_na; k++) {
    if (sortedSample(k, 1) <= Lxyui) {
      Pxyui += weights[sortedSample(k, 0)];  // weights[k];
      continue;
    }
    Lxyui = sortedSample(k, 1);
    int Nxyui = 0;
    if (n_samples_eff != n_samples) {
      if (randomrescaling) {
        Nxyui = (int)floor(Pxyui);
        double r = Pxyui - Nxyui;
        double rr = R::runif(0, 1);
        if (r > rr) ++Nxyui;
      }
    } else {
      Nxyui = (int)Pxyui;
    }

    if (Nxyui > 0) {
      Nui += Nxyui;
      Nyui += Nxyui;
      Nxui[X] += Nxyui;
      Nx[X] += Nxyui;
      Ny[Y] += Nxyui;
      Ntot += Nxyui;
      Nxyuis += Nxyui;

      double NlogN = Nxyui * cache->getLog(Nxyui);
      info_xui_y += NlogN;
      info_yui_x += NlogN;
    }
    if (sortedSample(k, 0) < n_samples)
      Pxyui = weights[sortedSample(k, 0)];  // weights[k];

    if (k < n_samples_non_na) X = sortedSample(k, 4);
    if (k < n_samples_non_na) Y = sortedSample(k, 5);

    if (sortedSample(k, 2) <= Lyui) continue;
    Lyui = sortedSample(k, 2);

    if (Nyui > 0) {
      double NlogN = Nyui * cache->getLog(Nyui);
      info_yui_x -= NlogN;
      info_ui_y += NlogN;
      if (cplx != MDL) {
        logC_yui_x += cache->getLogC(Nyui, r_list[0]);
      }
      Nyuis += Nyui;
      Nyui = 0;
    }

    if (sortedSample(k, 3) <= Lui) continue;
    Lui = sortedSample(k, 3);

    for (int j = 0; j < r_list[0]; j++) {
      int Nxuij = Nxui[j];
      if (Nxuij > 0) {
        double NlogN = Nxuij * cache->getLog(Nxuij);
        info_xui_y -= NlogN;
        info_ui_x += NlogN;
        if (cplx != MDL) {
          logC_xui_y += cache->getLogC(Nxuij, r_list[1]);
        }

        Nxuis += Nxuij;
        Nxui[j] = 0;
      }
    }

    if (Nui > 0) {
      double NlogN = Nui * cache->getLog(Nui);
      info_ui_y -= NlogN;
      info_ui_x -= NlogN;
      if (cplx != MDL) {
        logC_ui_x += cache->getLogC(Nui, r_list[0]);
        logC_ui_y += cache->getLogC(Nui, r_list[1]);
      }
      Nuis += Nui;
      Nui = 0;
    }
  }
  // increment for info for Nx[X] and Ny[Y] contributions
  for (int j = 0; j < r_list[0]; j++) {
    int Nxj = Nx[j];
    if (Nxj > 0) {
      double NlogN = Nxj * log(Nxj / (1.0 * Ntot));
      info_yui_x -= NlogN;
      info_ui_x -= NlogN;
      Nxs += Nxj;
    }
  }
  for (int j = 0; j < r_list[1]; j++) {
    int Nyj = Ny[j];
    if (Nyj > 0) {
      double NlogN = Nyj * log(Nyj / (1.0 * Ntot));
      info_xui_y -= NlogN;
      info_ui_y -= NlogN;
      Nys += Nyj;
    }
  }

  // check maximum mutual infos - cplx terms
  double info3xy_ui = 0.5 * (info_xui_y + info_yui_x);
  double info2xy_ui = 0.5 * (info_ui_y + info_ui_x);

  if (cplx == MDL) {
    int Prui = 1;
    double logN = cache->getLog(Ntot);
    for (int j = 2; j < n_nodes_xyui; j++)
      Prui *= r_list[j];
    logC_xui_y = 0.5 * (r_list[1] - 1) * (r_list[0] * Prui - 1) * logN;
    logC_yui_x = 0.5 * (r_list[0] - 1) * (r_list[1] * Prui - 1) * logN;
    logC_ui_y = 0.5 * (r_list[1] - 1) * (Prui - 1) * logN;
    logC_ui_x = 0.5 * (r_list[0] - 1) * (Prui - 1) * logN;
  }

  double logC3xy_ui = 0.5 * (logC_xui_y + logC_yui_x);
  double logC2xy_ui = 0.5 * (logC_ui_y + logC_ui_x);

  // forbid negative 2-point Ik // change 20200206
  // if(info3xy_ui - logC3xy_ui - info2xy_ui + logC2xy_ui<0) {
  //  info3xy_ui=0;
  //  logC3xy_ui=0;
  //  info2xy_ui=0;
  //  logC2xy_ui=0;
  //}

  int N_xy_ui{0};
  double NIxy_ui = -1.0;
  double k_xy_ui = -1.0;
  N_xy_ui = Ntot;
  NIxy_ui = info3xy_ui - info2xy_ui;  // info to be returned if no z
  k_xy_ui = logC3xy_ui - logC2xy_ui;  // cplx to be returned if no z

  if (n_zi == 0) {
    ptrRetValues[0] = N_xy_ui;
    ptrRetValues[1] = NIxy_ui;
    ptrRetValues[2] = k_xy_ui;

    if (rx_reduced == 1 || ry_reduced == 1 || rx_reduced == N_xy_ui ||
        ry_reduced == N_xy_ui) {
      ptrRetValues[1] = 0;
    }
    return ptrRetValues;
  }
  // n_zi > 0
  // pick next z and compute score
  int N_xyuiz_top = 0;
  double NIxy_ui_top = -1;
  double k_xy_ui_top = -1;

  int z_top = -1;
  double NIxy_uiz_top = -1;
  double k_xy_uiz_top = -1;

  double R_top = std::numeric_limits<double>::lowest();
  double NIxyz_ui_top = -1;
  double k_xyz_ui_top = -1;

  int z = n_nodes_xyui;
  for (int zi = 0; zi < n_zi; zi++) {
    TempAllocatorScope scope;

    int idxzi = zi_list[zi];

    int kz0{0};
    bool isGoodCandidate = false;
    for (int k = 0; k < n_samples_non_na; k++) {
      int i = sortedSample(k, 0);
      // find the first sample for which zi does not contain NA
      if (data[i][idxzi] > -1) {
        kz0 = k;
        isGoodCandidate = true;
        break;
      }
    }
    if (!isGoodCandidate) continue;

    TempVector<int> sample_non_na_with_z;
    sample_non_na_with_z.reserve(n_samples);
    for (int i = 0; i < n_samples; i++) {
      bool ok = true;
      // check if X, Y, Ui contain NA(-1)
      for (const auto node : var_idx) {
        if (data[i][node] == -1) {
          ok = false;
          break;
        }
      }
      if (data[i][idxzi] == -1) ok = false;

      if (ok) sample_non_na_with_z.push_back(i);
    }
    sample_non_na_with_z.shrink_to_fit();
    int n_samples_with_z = sample_non_na_with_z.size();

    isGoodCandidate = true;
    int rx_reduced = all_levels[id_x];
    int ry_reduced = all_levels[id_y];
    if (n_samples_with_z < n_samples) {
      std::unordered_set<int> sx, sy;
      for (const auto i : sample_non_na_with_z) {
        sx.insert(data[i][id_x]);
        sy.insert(data[i][id_y]);
      }
      rx_reduced = sx.size();
      ry_reduced = sy.size();

      if (test_mar) {
        TempAllocatorScope scope;

        TempGrid2d<int> counts2(all_levels[id_x], all_levels[id_y], 0);
        // fill table
        for (int k = 0; k < n_samples_with_z; k++) {
          int i = sample_non_na_with_z[k];
          ++counts2(data[i][id_x], data[i][id_y]);
        }
        double cplxMdl = cache->getLog(n_samples_with_z);
        double kldiv = exp(-n_samples_with_z * kl(counts2, freqs1) + cplxMdl);

        if (kldiv < 1) isGoodCandidate = false;
      }
    }

    if (!isGoodCandidate) continue;

    // Initialze variables
    int Lxyui = sortedSample(kz0, 1);  // min Lxyui
    int Lyui = sortedSample(kz0, 2);   // min Lyui
    int Lui = sortedSample(kz0, 3);    // min Lui

    int X = sortedSample(kz0, 4);
    int Y = sortedSample(kz0, 5);
    int i_sample = sortedSample(kz0, 0);
    int Z = data[i_sample][idxzi];  // first Z
    // to terminate loop properly below
    sortedSample(n_samples_non_na, 0) = i_sample;

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
    TempVector<int> Nxui(r_list[0], 0);
    TempVector<int> Nx(r_list[0], 0);
    TempVector<int> Ny(r_list[1], 0);
    r_list[z] = all_levels[idxzi];
    TempVector<double> Pxyuiz(r_list[z], 0);
    TempVector<int> Nyuiz(r_list[z], 0);
    TempVector<int> Nuiz(r_list[z], 0);
    TempVector<int> Nz(r_list[z], 0);
    Pxyuiz[Z] = weights[sortedSample(kz0, 0)];
    TempGrid2d<int> Nxuiz(r_list[0], r_list[z], 0);  // [X][Z]
    // make the counts and compute mutual infos & logCs
    for (int k = kz0 + 1; k <= n_samples_non_na; k++) {
      // check whether zi contains NA
      Z = data[sortedSample(k, 0)][idxzi];
      if (Z == -1) continue;

      if (sortedSample(k, 1) <= Lxyui) {
        Pxyuiz[Z] += weights[sortedSample(k, 0)];
        continue;
      }
      Lxyui = sortedSample(k, 1);

      int Nxyuiz = 0;
      for (int l = 0; l < r_list[z]; l++) {
        double Pxyuizl = Pxyuiz[l];
        if (Pxyuizl <= 0) continue;

        int Nxyuizl = (int)Pxyuizl;
        if (n_samples_eff != n_samples) {
          if (randomrescaling) {
            Nxyuizl = (int)floor(Pxyuizl);
            double r = Pxyuizl - Nxyuiz;
            double rr = R::runif(0, 1);
            if (r > rr) ++Nxyuizl;
          }
        }

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
      Pxyuiz[Z] = weights[sortedSample(k, 0)];

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

      if (k < n_samples_non_na) X = sortedSample(k, 4);
      if (k < n_samples_non_na) Y = sortedSample(k, 5);

      if (sortedSample(k, 2) == Lyui) continue;
      Lyui = sortedSample(k, 2);

      if (Nyui > 0) {
        double NlogN = Nyui * cache->getLog(Nyui);
        info_yui_x -= NlogN;
        info_yui_z -= NlogN;
        info_ui_y += NlogN;
        if (cplx != MDL) {
          logC_yui_x += cache->getLogC(Nyui, r_list[0]);
          logC_yui_z += cache->getLogC(Nyui, r_list[z]);
        }

        for (int l = 0; l < r_list[z]; l++) {
          int Nyuizl = Nyuiz[l];
          if (Nyuizl > 0) {
            double NlogN = Nyuizl * cache->getLog(Nyuizl);
            info_yui_z += NlogN;
            info_uiz_y += NlogN;
            info_yuiz_x -= NlogN;
            if (cplx != MDL) {
              logC_yuiz_x += cache->getLogC(Nyuizl, r_list[0]);
            }
            Nyuizs += Nyuizl;
            Nyuiz[l] = 0;
          }
        }
        Nyuis += Nyui;
        Nyui = 0;
      }

      if (sortedSample(k, 3) == Lui) continue;
      Lui = sortedSample(k, 3);

      if (Nui <= 0) continue;

      double NlogN = Nui * cache->getLog(Nui);
      info_ui_x -= NlogN;
      info_ui_y -= NlogN;
      info_ui_z -= NlogN;
      if (cplx != MDL) {
        logC_ui_x += cache->getLogC(Nui, r_list[0]);
        logC_ui_y += cache->getLogC(Nui, r_list[1]);
        logC_ui_z += cache->getLogC(Nui, r_list[z]);
      }
      Nuis += Nui;
      Nui = 0;

      for (int l = 0; l < r_list[z]; l++) {
        int Nuizl = Nuiz[l];
        if (Nuizl > 0) {
          double NlogN = Nuizl * cache->getLog(Nuizl);
          info_ui_z += NlogN;
          info_uiz_x -= NlogN;
          info_uiz_y -= NlogN;
          if (cplx != MDL) {
            logC_uiz_x += cache->getLogC(Nuizl, r_list[0]);
            logC_uiz_y += cache->getLogC(Nuizl, r_list[1]);
          }
          Nuizs += Nuizl;
          Nuiz[l] = 0;
        }
      }

      for (int j = 0; j < r_list[0]; j++) {
        int Nxuij = Nxui[j];
        if (Nxuij == 0) continue;

        double NlogN = Nxuij * cache->getLog(Nxuij);
        info_xui_y -= NlogN;
        info_xui_z -= NlogN;
        info_ui_x += NlogN;
        if (cplx != MDL) {
          logC_xui_y += cache->getLogC(Nxuij, r_list[1]);
          logC_xui_z += cache->getLogC(Nxuij, r_list[z]);
        }
        Nxuis += Nxuij;
        Nxui[j] = 0;

        for (int l = 0; l < r_list[z]; l++) {
          int Nxuizjl = Nxuiz(j, l);
          if (Nxuizjl == 0) continue;

          double NlogN = Nxuizjl * cache->getLog(Nxuizjl);
          info_xui_z += NlogN;
          info_uiz_x += NlogN;
          info_xuiz_y -= NlogN;
          if (cplx != MDL) {
            logC_xuiz_y += cache->getLogC(Nxuizjl, r_list[1]);
          }
          Nxuizs += Nxuizjl;
          Nxuiz(j, l) = 0;
        }
      }
    }
    // increment info with Nx[X], Ny[Y] and Nz[Z] contributions
    for (int j = 0; j < r_list[0]; j++) {
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
    for (int j = 0; j < r_list[1]; j++) {
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
    for (int l = 0; l < r_list[z]; l++) {
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
        Prui *= r_list[j];
      logC_xui_y = 0.5 * (r_list[1] - 1) * (r_list[0] * Prui - 1) * logN;
      logC_yui_x = 0.5 * (r_list[0] - 1) * (r_list[1] * Prui - 1) * logN;
      logC_ui_y = 0.5 * (r_list[1] - 1) * (Prui - 1) * logN;
      logC_ui_x = 0.5 * (r_list[0] - 1) * (Prui - 1) * logN;
      logC_uiz_y = 0.5 * (r_list[1] - 1) * (r_list[z] * Prui - 1) * logN;
      logC_yui_z = 0.5 * (r_list[z] - 1) * (r_list[1] * Prui - 1) * logN;
      logC_ui_z = 0.5 * (r_list[z] - 1) * (Prui - 1) * logN;
    }
    double logC3xy_ui = 0.5 * (logC_xui_y + logC_yui_x);
    double logC2xy_ui = 0.5 * (logC_ui_y + logC_ui_x);

    // forbid negative 2-point Ik // change 20200206
    // if(info3xy_ui - logC3xy_ui - info2xy_ui + logC2xy_ui<0) {
    //  info3xy_ui=0;
    //  logC3xy_ui=0;
    //  info2xy_ui=0;
    //  logC2xy_ui=0;
    //}

    // NI(yz|ui)
    double info3yz_ui = 0.5 * (info_uiz_y + info_yui_z);
    double info2yz_ui = 0.5 * (info_ui_y + info_ui_z);

    double logC3yz_ui = 0.5 * (logC_uiz_y + logC_yui_z);
    double logC2yz_ui = 0.5 * (logC_ui_y + logC_ui_z);

    // forbid negative 2-point Ik // change 20200206
    // if(info3yz_ui - logC3yz_ui - info2yz_ui + logC2yz_ui<0) {
    //  info3yz_ui=0;
    //  logC3yz_ui=0;
    //  info2yz_ui=0;
    //  logC2yz_ui=0;
    //}

    // NI(xz|ui)
    double info3xz_ui = 0.5 * (info_xui_z + info_uiz_x);
    double info2xz_ui = 0.5 * (info_ui_z + info_ui_x);

    if (cplx == MDL) {
      logC_uiz_x = 0.5 * (r_list[0] - 1) * (r_list[z] * Prui - 1) * logN;
      logC_xui_z = 0.5 * (r_list[z] - 1) * (r_list[0] * Prui - 1) * logN;
    }
    double logC3xz_ui = 0.5 * (logC_uiz_x + logC_xui_z);
    double logC2xz_ui = 0.5 * (logC_ui_x + logC_ui_z);

    // forbid negative 2-point Ik // change 20200206
    // if(info3xz_ui - logC3xz_ui - info2xz_ui + logC2xz_ui<0) {
    //  info3xz_ui=0;
    //  logC3xz_ui=0;
    //  info2xz_ui=0;
    //  logC2xz_ui=0;
    //}

    // NI(xy|uiz)
    double info3xy_uiz = 0.5 * (info_xuiz_y + info_yuiz_x);
    double info2xy_uiz = 0.5 * (info_uiz_y + info_uiz_x);

    if (cplx == MDL) {
      Prui *= r_list[z];
      logC_xuiz_y = 0.5 * (r_list[1] - 1) * (r_list[0] * Prui - 1) * logN;
      logC_yuiz_x = 0.5 * (r_list[0] - 1) * (r_list[1] * Prui - 1) * logN;
      logC_uiz_y = 0.5 * (r_list[1] - 1) * (Prui - 1) * logN;
      logC_uiz_x = 0.5 * (r_list[0] - 1) * (Prui - 1) * logN;
    }
    double logC3xy_uiz = 0.5 * (logC_xuiz_y + logC_yuiz_x);
    double logC2xy_uiz = 0.5 * (logC_uiz_y + logC_uiz_x);

    // forbid negative 2-point Ik // change 20200206
    // if(info3xy_uiz - logC3xy_uiz - info2xy_uiz + logC2xy_uiz<0) {
    //  info3xy_uiz=0;
    //  logC3xy_uiz=0;
    //  info2xy_uiz=0;
    //  logC2xy_uiz=0;
    //}

    int N_xyuiz = n_samples_zi;
    info_xy_ui = info3xy_ui - info2xy_ui;     // info NI(xy|ui)
    logC_xy_ui = logC3xy_ui - logC2xy_ui;     // cplx k_(xy|ui)
    info_yz_ui = info3yz_ui - info2yz_ui;     // info NI(yz|ui)
    logC_yz_ui = logC3yz_ui - logC2yz_ui;     // cplx k_(yz|ui)
    info_xz_ui = info3xz_ui - info2xz_ui;     // info NI(xz|ui)
    logC_xz_ui = logC3xz_ui - logC2xz_ui;     // cplx k_(xz|ui)
    info_xy_uiz = info3xy_uiz - info2xy_uiz;  // info NI(xy|uiz)
    logC_xy_uiz = logC3xy_uiz - logC2xy_uiz;  // cplx k_(xy|uiz)

    // compute score and store z with max score
    double xz = info_xz_ui - info_xy_ui;
    double yz = info_yz_ui - info_xy_ui;
    double xyz = info_xy_ui - info_xy_uiz;

    if (k23 && rx_reduced > 1 && ry_reduced > 1) {
      xz -= logC_xz_ui - logC_xy_ui;
      yz -= logC_yz_ui - logC_xy_ui;
      xyz -= logC_xy_ui - logC_xy_uiz;
    }
    double first{0}, second{0};
    if (xz < yz) {
      first = xz;
      second = yz;
    } else {
      first = yz;
      second = xz;
    }
    double dpi = first - log1p(exp(first - second));

    double Rzi{0};
    if (xyz < dpi) {
      Rzi = xyz;
    } else {
      Rzi = dpi;  // final score for each zi ;)
    }

    if (Rzi > R_top) {
      N_xyuiz_top = N_xyuiz;
      NIxy_ui_top = info_xy_ui;
      k_xy_ui_top = logC_xy_ui;

      z_top = zi;
      NIxy_uiz_top = info_xy_uiz;
      k_xy_uiz_top = logC_xy_uiz;

      R_top = Rzi;
      NIxyz_ui_top = info_xy_ui - info_xy_uiz;
      // to fit eq(20) and eq(22) in BMC Bioinfo 2016
      k_xyz_ui_top = -logC_xy_ui + logC_xy_uiz;
    }
  }  // end loop on all z

  ptrRetValues[0] = N_xyuiz_top;
  ptrRetValues[1] = NIxy_ui_top;
  ptrRetValues[2] = k_xy_ui_top;

  ptrRetValues[3] = z_top;
  ptrRetValues[4] = NIxy_uiz_top;
  ptrRetValues[5] = k_xy_uiz_top;

  ptrRetValues[6] = R_top;
  ptrRetValues[7] = NIxyz_ui_top;
  ptrRetValues[8] = k_xyz_ui_top;
  return ptrRetValues;
}

}  // namespace computation
}  // namespace miic
