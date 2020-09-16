#ifndef MIIC_MUTUAL_INFORMATION_H_
#define MIIC_MUTUAL_INFORMATION_H_

#define _USE_MATH_DEFINES
#include <cmath>    // std::log
#include <vector>

#include "computation_cache.h"
#include "structure.h"

namespace miic {
namespace computation {
namespace detail {

using namespace structure;
using namespace utility;
using std::log;
using std::vector;

// rux -> 0:x,1;u,2:ux
template <typename Cx, typename Cu, typename Cux, typename Crux,
    typename = void_t<IsIntContainer<Cx>, IsIntContainer<Cu>,
        IsIntContainer<Cux>, IsIntContainer<Crux>>>
InfoBlock computeMI_knml(const Cx& xfactors, const Cu& ufactors,
    const Cux& uxfactors, const Crux& rux, int n_eff,
    const TempVector<double>& sample_weights, std::shared_ptr<CtermCache> cache,
    int flag) {
  TempAllocatorScope scope;

  int n_samples = ufactors.size();
  TempVector<double> nx(rux[0]);
  TempVector<double> nu(rux[1]);
  TempVector<double> nux(rux[2]);
  for (int i = 0; i < n_samples; i++) {
    nx[xfactors[i]] += sample_weights[i];
    nu[ufactors[i]] += sample_weights[i];
    nux[uxfactors[i]] += sample_weights[i];
  }

  double Hux = 0, Hu = 0, Hx = 0, SC = 0;
  for (const auto x : nx) {
    if (x <= 0) continue;

    Hx -= x * log(x);
    if (flag == 0)
      SC += cache->getLogC(std::max(1, static_cast<int>(x + 0.5)), rux[1]);
  }
  for (const auto u : nu) {
    if (u <= 0) continue;

    Hu -= u * log(u);
    if (flag == 0 || flag == 1)
      SC += cache->getLogC(std::max(1, static_cast<int>(u + 0.5)), rux[0]);
  }
  for (const auto ux : nux) {
    if (ux <= 0) continue;

    Hux -= ux * log(ux);
  }

  if (flag == 0) {
    SC -= cache->getLogC(n_eff, rux[0]);
    SC -= cache->getLogC(n_eff, rux[1]);
    SC *= 0.5;
  }

  double Iux = n_eff * cache->getLog(n_eff) + (Hu + Hx - Hux);

  return InfoBlock{n_eff, Iux, SC};
}

// rux -> 0:x,1;u,2:ux
template <typename Cx, typename Cu, typename Cux, typename Crux,
    typename = void_t<IsIntContainer<Cx>, IsIntContainer<Cu>,
        IsIntContainer<Cux>, IsIntContainer<Crux>>>
InfoBlock computeMI_kmdl(const Cx& xfactors, const Cu& ufactors,
    const Cux& uxfactors, const Crux& rux,
    std::shared_ptr<CtermCache> cache, int flag = 0) {
  TempAllocatorScope scope;

  int n_samples = xfactors.size();
  TempVector<int> nx(rux[0]);
  TempVector<int> nu(rux[1]);
  TempVector<int> nux(rux[2]);
  for (int i = 0; i < n_samples; i++) {
    ++nx[xfactors[i]];
    ++nu[ufactors[i]];
    ++nux[uxfactors[i]];
  }

  double Hux = 0, Hu = 0, Hx = 0, SC = 0;
  for (const auto x : nx)
    if (x > 0) Hx -= x * cache->getLog(x);
  for (const auto u : nu)
    if (u > 0) Hu -= u * cache->getLog(u);
  for (const auto ux : nux)
    if (ux > 0) Hux -= ux * cache->getLog(ux);

  SC = 0.5 * cache->getLog(n_samples);
  if (flag == 0 || flag == 1) SC *= (rux[0] - 1);
  if (flag == 0 || flag == 2) SC *= (rux[1] - 1);

  double Iux = n_samples * cache->getLog(n_samples) + (Hu + Hx - Hux);

  return InfoBlock{n_samples, Iux, SC};
}

}  // namespace detail
using detail::computeMI_knml;
using detail::computeMI_kmdl;

void jointfactors_uiyx(const structure::TempGrid2d<int>& datafactors,
    const structure::TempVector<int>& r_list, int exclude,
    structure::TempGrid2d<int>& uiyxfactors,
    structure::TempVector<int>& r_joint_list);
void jointfactors_u(const structure::TempGrid2d<int>& datafactors,
    const structure::TempVector<int>& r_list,
    structure::TempVector<int>& ufactors, int& r_joint);

inline int getHashU(const structure::TempGrid2d<int>& data,
    const structure::TempVector<int>& r_list,
    const structure::TempVector<int>& ui_list, int i) {
  int rui_joint{1};
  int hash_u{0};
  for (const auto u : ui_list) {
    hash_u += data(u, i) * rui_joint;
    rui_joint *= r_list[u];
  }
  return hash_u;
};

structure::TempVector<int> getDataOrder(const structure::TempGrid2d<int>& data,
    const structure::TempVector<int>& r_list,
    const structure::TempVector<int>& var_idx);
}  // namespace computation
}  // namespace miic

#endif  // MIIC_MUTUAL_INFORMATION_H_
