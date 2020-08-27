#ifndef MIIC_MUTUAL_INFORMATION_H_
#define MIIC_MUTUAL_INFORMATION_H_

#include <vector>

#include "computation_cache.h"
#include "structure.h"

namespace miic {
namespace computation {
namespace detail {

using namespace structure;
using namespace utility;
using std::vector;

// update datafactors of a variable from the cut positions vector <cut>
// INPUT:
// d: index of variable in datafactors
// varidx: index of variable in sortidx
template <typename Cdf, typename = IsIntContainer<Cdf>, typename Ccut,
    typename = IsIntContainer<Ccut>>
void update_datafactors(
    const vector<int>& sortidx, Cdf& datafactor, const Ccut& cut) {
  int index = 0;
  for (size_t j = 0; j < sortidx.size(); ++j) {
    int level = sortidx[j];
    if (static_cast<int>(j) > cut[index]) index++;
    datafactor[level] = index;
  }
  return;
}

// rux -> 0:x,1;u,2:ux
template <typename Cx, typename = IsIntContainer<Cx>, typename Cu,
    typename = IsIntContainer<Cu>, typename Cux, typename = IsIntContainer<Cux>>
vector<double> computeMI_knml(const Cx& xfactors, const Cu& ufactors,
    const Cux& uxfactors, int* rux, int n, int n_eff,
    std::vector<double> sample_weights, std::shared_ptr<CtermCache> cache,
    int flag = 0) {
  TempAllocatorScope scope;

  TempVector<double> nx(rux[0]);
  TempVector<double> nu(rux[1]);
  TempVector<double> nux(rux[2]);
  for (int i = 0; i < n; i++) {
    nx[xfactors[i]] += sample_weights[i];
    nu[ufactors[i]] += sample_weights[i];
    nux[uxfactors[i]] += sample_weights[i];
  }

  double Hux = 0, Hu = 0, Hx = 0, SC = 0;
  for (const auto x : nx) {
    if (x <= 0) continue;

    Hx -= x * log(x);
    if (flag == 0 || flag == 2)
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
  }

  vector<double> I(2);
  I[0] = cache->getLog(n_eff) + (Hu + Hx - Hux) / n_eff;
  if (flag == 0)
    I[1] = I[0] - 0.5 * SC / n_eff;
  else
    I[1] = I[0] - SC / n_eff;

  return I;
}

// rux -> 0:x,1;u,2:ux
template <typename Cx, typename = IsIntContainer<Cx>, typename Cu,
    typename = IsIntContainer<Cu>, typename Cux, typename = IsIntContainer<Cux>>
vector<double> computeMI_kmdl(const Cx& xfactors, const Cu& ufactors,
    const Cux& uxfactors, int* rux, int n, std::shared_ptr<CtermCache> cache,
    int flag = 0) {
  TempAllocatorScope scope;

  TempVector<int> nx(rux[0]);
  TempVector<int> nu(rux[1]);
  TempVector<int> nux(rux[2]);
  for (int i = 0; i < n; i++) {
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

  SC = 0.5 * cache->getLog(n);
  if (flag == 0 || flag == 1) SC *= (rux[0] - 1);
  if (flag == 0 || flag == 2) SC *= (rux[1] - 1);

  vector<double> I(2);
  I[0] = cache->getLog(n) + (Hu + Hx - Hux) / n;
  I[1] = I[0] - SC / n;

  return I;
}

}  // namespace detail
using detail::computeMI_knml;
using detail::computeMI_kmdl;
using detail::update_datafactors;

void jointfactors_uiyx(int** datafactors, int exclude, int n, int n_ui,
    const structure::TempVector<int>& r, int** uiyxfactors, int* ruiyx);
void jointfactors_u(int** datafactors, int* ptrIdx, int n, int n_ui,
    const structure::TempVector<int>& r, int* ufactors, int* ru);
}  // namespace computation
}  // namespace miic

#endif  // MIIC_MUTUAL_INFORMATION_H_
