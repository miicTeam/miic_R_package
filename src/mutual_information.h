#ifndef MIIC_MUTUAL_INFORMATION_H_
#define MIIC_MUTUAL_INFORMATION_H_

#define _USE_MATH_DEFINES
#include <cmath>    // std::log

#include "computation_cache.h"
#include "structure.h"

namespace miic {
namespace computation {

void setUyxJointFactors(const structure::TempGrid2d<int>& datafactors,
    const structure::TempVector<int>& r_list, int exclude,
    structure::TempGrid2d<int>& uiyxfactors,
    structure::TempVector<int>& r_joint_list);

structure::TempVector<int> getDataOrder(const structure::TempGrid2d<int>& data,
    const structure::TempVector<int>& r_list,
    const structure::TempVector<int>& var_idx);

int fillHashList(const structure::TempGrid2d<int>& data,
    const structure::TempVector<int>& r_list,
    const structure::TempVector<int>& ui_list,
    structure::TempVector<int>& hash_list);

namespace detail {

using namespace structure;
using namespace utility;
using std::log;
using std::lround;
constexpr double kPrecision = 1.e-10;

// rux: number of levels of each (joint) variable [x, u, ux]
// cplx 0: BIC, 1: NML
// flag (for cplx == 1 only) 0: mutual info, 1: conditional mutual info
// When flag == 1 && cplx == 1, x and u are not symmetrical, x represents single
// variable, whereas u represents joint variable (see def of cond mutual info)
template <typename Cx, typename Cu, typename Cux, typename Crux,
    typename = void_t<IsIntContainer<Cx>, IsIntContainer<Cu>,
        IsIntContainer<Cux>, IsIntContainer<Crux>>>
InfoBlock computeMI(const Cx& xfactors, const Cu& ufactors,
    const Cux& uxfactors, const Crux& rux, double n_eff,
    const TempVector<double>& sample_weights, std::shared_ptr<CtermCache> cache,
    int cplx, int flag) {
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

  double Hux{0}, Hu{0}, Hx{0}, sc{0};
  for (const auto x : nx) {
    if (x <= 0) continue;

    Hx -= x * log(x);
    if (cplx == 1 && flag == 0)
      sc += cache->getLogC(std::max((long)1, lround(x)), rux[1]);
  }
  for (const auto u : nu) {
    if (u <= 0) continue;

    Hu -= u * log(u);
    if (cplx == 1)
      sc += cache->getLogC(std::max((long)1, lround(u)), rux[0]);
  }
  for (const auto ux : nux) {
    if (ux <= 0) continue;

    Hux -= ux * log(ux);
  }

  if (cplx == 1) {
    if (flag == 0) {
      auto n_eff_long = lround(n_eff);
      sc -= cache->getLogC(n_eff_long, rux[0]);
      sc -= cache->getLogC(n_eff_long, rux[1]);
      sc *= 0.5;
    }
  } else {
    sc = 0.5 * log(n_eff) * (rux[0] - 1) * (rux[1] - 1);
  }

  double Iux = n_eff * log(n_eff) + (Hu + Hx - Hux);
  if(rux[0] == 1 || rux[1] == 1 || Iux < kPrecision) Iux = 0; // The summed Iux may not be 0 due to precision issues

  return InfoBlock{n_eff, Iux, sc};
}

template <typename Cjf, typename = IsIntContainer<Cjf>>
int setJointFactors(const TempGrid2d<int>& factors,
    const TempVector<int>& r_list, const TempVector<int>& var_idx,
    Cjf&& joint_factors) {
  if (var_idx.size() == 1) {
    const auto row = factors.getConstRow(var_idx[0]);
    std::copy(std::begin(row), std::end(row), std::begin(joint_factors));
    return r_list[var_idx[0]];
  }
  int n_samples = factors.n_cols();
  TempAllocatorScope scope;
  // Compute unique hash value for each sample in the joint space
  TempVector<int> hash_u(n_samples, 0);
  int level_product = fillHashList(factors, r_list, var_idx, hash_u);

  int r_joint{0};  // get ready to count
  if (level_product <= 8 * n_samples) {
    // Use large sparse vectors, no sort
    TempVector<int> counts(level_product);
    for (const auto h : hash_u)
      counts[h] = 1;
    // Order of the levels follow the order of the hash values,
    // which are sorted automatically (as indices) with the sparse vector.
    for (auto& l : counts)
      if (l == 1) l = r_joint++;

    for (int i = 0; i < n_samples; ++i)
      joint_factors[i] = counts[hash_u[i]];
  } else {
    // Fall back to radix sort
    TempVector<int> order = getDataOrder(factors, r_list, var_idx);

    auto hash_u_prev = hash_u[order[0]];
    for (const auto index : order) {
      auto hash_u_current = hash_u[index];
      if (hash_u_current > hash_u_prev) ++r_joint;

      joint_factors[index] = r_joint;
      hash_u_prev = hash_u_current;
    }
    ++r_joint;
  }
  return r_joint;
}

}  // namespace detail
using detail::computeMI;
using detail::setJointFactors;
}  // namespace computation
}  // namespace miic

#endif  // MIIC_MUTUAL_INFORMATION_H_
