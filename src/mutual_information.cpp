#include "mutual_information.h"

#include <algorithm>  // std::max, std::sort
#include <numeric>  // std::iota

namespace miic {
namespace computation {

using structure::TempGrid2d;
using structure::TempVector;
using utility::TempAllocatorScope;

// INPUT
// datafactors: [0, ]: x, [1, ]: y, [2 ... ]: {ui}
// exclude: ignore index <exclude> in the {ui}
// n: Number of samples
// n_ui: number of {ui}
// r: number of levels of each variable [x, y, {ui}]
// OUTPUT
// uiyxfactors: Joint datafactors. [ui, uiy, uix, uixy]
// ruiyx[0,1,2,3]: Number of joint levels. 0: u, 1: uy, 2: ux, 3: uyx
void jointfactors_uiyx(const TempGrid2d<int>& datafactors,
    const TempVector<int>& r_list, int exclude, TempGrid2d<int>& uiyxfactors,
    TempVector<int>& r_joint_list) {
  TempAllocatorScope scope;

  int n_ui = datafactors.n_rows() - 2;
  int n = datafactors.n_cols();
  // Compute unique hash values for each sample in each of the joint spaces
  TempVector<int> hash_ui(n);
  TempVector<int> hash_uiy(n);
  TempVector<int> hash_uix(n);
  TempVector<int> hash_uiyx(n);
  for (int i = 0; i < n; ++i) {
    hash_ui[i] = 0;
    hash_uiy[i] = datafactors(1, i);
    hash_uix[i] = datafactors(0, i);
    hash_uiyx[i] = hash_uix[i] + hash_uiy[i] * r_list[0];

    int Pbin_ui = 1;
    for (int l = n_ui + 1; l >= 2; --l) {
      if (l == exclude) continue;

      int df = datafactors(l, i) * Pbin_ui;
      hash_uiyx[i] += df * r_list[1] * r_list[0];
      hash_uix[i] += df * r_list[0];
      hash_uiy[i] += df * r_list[1];
      hash_ui[i] += df;
      Pbin_ui *= r_list[l];
    }
  }

  bool too_many_levels = false;
  int n_joint_levels = 1;
  for (int i = 0; i < n_ui + 2 && !too_many_levels; ++i) {
    n_joint_levels *= r_list[i];
    if (n_joint_levels > 8 * n) {
      too_many_levels = true;  // Too large for the sparse vectors
    }
  }

  std::fill(begin(r_joint_list), end(r_joint_list), 0);  // get ready to count
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
      if (levels_ui[i] == 1) levels_ui[i] = r_joint_list[0]++;
      if (levels_uiy[i] == 1) levels_uiy[i] = r_joint_list[1]++;
      if (levels_uix[i] == 1) levels_uix[i] = r_joint_list[2]++;
      if (levels_uiyx[i] == 1) levels_uiyx[i] = r_joint_list[3]++;
    }

    for (int i = 0; i < n; ++i) {
      uiyxfactors(0, i) = levels_ui[hash_ui[i]];      // ui
      uiyxfactors(1, i) = levels_uiy[hash_uiy[i]];    // uiy
      uiyxfactors(2, i) = levels_uix[hash_uix[i]];    // uix
      uiyxfactors(3, i) = levels_uiyx[hash_uiyx[i]];  // uiyx
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
      if (hash_ui_current > hash_ui_prev) ++r_joint_list[0];
      if (hash_uix_current > hash_uix_prev) ++r_joint_list[2];

      uiyxfactors(0, index) = r_joint_list[0];  // ui
      uiyxfactors(2, index) = r_joint_list[2];  // uix
      hash_ui_prev = hash_ui_current;
      hash_uix_prev = hash_uix_current;
    }
    // hash_uixy[a] < hash_uixy[b] -> hash_uiy[a] <= hash_uiy[b]
    int hash_uiy_prev = hash_uiy[orderSample_uiyx[0]];
    int hash_uiyx_prev = hash_uiyx[orderSample_uiyx[0]];
    for (const auto index : orderSample_uiyx) {
      auto hash_uiy_current = hash_uiy[index];
      auto hash_uiyx_current = hash_uiyx[index];
      if (hash_uiy_current > hash_uiy_prev) ++r_joint_list[1];
      if (hash_uiyx_current > hash_uiyx_prev) ++r_joint_list[3];

      uiyxfactors(1, index) = r_joint_list[1];  // uiy
      uiyxfactors(3, index) = r_joint_list[3];  // uiyx
      hash_uiy_prev = hash_uiy_current;
      hash_uiyx_prev = hash_uiyx_current;
    }
    // number of joint levels
    ++r_joint_list[0];  // ui
    ++r_joint_list[1];  // uiy
    ++r_joint_list[2];  // uix
    ++r_joint_list[3];  // uiyx
  }
  return;
}

// INPUT:
// datarank, datafactors, cut
// OUTPUT
// return joint datafactors ui , with number of levels rui
void jointfactors_u(const TempGrid2d<int>& datafactors,
    const TempVector<int>& r_list, TempVector<int>& ufactors, int& r_joint) {
  int n_ui = datafactors.n_rows();
  int n = datafactors.n_cols();
  if (n_ui == 1) {
    const auto row = datafactors.getConstRow(0);
    std::copy(std::begin(row), std::end(row), begin(ufactors));
    r_joint = r_list[0];
    return;
  }
  TempAllocatorScope scope;
  // Compute unique hash value for each sample in the joint space
  TempVector<int> hash_u(n, 0);
  for (int i = 0; i < n; ++i) {
    int Pbin_ui = 1;
    for (int ui = n_ui - 1; ui >= 0; --ui) {
      hash_u[i] += datafactors(ui, i) * Pbin_ui;
      Pbin_ui *= r_list[ui];
    }
  }
  bool too_many_levels = false;
  int n_joint_levels = 1;
  for (int ui = 0; ui < n_ui && !too_many_levels; ++ui) {
    n_joint_levels *= r_list[ui];
    if (n_joint_levels > 8 * n) {
      too_many_levels = true;  // Too large for the sparse vectors
    }
  }

  r_joint = 0;  // get ready to count
  if (!too_many_levels) {
    // Use large sparse vectors to have O(n) time complexity (no sort)
    TempVector<int> levels_ui(n_joint_levels);
    for (const auto h : hash_u)
      levels_ui[h] = 1;
    // Use ru as level indices, whose final value is the total numbers of
    // joint levels. Order of the levels follow the order of the hash values,
    // which are sorted automatically (as indices) with the sparse vector.
    for (auto& l : levels_ui)
      if (l == 1) l = r_joint++;

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
      if (hash_u_current > hash_u_prev) ++r_joint;

      ufactors[index] = r_joint;
      hash_u_prev = hash_u_current;
    }
    ++r_joint;
  }
  return;
}

// Sort data w.r.t. each of the variables in var_idx from begin to end, in bulk
// of variables with limited number of joint levels
TempVector<int> getDataOrder(const TempGrid2d<int>& data,
    const TempVector<int>& r_list, const TempVector<int>& var_idx) {
  int n_samples = data.n_cols();
  int n_vars = var_idx.size();

  TempVector<int> order(n_samples);
  std::iota(begin(order), end(order), 0);
  TempVector<int> new_order(order);

  TempAllocatorScope scope;

  TempVector<int> temp_var_idx;
  temp_var_idx.reserve(n_vars);
  int n_vars_done{0};
  while (n_vars_done < n_vars) {
    temp_var_idx.clear();
    int r_joint = 1;
    for (auto it = begin(var_idx) + n_vars_done; it < end(var_idx); ++it) {
      int var = *it;
      if (r_joint != 1 && r_joint * r_list[var] > 8 * n_samples) {
        break;
      } else {
        temp_var_idx.push_back(var);
        r_joint *= r_list[var];
      }
    }
    // Counting sort
    TempAllocatorScope scope;
    TempVector<int> temp_hash_list(n_samples);
    for (int i = 0; i < n_samples; ++i)
      temp_hash_list[i] = getHashU(data, r_list, temp_var_idx, i);

    TempVector<int> counts(r_joint);
    for (const auto index : order)
      ++counts[temp_hash_list[index]];
    int sum{0};
    for (auto& c : counts) {
      int temp = c;
      c = sum;
      sum += temp;
    }
    for (const auto index : order)
      new_order[counts[temp_hash_list[index]]++] = index;

    order.swap(new_order);
    n_vars_done += temp_var_idx.size();
  }
  return order;
}

}  // namespace computation
}  // namespace miic
