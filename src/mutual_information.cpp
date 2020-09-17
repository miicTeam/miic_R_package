#include "mutual_information.h"

#include <algorithm>  // std::copy
#include <numeric>    // std::iota

namespace miic {
namespace computation {

using structure::TempGrid2d;
using structure::TempVector;
using utility::TempAllocatorScope;

// INPUT
// datafactors: [0, ]: x, [1, ]: y, [2 ... ]: {ui}
// exclude: ignore index <exclude>
// r_list: number of levels of each variable [x, y, {ui}]
// OUTPUT
// uiyxfactors: Joint datafactors. [ui, uiy, uix, uixy]
// r_joint_list: Number of joint levels. 0: u, 1: uy, 2: ux, 3: uyx
void setUyxJointFactors(const TempGrid2d<int>& factors,
    const TempVector<int>& r_list, int exclude, TempGrid2d<int>& uyxfactors,
    TempVector<int>& ruyx) {
  TempAllocatorScope scope;

  int n_ui = factors.n_rows() - 2;
  TempVector<int> ui_list;
  ui_list.reserve(n_ui);
  for (int u = 2; u < n_ui + 2; ++u) {
    if (u != exclude)
      ui_list.push_back(u);
  }
  // Set u joint factors
  ruyx[0] = setJointFactors(factors, r_list, ui_list, uyxfactors.getRow(0));
  // Copy x and y factors and n_levels
  std::copy(factors.row_begin(0), factors.row_end(0), uyxfactors.row_begin(2));
  ruyx[2] = r_list[0];
  std::copy(factors.row_begin(1), factors.row_end(1), uyxfactors.row_begin(1));
  ruyx[1] = r_list[1];
  // Combine y and u to get uy, overwriting column y
  TempVector<int> var_idx{1, 0};
  ruyx[1] = setJointFactors(uyxfactors, ruyx, var_idx, uyxfactors.getRow(1));
  // Combine x and uiy to get uiyx
  var_idx.assign({2, 1});
  ruyx[3] = setJointFactors(uyxfactors, ruyx, var_idx, uyxfactors.getRow(3));
  // Combine x and ui to get uix, overwriting column x
  var_idx.assign({2, 0});
  ruyx[2] = setJointFactors(uyxfactors, ruyx, var_idx, uyxfactors.getRow(2));
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
      // The value of each element in r_list can't be larger than n_samples
      if (r_joint * r_list[var] > 8 * n_samples) {
        break;
      } else {
        temp_var_idx.push_back(var);
        r_joint *= r_list[var];
      }
    }
    // Counting sort
    TempAllocatorScope scope;
    TempVector<int> temp_hash_list(n_samples, 0);
    fillHashList(data, r_list, temp_var_idx, temp_hash_list);

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

int fillHashList(const structure::TempGrid2d<int>& data,
    const structure::TempVector<int>& r_list,
    const structure::TempVector<int>& ui_list,
    structure::TempVector<int>& hash_list) {
  int n_ui = ui_list.size();
  if (n_ui == 1) {
    int u = ui_list[0];
    std::copy(data.row_begin(u), data.row_end(u), begin(hash_list));
    return r_list[u];
  }
  int n_samples = data.n_cols();
  if (n_ui == 2) {
    int u0 = ui_list[0], u1 = ui_list[1];
    int r0 = r_list[u0];
    for (int i = 0; i < n_samples; ++i) {
      hash_list[i] += data(u0, i) + data(u1, i) * r0;
    }
    return r0 * r_list[u1];
  }
  utility::TempAllocatorScope scope;

  structure::TempVector<int> r_joint_list(n_ui);
  int n_levels_product{1};
  for (const auto u : ui_list) {
    r_joint_list[u] = n_levels_product;
    n_levels_product *= r_list[u];
  }
  for (int i = 0; i < n_samples; ++i) {
    for (const auto u : ui_list)
      hash_list[i] += data(u, i) * r_joint_list[u];
  }
  return n_levels_product;
}

}  // namespace computation
}  // namespace miic
