#include "computation_cache.h"

namespace miic {
namespace computation {
namespace computation_impl {

double CtermCache::getC(int n, int level) {
  if (n == 0 || level == 0) return 1;
  double res{0};
  if (level <= kMaxLevel) {
    res = c_(n - 1, level - 1);
    if (res != -1) return res;
  }
  if (level == 1) {
    c_(n - 1, level - 1) = 1;
    return 1;
  } else if (level == 2) {
    if (n > kApproxLimit) {
      res = sqrt(n * M_PI_2) *
            exp(sqrt(8 / (9 * n * M_PI)) + (3 * M_PI - 16) / (36 * n * M_PI));
    } else {
      res = 0;
      for (int h = 0; h <= n; ++h)
        res += exp(
            getLogChoose(n, h) + n_log_n_[h] + n_log_n_[n - h] - n_log_n_[n]);
    }
    c_(n - 1, level - 1) = res;
  } else if (level <= kMaxLevel) {
    res = getC(n, level - 1) + getC(n, level - 2) / (level - 2) * n;
    c_(n - 1, level - 1) = res;
  } else {
    double r1{getC(n, kMaxLevel)}, r2{getC(n, kMaxLevel - 1)};
    for (int i = kMaxLevel + 1; i <= level; ++i) {
      res = r1 + r2 / (i - 2) * n;
      r2 = r1;
      r1 = res;
    }
  }
  return res;
}

}  // namespace computation_impl
}  // namespace computation
}  // namespace miic
