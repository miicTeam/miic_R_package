#include "computation_cache.h"

namespace miic {
namespace computation {
namespace computation_impl {

double CtermCache::getLogC(int n, int level) {
  if (n == 0 || level == 0) return 0;
  if (level <= kLevelLimit) {
    double res = log_c_(n - 1, level - 1);
    if (res != -1) return res;
  }
  if (level == 1) {
    log_c_(n - 1, level - 1) = 0;
    return 0;
  } else if (level == 2) {
    double c2{0};
    if (n > kApproxLimit) {
      c2 = sqrt(n * M_PI_2) *
            exp(sqrt(8 / (9 * n * M_PI)) + (3 * M_PI - 16) / (36 * n * M_PI));
    } else {
      for (int h = 0; h <= n; ++h)
        c2 += exp(
            getLogChoose(n, h) + n_log_n_[h] + n_log_n_[n - h] - n_log_n_[n]);
    }
    double res = log(c2);
    log_c_(n - 1, level - 1) = res;
    return res;
  } else {
    // Use the recurrence C(n, r) = C(n, r - 1) + n / (r - 2) * C(n, r - 2)
    // When n and r are relatively large (e.g., n = 20000, r = 100), C(n, r)
    // will cause overflow of double (of the order pow(n, r)), thus we need the
    // log version of the above recurrence:
    // C(n, r) / C(n, r - 1) = 1 + n / ((r - 2) * C(n, r - 1) / C(n, r - 2))
    // log(C(n, r)) = sum_(i=r_0)^(i=r)[log(C(n, i) / C(n, i - 1))]
    int r_0 = level <= kLevelLimit ? level : kLevelLimit + 1;
    double res = getLogC(n, r_0 - 1);
    double c_ratio = exp(res - getLogC(n, r_0 - 2));
    for (int i = r_0; i <= level; ++i) {
      c_ratio = 1 + static_cast<double>(n) / ((i - 2) * c_ratio);
      res += log(c_ratio);
    }
    if (level <= kLevelLimit)
      log_c_(n - 1, level - 1) = res;
    return res;
  }
}

}  // namespace computation_impl
}  // namespace computation
}  // namespace miic
