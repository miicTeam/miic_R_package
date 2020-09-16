#include "computation_cache.h"

namespace miic {
namespace computation {
namespace detail {

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
    // log version of the above recursion:
    // C(n, r) / C(n, r - 1) = 1 + n / ((r - 2) * C(n, r - 1) / C(n, r - 2))
    // log(C(n, r)) = sum_(i=r_0)^(i=r)[log(C(n, i) / C(n, i - 1))]

    int r = level <= kLevelLimit ? level : kLevelLimit + 1;
    // res stores log(C(n, r)), res_aux stores log(C(n, r - 1))
    double res{-1}, res_aux{-1};
    while (res == -1 || res_aux == -1) {
      --r;  // Backtrack to find nearest cached results, avoid recursive call
      res = log_c_(n - 1, r - 1);      // log(C(n, r))
      res_aux = log_c_(n - 1, r - 2);  // log(C(n, r - 1))
      // Will break for sure for r == 2, since log(C(n, 2)) and log(C(n, 1)) are
      // already computed during the construction of cache.
    }
    double c_ratio = exp(res - res_aux);  // C(n, r) / C(n, r - 1)
    for (int l = r + 1; l <= level; ++l) {
      c_ratio = 1 + static_cast<double>(n) / ((l - 2) * c_ratio);  // recursion
      res += log(c_ratio);  // res == log(C(n, l))
      if (l <= kLevelLimit)
        log_c_(n - 1, l - 1) = res;
    }
    return res;
  }
}

}  // namespace detail
}  // namespace computation
}  // namespace miic
