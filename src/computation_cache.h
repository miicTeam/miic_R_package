#ifndef MIIC_COMPUTATION_CACHE
#define MIIC_COMPUTATION_CACHE

#include "structure.h"

namespace miic {
namespace computation {

namespace computation_impl {

using std::vector;
using structure::Grid2d;

class CtermCache {
 private:
  static constexpr int kLevelLimit = 50;
  // Limit beyond which cterm is calculated by approximation instead of log
  static constexpr int kApproxLimit = 1000;
  int size_n_{1 + kApproxLimit};
  // Of size size_n_, log_n_[i] = log(i), log_n_[0] = 0 (won't be called alone)
  vector<double> log_n_;
  // Of size size_n_, n_log_n_[i] = i * log(i)
  vector<double> n_log_n_;
  // Of size size_n_, log_factorial_[i] = log(i!)
  vector<double> log_factorial_;
  // Hold <n_samples> * <kMaxLevel> C(omplexity)_n^level terms
  // with n in [1, n_samples] and level in [1, kMaxLevel]
  Grid2d<double> log_c_;

 public:
  CtermCache(int n_samples)
      : size_n_(1 + n_samples),
        log_n_(size_n_, 0),
        n_log_n_(size_n_, 0),
        log_factorial_(size_n_, 0),
        log_c_(n_samples, kLevelLimit, -1) {
    for (int i = 2; i < size_n_; ++i) {  // first two terms are zero
      double logi = log(i);
      log_n_[i] = logi;
      n_log_n_[i] = i * logi;
      log_factorial_[i] = log_factorial_[i - 1] + logi;
    }
  }
  CtermCache() = default;

  double getH(int n) const { return n_log_n_[n]; }
  double getLog(int n) const { return log_n_[n]; }
  double getLogC(int n, int level);
  double getLogChoose(int n, int k) const {
    if (k == n || k == 0) return 0;
    return log_factorial_[n] - log_factorial_[k] -
           log_factorial_[n - k];
  }
};

struct CompCache {
  std::shared_ptr<CtermCache> cterm;

  CompCache(int n_samples) : cterm(std::make_shared<CtermCache>(n_samples)) {}
  CompCache() = default;
};
}  // namespace computation_impl
using computation_impl::CompCache;
using computation_impl::CtermCache;
}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTATION_CACHE
