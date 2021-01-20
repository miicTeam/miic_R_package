#ifndef MIIC_COMPUTATION_CACHE
#define MIIC_COMPUTATION_CACHE

#ifdef _OPENMP
#include <omp.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <map>
#include <set>

#include "structure.h"

namespace miic {
namespace computation {

namespace detail {

using std::pair;
using std::set;
using std::vector;
using structure::Grid2d;
using structure::InfoBlock;

class CtermCache {
 public:
  CtermCache(int n_samples)
      : size_n_(1 + n_samples),
        log_n_(size_n_, 0),
        n_log_n_(size_n_, 0),
        log_factorial_(size_n_, 0),
        log_c_(n_samples, kLevelLimit, -1) {
    for (int i = 2; i < size_n_; ++i) {  // first two terms are zero
      double logi = log(static_cast<double>(i));
      log_n_[i] = logi;
      n_log_n_[i] = i * logi;
      log_factorial_[i] = log_factorial_[i - 1] + logi;
    }
    for (int n = 1; n < size_n_; ++n) {
      getLogC(n, 1);
      getLogC(n, 2);
    }
  }
  CtermCache() = default;

  double getH(int n) const { return n_log_n_[n]; }
  double getLog(int n) const { return log_n_[n]; }
  double getLogC(int n, int level);
  double getLogChoose(int n, int k) const {
    if (k == n || k == 0) return 0;
    return log_factorial_[n] - log_factorial_[k] - log_factorial_[n - k];
  }

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
};

struct MutualInfoKey {
  set<int> xy;
  set<int> ui;

  MutualInfoKey(int X, int Y, const vector<int>& ui)
      : xy({X, Y}), ui(begin(ui), end(ui)) {}

  bool operator<(const MutualInfoKey& other) const {
    if (xy == other.xy) {
      return ui < other.ui;
    }
    return xy < other.xy;
  }
};

// value: Shifted 3-point information I(X;Y;Z|ui) - k(X;Y;Z|ui)
struct Info3PointKey {
  set<int> xyz;
  set<int> ui;

  Info3PointKey(int X, int Y, int Z, const vector<int>& ui)
      : xyz({X, Y, Z}), ui(begin(ui), end(ui)) {}

  bool operator<(const Info3PointKey& other) const {
    if (xyz == other.xyz) {
      return ui < other.ui;
    }
    return xyz < other.xyz;
  }
};

// value: Contributing score R(X,Y;Z|ui)
struct ScoreKey {
  set<int> XY;
  int Z;
  set<int> ui;

  ScoreKey(int X, int Y, int Z, const vector<int>& ui)
      : XY({X, Y}), Z(Z), ui(begin(ui), end(ui)) {}

  bool operator<(const ScoreKey& other) const {
    if (XY == other.XY) {
      if (Z == other.Z) {
        return ui < other.ui;
      }
      return Z < other.Z;
    }
    return XY < other.XY;
  }
};

using MutualInfoMap = std::map<MutualInfoKey, InfoBlock>;
using Info3PointMap = std::map<Info3PointKey, double>;
using ScoreMap = std::map<ScoreKey, double>;
using EntropyMap = std::map<ScoreKey, double>;

class InfoScoreCache {
 public:
  InfoScoreCache() = default;

  pair<InfoBlock, bool> getMutualInfo(int X, int Y, const vector<int>& ui) {
    auto it = mi_map_.find(MutualInfoKey(X, Y, ui));
    bool found = it != mi_map_.end();
    return std::make_pair(found ? it->second : InfoBlock{0, 0, 0}, found);
  }

  void saveMutualInfo(int X, int Y, const vector<int>& ui, InfoBlock block) {
#ifdef _OPENMP
#pragma omp critical
#endif
    mi_map_.insert({MutualInfoKey(X, Y, ui), std::move(block)});
  }

  pair<double, bool> getInfo3Point(int X, int Y, int Z, const vector<int>& ui) {
    auto it = i3_map_.find(Info3PointKey(X, Y, Z, ui));
    bool found = it != i3_map_.end();
    return std::make_pair(found ? it->second : 0, found);
  }

  void saveInfo3Point(int X, int Y, int Z, const vector<int>& ui, double I3) {
#ifdef _OPENMP
#pragma omp critical
#endif
    i3_map_.insert({Info3PointKey(X, Y, Z, ui), I3});
  }

  pair<double, bool> getScore(int X, int Y, int Z, const vector<int>& ui) {
    auto it = score_map_.find(ScoreKey(X, Y, Z, ui));
    bool found = it != score_map_.end();
    return std::make_pair(
        found ? it->second : std::numeric_limits<double>::lowest(), found);
  }

  void saveScore(int X, int Y, int Z, const vector<int>& ui, double score) {
#ifdef _OPENMP
#pragma omp critical
#endif
    score_map_.insert({ScoreKey(X, Y, Z, ui), score});
  }

  pair<double, bool> getEntropy(int X, int Y, int Z) {
    auto it = entropy_map_.find(ScoreKey(X, Y, Z, vector<int>()));
    bool found = it != entropy_map_.end();
    return std::make_pair(found ? it->second : 0, found);
  }

  void saveEntropy(int X, int Y, int Z, double H) {
    // Already in critical block
    entropy_map_.insert({ScoreKey(X, Y, Z, vector<int>()), H});
  }


 private:
  MutualInfoMap mi_map_;
  Info3PointMap i3_map_;
  ScoreMap score_map_;
  EntropyMap entropy_map_;
};

struct CompCache {
  std::shared_ptr<CtermCache> cterm;
  std::shared_ptr<InfoScoreCache> info_score;

  CompCache(int n_samples)
      : cterm(std::make_shared<CtermCache>(n_samples)),
        info_score(std::make_shared<InfoScoreCache>()) {}
  CompCache() = default;
};
}  // namespace detail
using detail::CompCache;
using detail::CtermCache;
}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTATION_CACHE
