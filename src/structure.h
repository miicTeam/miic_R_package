#ifndef MIIC_STRUCTURE_H_
#define MIIC_STRUCTURE_H_

#include <Rcpp.h>

#include <array>
#include <memory>  // std::shared_ptr
#include <set>
#include <string>
#include <vector>

namespace miic {
namespace structure {

namespace structure_impl {

using std::string;
using std::vector;

template <typename T>
struct Grid2d {
 private:
  size_t rows_, cols_;
  vector<T> data_;

 public:
  Grid2d() = default;
  Grid2d(size_t rows, size_t cols)
      : rows_(rows), cols_(cols), data_(rows * cols) {}

  Grid2d(size_t rows, size_t cols, T&& init)
      : rows_(rows), cols_(cols), data_(rows * cols, init) {}

  Grid2d(const Grid2d&) = default;
  Grid2d(Grid2d&&) = default;
  Grid2d& operator=(const Grid2d&) = default;
  Grid2d& operator=(Grid2d&&) = default;

  T& operator()(size_t row, size_t col) { return data_[row * cols_ + col]; }
  const T& operator()(size_t row, size_t col) const {
    return data_[row * cols_ + col];
  }

  size_t n_rows() { return rows_; }
  size_t n_cols() { return cols_; }
  size_t size() { return data_.size(); }

  auto begin() { return data_.begin(); }
  auto end() { return data_.end(); }
  auto cbegin() const { return data_.cbegin(); }
  auto cend() const { return data_.cend(); }

  auto row_begin(size_t row) { return data_.begin() + row * cols_; }
  auto row_end(size_t row) { return data_.end() + (row + 1) * cols_; }
  auto row_cbegin(size_t row) const { return data_.cbegin() + row * cols_; }
  auto row_cend(size_t row) const { return data_.cend() + (row + 1) * cols_; }
};

struct EdgeSharedInfo {
  // {ui}: indices of separating nodes
  vector<int> ui_list;
  // {zi}: indices of candidate conditioning nodes
  vector<int> zi_list;
  // Index of the last best contributor in zi_list
  int z_name_idx = -1;
  // Score of the best contributor
  double Rxyz_ui = 0;
  // Conditional mutual information
  double Ixy_ui = 0;
  // Complexity with conditioning
  double cplx = 0;
  // Count of joint factors without NA
  int Nxy_ui = -1;
  // 1 or 0. An edge is by default connected.
  short int connected = 1;
  // Mutual information without conditioning
  double Ixy = 0;
  // Complexity without conditioning
  double cplx_no_u = 0;
  // Count of joint factors without NA
  int Nxy = -1;
  // if doing shuffling, exp(-I_shuffle)
  double exp_shuffle = -1;

  EdgeSharedInfo() = default;
  // Remove knowledge about all contributing nodes.
  void reset() {
    zi_list.clear();
    ui_list.clear();
    z_name_idx = -1;
    Rxyz_ui = 0;
    Ixy_ui = Ixy;
    cplx = cplx_no_u;
    Nxy_ui = Nxy;
    connected = 1;
  }

  void setUndirected() {
    ui_list.clear();
    z_name_idx = -1;
    Rxyz_ui = 0;
    Ixy_ui = Ixy;
    cplx = cplx_no_u;
    Nxy_ui = Nxy;
    connected = 1;
  }
};

struct Node {
  string name;
  Node(string name) : name(std::move(name)) {}
};

struct Edge {
  // Edge is stored in Edge** edges
  // Status code (suppose edges[X][Y]):
  // 0: not connected;
  // 1: connected and undirected;
  // 2: connected directed X -> Y;
  // -2: connected directed X <- Y;
  // 6: connected bidirected X <-> Y;
  short int status;       // Current status.
  short int status_init;  // Status after initialization.
  short int status_prev;  // Status in the previous iteration.
  std::shared_ptr<EdgeSharedInfo> shared_info;
};

// Observer of Edge
class EdgeID {
 private:
  std::reference_wrapper<const Edge> edge_;

 public:
  int i, j;
  EdgeID() = delete;
  EdgeID(int i, int j, const Edge& edge) : edge_(edge), i(i), j(j) {}
  EdgeID(int i, int j, const Edge&&) = delete;  // forbid rvalue

  bool operator<(const EdgeID& rhs) const {
    const auto info1 = this->edge_.get().shared_info;
    const auto info2 = rhs.edge_.get().shared_info;
    //  connected can be 0 or 1, prefer connected over non-connected
    if (info1->connected != info2->connected)
      return info1->connected > info2->connected;
    if (info1->connected) {
      return info1->Ixy_ui > info2->Ixy_ui;
    } else {
      return info1->Rxyz_ui > info2->Rxyz_ui;
    }
  }
};

struct CacheInfoKey{
  std::set<int> xyz;
  std::set<int> Ui;

  //Two point information constructor : I(X;Y|Ui,Z)
  CacheInfoKey(int x, int y, const std::set<int>& Ui_) {
    xyz.insert({x,y});
    Ui = Ui_;
  }
  //Three point information constructor : I(X;Y;Z|Ui)
  CacheInfoKey(int x, int y, int z, const std::set<int>& Ui_) {
    xyz.insert({x,y,z});
    Ui = Ui_;
  }

  bool operator<(const CacheInfoKey& other) const {
    if (xyz == other.xyz) {
      return Ui < other.Ui;
    }
    return xyz < other.xyz;
  }
};

struct CacheScoreValue {
  int    n_samples;
  double I_xyzUi;
  double cplx;
};

class MemorySpace {
 public:
  Grid2d<int> sample;
  Grid2d<int> sortedSample;
  Grid2d<int> Opt_sortedSample;
  vector<int> orderSample;
  vector<int> sampleKey;
  vector<int> Nyuiz;
  vector<int> Nuiz;
  vector<int> Nz;
  vector<int> Ny;
  vector<int> Nxui;
  vector<int> Nx;
  Grid2d<int> Nxuiz;
  vector<int> bridge;
  vector<double> Pxyuiz;

 public:
  MemorySpace(int n_nodes, int n_samples, int max_level)
      : sample(n_samples + 1, 7),
        sortedSample(n_samples + 1, 7),
        Opt_sortedSample(n_samples + 1, 7),
        orderSample(n_samples + 2),
        sampleKey(n_samples + 2),
        Nyuiz(max_level + 1),
        Nuiz(max_level + 1),
        Nz(max_level + 1),
        Ny(max_level + 1),
        Nxui(max_level + 1),
        Nx(max_level + 1),
        Nxuiz(max_level + 1, max_level + 1),
        bridge(n_samples + 2),
        Pxyuiz(max_level + 1) {}

  MemorySpace() = default;
  MemorySpace(const MemorySpace&) = default;
  MemorySpace(MemorySpace&&) = default;
  MemorySpace& operator=(const MemorySpace&) = default;
  MemorySpace& operator=(MemorySpace&&) = default;
};

struct ExecutionTime {
  double init{0};  // skeletonInitialization
  double iter{0};  // firstStepIteration + skeletonIteration
  double cut{0};   // setConfidence + confidenceCut
  double ori{0};   // orientationProbability

  double getTotal() const { return init + iter + cut + ori; }
};

}  // namespace structure_impl
using structure_impl::CacheInfoKey;
using structure_impl::CacheScoreValue;
using structure_impl::Edge;
using structure_impl::EdgeID;
using structure_impl::EdgeSharedInfo;
using structure_impl::ExecutionTime;
using structure_impl::Grid2d;
using structure_impl::MemorySpace;
using structure_impl::Node;
}  // namespace structure
}  // namespace miic

#endif  // MIIC_STRUCTURE_H_
