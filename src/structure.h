#ifndef MIIC_STRUCTURE_H_
#define MIIC_STRUCTURE_H_

#include <Rcpp.h>

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
  int z_name_idx = -1;  // Index of the last best contributor in zi_list
  double Rxyz_ui = 0;   // Score of the best contributor
  double Ixy_ui = 0;
  double cplx = 0;
  int Nxy_ui = -1;
  short int connected = 1;  // 1 or 0. An edge is by default connected.
  double Ixy = 0;           // Mutual information without conditioning
  double cplx_no_u = 0;     // Complexity without conditioning
  int Nxy = -1;             // Count of joint factors without NA

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

struct MemorySpace {
  int max_level;
  int** sample;
  int** sortedSample;
  int** Opt_sortedSample;
  int* orderSample;
  int* sampleKey;
  int* Nxyuiz;
  int* Nyuiz;
  int* Nuiz;
  int* Nz;
  int* Ny;
  int* Nxui;
  int* Nx;
  int** Nxuiz;
  int* bridge;
  double* Pxyuiz;
  // continuous data
  int* sample_is_not_NA;
  int* NAs_count;

  int** dataNumeric_red;
  int** dataNumericIdx_red;

  int* AllLevels_red;
  int* cnt_red;
  int* posArray_red;
};

struct ExecutionTime {
  double start_time_init;
  double start_time_iter;
  long double init;
  long double iter;
  long double init_iter;
  long double ort;
  long double cut;
  long double ort_after_cut;
  long double total;
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
