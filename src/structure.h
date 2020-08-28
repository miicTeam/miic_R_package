#ifndef MIIC_STRUCTURE_H_
#define MIIC_STRUCTURE_H_

#include <Rcpp.h>

#include <array>
#include <memory>  // std::shared_ptr
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "linear_allocator.h"

namespace miic {
namespace structure {

namespace detail {

using std::size_t;
using std::string;
using std::vector;

// In absence of c++17 and to accomodate older compiler (ref: CWG 1558)
template <typename... Ts>
struct make_void {
  typedef void type;
};
template <typename... Ts>
using void_t = typename make_void<Ts...>::type;

// SFINAE classes
template <class T, class Index = size_t, class = void>
struct has_subscript_operator : std::false_type {};
template <class T>
struct has_subscript_operator<T, size_t,
    void_t<decltype(std::declval<T>()[std::declval<size_t>()])>>
    : std::true_type {};

template <class T, class Reduced = std::remove_reference_t<T>, class = void>
struct is_int_container : std::false_type {};
template <class T>
struct is_int_container<T, std::remove_reference_t<T>,
    void_t<std::enable_if_t<
               has_subscript_operator<std::remove_reference_t<T>>::value>,
        std::enable_if_t<std::is_same<
            typename std::remove_reference_t<T>::value_type, int>::value>>>
    : std::true_type {};

template <typename T, typename Allocator = std::allocator<T>>
struct Grid2d {
 public:
  typedef T value_type;
  // Wrapper class only for non-const instance
  struct Row {
   public:
    typedef T value_type;
    Row() = delete;
    Row(Grid2d& parent, size_t row) : parent_(parent), row_(row) {}

    T& operator()(size_t col) { return parent_(row_, col); }
    const T& operator()(size_t col) const { return parent_(row_, col); }
    T& operator[](size_t col) { return parent_(row_, col); }
    const T& operator[](size_t col) const { return parent_(row_, col); }
    size_t size() { return parent_.n_cols(); }

    auto begin() { return parent_.row_begin(row_); }
    auto end() { return parent_.row_end(row_); }
    auto begin() const { return parent_.row_begin(row_); }
    auto end() const { return parent_.row_end(row_); }

   private:
    Grid2d& parent_;
    size_t row_;
  };
  // Wrapper class for const instance
  struct ConstRow {
   public:
    typedef T value_type;
    ConstRow() = delete;
    ConstRow(const Grid2d& parent, size_t row) : parent_(parent), row_(row) {}

    const T& operator()(size_t col) const { return parent_(row_, col); }
    const T& operator[](size_t col) const { return parent_(row_, col); }
    size_t size() { return parent_.n_cols(); }

    auto begin() const { return parent_.row_begin(row_); }
    auto end() const { return parent_.row_end(row_); }

   private:
    const Grid2d& parent_;
    size_t row_;
  };

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
  Row getRow(size_t row) { return Row(*this, row); }
  ConstRow getConstRow(size_t row) const { return ConstRow(*this, row); }

  size_t n_rows() const { return rows_; }
  size_t n_cols() const { return cols_; }
  size_t size() { return data_.size(); }

  auto begin() { return data_.begin(); }
  auto end() { return data_.end(); }
  auto begin() const { return data_.cbegin(); }
  auto end() const { return data_.cend(); }

  auto row_begin(size_t row) { return data_.begin() + row * cols_; }
  auto row_end(size_t row) { return data_.begin() + (row + 1) * cols_; }
  auto row_begin(size_t row) const { return data_.cbegin() + row * cols_; }
  auto row_end(size_t row) const { return data_.cbegin() + (row + 1) * cols_; }

 private:
  size_t rows_, cols_;
  vector<T, Allocator> data_;
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

struct ExecutionTime {
  double init{0};  // skeletonInitialization
  double iter{0};  // firstStepIteration + skeletonIteration
  double cut{0};   // setConfidence + confidenceCut
  double ori{0};   // orientationProbability

  double getTotal() const { return init + iter + cut + ori; }
};

}  // namespace detail
using detail::CacheInfoKey;
using detail::CacheScoreValue;
using detail::Edge;
using detail::EdgeID;
using detail::EdgeSharedInfo;
using detail::ExecutionTime;
using detail::Grid2d;
using detail::Node;
using detail::void_t;

template <class T>
using IsIntContainer = std::enable_if_t<detail::is_int_container<T>::value>;
// types using linear allocator
using TempString = std::basic_string<char, std::char_traits<char>,
    utility::TempStdAllocator<char>>;

template <class T>
using TempVector = std::vector<T, utility::TempStdAllocator<T>>;

template <class T>
using TempGrid2d = Grid2d<T, utility::TempStdAllocator<T>>;
}  // namespace structure
}  // namespace miic

#endif  // MIIC_STRUCTURE_H_
