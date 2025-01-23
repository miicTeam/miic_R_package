#ifndef MIIC_STRUCTURE_H_
#define MIIC_STRUCTURE_H_

#include <functional>  // std::reference_wrapper
#include <limits>      // std::numeric_limits
#include <memory>      // std::shared_ptr
#include <string>
#include <type_traits>
#include <vector>

#include "linear_allocator.h"

namespace miic {
namespace structure {

namespace detail {

using std::size_t;
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
    size_t size() const { return parent_.n_cols(); }

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
    size_t size() const { return parent_.n_cols(); }

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

  Grid2d(size_t rows, size_t cols, vector<T, Allocator> data)
      : rows_(rows), cols_(cols), data_(std::move(data)) {}

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
  size_t size() const { return data_.size(); }
  bool empty() const { return data_.empty(); }

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

// Shifted conditional mutual information Nxy_ui * I(X;Y|ui) - k(X;Y|ui)
struct InfoBlock {
  double n_samples;  // double to include weighted case
  double I;
  double k;

  constexpr InfoBlock() : n_samples(0), I(0), k(0) {}
  constexpr InfoBlock(double N, double I, double k)
      : n_samples(N), I(I), k(k) {}
};

struct Info3PointBlock {
  double score;
  double Ixyz_ui;
  double kxyz_ui;
  // The two terms below are kept to facilitate the computation of contribution
  // of each 3-point information in terms of reduction in mutual information.
  // Ixyz_ui = Ixy_ui - Ixy_uiz
  double Ixy_ui;
  double kxy_ui;

  constexpr Info3PointBlock()
      : score(std::numeric_limits<double>::lowest()),
        Ixyz_ui(0),
        kxyz_ui(0),
        Ixy_ui(0),
        kxy_ui(0) {}
  constexpr Info3PointBlock(
      double R, double I3, double k3, double I2, double k2)
      : score(R), Ixyz_ui(I3), kxyz_ui(k3), Ixy_ui(I2), kxy_ui(k2) {}
};

struct EdgeSharedInfo {
  // {ui}: indices of separating nodes
  vector<int> ui_list;
  // The raw contribution of each ui to the conditional independence, measured
  // by I'(X;Y;ui|{uj}) / I'(X;Y)
  vector<double> raw_contributions;
  // The contribution of each ui to the reduction of conditional mutual
  // information, measured by I'(X;Y;ui|{uj}) / I'(X;Y|{uj})
  vector<double> contributions;
  // {zi}: indices of candidate conditioning nodes
  vector<int> zi_list;
  // Best candidate separating node
  int top_z = -1;
  // Score of the best contributor, this is the exponential part of the full
  // score as defined in Verny et al., 2017 (Supplementary Text)
  double Rxyz_ui = 0;
  // The raw contribution of top_z to the conditional independence, measured by
  // I'(X;Y;top_z|{ui}) / I'(X;Y)
  double top_raw_contribution = 0;
  // The contribution of top_z to the reduction of conditional mutual
  // information, measured by I'(X;Y;top_z|{ui}) / I'(X;Y|{ui})
  double top_contribution = 0;
  // Conditional mutual information
  double Ixy_ui = 0;
  // Complexity with conditioning
  double kxy_ui = 0;
  // Count of joint factors without NA
  int Nxy_ui = -1;
  // 1 or 0. An edge is by default connected.
  short int connected = 1;
  // Mutual information without conditioning
  double Ixy = 0;
  // Complexity without conditioning
  double kxy = 0;
  // Count of joint factors without NA
  int Nxy = -1;
  // if doing shuffling, exp(-I_shuffle)
  double exp_shuffle = -1;

  EdgeSharedInfo() = default;
  // Remove knowledge about all contributing nodes.
  void reset() {
    zi_list.clear();
    ui_list.clear();
    raw_contributions.clear();
    contributions.clear();
    top_z = -1;
    top_raw_contribution = 0;
    top_contribution = 0;
    Rxyz_ui = 0;
    Ixy_ui = Ixy;
    kxy_ui = kxy;
    Nxy_ui = Nxy;
    connected = 1;
  }

  void setUndirected() {
    ui_list.clear();
    top_z = -1;
    Rxyz_ui = 0;
    Ixy_ui = Ixy;
    kxy_ui = kxy;
    Nxy_ui = Nxy;
    connected = 1;
  }
};

struct Node {
  std::string name;
  Node(std::string name) : name(std::move(name)) {}
};

struct Edge {
  // For a pair of nodes (X, Y), status describes the arrow tip of the edge
  // (X *-* Y) at the Y side ('*' means either head '>' or tail '-')
  // Status code
  // 0: not connected;
  // 1: tail (X *-- Y);
  // 2: head (X *-> Y);
  short int status{1};       // Current status.
  short int status_init{1};  // Status after initialization.
  short int status_prev{1};  // Status in the previous iteration.
  double proba_head{0.5};    // Probability that the arrow tip is head (X *-> Y)
  std::shared_ptr<EdgeSharedInfo> shared_info;
};

struct interactEdge {
  // For a pair of nodes (X, Y), status describes the arrow tip of the edge
  // (X *-* Y) at the Y side ('*' means either head '>' or tail '-')
  // Status code
  // 0: not connected;
  // 1: connected;
  short int status{0};       // Current status
  double weight{0};          // Weight given by the Cell-Cell Communication (CCC) method
  double proba_head{0.5};    // Probability that the arrow tip is head (X *-> Y)
};

// Observer of Edge
class EdgeID {
 private:
  std::reference_wrapper<const Edge> edge_;

 public:
  int X, Y;
  EdgeID() = delete;
  EdgeID(int i, int j, const Edge& edge) : edge_(edge), X(i), Y(j) {}
  EdgeID(int i, int j, const Edge&&) = delete;  // forbid rvalue

  const Edge& getEdge() const { return edge_.get(); }

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

struct ExecutionTime {
  double init{0};  // skeletonInitialization
  double iter{0};  // firstStepIteration + skeletonIteration
  double cut{0};   // setConfidence + confidenceCut
  double ori{0};   // orientationProbability

  double getTotal() const { return init + iter + cut + ori; }
};

struct CutPointsInfo {
  Grid2d<int> cutpoints;
  double I{0};
  double Ik{0};
  double I_equal_freq_max{0};
  int n_iterations{0};

  CutPointsInfo() = default;
  CutPointsInfo(size_t n_rows, size_t n_cols) : cutpoints(n_rows, n_cols, -1) {}
};

}  // namespace detail
using detail::CutPointsInfo;
using detail::Edge;
using detail::interactEdge;
using detail::EdgeID;
using detail::EdgeSharedInfo;
using detail::ExecutionTime;
using detail::Grid2d;
using detail::Info3PointBlock;
using detail::InfoBlock;
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
