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

struct EdgeID {
  int i, j;
  EdgeID() = delete;
  EdgeID(int i, int j) : i(i), j(j) {}
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
using structure_impl::MemorySpace;
using structure_impl::Node;
}  // namespace structure
}  // namespace miic

#endif  // MIIC_STRUCTURE_H_
