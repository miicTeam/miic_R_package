#ifndef MIIC_STRUCTURE_H_
#define MIIC_STRUCTURE_H_

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>

namespace miic {
namespace structure {

namespace structure_impl {

using uint = unsigned int;
using std::string;
using std::vector;

struct EdgeSharedInfo {
  // Indices of separating nodes
  vector<int> ui_vect_idx;
  // Indices of candidate nodes contributing to the conditional independence
  vector<int> zi_vect_idx;
  int z_name_idx = -1;  // Index of the last best contributor
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
    zi_vect_idx.clear();
    ui_vect_idx.clear();
    z_name_idx = -1;
    Rxyz_ui = 0;
    Ixy_ui = Ixy;
    cplx = cplx_no_u;
    Nxy_ui = Nxy;
    connected = 1;
  }

  void setUndirected() {
    ui_vect_idx.clear();
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
  Node(string name) : name(name) {}
};

struct EdgeID {
  uint i, j;
  EdgeID() = delete;
  EdgeID(uint i, uint j) : i(i), j(j) {}
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

struct EdgeKey {
  int           x;
  int           y;
  std::set<int> Ui;
  
  bool operator<(const EdgeKey& other) const {
      if ( std::min(x, y) == std::min(other.x, other.y) ) {
        if ( std::max(x, y) == std::max(other.x, other.y) ) {
          return Ui < other.Ui;
        }
        return std::max(x, y) < std::max(other.x, other.y);
      }
      return std::min(x, y) < std::min(other.x, other.y);
  }
};

struct ScoreValue {
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

// Structure for all the needed parameters (input plus state variables)
struct Environment {
  ExecutionTime exec_time;
  // level of consistency required for the graph
  // 0: no consistency requirement
  // 1: skeleton consistent
  // 2: orientation consistent
  int consistent;
  // when consistent > 0, the maximum number of iterations allowed when trying
  // to find a consistent graph.
  int max_iteration;
  double** data_double;
  int** iterative_cuts;
  vector<double> sample_weights;
  bool flag_sample_weights;
  // whether or not do MAR (Missing at random) test using KL-divergence
  bool test_mar;
  uint n_threads;
  MemorySpace m;
  MemorySpace* memoryThreads;

  double* c2terms;
  double** cterms;
  double** lookchoose;
  vector<int> is_continuous;
  int* oneLineMatrix;

  uint n_nodes;
  uint n_samples;
  // if firstStepIteration is done
  bool first_iter_done = false;
  // List of ids of edge whose status is not yet determined
  vector<EdgeID*> unsettled_list;
  // List of ids of edge whose status is sure to be connected
  vector<EdgeID*> connected_list;
  int numSearchMore;
  int numNoMore;

  vector<Node> nodes;
  Edge** edges;
  vector<vector<string>> data;
  int** data_numeric;
  int** data_numeric_idx;
  uint* levels;

  double log_eta = 0;

  bool degenerate;
  bool verbose;
  bool latent;
  bool latent_orientation;
  bool no_init_eta = false;
  bool is_k23;
  bool propagation;

  int cplx;
  int half_v_structure;

  int n_shuffles;
  double conf_threshold;

  int n_eff;
  int thresPc;

  int maxbins;
  int initbins;
  double* looklog;
  double* lookH;
  std::map<EdgeKey, double> look_scores;
  std::map<EdgeKey, ScoreValue> look_scores_orientation;

  double* noise_vec;
};

}  // namespace structure_impl
using structure_impl::Edge;
using structure_impl::EdgeID;
using structure_impl::EdgeKey;
using structure_impl::EdgeSharedInfo;
using structure_impl::ScoreValue;
using structure_impl::Environment;
using structure_impl::ExecutionTime;
using structure_impl::MemorySpace;
using structure_impl::Node;
}  // namespace structure
}  // namespace miic

#endif  // MIIC_STRUCTURE_H_
