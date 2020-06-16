#ifndef MIIC_STRUCTURE_H_
#define MIIC_STRUCTURE_H_

#include <memory>
#include <string>
#include <vector>

namespace miic {
namespace structure {

namespace structure_impl {

using uint = unsigned int;
using std::string;
using std::vector;

struct EdgeSharedInfo {
  vector<int> ui_vect_idx;  // Indice of separating nodes
  // Indice of candidate nodes contributing to the conditional independence
  vector<int> zi_vect_idx;
  int z_name_idx = -1;  // Index of the last best contributor
  double Rxyz_ui = 0;   // Score of the best contributor
  double Ixy_ui = 0;
  double cplx = 0;
  int Nxy_ui = -1;
  short int connected = 1;  // 1 or 0. An edge is by default connceted.
  double mutInfo = 0;       // Mutual information without conditioning
  double cplx_noU = 0;      // Complexity without conditioning
  int Nxy = -1;             // Count of joint factors without NA

  EdgeSharedInfo() = default;
  // Remove knowledge about all contributing nodes.
  void reset() {
    zi_vect_idx.clear();
    ui_vect_idx.clear();
    z_name_idx = -1;
    Rxyz_ui = 0;
    Ixy_ui = mutInfo;
    cplx = cplx_noU;
    Nxy_ui = Nxy;
    connected = 1;
  }

  void setUndirected() {
    ui_vect_idx.clear();
    z_name_idx = -1;
    Rxyz_ui = 0;
    Ixy_ui = mutInfo;
    cplx = cplx_noU;
    Nxy_ui = Nxy;
    connected = 1;
  }
};

struct Node {
  string name;
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

struct MemorySpace {
  int maxlevel;
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
  double startTimeInit;
  double startTimeIter;
  long double init;
  long double iter;
  long double initIter;
  long double ort;
  long double cut;
  long double ort_after_cut;
  long double total;
};

// Structure for all the needed parameters (input plus state variables)
struct Environment {
  ExecutionTime execTime;
  int consistentPhase;
  // for gaussian case
  double** dataDouble;
  int** iterative_cuts;
  double* sampleWeights;
  vector<double> sampleWeightsVec;
  bool flag_sample_weights;

  bool testDistribution;
  uint nThreads;
  MemorySpace m;
  MemorySpace* memoryThreads;

  double* c2terms;
  double** cterms;
  double** lookchoose;
  int* columnAsContinuous;
  vector<int> cntVarVec;
  int* oneLineMatrix;

  string myVersion;  // -i parameter

  string outDir;    // -i parameter
  string inData;    // -o parameter
  string cplxType;  // -c parameter
  string blackbox_name;
  string edgeFile;
  string dataTypeFile;
  // string sampleWeightsFile;// -w parameter

  uint numNodes;
  uint numSamples;
  bool firstIterationDone;

  vector<EdgeID*> searchMoreAddress;
  vector<EdgeID*> noMoreAddress;
  int numSearchMore;
  int numNoMore;

  int typeOfData;
  int atLeastOneContinuous;
  int atLeastTwoDiscrete;

  // Keep a trace of the number of edges for every state
  int phantomEdgenNum;

  Node* nodes;
  Edge** edges;
  vector<string> vectorData;
  vector<vector<string> > data;
  int** dataNumeric;
  int** dataNumericIdx;
  uint* allLevels;

  double** proportions;

  int** cut;

  double logEta;

  bool myTest;
  bool isDegeneracy;             // -d parameter
  bool isVerbose;                // -v parameter
  bool isLatent;                 // -l parameter
  bool isLatentOnlyOrientation;  // -l parameter
  bool isNoInitEta;              // -f parameter
  bool propag;
  bool isK23;  // -k parameter
  bool isPropagation;

  int isTplReuse;  // -r parameter
  int cplx;
  int halfVStructures;

  int numberShuffles;          // -s parameter
  double confidenceThreshold;  // -e parameter

  int effN;  // -n parameter
  int minN;
  int thresPc;

  int maxbins;
  int initbins;
  double* looklog;
  double* lookH;

  double* noiseVec;
};

}  // namespace structure_impl
using structure_impl::Edge;
using structure_impl::EdgeID;
using structure_impl::EdgeSharedInfo;
using structure_impl::Environment;
using structure_impl::ExecutionTime;
using structure_impl::MemorySpace;
using structure_impl::Node;
}  // namespace structure
}  // namespace miic

#endif  // MIIC_STRUCTURE_H_
