#include "utilities.h"

#include <Rcpp.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "compute_info.h"
#include "nanoflann.hpp"

// for memory space on continuous data
#define MAX_NBRUI 10
#define N_COL_NML 1000
#define MAGNITUDE_TIES 0.00005f
#define KNN_K 5

namespace miic {
namespace utility {

using uint = unsigned int;
using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;
using namespace miic::computation;
using namespace miic::structure;

double kl(double** freqs1, double** freqs2, int nrows, int ncols) {
  double** lr;
  double kl_div = 0;

  lr = new double*[nrows];
  for (int i = 0; i < nrows; i++) lr[i] = new double[ncols];

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (freqs1[i][j] != 0) {
        lr[i][j] = log(freqs1[i][j] / freqs2[i][j]);
      } else
        lr[i][j] = 0;
    }
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      kl_div += freqs1[i][j] * lr[i][j];
    }
  }

  for (int i = 0; i < nrows; i++) delete[] lr[i];
  delete[] lr;

  return kl_div;
}

double kl(int** counts1, double** freqs2, int nrows, int ncols) {
  double** freqs1 = new double*[nrows];
  int n2 = 0;
  for (int i = 0; i < nrows; i++) {
    freqs1[i] = new double[ncols];
    for (int j = 0; j < ncols; j++) {
      n2 += counts1[i][j];
    }
  }
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      freqs1[i][j] = 1.0 * counts1[i][j] / n2;
    }
  }

  double kl_div = kl(freqs1, freqs2, nrows, ncols);

  for (int i = 0; i < nrows; i++) delete[] freqs1[i];
  delete[] freqs1;

  return (kl_div);
}

double kl(double* freqs1, double* freqs2, int nlevels) {
  double kl_div = 0;
  for (int level = 0; level < nlevels; level++) {
    if (freqs1[level] != 0) {
      kl_div += freqs1[level] * log(freqs1[level] / freqs2[level]);
    }
  }
  return (kl_div);
}

double kl(int* count1, int* count2, int n1, int n2, int nlevels) {
  double* freqs1 = new double[nlevels];
  double* freqs2 = new double[nlevels];

  for (int level = 0; level < nlevels; level++) {
    freqs1[level] = 1.0 * count1[level] / n1;
    freqs2[level] = 1.0 * count2[level] / n2;
  }
  double kl_div = kl(freqs1, freqs2, nlevels);
  delete[] freqs1;
  delete[] freqs2;

  return (kl_div);
}

bool allVariablesDiscrete(int* array, int* posArray, int num) {
  for (int i = 0; i < num; i++) {
    if (array[posArray[i]] == 1) return false;
  }
  return true;
}

class sort_indices {
 private:
  int* mparr;

 public:
  sort_indices(int* parr) : mparr(parr) {}
  bool operator()(int i, int j) const { return mparr[i] < mparr[j]; }
};

void sort2arrays(int len, int a[], int brr[], int bridge[]) {
  int i;

  int* pArray = &a[1];
  int* pArray2 = &brr[1];

  std::sort(pArray2, pArray2 + len, sort_indices(pArray));

  for (i = 1; i < len + 1; i++) {
    bridge[i] = pArray2[i - 1];
  }

  brr = bridge;
  brr[0] = 0;
}

double ramanujan(int n) {
  // Returns log(fac(n)) with Ramanujan's approximation.
  if (n == 0) {
    return (0);
  }
  double N = n * log(n) - n + log(1.0 * n * (1 + 4 * n * (1 + 2 * n))) / 6 +
             log(M_PI) / 2L;
  return (N);
}

double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    // Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void sort2arraysConfidence(int len, int a[], int brr[]) {
  std::sort(brr, brr + len, sort_indices(a));
}

vector<vector<int>> getAdjMatrix(const Environment& env) {
  vector<vector<int>> adj_matrix(env.numNodes, vector<int>(env.numNodes, 0));
  for (uint i = 1; i < env.numNodes; i++) {
    for (uint j = 0; j < i; j++) {
      adj_matrix[i][j] = env.edges[i][j].status;
      adj_matrix[j][i] = env.edges[j][i].status;
    }
  }
  return adj_matrix;
}

void createMemorySpace(Environment& environment, MemorySpace& m) {
  if (environment.atLeastTwoDiscrete) {
    uint maxLevel = 0;
    for (uint i = 0; i < environment.numNodes; i++) {
      //if (environment.columnAsContinuous[i] == 0 &&
      if (environment.allLevels[i] > maxLevel)
        maxLevel = environment.allLevels[i];
    }
    m.maxlevel = maxLevel;
    // cout<< "samples" << environment.numSamples << endl;
    int nrow = environment.numSamples + 1;
    int sampleSize = environment.numSamples;
    int ncol = 7;
    int bin_max = maxLevel;
    int iii;

    m.sample = (int**)calloc(nrow, sizeof(int*));
    for (iii = 0; iii < nrow; iii++)
      m.sample[iii] = (int*)calloc(ncol, sizeof(int));

    m.sortedSample = (int**)calloc(nrow, sizeof(int*));
    for (iii = 0; iii < nrow; iii++)
      m.sortedSample[iii] = (int*)calloc(ncol, sizeof(int));

    m.Opt_sortedSample = (int**)calloc(nrow, sizeof(int*));
    for (iii = 0; iii < nrow; iii++)
      m.Opt_sortedSample[iii] = (int*)calloc(ncol, sizeof(int));

    m.Nxuiz = (int**)calloc(bin_max + 1, sizeof(int*));
    for (iii = 0; iii < bin_max + 1; iii++)
      m.Nxuiz[iii] = (int*)calloc(bin_max + 1, sizeof(int));

    m.orderSample = (int*)calloc((sampleSize + 2), sizeof(int));
    m.sampleKey = (int*)calloc((sampleSize + 2), sizeof(int));
    m.Pxyuiz = (double*)calloc((bin_max + 1), sizeof(double));
    m.Nyuiz = (int*)calloc((bin_max + 1), sizeof(int));
    m.Nuiz = (int*)calloc((bin_max + 1), sizeof(int));
    m.Nz = (int*)calloc((bin_max + 1), sizeof(int));
    m.Ny = (int*)calloc((bin_max + 1), sizeof(int));
    m.Nxui = (int*)calloc((bin_max + 1), sizeof(int));
    m.Nx = (int*)calloc((bin_max + 1), sizeof(int));
    m.bridge = (int*)calloc(sampleSize + 2, sizeof(int));
  }
  // continuous part
  if (environment.atLeastOneContinuous) {
    m.sample_is_not_NA = (int*)new int[environment.numSamples];
    m.NAs_count = (int*)new int[environment.numSamples];
    m.dataNumericIdx_red = (int**)new int*[(MAX_NBRUI + 3)];
    m.dataNumeric_red = (int**)new int*[(MAX_NBRUI + 3)];

    for (int j = 0; (j < MAX_NBRUI + 3); j++) {
      m.dataNumericIdx_red[j] = (int*)new int[environment.numSamples];
      m.dataNumeric_red[j] = (int*)new int[environment.numSamples];
    }
    m.AllLevels_red = (int*)new int[(MAX_NBRUI + 3)];
    m.cnt_red = (int*)new int[(MAX_NBRUI + 3)];
    m.posArray_red = (int*)new int[(MAX_NBRUI + 3)];
  }
}

void deleteMemorySpace(Environment& environment, MemorySpace& m) {
  if (environment.atLeastTwoDiscrete) {
    uint maxLevel = 0;
    for (uint i = 0; i < environment.numNodes; i++) {
      if (environment.columnAsContinuous[i] == 0 &&
          environment.allLevels[i] > maxLevel)
        maxLevel = environment.allLevels[i];
    }

    int nrow = environment.numSamples + 1;
    int bin_max = maxLevel;
    int i;
    for (i = 0; i < nrow; i++) free(m.sample[i]);
    free(m.sample);

    for (i = 0; i < nrow; i++) free(m.sortedSample[i]);
    free(m.sortedSample);

    for (i = 0; i < bin_max + 1; i++) free(m.Nxuiz[i]);

    for (i = 0; i < bin_max + 1; i++) free(m.Opt_sortedSample[i]);

    free(m.Opt_sortedSample);
    free(m.Nxuiz);
    free(m.orderSample);
    free(m.sampleKey);
    free(m.Pxyuiz);
    free(m.Nyuiz);
    free(m.Nuiz);
    free(m.Nz);
    free(m.Ny);
    free(m.Nxui);
    free(m.Nx);
  }
  // cotinuous part
  if (environment.atLeastOneContinuous) {
    delete[] m.sample_is_not_NA;
    delete[] m.NAs_count;
    delete[] m.AllLevels_red;
    delete[] m.cnt_red;
    delete[] m.posArray_red;
    for (int i = 0; i < MAX_NBRUI + 3; i++) {
      delete[] m.dataNumericIdx_red[i];
      delete[] m.dataNumeric_red[i];
    }
    delete[] m.dataNumericIdx_red;
    delete[] m.dataNumeric_red;
  }
}

void deleteStruct(Environment& environment) {
  for (auto& address : environment.noMoreAddress) delete address;
  delete[] environment.oneLineMatrix;
  delete[] environment.allLevels;
  for (uint i = 0; i < environment.numSamples; i++)
    delete[] environment.dataNumeric[i];
  delete[] environment.dataNumeric;
  delete[] environment.c2terms;
  delete[] environment.nodes;
  for (uint i = 0; i < environment.numNodes; i++) delete[] environment.edges[i];
  delete[] environment.edges;
}

string toNameString(const Environment& env, const vector<int>& vec) {
  if (vec.empty()) {
    return "NA";
  } else {
    stringstream ss;
    std::transform(vec.begin(), vec.end() - 1,
        std::ostream_iterator<std::string>(ss, ","),
        [&env](int i) { return env.nodes[i].name; });
    ss << env.nodes[vec.back()].name;
    return ss.str();
  }
}

bool readBlackbox(vector<string> v, Environment& environment) {
  string s1, s2;
  int posX, posY;
  for (uint pos = 0; pos < v.size(); pos++) {
    posX = -1;
    posY = -1;

    s1 = v[pos];
    for (uint i = 0; i < environment.numNodes; i++) {
      if (environment.nodes[i].name.compare(s1) == 0) posX = i;
    }

    pos++;
    s2 = v[pos];
    for (uint i = 0; i < environment.numNodes; i++) {
      if (environment.nodes[i].name.compare(s2) == 0) posY = i;
    }

    if (posX != -1 && posY != -1) {
      environment.edges[posX][posY].status = 0;
      environment.edges[posY][posX].status = 0;
    }
  }

  return true;
}

vector<vector<string>> getEdgesInfoTable(Environment& env) {
  vector<vector<string>> table;

  vector<EdgeID> edge_list;
  for (uint i = 0; i < env.numNodes - 1; i++) {
    for (uint j = i + 1; j < env.numNodes; j++) {
      edge_list.emplace_back(EdgeID(i, j));
    }
  }
  std::sort(edge_list.begin(), edge_list.end(), EdgeSorter(env));

  table.emplace_back(std::initializer_list<std::string>{"x", "y", "z.name",
      "ai.vect", "zi.vect", "Ixy", "Ixy_ai", "cplx", "Rxyz_ai", "category",
      "Nxy_ai"});
  for (const auto& edge : edge_list) {
    auto i = edge.i, j = edge.j;
    auto info = env.edges[i][j].shared_info;

    table.emplace_back(std::initializer_list<std::string>{
        env.nodes[i].name,
        env.nodes[j].name,
        info->z_name_idx == -1
          ? "NA"
          : env.nodes[info->zi_vect_idx[info->z_name_idx]].name,
        toNameString(env, info->ui_vect_idx),
        toNameString(env, info->zi_vect_idx),
        std::to_string(info->mutInfo),
        std::to_string(info->Ixy_ui),
        std::to_string(info->cplx),
        std::to_string(info->Rxyz_ui),
        std::to_string(info->connected),
        std::to_string(info->Nxy_ui)
    });
  }

  return table;
}

// Initialize all the elements of the array to the given value
bool setArrayValuesInt(int* array, int length, int value) {
  for (int i = 0; i < length; i++) {
    array[i] = value;
  }
  return true;
}

bool readData(Environment& environment, bool& isNA) {
  vector<string> vec;
  environment.nodes = new Node[environment.numNodes];

  // convert input data
  for (uint i = 0; i < environment.vectorData.size(); i++) {
    if (i < environment.numNodes) {
      environment.nodes[i].name = environment.vectorData[i];
    } else {
      if (i % environment.numNodes == 0) {
        if (i != environment.numNodes) {
          environment.data.push_back(vec);
          vec.clear();
        }
      }
      vec.push_back(environment.vectorData[i]);
    }
  }

  environment.data.push_back(vec);

  environment.numSamples = environment.data.size();

  // set effN if not set before to the number of rows from the input
  // environment.data
  if (environment.effN == -1) environment.effN = environment.numSamples;

  // set maxbin coarse
  environment.maxbins = 50;

  return true;
}

bool removeRowsAllNA(Environment& environment) {
  if (environment.isVerbose) cout << "# Removing NA rows\n";

  int* indexNA = new int[environment.numSamples];
  setArrayValuesInt(indexNA, environment.numSamples, -1);

  int pos = 0;
  for (uint i = 0; i < environment.numSamples; i++) {
    bool isNA = true;
    for (uint j = 0; j < environment.numNodes && isNA; j++) {
      if ((environment.data[i][j].compare("NA") != 0) &&
          (environment.data[i][j].compare("") != 0)) {
        isNA = false;
      }
    }
    if (!isNA) {
      indexNA[pos] = i;
      pos++;
    }
  }

  if (environment.isVerbose)
    cout << "-------->Number of NA rows: " << pos << "\n";
  cout << "environment.numSamples :" << environment.numSamples << endl;

  // if there are rows of NA value
  if (pos != 0) {
    // correct variable numSamples

    // save the values
    int pos = 0;
    for (uint i = 0; i < environment.numSamples; i++) {
      if (indexNA[i] != -1) {
        for (uint j = 0; j < environment.numNodes; j++) {
          environment.data[pos][j] = environment.data[indexNA[i]][j];
        }
        pos++;
      }
    }
    environment.numSamples = pos;
    if (environment.effN > pos) {
      environment.effN = pos;
    }
  }
  // cout << "environment.numSamples :" << environment.numSamples << endl;
  return true;
}

void transformToFactors(Environment& environment, int i) {
  if (environment.isVerbose) cout << "# Transforming matrix to factors\n";

  // create a dictionary to store the factors of the strings
  std::map<string, int> myMap;

  // clean the dictionary since it is used column by column
  myMap.clear();
  myMap["NA"] = -1;
  myMap[""] = -1;
  int factor = 0;

  for (uint j = 0; j < environment.numSamples; j++) {
    std::map<string, int>::iterator it = myMap.find(environment.data[j][i]);
    if (it != myMap.end()) {
      environment.dataNumeric[j][i] = it->second;
    } else {
      myMap[environment.data[j][i]] = factor;
      environment.dataNumeric[j][i] = factor;
      factor++;
    }
  }
}

void transformToFactorsContinuous(Environment& environment, int i) {
  if (environment.isVerbose)
    cout << "# Transforming matrix to factors continuous\n";

  std::multimap<double, int> myMap;

  // clean the dictionary since it is used column by column
  myMap.clear();

  vector<double> clmn;
  for (uint j = 0; j < environment.numSamples; j++) {
    string entry = environment.data[j][i];
    if (entry.compare("NA") != 0 && entry.compare("") != 0)
      clmn.push_back(atof(entry.c_str()));
  }

  sort(clmn.begin(), clmn.end());

  for (uint j = 0; j < clmn.size(); j++) {
    myMap.insert(pair<double, int>(clmn[j], j));
  }

  for (uint j = 0; j < environment.numSamples; j++) {
    string entry = environment.data[j][i];
    if (entry.compare("NA") != 0 && entry.compare("") != 0) {
      environment.dataNumeric[j][i] = myMap.find(atof(entry.c_str()))->second;

      typedef std::multimap<double, int>::iterator iterator;
      std::pair<iterator, iterator> iterpair =
          myMap.equal_range(atof(entry.c_str()));
      iterator it = iterpair.first;
      for (; it != iterpair.second; ++it) {
        if (it->second == environment.dataNumeric[j][i]) {
          myMap.erase(it);
          break;
        }
      }
    } else {
      environment.dataNumeric[j][i] = -1;
    }
  }
}

void transformToFactorsContinuousIdx(Environment& environment, int i) {
  if (environment.isVerbose)
    cout << "# Transforming matrix to factors continuous\n";

  std::map<int, int> myMap;

  // clean the dictionary since it is used column by column
  myMap.clear();

  int entry;
  for (uint j = 0; j < environment.numSamples; j++) {
    entry = environment.dataNumeric[j][i];
    if (entry != -1) myMap[entry] = j;
  }

  int j = 0;
  for (std::map<int, int>::iterator it = myMap.begin(); it != myMap.end();
       ++it) {
    environment.dataNumericIdx[i][j] = it->second;
    j++;
  }
}

// Set the number of levels for each node (the maximum level of each column)
void setNumberLevels(Environment& environment) {
  environment.allLevels = new uint[environment.numNodes];
  int max;
  for (uint i = 0; i < environment.numNodes; i++) {
    max = 0;
    for (uint j = 0; j < environment.numSamples; j++) {
      if (environment.dataNumeric[j][i] > max)
        max = environment.dataNumeric[j][i];
    }
    environment.allLevels[i] = max + 1;
  }
}

int getNumSamples_nonNA(Environment& environment, int i, int j) {
  bool sampleOk;
  int numSamples_nonNA = 0;
  for (uint k = 0; k < environment.numSamples; k++) {
    sampleOk = true;
    if (environment.dataNumeric[k][i] != -1 &&
        environment.dataNumeric[k][j] != -1) {
      for (uint u = 0;
           u < environment.edges[i][j].shared_info->ui_vect_idx.size(); u++) {
        if (environment.dataNumeric[k][environment.edges[i][j]
                                           .shared_info->ui_vect_idx[u]] ==
            -1) {
          sampleOk = false;
        }
      }
      if (sampleOk) {
        numSamples_nonNA++;
      }
    }
  }
  return (numSamples_nonNA);
}

void getJointSpace(Environment& environment, int i, int j,
    vector<vector<double> >& jointSpace, int* curr_sample_is_not_NA) {
  int numSamples_nonNA = 0;
  bool sampleOk;
  for (uint k = 0; k < environment.numSamples; k++) {
    sampleOk = true;
    curr_sample_is_not_NA[k] = 0;
    if (environment.dataNumeric[k][i] != -1 &&
        environment.dataNumeric[k][j] != -1) {
      for (uint u = 0;
           u < environment.edges[i][j].shared_info->ui_vect_idx.size(); u++) {
        if (environment.dataNumeric[k][environment.edges[i][j]
                                           .shared_info->ui_vect_idx[u]] ==
            -1) {
          sampleOk = false;
        }
      }
      if (sampleOk) {
        curr_sample_is_not_NA[k] = 1;
        jointSpace[numSamples_nonNA][0] = environment.dataDouble[k][i];
        jointSpace[numSamples_nonNA][1] = environment.dataDouble[k][j];
        numSamples_nonNA++;
      }
    }
  }
}

double** getJointFreqs(
    Environment& environment, int i, int j,
    const vector<int> &sample_is_not_NA) {
  double** jointFreqs = new double*[environment.allLevels[i]];
  for (uint k = 0; k < environment.allLevels[i]; k++) {
    jointFreqs[k] = new double[environment.allLevels[j]];
    for (uint l = 0; l < environment.allLevels[j]; l++) {
      jointFreqs[k][l] = 0;
    }
  }

  int numSamples_nonNA(0);
  bool sampleOk;
  for (uint k = 0; k < environment.numSamples; k++) {
    sampleOk = true;
    if(sample_is_not_NA.empty()){
      // Check if any value is NA in X,Y or the Us
      if (environment.dataNumeric[k][i] != -1 &&
          environment.dataNumeric[k][j] != -1) {
        for (uint u = 0;
            u < environment.edges[i][j].shared_info->ui_vect_idx.size(); u++) {
          if (environment.dataNumeric[k][environment.edges[i][j]
                                            .shared_info->ui_vect_idx[u]] == -1)
            sampleOk = false;
        }
      }
      else sampleOk = false;
    }
    else {
      // Use the sample_is_not_NA vector
      sampleOk = sample_is_not_NA[k];
    }
    if (sampleOk) {
      jointFreqs[environment.dataNumeric[k][i]]
                [environment.dataNumeric[k][j]]++;
      numSamples_nonNA++;
    }
  }

  for (uint k = 0; k < environment.allLevels[i]; k++)
    for (uint l = 0; l < environment.allLevels[j]; l++)
      jointFreqs[k][l] /= numSamples_nonNA;

  return (jointFreqs);
}

void getJointMixed(Environment& environment, int i, int j, int* mixedDiscrete,
    double* mixedContinuous, int* curr_sample_is_not_NA) {
  int discrete_pos = environment.columnAsContinuous[i] == 0 ? i : j;
  int continuous_pos = environment.columnAsContinuous[i] == 1 ? i : j;

  // Fill marginal distributions
  int numSamples_nonNA = 0;
  bool sampleOk;
  for (uint k = 0; k < environment.numSamples; k++) {
    sampleOk = true;
    curr_sample_is_not_NA[k] = 0;
    if (environment.dataNumeric[k][i] != -1 &&
        environment.dataNumeric[k][j] != -1) {
      for (uint u = 0;
           u < environment.edges[i][j].shared_info->ui_vect_idx.size(); u++) {
        if (environment.dataNumeric[k][environment.edges[i][j]
                                           .shared_info->ui_vect_idx[u]] ==
            -1) {
          sampleOk = false;
        }
      }
      if (sampleOk) {
        curr_sample_is_not_NA[k] = 1;
        mixedContinuous[numSamples_nonNA] =
            environment.dataDouble[k][continuous_pos];
        mixedDiscrete[numSamples_nonNA] =
            environment.dataNumeric[k][discrete_pos];
        numSamples_nonNA++;
      }
    }
  }
}

void readFileType(Environment& environment) {
  for (uint pos = 0; pos < environment.numNodes; pos++) {
    environment.columnAsContinuous[pos] = environment.cntVarVec[pos];
  }
  if (environment.typeOfData != 0) {
    environment.dataDouble = new double*[environment.numSamples];

    for (uint i = 0; i < environment.numSamples; i++) {
      environment.dataDouble[i] = new double[environment.numNodes];
      for (uint j = 0; j < environment.numNodes; j++) {
        if ((environment.columnAsContinuous[j] == 1) ||
            (environment.columnAsContinuous[j] == 1)) {
          if ((environment.data[i][j].compare("NA") == 0) ||
              (environment.data[i][j].compare("") == 0)) {
            environment.dataDouble[i][j] =
                std::numeric_limits<double>::quiet_NaN();
          } else {
            environment.dataDouble[i][j] = atof(environment.data[i][j].c_str());
          }
        }
      }
    }
  }
}

void setEnvironment(Environment& environment) {
  environment.noMoreAddress.clear();
  environment.numNoMore = 0;
  environment.searchMoreAddress.clear();
  environment.numSearchMore = 0;

  environment.atLeastTwoDiscrete = 0;
  environment.atLeastOneContinuous = 0;

  bool isNA = false;

  readData(environment, isNA);

  if (isNA) removeRowsAllNA(environment);

  environment.columnAsContinuous = new int[environment.numNodes];
  if (environment.typeOfData == 0) {
    for (uint i = 0; i < environment.numNodes; i++) {
      environment.columnAsContinuous[i] = 0;
    }
  } else {
    readFileType(environment);
  }

  // Set the effN if not already done
  if (environment.effN == -1 || environment.effN > (int)environment.numSamples)
    environment.effN = environment.numSamples;

  environment.sampleWeights = new double[environment.numSamples];
  if (environment.sampleWeightsVec[0] != -1) {
    for (uint i = 0; i < environment.numSamples; i++) {
      environment.sampleWeights[i] = environment.sampleWeightsVec[i];
    }
  } else {
    for (uint i = 0; i < environment.numSamples; i++) {
      if (environment.effN == (int)environment.numSamples)
        environment.sampleWeights[i] = 1;
      else
        environment.sampleWeights[i] =
            (environment.effN * 1.0) / environment.numSamples;
    }
  }

  int count = 0;
  for (uint i = 0; i < environment.numNodes; i++) {
    if (environment.columnAsContinuous[i] == 0) count++;
  }
  if (count > 1) environment.atLeastTwoDiscrete = 1;

  for (uint i = 0; i < environment.numNodes; i++) {
    if (environment.columnAsContinuous[i] == 1)
      environment.atLeastOneContinuous = 1;
  }

  // create the data matrix for factors
  environment.dataNumeric = new int*[environment.numSamples];
  for (uint i = 0; i < environment.numSamples; i++) {
    environment.dataNumeric[i] = new int[environment.numNodes];
  }

  // for continuous non all gaussians
  if (environment.atLeastOneContinuous) {
    // create the data matrix for factors indexes
    environment.dataNumericIdx = new int*[environment.numNodes];
    for (uint i = 0; i < environment.numNodes; i++) {
      environment.dataNumericIdx[i] = new int[environment.numSamples];
      for (uint j = 0; j < environment.numSamples; j++)
        environment.dataNumericIdx[i][j] = -1;
    }
  }

  // transform to factors
  for (uint i = 0; i < environment.numNodes; i++) {
    if (environment.columnAsContinuous[i] == 0)
      transformToFactors(environment, i);
    else
      transformToFactorsContinuous(environment,
          i);  // update environment.dataNumeric not
               // taking into account repetition
  }

  // for continuous non gaussian
  if (environment.atLeastOneContinuous) {
    for (uint j = 0; j < environment.numNodes; j++) {
      if (environment.columnAsContinuous[j] != 0) {
        transformToFactorsContinuousIdx(environment, j);
        transformToFactors(environment, j);  // update environment.dataNumeric
                                             // taking into account repetition
      }
    }
  }

  // Set a variables with all properties name and levels
  setNumberLevels(environment);
  // create the 1000 entries to store c2 values
  environment.c2terms = new double[environment.numSamples + 1];
  for (uint i = 0; i < environment.numSamples + 1; i++) {
    environment.c2terms[i] = -1;
  }

  environment.initbins = std::min(30, int(0.5 + cbrt(environment.numSamples)));

  // for mixed
  // if(environment.atLeastOneContinuous){
  // create the log(j) lookup table with j=0..numSamples;
  environment.looklog = new double[environment.numSamples + 2];
  environment.looklog[0] = 0.0;
  for (uint i = 1; i < environment.numSamples + 2; i++) {
    environment.looklog[i] = log(1.0 * i);
  }

  environment.lookH = new double[environment.numSamples + 2];
  environment.lookH[0] = 0.0;
  for (uint i = 1; i < environment.numSamples + 2; i++) {
    // environment.lookH[i] =
    // i*environment.looklog[i]-(i+1)*environment.looklog[(i+1)];
    environment.lookH[i] = i * environment.looklog[i];
  }

  environment.logEta = 0;
  environment.isNoInitEta = 0;
  environment.firstIterationDone = false;

  int ncol = N_COL_NML;  // Number of levels r for which we want to store the
                         // stochastic NML complexity LogC(n,r) for n in [1,N].
                         // For r>N_COL_NML LogC() is computed with the normal
                         // recurrence and the result is not stored.
  environment.cterms = new double*[ncol];
  for (int K = 0; K < (ncol); K++) {
    environment.cterms[K] = new double[environment.numSamples + 1];
    for (uint i = 0; i < (environment.numSamples + 1); i++) {
      if (K == 1)
        environment.cterms[K][i] = 0;
      else if (i == 0)
        environment.cterms[K][i] = 0;
      else
        environment.cterms[K][i] = -1;
    }
  }
  for (uint i = 0; i < (environment.numSamples + 1); i++) {
    computeLogC(i, 2, environment.looklog,
        environment.cterms);  // Initialize the c2 terms
  }

  environment.lookchoose = new double*[ncol];
  for (int K = 0; K < (ncol); K++) {
    environment.lookchoose[K] = new double[environment.numSamples + 1];
    for (uint i = 0; i < (environment.numSamples + 1); i++) {
      environment.lookchoose[K][i] = -1;
    }
  }

  // Set the probability threshold for the rank
  environment.thresPc = 0;  // if the contribution probability is the min value

  // Stats test correction
  // environment.logEta = log( environment.eta);

  // create the edge structure and keep track of how many searchMore we have
  environment.edges = new Edge*[environment.numNodes];

  for (uint i = 0; i < environment.numNodes; i++)
    environment.edges[i] = new Edge[environment.numNodes];

  for (uint i = 0; i < environment.numNodes; i++) {
    for (uint j = 0; j < environment.numNodes; j++) {
      if ((environment.columnAsContinuous[i] == 0 &&
              environment.allLevels[i] == environment.numSamples) ||
          (environment.columnAsContinuous[j] == 0 &&
              environment.allLevels[j] == environment.numSamples)) {
        // If a node is discrete with as many levels as there are samples, its
        // information with other nodes is null.
        environment.edges[i][j].status = 0;
        environment.edges[i][j].status_prev = 0;
      } else {
        // Initialise all other edges.
        environment.edges[i][j].status = 1;
        environment.edges[i][j].status_prev = 1;
      }
    }
  }

  for (uint i = 0; i < environment.numNodes; i++) {
    environment.edges[i][i].status = 0;
    environment.edges[i][i].status_prev = 0;
  }

  environment.noiseVec = new double[2 * environment.numSamples];
  for (uint i = 0; i < 2 * environment.numSamples; i++) {
    environment.noiseVec[i] =
        std::rand() / ((RAND_MAX + 1u) / MAGNITUDE_TIES) - MAGNITUDE_TIES / 2;
  }
}

static void chkIntFn(void* dummy) { R_CheckUserInterrupt(); }

bool checkInterrupt(bool check /*=true*/) {
  if (check)
    return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
  else
    return false;
}

int printProgress(double percentage, double startTime, int prg_numSearchMore) {
  int pbwidth(40);
  string pbstr = string(pbwidth, '|');
  if (std::isnan(percentage) || std::isinf(percentage)) return 0;
  int val = (int)(percentage * 100);
  if (val != prg_numSearchMore) {
    int lpad = (int)(percentage * pbwidth);
    int rpad = pbwidth - lpad;
    double remaining_time =
        (get_wall_time() - startTime) / percentage * (1 - percentage);
    stringstream sremaining_time;
    if ((remaining_time > 60) & (!std::isinf(remaining_time))) {
      int minutes = remaining_time / 60;
      if (minutes > 60) {
        int hours = minutes / 60;
        sremaining_time << hours << "h";
      }
      sremaining_time << minutes % 60 << "m";
    }
    sremaining_time << int(remaining_time) % 60 << "s";
    printf("\r\t %3d%% [%.*s%*s] est. remaining time : %10s", val, lpad,
        pbstr.c_str(), rpad, "", sremaining_time.str().c_str());
    fflush(stdout);
  }
  return val;
}

// Compute the distance to the kth nearest neigbhour of the given point in the
// given space. The point object is the coordinates of the point with has as
// many dimensions as the space.
//
// Args
// 	-point		: a ndims-dimensional vector containing the coordinates
// of the given point. 	-kdTree		:
// 	-k			: the rank of the nearest neigbhour's distance
// to return
using my_kd_tree_t =
    KDTreeVectorOfVectorsAdaptor<vector<vector<double> >, double>;
double compute_k_nearest_distance(
    vector<double> point, my_kd_tree_t::index_t* index, int k) {
  vector<size_t> ret_indexes(k);
  vector<double> out_dists_sqr(k);
  nanoflann::KNNResultSet<double> resultSet(k);
  resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

  index->findNeighbors(resultSet, &point[0], nanoflann::SearchParams(10));

  return (sqrt(out_dists_sqr[k - 1]));
}

// Computes the Kullback-Leibler divergence between two joint (2D) distributions
// of real values based on the KNN estimation (F. Perez-Cruz 2004).
//
// <space 1> is the subsampling of <space 2> after removing NAs.
double compute_kl_divergence_continuous(vector<vector<double> > space1,
    vector<vector<double> > space2, int n1, int n2, int ndims, int k,
    bool* flag_break_ties, int* map_samples, double* noiseVec) {
  double D;
  double sumlog = 0;
  double noise;
  int i_map;
  for (int j = 0; j < ndims; j++) {
    i_map = 0;
    for (int i = 0; i < n2; i++) {
      if (flag_break_ties[j]) {
        noise = noiseVec[(j * n2) + i];
        space2[i][j] += noise;
        if (map_samples[i] == 1) {
          space1[i_map][j] += noise;
          i_map++;
        }
      }
    }
  }
  // construct a kd-tree index:
  // Dimensionality set at run-time (default: L2)
  my_kd_tree_t mat_index1(ndims, space1, 10);
  mat_index1.index->buildIndex();
  my_kd_tree_t mat_index2(ndims, space2, 10);
  mat_index2.index->buildIndex();

  for (int i = 0; i < n1; i++) {
    vector<double> point(ndims);
    for (int j = 0; j < ndims; j++) {
      point[j] = space1[i][j];
    }

    sumlog += log(compute_k_nearest_distance(point, mat_index2.index, k) /
                  compute_k_nearest_distance(point, mat_index1.index, k));
  }
  D = ndims * (sumlog / n1) + log(1.0 * (n2 - 1) / (n1 - 1));

  return (D);
}

double compute_kl_divergence(int* posArray, Environment& environment,
    int samplesNotNA, const vector<int> &AllLevels_red,
    const vector<int> &sample_is_not_NA) {
  int current_samplesNotNA =
      getNumSamples_nonNA(environment, posArray[0], posArray[1]);
  double kldiv = 0;

  // 1 - XY discrete
  if (environment.columnAsContinuous[posArray[0]] == 0 &&
      environment.columnAsContinuous[posArray[1]] == 0) {

    // Joint freqs X,Y before adding the new contributor (with the current
    // conditioning Us)
    double** jointFreqs = getJointFreqs(
        environment, posArray[0], posArray[1]);

    // Joint freqs X,Y after adding the new contributor (Z)
    double** freqs2 = getJointFreqs(
        environment, posArray[0], posArray[1], sample_is_not_NA);

    kldiv = samplesNotNA * kl(freqs2, jointFreqs,
                               environment.allLevels[posArray[0]],
                               environment.allLevels[posArray[1]]);

    for (uint j = 0; j < environment.allLevels[posArray[0]]; j++) {
      delete[] freqs2[j];
      delete[] jointFreqs[j];
    }
    delete[] freqs2;
    delete[] jointFreqs;
  } else if (environment.columnAsContinuous[posArray[0]] == 1 &&
             environment.columnAsContinuous[posArray[1]] == 1) {
    // 2 - XY continuous
    // Retrieve marginal distibutions with the current conditioning Us
    vector<vector<double> > joint_base(current_samplesNotNA, vector<double>(2));
    int* curr_sample_is_not_NA = new int[environment.numSamples];
    getJointSpace(environment, posArray[0], posArray[1], joint_base,
        curr_sample_is_not_NA);

    int* map_samples = new int[current_samplesNotNA];
    int i_map = 0;
    for (uint i = 0; i < environment.numSamples; i++) {
      if (curr_sample_is_not_NA[i] == 1) {  // sample i is present in X;Y|U
        map_samples[i_map] = 0;
        if (sample_is_not_NA[i] == 1) {  // sample i is also present in X;Y|U,Z
          map_samples[i_map] = 1;
        }
        i_map++;
      }
    }
    vector<vector<double> > joint_nonNA(samplesNotNA, vector<double>(2));
    int i_nonNA = 0;
    for (uint i = 0; i < environment.numSamples; i++) {
      if (sample_is_not_NA[i] == 1) {
        for (int k = 0; k < 2; k++) {
          joint_nonNA[i_nonNA][k] = environment.dataDouble[i][posArray[k]];
        }
        i_nonNA++;
      }
    }
    bool* flag_break_ties = new bool[2];
    for (int k = 0; k < 2; k++) {
      flag_break_ties[k] = false || (AllLevels_red[k] != samplesNotNA) ||
                           (AllLevels_red[k] != current_samplesNotNA);
    }

    kldiv =
        samplesNotNA * compute_kl_divergence_continuous(joint_nonNA, joint_base,
                           samplesNotNA, current_samplesNotNA, 2, KNN_K,
                           flag_break_ties, map_samples, environment.noiseVec);

    delete[] curr_sample_is_not_NA;
    delete[] map_samples;
    delete[] flag_break_ties;
  } else {
    // 3 - One discrete and one continuous
    int discrete_pos, continuous_pos, discrete_pos_binary, continuous_pos_binary;
    if (environment.columnAsContinuous[posArray[0]] == 0) {
      discrete_pos_binary = 0;
      continuous_pos_binary = 1;
    } else {
      discrete_pos_binary = 1;
      continuous_pos_binary = 0;
    }
    discrete_pos = posArray[discrete_pos_binary];
    continuous_pos = posArray[continuous_pos_binary];
    int n_discrete_levels = environment.allLevels[discrete_pos];
    // Full and reduced data may not have the same number of unique levels

    // Retrieve marginal distibutions with the current conditioning Us
    int current_samplesNotNA =
        getNumSamples_nonNA(environment, posArray[0], posArray[1]);
    int* mixedDiscrete = new int[current_samplesNotNA];
    double* mixedContinuous = new double[current_samplesNotNA];
    int* curr_sample_is_not_NA = new int[environment.numSamples];
    getJointMixed(environment, posArray[0], posArray[1], mixedDiscrete,
        mixedContinuous, curr_sample_is_not_NA);

    // Create count vectors for the discrete variable
    int* count_base = new int[n_discrete_levels];
    int* count_nonNA = new int[n_discrete_levels];
    for (int level = 0; level < n_discrete_levels; level++) {
      count_base[level] = 0;
      count_nonNA[level] = 0;
    }
    for (int i = 0; i < current_samplesNotNA; i++) {
      count_base[mixedDiscrete[i]]++;
    }
    for (uint i = 0; i < environment.numSamples; i++) {
      // Make sure to use environment data so that the levels match (may be
      // recoded in reduced data)
      if(sample_is_not_NA[i])
        count_nonNA[environment.dataNumeric[i][discrete_pos]]++;
    }

    // Compute the sum count(y) * KL(X_nonNA|y || X|y) over all values of Y y
    for (int level = 0; level < n_discrete_levels; level++) {
      int* map_level = new int[count_base[level]];
      int i_level = 0;
      for (uint i = 0; i < environment.numSamples; i++) {
        if (environment.dataNumeric[i][discrete_pos] == level) {
          if (curr_sample_is_not_NA[i]) {
            map_level[i_level] = 0;
            if (sample_is_not_NA[i] == 1) map_level[i_level] = 1;
            i_level++;
          }
        }
      }

      vector<vector<double> > continuous_base(
          count_base[level], vector<double>(1));
      i_level = 0;
      for (int i = 0; i < current_samplesNotNA; i++) {
        if (mixedDiscrete[i] == level) {
          continuous_base[i_level][0] = mixedContinuous[i];
          i_level++;
        }
      }

      vector<vector<double> > continuous_nonNA(
          count_nonNA[level], vector<double>(1));
      int i_level_nonNA = 0;
      for (uint i = 0; i < environment.numSamples; i++) {
        if (sample_is_not_NA[i] == 1 &&
            environment.dataNumeric[i][discrete_pos] == level) {
          continuous_nonNA[i_level_nonNA][0] =
              environment.dataDouble[i][continuous_pos];
          i_level_nonNA++;
        }
      }
      bool flag_break_ties[1];
      flag_break_ties[0] =
          false || (AllLevels_red[continuous_pos_binary] != samplesNotNA) ||
          (AllLevels_red[continuous_pos_binary] != current_samplesNotNA);

      if (count_nonNA[level] > KNN_K) {
        kldiv += fmax(0,
            count_nonNA[level] *
                compute_kl_divergence_continuous(continuous_nonNA,
                    continuous_base, count_nonNA[level], count_base[level], 1,
                    KNN_K, flag_break_ties, map_level, environment.noiseVec));
      }
      delete[] map_level;
    }  // level loop

    // add KL(Y(!NA) || Y)
    kldiv += samplesNotNA * kl(count_nonNA, count_base, samplesNotNA,
                                current_samplesNotNA, n_discrete_levels);

    delete[] count_nonNA;
    delete[] count_base;
    delete[] mixedDiscrete;
    delete[] mixedContinuous;
    delete[] curr_sample_is_not_NA;
  }
  return (kldiv);
}

int sign(double val) {
  if (val < 0)
    return -1;
  else if (val > 0)
    return 1;
  else
    return 0;
}

/**
 * Counts and marks the rows that contain "NA"s for the edge given by
 * \a posArray ([0] is X, [1] is Y and above are Us) and an optional \a z.
 *
 * \return The number of non NA samples and modifies the vectors
 * sample_is_not_NA and NAs_count
 */
uint count_non_NAs(int nbrUi, vector<int> &sample_is_not_NA,
    vector<int> &NAs_count, int* posArray,
    Environment& environment, int z){

  uint samplesNotNA = 0;
  bool is_NA;

  for (uint i = 0; i < environment.numSamples; i++) {
    sample_is_not_NA[i] = 1;
    if (i != 0)
      NAs_count[i] = NAs_count[i - 1];
    else
      NAs_count[i] = 0;

    is_NA = false;
    for (int j = 0; (j < nbrUi + 2) && (!is_NA); j++) {
      is_NA = environment.dataNumeric[i][posArray[j]] == -1;
    }
    if (z != -1) is_NA = is_NA || environment.dataNumeric[i][z] == -1;
    if (is_NA) {
      sample_is_not_NA[i] = 0;
      NAs_count[i]++;
    } else
      samplesNotNA++;
  }

  return(samplesNotNA);
}

/**
 * Reduces the full data to its "reduced" form without NA rows given the edge
 * X-Y | U,(Z). In the case of discrete variables, levels are remapped so that
 * the reduced data still has levels that start at 0. For continuous variables,
 * the ranks are updated (\a dataNumericIdx).
 *
 * \return A boolean telling whether or not sample weights are needed on the
 * reduced data. More importantly, modifies the vectors passed as argument :
 * \a AllLevels contains the number of unique levels of the reduced variables
 * \a cnt tells whether the reduced variables are continuous or not
 * \a posArray_red contains the positions of variables in the reduced data
 * \a dataNumeric contains the discrete levels of the reduced data
 * \a dataNumericIdx contains the continuous ranks of the reduced data
 * \a sample_weights contains the individual sample weights of the reduced data
 */
bool filter_NAs(int nbrUi, vector<int> &AllLevels, vector<int> &cnt,
    vector<int> &posArray_red, int* posArray, vector<vector<int> > &dataNumeric,
    vector<vector<int> > &dataNumericIdx, vector<double> &sample_weights,
    const vector<int> &sample_is_not_NA,
    const vector<int> &NAs_count,
    Environment& environment, int z) {

  int k1, k2, nnr, prev_val, si, old_val, new_val, updated_discrete_level;
  int column;
  bool flag_sample_weights(false);

  // Map to make sure that the levels of the reduced data start at zero
  std::map<int, int> new_levels;

  for (int j = 0; j < (nbrUi + 2 + int(z!=-1)); j++) {
    if(j < nbrUi+2) {
      column = posArray[j];
    } else {
      column = z;
    }
    posArray_red[j] = j;
    cnt[j] = environment.columnAsContinuous[column];

    k1 = 0; // position in the full data
    k2 = 0; // position in the reduced data
    nnr = 0; // number of non repeated values
    prev_val = -1; // previous value
    new_levels.clear();
    new_val = 0;
    updated_discrete_level = 0;

    for (uint i = 0; i < environment.numSamples; i++) {

      if (sample_is_not_NA[i] == 1) {
        // Row at index i does not contain any NA
        old_val = environment.dataNumeric[i][column];
        if (new_levels.count(old_val)==0){
          // If level has not already been seen add it to the map
          // and increment the number of unique levels in reduced data
          new_levels.insert({old_val, updated_discrete_level});
          updated_discrete_level ++;
        }
        new_val = new_levels[old_val];

        dataNumeric[j][k1] = new_val;
        if(j==0) {
          sample_weights[k1] = environment.sampleWeights[i];
          if(sample_weights[k1] != 1.0) flag_sample_weights = true;
        }
        k1++;
      }
      if (cnt[j] != 0) {
        // Variable j is continuous
        si = environment.dataNumericIdx[column][i]; //position of ith sample (order)
        if (si != -1 && sample_is_not_NA[si] == 1) {
          // Row at position si does not contain any NA, rank is updated taking
          // into account the number of NAs up to si.
          dataNumericIdx[j][k2] = si - NAs_count[si];
          k2++;
          // check whether is a different values or repeated
          if (environment.dataNumeric[si][column] != prev_val) {
            nnr++;
            prev_val = environment.dataNumeric[si][column];
          }
        }
      }
    }
    // Update with the effective number of levels
    if (cnt[j] == 1) {
      if (nnr < 3) cnt[j] = 0;
      AllLevels[j] = nnr;
    }
    else AllLevels[j] = updated_discrete_level;
  }
  return(flag_sample_weights);
}

}  // namespace utility
}  // namespace miic
