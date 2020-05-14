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

#define _DEBUG 0

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

vector<vector<string>> getAdjMatrix(const Environment& environment) {
  stringstream ss;
  vector<vector<string>> adjMatrix;
  vector<string> vec;
  for (uint i = 0; i < environment.numNodes; i++) {
    vec.push_back(environment.nodes[i].name);
  }

  adjMatrix.push_back(vec);

  for (uint i = 0; i < environment.numNodes; i++) {
    vec.clear();
    vec.push_back(environment.nodes[i].name);
    for (uint j = 0; j < environment.numNodes; j++) {
      ss.str("");
      ss << environment.edges[i][j].status;
      vec.push_back(ss.str());
    }
    adjMatrix.push_back(vec);
  }
  return adjMatrix;
}

void createMemorySpace(Environment& environment, MemorySpace& m) {
  if (environment.atLeastTwoDiscrete) {
    uint maxLevel = 0;
    for (uint i = 0; i < environment.numNodes; i++) {
      if (environment.columnAsContinuous[i] == 0 &&
          environment.allLevels[i] > maxLevel)
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
    m.samplesToEvaluate = (int*)new int[environment.numSamples];
    m.samplesToEvaluateTemplate = (int*)new int[environment.numSamples];
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
    free(m.samplesToEvaluate);
    free(m.samplesToEvaluateTemplate);
    free(m.AllLevels_red);
    free(m.cnt_red);
    free(m.posArray_red);
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

void computeMeansandStandardDeviations(Environment& environment) {
  for (uint c = 0; c < environment.numNodes; c++) {
    if (environment.columnAsGaussian[c] == 1) {
      int nsamples = 0;
      for (uint r = 0; r < environment.numSamples; r++) {
        if (!std::isnan(environment.dataDouble[r][c])) {
          environment.means[c] += environment.dataDouble[r][c];
          nsamples++;
        }
      }

      environment.means[c] /= (double)nsamples;

      for (uint r = 0; r < environment.numSamples; r++) {
        if (!std::isnan(environment.dataDouble[r][c]))
          environment.standardDeviations[c] +=
              pow((environment.dataDouble[r][c] - environment.means[c]), 2);
      }

      environment.standardDeviations[c] /= (double)nsamples;
      environment.standardDeviations[c] =
          sqrt(environment.standardDeviations[c]);
    }
  }
}

void computeCorrelations(Environment& environment) {
  double covariance = 0.0;
  for (uint i = 0; i < environment.numNodes; i++) {
    for (uint j = 0; j < environment.numNodes; j++) {
      if (environment.columnAsGaussian[i] == 1 &&
          environment.columnAsGaussian[j] == 1) {
        // cout << "i: " << i << "j: " << j <<  endl << flush;
        if (i != j) {
          covariance = 0.0;
          int nSamples = 0;
          // cout << "environment.numSamples: " << environment.numSamples  <<
          // endl << flush;
          for (uint k = 0; k < environment.numSamples; k++) {
            // cout << "environment.dataDouble[k][i]: " <<
            // environment.dataDouble[k][i]  <<  endl << flush;

            if (!std::isnan(environment.dataDouble[k][i]) &&
                !std::isnan(environment.dataDouble[k][j])) {
              covariance +=
                  (environment.dataDouble[k][i] - environment.means[i]) *
                  (environment.dataDouble[k][j] - environment.means[j]);
              nSamples++;
            }
          }

          // divide covariance by nCols
          covariance /= nSamples;
          environment.nSamples[i][j] = nSamples;
          environment.rho[i][j] =
              covariance / (environment.standardDeviations[i] *
                               environment.standardDeviations[j]);
        } else {
          environment.rho[i][j] = 1;
        }
      }
    }
  }
}

bool isOnlyDouble(const char* str) {
  char* endptr = 0;
  strtod(str, &endptr);

  if (*endptr != '\0' || endptr == str) return false;
  return true;
}

bool comparatorPairs(const pair<double, int>& l, const pair<double, int>& r) {
  return l.first < r.first;
}

bool SortFunctionNoMore1(
    const EdgeID* a, const EdgeID* b, const Environment& environment) {
  return environment.edges[a->i][a->j].shared_info->Ixy_ui >
         environment.edges[b->i][b->j].shared_info->Ixy_ui;
}

class sorterNoMore {
  Environment& environment;

 public:
  sorterNoMore(Environment& env) : environment(env) {}
  bool operator()(EdgeID const* o1, EdgeID const* o2) const {
    return SortFunctionNoMore1(o1, o2, environment);
  }
};

bool SortFunction(
    const EdgeID* a, const EdgeID* b, const Environment& environment) {
  if (environment.edges[a->i][a->j].shared_info->connected >
      environment.edges[b->i][b->j].shared_info->connected)
    return true;
  else if (environment.edges[a->i][a->j].shared_info->connected <
           environment.edges[b->i][b->j].shared_info->connected)
    return false;

  if (environment.edges[a->i][a->j].shared_info->connected == 0 &&
      environment.edges[b->i][b->j].shared_info->connected == 0) {
    if (environment.edges[a->i][a->j].shared_info->Rxyz_ui == 0 &&
        environment.edges[b->i][b->j].shared_info->Rxyz_ui != 0)
      return true;
    else if (environment.edges[a->i][a->j].shared_info->Rxyz_ui != 0 &&
             environment.edges[b->i][b->j].shared_info->Rxyz_ui == 0)
      return false;

    if (environment.edges[a->i][a->j].shared_info->Rxyz_ui >
        environment.edges[b->i][b->j].shared_info->Rxyz_ui)
      return true;
    else if (environment.edges[a->i][a->j].shared_info->Rxyz_ui <
             environment.edges[b->i][b->j].shared_info->Rxyz_ui)
      return false;
  }

  if (environment.edges[a->i][a->j].shared_info->connected == 1 &&
      environment.edges[b->i][b->j].shared_info->connected == 1) {
    if (environment.edges[a->i][a->j].shared_info->Ixy_ui >
        environment.edges[b->i][b->j].shared_info->Ixy_ui)
      return true;
    else if (environment.edges[a->i][a->j].shared_info->Ixy_ui <
             environment.edges[b->i][b->j].shared_info->Ixy_ui)
      return false;
  }
  return false;
}

class sorter {
  const Environment& environment;

 public:
  sorter(const Environment& env) : environment(env) {}
  bool operator()(EdgeID const* o1, EdgeID const* o2) const {
    return SortFunction(o1, o2, environment);
  }
};

void readTime(Environment& environment, string name) {
  const char* c = name.c_str();
  std::ifstream input(c);
  string lineData;
  string s;
  int row = 0;
  int col = 0;
  while (getline(input, lineData)) {
    if (row == 1) {
      std::istringstream f(lineData);
      while (getline(f, s, '\t')) {
        if (col == 0)
          environment.execTime.init = atof(s.c_str());
        else if (col == 1)
          environment.execTime.iter = atof(s.c_str());
        else if (col == 2)
          environment.execTime.initIter = atof(s.c_str());
        else if (col == 3)
          environment.execTime.ort = atof(s.c_str());
        else if (col == 4)
          environment.execTime.cut = atof(s.c_str());
        else if (col == 5)
          environment.execTime.ort_after_cut = atof(s.c_str());
        else if (col == 6)
          environment.execTime.total = atof(s.c_str());
        col++;
      }
    }
    row++;
  }
}

void saveAdjMatrix(const Environment& environment, const string filename) {
  if (environment.isVerbose) cout << "Saving adjacency matrix\n";
  std::ofstream output;
  output.open(filename.c_str());
  for (uint i = 0; i < environment.numNodes; i++) {
    output << environment.nodes[i].name;
    if (i + 1 < environment.numNodes) output << "\t";
  }
  output << endl;

  for (uint i = 0; i < environment.numNodes; i++) {
    output << environment.nodes[i].name << "\t";
    for (uint j = 0; j < environment.numNodes; j++) {
      output << environment.edges[i][j].status;
      if (j + 1 < environment.numNodes) output << "\t";
    }
    output << endl;
  }
  output.close();
}

void saveAdjMatrixState(const Environment& environment, const string filename) {
  if (environment.isVerbose) cout << "Saving adjacency matrix\n";
  std::ofstream output;
  output.open(filename.c_str());
  output << "\t";
  for (uint i = 0; i < environment.numNodes; i++) {
    output << environment.nodes[i].name;
    if (i + 1 < environment.numNodes) output << "\t";
  }
  output << endl;

  for (uint i = 0; i < environment.numNodes; i++) {
    output << environment.nodes[i].name << "\t";
    for (uint j = 0; j < environment.numNodes; j++) {
      if (j > i) {
        if (environment.edges[i][j].shared_info->connected == 1)
          output << "1";
        else
          output << "0";
        if (j + 1 < environment.numNodes) output << "\t";
      } else if (i > j) {
        if (environment.edges[j][i].shared_info->connected == 1)
          output << "1";
        else
          output << "0";
        if (j + 1 < environment.numNodes) output << "\t";
      } else {
        output << "0";
        if (j + 1 < environment.numNodes) output << "\t";
      }
    }
    output << endl;
  }
}

string vectorToStringNodeName(
    const Environment& environment, const vector<int>& vec) {
  stringstream ss;
  int length = vec.size();
  if (length > 0) {
    for (int temp = 0; temp < length; temp++) {
      if (vec[temp] != -1) ss << environment.nodes[vec[temp]].name;
      if (temp + 1 < length) ss << ",";
    }
  } else {
    ss << "NA";
  }
  return ss.str();
}

string vectorToString(const vector<int>& vec) {
  stringstream ss;
  int length = vec.size();
  if (length > 0) {
    for (int temp = 0; temp < length; temp++) {
      ss << vec[temp];
      if (temp + 1 < length) ss << ",";
    }
  }
  return ss.str();
}

string arrayToString(const double* int_array, const int length) {
  stringstream ss;
  if (length > 0) {
    for (int temp = 0; temp < length; temp++) {
      if (int_array[temp] != -1) ss << int_array[temp] << ", ";
    }
  } else {
    ss << "NA";
  }
  return ss.str();
}

string zNameToString(const Environment& environment, vector<int> vec, int pos) {
  stringstream ss;
  if (pos != -1)
    ss << environment.nodes[vec[pos]].name;
  else
    ss << "NA";
  return ss.str();
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

vector<vector<string>> saveEdgesListAsTable(Environment& environment) {
  vector<vector<string>> data;

  vector<EdgeID*> allEdges;

  for (uint i = 0; i < environment.numNodes - 1; i++) {
    for (uint j = i + 1; j < environment.numNodes; j++) {
      allEdges.emplace_back(new EdgeID(i, j));
    }
  }

  vector<string> row;

  std::sort(allEdges.begin(), allEdges.end(), sorter(environment));

  row.push_back("x");
  row.push_back("y");
  row.push_back("z.name");
  row.push_back("ai.vect");
  row.push_back("zi.vect");
  row.push_back("Ixy");
  row.push_back("Ixy_ai");
  row.push_back("cplx");
  row.push_back("Rxyz_ai");
  row.push_back("category");
  row.push_back("Nxy_ai");

  data.push_back(row);

  for (uint i = 0; i < allEdges.size(); i++) {
    stringstream output;
    row.clear();
    for (uint j = 0; j < environment.numNodes; j++) {
      if (j == allEdges[i]->i || j == allEdges[i]->j)
        output << "1";
      else
        output << "0";
    }
    row.push_back(output.str());
    row.push_back(environment.nodes[allEdges[i]->i].name);
    row.push_back(environment.nodes[allEdges[i]->j].name);
    row.push_back(zNameToString(environment,
        environment.edges[allEdges[i]->i][allEdges[i]->j]
            .shared_info->zi_vect_idx,
        environment.edges[allEdges[i]->i][allEdges[i]->j]
            .shared_info->z_name_idx));
    row.push_back(vectorToStringNodeName(
        environment, environment.edges[allEdges[i]->i][allEdges[i]->j]
                         .shared_info->ui_vect_idx));
    row.push_back(vectorToStringNodeName(
        environment, environment.edges[allEdges[i]->i][allEdges[i]->j]
                         .shared_info->zi_vect_idx));

    output.str("");
    output << environment.edges[allEdges[i]->i][allEdges[i]->j]
                  .shared_info->mutInfo;
    row.push_back(output.str());

    output.str("");
    output << environment.edges[allEdges[i]->i][allEdges[i]->j]
                  .shared_info->Ixy_ui;
    row.push_back(output.str());

    output.str("");
    output
        << environment.edges[allEdges[i]->i][allEdges[i]->j].shared_info->cplx;
    row.push_back(output.str());

    output.str("");
    output << environment.edges[allEdges[i]->i][allEdges[i]->j]
                  .shared_info->Rxyz_ui;
    row.push_back(output.str());

    output.str("");
    output << environment.edges[allEdges[i]->i][allEdges[i]->j]
                  .shared_info->connected;
    row.push_back(output.str());

    output.str("");
    output << environment.edges[allEdges[i]->i][allEdges[i]->j]
                  .shared_info->Nxy_ui;
    row.push_back(output.str());

    data.push_back(row);
  }

  for (uint i = 0; i < allEdges.size(); i++) {
    delete allEdges[i];
  }

  return data;
}

void saveExecTime(const Environment& environment, const string filename) {
  if (environment.isVerbose) cout << "Saving execution time\n";
  std::ofstream output;
  output.open(filename.c_str());

  output << "init"
         << "\t"
         << "iter"
         << "\t"
         << "initIter"
         << "\t"
         << "ort"
         << "\t"
         << "cut"
         << "\t"
         << "ort_after_cut"
         << "\t"
         << "total"
         << "\n";
  output << environment.execTime.init << "\t" << environment.execTime.iter
         << "\t" << environment.execTime.initIter << "\t"
         << environment.execTime.ort << "\t" << environment.execTime.cut << "\t"
         << environment.execTime.ort_after_cut << "\t"
         << environment.execTime.total;
}

bool existsTest(const string& name) {
  std::ifstream f(name.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }
}

int** copyMatrix(int** oldmatrix, int numRows, int numColumns) {
  int** newMatrix;
  newMatrix = new int*[numRows];
  for (int i = 0; i < numRows; i++) newMatrix[i] = new int[numColumns];

  for (int i = 0; i < numRows; i++) {
    for (int j = 0; j < numColumns; j++) {
      newMatrix[i][j] = oldmatrix[i][j];
    }
  }
  return newMatrix;
}

bool checkNA(int** data, int numRows, int numColumns) {
  for (int i = 0; i < numRows; i++)
    for (int j = 0; j < numColumns; j++)
      if (data[i][j] == -1) return true;

  return false;
}

string printNodesName(const Environment& environment) {
  string s = "";
  for (uint i = 0; i < environment.numNodes; i++) {
    cout << environment.nodes[i].name;
    if (i + 1 < environment.numNodes) cout << " ";
  }
  cout << "\n";
  return s;
}

//-----------------------------------------------------------------------------
// printEdges
//-----------------------------------------------------------------------------
// Description: print the list of edges in environment.edges 
// (only half of them as the matrix is symetrical)
//
// Params: 
// - Environment&: the environment structure
// - filter_status: boolean. Optional, true by default. 
//   When true, displays only edges with status != 0.
//   When false, displays all edges.
// - half_only: boolean. Optional, true by default. 
//   When true, displays only edges with col > row (matrix is symetrical).
//   When false, displays all edges.
//
// Returns: None
//-----------------------------------------------------------------------------
void printEdges (Environment& environment, bool filter_status, bool half_only) 
  {
  string half_or_full_str = "";
  if (half_only)
    half_or_full_str = "half list";
  else
    half_or_full_str = "full list";
  if (filter_status)
    cout << "List of edges in environment.edges (" << half_or_full_str << " filtered on status != 0):\n";
  else
    cout << "List of edges in environment.edges (" << half_or_full_str << "):\n";
    
  cout << "node\tnode\tstat\tstat\tstat\tconnect.  Nxy      mutInfo      cplx_noU   Nxy_ui cplx Ixy_ui Rxyz_ui z_name_idx zi_vect_idx  ui_vect_idx" << endl;
  cout << "  1 \t  2 \t    \tinit\tprev\t        nb joint mutual info   Complexity                      Score  Index last candid.nodes  Indice of" << endl;
  cout << "    \t    \t    \t   \t    \t         factors    without       without                       best     best     contributing  separating" << endl;
  cout << "    \t    \t    \t   \t    \t          not NA  conditioning conditioning                    contrib  contrib   cond. indep.   nodes" << endl;

  for (uint i = 0; i < environment.numNodes; i++) 
    {
    uint j;
    if (half_only)
      j = i + 1;
    else
      j = 0;
    for (; j < environment.numNodes; j++) 
      {
      const Edge& one_edge = environment.edges[i][j];
      std::shared_ptr<EdgeSharedInfo> shared_info_ptr = one_edge.shared_info;
      
      if ( (!filter_status) || (environment.edges[i][j].status) )
        {
        cout << environment.nodes[i].name << "\t"
             << environment.nodes[j].name << "\t"
             << environment.edges[i][j].status << "\t"
             << environment.edges[i][j].status_init << "\t"
             << environment.edges[i][j].status_prev << "\t";
          
        if (shared_info_ptr != NULL)
          cout << environment.edges[i][j].shared_info->connected << "\t"
               << environment.edges[i][j].shared_info->Nxy << "\t"
               << environment.edges[i][j].shared_info->mutInfo << "\t"
               << environment.edges[i][j].shared_info->cplx_noU << "\t"
               << environment.edges[i][j].shared_info->Nxy_ui << "\t"
               << environment.edges[i][j].shared_info->cplx << "\t"
               << environment.edges[i][j].shared_info->Ixy_ui << "\t"
               << environment.edges[i][j].shared_info->Rxyz_ui << "\t"
               << environment.edges[i][j].shared_info->z_name_idx << "\t"
               << vectorToStringNodeName (environment, environment.edges[i][j].shared_info->zi_vect_idx) << "\t"
               << vectorToStringNodeName (environment, environment.edges[i][j].shared_info->ui_vect_idx) << "\t";
        cout << endl; 
        }
      }
    }
  }

//-----------------------------------------------------------------------------
// printNoMoreAdress
//-----------------------------------------------------------------------------
// Desciption: print the list of edges in environment.noMoreAdress
//
// Params: 
// - Environment&: the environment structure
//
// Returns: None
//-----------------------------------------------------------------------------
void printNoMoreAdress (Environment& environment) 
  {
  cout << "List of edges in environment.noMoreAdress:\n";
  cout << "Node 1 idx + (name) - Node 2 idx + (name)" << endl;
  
  for (int edge_idx=0; edge_idx < environment.noMoreAddress.size(); edge_idx++)
    {
    int posX = environment.noMoreAddress[edge_idx]->i;
    int posY = environment.noMoreAddress[edge_idx]->j;
    cout << posX << " (" << environment.nodes[posX].name << ") - "
         << posY << " (" << environment.nodes[posY].name << ")" << endl;
    }
  cout << endl;
  }

//-----------------------------------------------------------------------------
// printAdjacencyMatrix
//-----------------------------------------------------------------------------
// Desciption: 
// print an adjacency matrix from the list of edges environment.edges
//
// Params: 
// - Environment&: the environment structure
// - status_field: string. Optional, "status" by default.
//   Control what information is used to print the adjacency matrix.
//   Possible values are "status", "status_init" and "status_prev"  
//
// Returns: None
//-----------------------------------------------------------------------------
void printAdjacencyMatrix (Environment& environment, std::string status_field) 
  {
  cout << "Adjacency matrix using environment.edges on col " << status_field << endl;
  
  for (uint i = 0; i < environment.numNodes; i++) 
    cout <<  "\t" << environment.nodes[i].name;
  cout << endl;
  for (uint i = 0; i < environment.numNodes; i++) 
    {
    cout << environment.nodes[i].name;
    for (uint j = 0; j < environment.numNodes; j++) 
      {
      if (status_field == "status_init")      
        cout << "\t" << environment.edges[i][j].status_init;
      else if (status_field == "status_prev")      
        cout << "\t" << environment.edges[i][j].status_prev;
      else
        cout << "\t" << environment.edges[i][j].status;
      }
    cout << endl;
    }
  }


void printMatrix(const Environment& environment, string type) {
  if (type.compare("string") == 0) {
    cout << "Data matrix of strings\n";
    printNodesName(environment);
    for (uint i = 0; i < environment.numSamples; i++) {
      for (uint j = 0; j < environment.numNodes; j++) {
        cout << environment.data[i][j] << " ";
      }
      cout << endl;
    }
  } else if (type.compare("factors") == 0) {
    cout << "Data matrix of factors\n";
    printNodesName(environment);
    for (uint i = 0; i < environment.numSamples; i++) {
      for (uint j = 0; j < environment.numNodes; j++) {
        cout << environment.dataNumeric[i][j] << " ";
      }
      cout << endl;
    }
  }
}

// Initialize all the elements of the array to the given value
bool setArrayValuesInt(int* array, int length, int value) {
  for (int i = 0; i < length; i++) {
    array[i] = value;
  }
  return true;
}

bool isInteger(const string& s) {
  if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+')))
    return false;

  char* p;
  strtol(s.c_str(), &p, 10);

  return (*p == 0);
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

void copyValue(Environment& environment, int i) {
  for (uint j = 0; j < environment.numSamples; j++) {
    environment.dataNumeric[j][i] = atof(environment.data[j][i].c_str());
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

  // vector <int> clmn;
  for (uint j = 0; j < environment.numSamples; j++) {
    int entry = environment.dataNumeric[j][i];
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

void setProportions(Environment& environment) {
  environment.proportions = new double*[environment.numNodes];
  int count = 0;
  for (uint i = 0; i < environment.numNodes; i++) {
    count = 0;
    if (environment.columnAsContinuous[i] == 0) {
      environment.proportions[i] = new double[environment.allLevels[i]];
      for (uint j = 0; j < environment.allLevels[i]; j++) {
        environment.proportions[i][j] = 0;
      }
      for (uint j = 0; j < environment.numSamples; j++) {
        if (environment.dataNumeric[j][i] > -1) {
          environment.proportions[i][environment.dataNumeric[j][i]]++;
          count++;
        }
      }

      for (uint j = 0; j < environment.allLevels[i]; j++) {
        environment.proportions[i][j] /= double(count);
      }
    }
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
    vector<vector<double>>& jointSpace, int* curr_samplesToEvaluate) {
  int numSamples_nonNA = 0;
  bool sampleOk;
  for (uint k = 0; k < environment.numSamples; k++) {
    sampleOk = true;
    curr_samplesToEvaluate[k] = 0;
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
        curr_samplesToEvaluate[k] = 1;
        jointSpace[numSamples_nonNA][0] = environment.dataDouble[k][i];
        jointSpace[numSamples_nonNA][1] = environment.dataDouble[k][j];
        numSamples_nonNA++;
      }
    }
  }
}

double** getJointFreqs(
    Environment& environment, int i, int j, int numSamples_nonNA) {
  double** jointFreqs = new double*[environment.allLevels[i]];
  for (uint k = 0; k < environment.allLevels[i]; k++) {
    jointFreqs[k] = new double[environment.allLevels[j]];
    for (uint l = 0; l < environment.allLevels[j]; l++) {
      jointFreqs[k][l] = 0;
    }
  }

  bool sampleOk;
  for (uint k = 0; k < environment.numSamples; k++) {
    sampleOk = true;
    if (environment.dataNumeric[k][i] != -1 &&
        environment.dataNumeric[k][j] != -1) {
      for (uint u = 0;
           u < environment.edges[i][j].shared_info->ui_vect_idx.size(); u++) {
        if (environment.dataNumeric[k][environment.edges[i][j]
                                           .shared_info->ui_vect_idx[u]] == -1)
          sampleOk = false;
      }
      if (sampleOk)
        jointFreqs[environment.dataNumeric[k][i]]
                  [environment.dataNumeric[k][j]]++;
    }
  }

  for (uint k = 0; k < environment.allLevels[i]; k++)
    for (uint l = 0; l < environment.allLevels[j]; l++)
      jointFreqs[k][l] /= numSamples_nonNA;

  return (jointFreqs);
}

void getJointMixed(Environment& environment, int i, int j, int* mixedDiscrete,
    double* mixedContinuous, int* curr_samplesToEvaluate) {
  int discrete_pos = environment.columnAsContinuous[i] == 0 ? i : j;
  int continuous_pos = environment.columnAsContinuous[i] == 1 ? i : j;

  // Fill marginal distributions
  int numSamples_nonNA = 0;
  bool sampleOk;
  for (uint k = 0; k < environment.numSamples; k++) {
    sampleOk = true;
    curr_samplesToEvaluate[k] = 0;
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
        curr_samplesToEvaluate[k] = 1;
        mixedContinuous[numSamples_nonNA] =
            environment.dataDouble[k][continuous_pos];
        mixedDiscrete[numSamples_nonNA] =
            environment.dataNumeric[k][discrete_pos];
        numSamples_nonNA++;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
// FRS 6 may 2010 parseCommandLine is not used anymore
//////////////////////////////////////////////////////////////////////////////////
// 
// bool parseCommandLine(Environment& environment, int argc, char** argv) {
//   int c;
//   environment.inData = "";
//   environment.outDir = "";
//   environment.blackbox_name = "";
//   environment.edgeFile = "";
//   environment.effN = -1;
//   environment.cplx = 1;
//   environment.isVerbose = false;
//   environment.numberShuffles = 0;
//   environment.isLatent = false;
//   environment.isLatentOnlyOrientation = false;
//   environment.isTplReuse = true;
//   environment.isK23 = true;
//   environment.isDegeneracy = false;
//   environment.isNoInitEta = false;
//   environment.isPropagation = true;
//   environment.halfVStructures = 0;
//   environment.typeOfData = 0;
//   environment.isAllGaussian = 0;
//   environment.atLeastTwoGaussian = 0;
//   environment.atLeastTwoDiscrete = 0;
//   environment.atLeastOneContinuous = 0;
//   environment.nThreads = 0;
//   environment.testDistribution = true;
//   environment.consistentPhase = 0;
// 
//   string s;
//   while ((c = getopt(argc, argv,
//               "j:i:o:b:d:c:e:s:r:q:k:n:p:a:h:m:t:u:z:x:l:gfv?")) != -1) {
//     switch (c) {
//       case 'i': {
//         environment.inData.append(optarg);
//         break;
//       }
//       case 'o': {
//         environment.outDir.append(optarg);
//         break;
//       }
//       case 'x': {
//         environment.seed = 0;
//         s = optarg;
//         if (!isInteger(s)) {
//           cout << "[ERR] Seed should be an integer!\n";
//         } else
//           environment.seed = atoi(optarg);
//         break;
//       }
//       case 'b': {
//         environment.blackbox_name.append(optarg);
//         break;
//       }
//       case 'd': {
//         s = optarg;
//         stringstream ss(s);  // Turn the string into a stream.
//         string tok;
//         char delimiter = ',';
//         while (getline(ss, tok, delimiter)) {
//           int ival = atoi(tok.c_str());
//           environment.steps.push_back(ival);
//         }
//         break;
//       }
//       case 'n': {
//         s = optarg;
//         if (!isInteger(s) && atoi(optarg) > 1) {
//           cout << "[ERR] EffeN should be a positive integer!\n";
//         } else
//           environment.effN = atoi(optarg);
//         break;
//       }
//       case 'h': {
//         environment.isAllGaussian = atoi(optarg);
//         break;
//       }
//       case 't': {
//         s = optarg;
//         if (!isInteger(s)) {
//           cout << "[ERR] Type of data should be an integer!\n";
//         } else {
//           if (atoi(optarg) == 0 || atoi(optarg) == 1 || atoi(optarg) == 2)
//             environment.typeOfData = atoi(optarg);
//           else
//             exit(1);
//         }
//         break;
//       }
//       case 'u': {
//         s = optarg;
//         environment.dataTypeFile.append(optarg);
//         break;
//       }
//       case 'm': {
//         environment.edgeFile.append(optarg);
//         break;
//       }
//       case 'z': {
//         environment.nThreads = atoi(optarg);
//         break;
//       }
//       case 'c': {
//         environment.cplxType.append(optarg);
// 
//         if (environment.cplxType.compare("mdl") != 0 &&
//             environment.cplxType.compare("nml") != 0) {
//           cout << "[ERR] Wrong complexity check option!\n";
//           exit(1);
//         } else if (environment.cplxType.compare("mdl") == 0) {
//           environment.cplx = 0;
//         }
//         break;
//       }
//       case 'e': {
//         s = optarg;
//         if (!isOnlyDouble(optarg)) {
//           cout << "[ERR] Confidence cut should be a double!\n";
//           exit(1);
//         }
// 
//         environment.confidenceThreshold = atof(optarg);
//         break;
//       }
//       case 's': {
//         s = optarg;
//         if (!isInteger(s)) {
//           cout << "[ERR] Shuffle should be an integer!\n";
//           exit(1);
//         } else
//           environment.numberShuffles = atoi(optarg);
//         break;
//       }
//       case 'r': {
//         s = optarg;
//         if (s.compare("1") != 0 && s.compare("0") != 0) {
//           cout << "[ERR] Wrong reuse/not reuse tpl argument!\n";
//           return false;
//         } else if (s.compare("0") == 0)
//           environment.isTplReuse = false;
//         break;
//       }
//       case 'k': {
//         s = optarg;
//         if (s.compare("1") != 0 && s.compare("0") != 0) {
//           cout << "[ERR] Case k: Wrong k23 argument!\n";
//           exit(1);
//         } else if (s.compare("0") == 0)
//           environment.isK23 = false;
//         break;
//       }
//       case 'p': {
//         s = optarg;
//         if (s.compare("0") != 0) {
//           cout << "[ERR] Case p: Wrong propagation argument!*" << s << "*\n";
//           exit(1);
//         } else if (s.compare("0") == 0)
//           environment.isPropagation = false;
//         break;
//       }
//       case 'j': {
//         s = optarg;
//         if (s.compare("1") != 0 && s.compare("0") != 0) {
//           cout << "[ERR] Case j: Wrong consistent argument!*" << s << "*\n";
//           exit(1);
//         } else if (s.compare("1") == 0)
//           environment.consistentPhase = true;
//         break;
//       }
//       case 'a': {
//         s = optarg;
//         if (s.compare("0") != 0 && s.compare("1") != 0) {
//           cout << "[ERR] Case a: Wrong half V-structures argument!*" << s
//                << "*\n";
//           exit(1);
//         } else
//           environment.halfVStructures = atoi(optarg);
//         break;
//       }
//       case 'l': {
//         s = optarg;
//         if (s.compare("0") != 0 && s.compare("1") != 0 && s.compare("2") != 0) {
//           cout << "[ERR] Wrong latent argument!*" << s << "*\n";
//           exit(1);
//         } else {
//           if (s.compare("1") == 0) {
//             environment.isLatent = true;
//           }
//           if (s.compare("2") == 0) {
//             environment.isLatentOnlyOrientation = true;
//           }
//         }
//         break;
//       }
//       case 'g': {
//         environment.isDegeneracy = true;
//         break;
//       }
//       case 'f': {
//         environment.isNoInitEta = true;
//         break;
//       }
//       case 'v': {
//         environment.isVerbose = true;
//         break;
//       }
// 
//       case '?': {
//         if (optopt == 'c')
//           fprintf(stderr, "Option -%c requires an argument.\n", optopt);
//         else if (isprint(optopt))
//           fprintf(stderr, "Unknown option `-%c'.\n", optopt);
//         else
//           fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
//         exit(1);
//       }
//     }
//   }
// 
//   if (!environment.inData.compare("")) {
//     cout << "The input data file is required (-i)\n";
//     exit(1);
//   }
// 
//   if (!existsTest(environment.inData)) {
//     cout << "The input file does not exist\n";
//     exit(1);
//   }
// 
//   if (!environment.outDir.compare("")) {
//     cout << "The output dir path is required (-o)\n";
//     exit(1);
//   }
// 
//   if (environment.steps.empty()) {
//     environment.steps.push_back(1);
//     environment.steps.push_back(2);
//   }
//   return true;
// }
//
//////////////////////////////////////////////////////////////////////////////////
// FRS 6 may 2010 END parseCommandLine is not used anymore
//////////////////////////////////////////////////////////////////////////////////

void printEnvironment(const Environment& environment) {
  stringstream s;
  s << "# --------\n# Inputs:\n# ----\n"
    << "# Input data file --> " << environment.inData << "\n"
    << "# Output directory --> " << environment.outDir << "\n"
    << "# All properties --> " << printNodesName(environment) << "\n"
    << "# Eff. N --> " << environment.effN << "\n"
    << "# Thres. Pc --> " << environment.thresPc << "\n"
    << "# N min --> " << environment.minN << "\n"
    << "# Clpx type --> " << environment.cplxType << "\n"
    << "# Clpx check --> " << environment.cplx << "\n"
    << "# Latent --> " << environment.isLatent << "\n"
    << "# Consistent --> " << environment.consistentPhase << "\n"
    << "# Reuse --> " << environment.isTplReuse << "\n"
    << "# K23 --> " << environment.isK23 << "\n"
    << "# propagation --> " << environment.isPropagation << "\n"
    << "# half V structures --> " << environment.halfVStructures << "\n"
    << "# Degeneracy --> " << environment.isDegeneracy << "\n"
    << "# No Init Eta --> " << environment.isNoInitEta << "\n"
    << "# Tau --> " << environment.tau << "\n"
    << "# VERSION --> " << environment.myVersion << "\n"
    << "# --------\n";
  cout << s.str();
}

void readFileType(Environment& environment) {
  uint numberGaussian = 0;
  for (uint pos = 0; pos < environment.numNodes; pos++) {
    environment.columnAsContinuous[pos] = environment.cntVarVec[pos];
    if (environment.columnAsContinuous[pos] == 2) {
      numberGaussian++;
      environment.columnAsContinuous[pos] = 1;
      environment.columnAsGaussian[pos] = 1;
    } else {
      environment.columnAsGaussian[pos] = 0;
    }
  }

  if (numberGaussian == environment.numNodes) environment.isAllGaussian = 1;

  if (numberGaussian >= 2) environment.atLeastTwoGaussian = 1;

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

  environment.isAllGaussian = 0;
  environment.atLeastTwoGaussian = 0;
  environment.atLeastTwoDiscrete = 0;
  environment.atLeastOneContinuous = 0;

  bool isNA = false;

  readData(environment, isNA);

  if (isNA) removeRowsAllNA(environment);

  environment.columnAsContinuous = new int[environment.numNodes];
  environment.columnAsGaussian = new int[environment.numNodes];
  if (environment.typeOfData == 0) {
    for (uint i = 0; i < environment.numNodes; i++) {
      environment.columnAsContinuous[i] = 0;
      environment.columnAsGaussian[i] = 0;
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

  //// Set a variables with all properties name and levels
  setNumberLevels(environment);

  // setProportions(environment);

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

  // Set the number of digits for the precision while using round( ..., digits =
  // ... ) Make sure the min levels for the data is 0
  environment.minN = 1;

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
  //
  // For temporal miic, only edges with a node of the last timestep are kept
  //
  if (environment.tau > 0)
    {
    for (uint i = 0; i < environment.numNodes; i++) 
      {
      bool first_node_not_lag0 = (environment.nodes[i].name.substr
        (environment.nodes[i].name.length() - 5, 5) != "_lag0");
      
      for (uint j = 0; j < environment.numNodes; j++) 
        {
        bool second_node_not_lag0 = (environment.nodes[j].name.substr
          (environment.nodes[j].name.length() - 5, 5) != "_lag0");
        
        if (first_node_not_lag0 && second_node_not_lag0)
          {
          environment.edges[i][j].status = 0;
          environment.edges[i][j].status_prev = 0;
          }
        }
      }

#if _DEBUG
    std::cout << "\nsetEnvironment:after status init\n\n";
    printAdjacencyMatrix (environment);
#endif
    }
  
  environment.noiseVec = new double[2 * environment.numSamples];
  for (uint i = 0; i < 2 * environment.numSamples; i++) {
    environment.noiseVec[i] =
        std::rand() / ((RAND_MAX + 1u) / MAGNITUDE_TIES) - MAGNITUDE_TIES / 2;
  }

  // for continuous gaussian
  if (environment.atLeastTwoGaussian == 1) {
    // set correlation part
    environment.means = new double[environment.numNodes];
    environment.standardDeviations = new double[environment.numNodes];

    // create rho matrix
    environment.rho = new double*[environment.numNodes];

    for (uint i = 0; i < environment.numNodes; i++) {
      environment.rho[i] = new double[environment.numNodes];
      for (uint j = 0; j < environment.numNodes; j++) {
        environment.rho[i][j] = 0;
      }
    }

    // create nsamples matrix
    environment.nSamples = new int*[environment.numNodes];

    for (uint i = 0; i < environment.numNodes; i++) {
      environment.nSamples[i] = new int[environment.numNodes];
    }

    for (uint i = 0; i < environment.numNodes; i++) {
      environment.means[i] = 0.0;
      environment.standardDeviations[i] = 0.0;
    }
    // compute the means and standard deviations
    computeMeansandStandardDeviations(environment);
    // compute the correlations coefficients
    computeCorrelations(environment);

    // Alloc the rho for correlations() (save time)
    environment.pMatrix = new double*[environment.numNodes];

    for (uint i = 0; i < environment.numNodes; i++) {
      environment.pMatrix[i] = new double[environment.numNodes];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
// FRS 6 may 2010 readFilesAndFillStructures is not used anymore
//////////////////////////////////////////////////////////////////////////////////
// 
// void readFilesAndFillStructures(
//     vector<string> edgesVectorOneLine, Environment& environment) {
//   setEnvironment(environment);
//   environment.oneLineMatrix =
//       new int[environment.numSamples * environment.numNodes];
//   for (uint i = 0; i < environment.numSamples; i++) {
//     for (uint j = 0; j < environment.numNodes; j++) {
//       environment.oneLineMatrix[j * environment.numSamples + i] =
//           environment.dataNumeric[i][j];
//     }
//   }
//   environment.edges = new Edge*[environment.numNodes];
// 
//   for (uint i = 0; i < environment.numNodes; i++)
//     environment.edges[i] = new Edge[environment.numNodes];
// 
//   for (uint i = 0; i < environment.numNodes; i++) {
//     environment.edges[i][i].status = 0;
//     environment.edges[i][i].status_init = 0;
//   }
// 
//   for (uint i = 0; i < environment.numNodes - 1; i++) {
//     for (uint j = i + 1; j < environment.numNodes; j++) {
//       environment.edges[i][j].shared_info = std::make_shared<EdgeSharedInfo>();
//       environment.edges[j][i].shared_info = environment.edges[i][j].shared_info;
//     }
//   }
// 
//   string lineData;
//   string s;
//   int posX = -1;
//   int posY = -1;
//   int numCols = 10;
// 
//   vector<vector<string>> vec;
//   vector<string> v;
// 
//   for (uint i = 0; i < edgesVectorOneLine.size(); i++) {
//     v.push_back(edgesVectorOneLine[i]);
//     if ((i + 1) % numCols == 0 && i != 0) {
//       vec.push_back(v);
//       v.clear();
//     }
//   }
// 
//   for (uint row = 0; row < vec.size(); row++) {
//     v = vec[row];
//     for (uint col = 0; col < v.size(); col++) {
//       string s = vec[row][col];
// 
//       if (col == 0) {
//         for (uint i = 0; i < environment.numNodes; i++)
//           if (environment.nodes[i].name.compare(s) == 0) posX = i;
// 
//       } else if (col == 1) {
//         for (uint i = 0; i < environment.numNodes; i++)
//           if (environment.nodes[i].name.compare(s) == 0) posY = i;
//       } else if (col == 2) {
//       } else if (col == 3) {
//         if (s.compare("NA") != 0) {
//           stringstream ss(s);  // Turn the string into a stream.
//           string tok;
//           char delimiter = ',';
//           while (getline(ss, tok, delimiter)) {
//             int ival;
//             for (uint i = 0; i < environment.numNodes; i++) {
//               if (environment.nodes[i].name.compare(tok) == 0) {
//                 ival = i;
//                 break;
//               }
//             }
//             environment.edges[posX][posY].shared_info->ui_vect_idx.push_back(
//                 ival);
//           }
//         }
//       } else if (col == 4) {
//         if (s.compare("NA") != 0) {
//           stringstream ss(s);  // Turn the string into a stream.
//           string tok;
//           char delimiter = ',';
//           while (getline(ss, tok, delimiter)) {
//             int ival;
//             for (uint i = 0; i < environment.numNodes; i++) {
//               if (environment.nodes[i].name.compare(tok) == 0) {
//                 ival = i;
//                 break;
//               }
//             }
//             environment.edges[posX][posY].shared_info->zi_vect_idx.push_back(
//                 ival);
//           }
//         }
//       } else if (col == 5) {
//         environment.edges[posX][posY].shared_info->Ixy_ui = atof(s.c_str());
//       } else if (col == 6) {
//         environment.edges[posX][posY].shared_info->cplx = atof(s.c_str());
//       } else if (col == 7) {
//         environment.edges[posX][posY].shared_info->Rxyz_ui = atof(s.c_str());
//       } else if (col == 8) {
//         int state = atoi(s.c_str());
//         environment.edges[posX][posY].shared_info->connected = state;
//         if (state == 1) {
//           environment.edges[posX][posY].status = 1;
//           environment.edges[posY][posX].status = 1;
//           // add the edge to Nomore
//           environment.noMoreAddress.emplace_back(new EdgeID(posX, posY));
//         } else {
//           environment.edges[posX][posY].status = 0;
//           environment.edges[posY][posX].status = 0;
//         }
//       } else if (col == 9) {
//         environment.edges[posX][posY].shared_info->Nxy_ui = atof(s.c_str());
//       }
//     }
//   }
//   environment.numNoMore = environment.noMoreAddress.size();
//   std::sort(environment.noMoreAddress.begin(), environment.noMoreAddress.end(),
//       sorterNoMore(environment));
// }
//////////////////////////////////////////////////////////////////////////////////
// FRS 6 may 2010 END readFilesAndFillStructures is not used anymore
//////////////////////////////////////////////////////////////////////////////////

static void chkIntFn(void* dummy) { R_CheckUserInterrupt(); }

bool checkInterrupt(bool check /*=true*/) {
  if (check)
    return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
  else
    return false;
}

int printProgress(
    double percentage, double startTime, string outdir, int prg_numSearchMore) {
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
    KDTreeVectorOfVectorsAdaptor<vector<vector<double>>, double>;
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
double compute_kl_divergence_continuous(vector<vector<double>> space1,
    vector<vector<double>> space2, int n1, int n2, int ndims, int k,
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
    int samplesNotNA, int** dataNumeric_red, int* AllLevels_red,
    int* samplesToEvaluate) {
  int current_samplesNotNA =
      getNumSamples_nonNA(environment, posArray[0], posArray[1]);
  double kldiv = 0;
  // 1 - XY discrete
  if (environment.columnAsContinuous[posArray[0]] == 0 &&
      environment.columnAsContinuous[posArray[1]] == 0) {
    // XY discrete
    // Retrieve marginal distibutions with the current conditioning Us
    int current_samplesNotNA =
        getNumSamples_nonNA(environment, posArray[0], posArray[1]);
    double** jointFreqs = getJointFreqs(
        environment, posArray[0], posArray[1], current_samplesNotNA);

    double** freqs2 = new double*[environment.allLevels[posArray[0]]];
    for (uint j = 0; j < environment.allLevels[posArray[0]]; j++) {
      freqs2[j] = new double[environment.allLevels[posArray[1]]];
      for (uint k = 0; k < environment.allLevels[posArray[1]]; k++) {
        freqs2[j][k] = 0;
      }
    }

    // fill table
    for (int k = 0; k < samplesNotNA; k++) {
      freqs2[dataNumeric_red[0][k]][dataNumeric_red[1][k]]++;
    }
    for (uint j = 0; j < environment.allLevels[posArray[0]]; j++) {
      for (uint k = 0; k < environment.allLevels[posArray[1]]; k++) {
        freqs2[j][k] = 1.0 * freqs2[j][k] / samplesNotNA;
      }
    }

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
    vector<vector<double>> joint_base(current_samplesNotNA, vector<double>(2));
    int* curr_samplesToEvaluate = new int[environment.numSamples];
    getJointSpace(environment, posArray[0], posArray[1], joint_base,
        curr_samplesToEvaluate);

    int* map_samples = new int[current_samplesNotNA];
    int i_map = 0;
    for (uint i = 0; i < environment.numSamples; i++) {
      if (curr_samplesToEvaluate[i] == 1) {  // sample i is present in X;Y|U
        map_samples[i_map] = 0;
        if (samplesToEvaluate[i] == 1) {  // sample i is also present in X;Y|U,Z
          map_samples[i_map] = 1;
        }
        i_map++;
      }
    }
    vector<vector<double>> joint_nonNA(samplesNotNA, vector<double>(2));
    int i_nonNA = 0;
    for (uint i = 0; i < environment.numSamples; i++) {
      if (samplesToEvaluate[i] == 1) {
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

    delete[] curr_samplesToEvaluate;
    delete[] map_samples;
    delete[] flag_break_ties;
  } else {
    // 3 - One discrete and one continuous
    int discrete_pos;
    int continuous_pos;
    int discrete_pos_binary;
    int continuous_pos_binary;
    if (environment.columnAsContinuous[posArray[0]] == 0) {
      discrete_pos = posArray[0];
      continuous_pos = posArray[1];
      discrete_pos_binary = 0;
      continuous_pos_binary = 1;
    } else {
      discrete_pos = posArray[1];
      continuous_pos = posArray[0];
      discrete_pos_binary = 1;
      continuous_pos_binary = 0;
    }
    int n_discrete_levels = environment.allLevels[discrete_pos];

    // Retrieve marginal distibutions with the current conditioning Us
    int current_samplesNotNA =
        getNumSamples_nonNA(environment, posArray[0], posArray[1]);
    int* mixedDiscrete = new int[current_samplesNotNA];
    double* mixedContinuous = new double[current_samplesNotNA];
    int* curr_samplesToEvaluate = new int[environment.numSamples];
    getJointMixed(environment, posArray[0], posArray[1], mixedDiscrete,
        mixedContinuous, curr_samplesToEvaluate);

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
    for (int i = 0; i < samplesNotNA; i++) {
      count_nonNA[dataNumeric_red[discrete_pos_binary][i]]++;
    }

    // Compute the sum count(y) * KL(X_nonNA|y || X|y) over all values of Y y
    for (int level = 0; level < n_discrete_levels; level++) {
      int* map_level = new int[count_base[level]];
      int i_level = 0;
      for (uint i = 0; i < environment.numSamples; i++) {
        if (environment.dataNumeric[i][discrete_pos] == level) {
          if (curr_samplesToEvaluate[i]) {
            map_level[i_level] = 0;
            if (samplesToEvaluate[i] == 1) map_level[i_level] = 1;
            i_level++;
          }
        }
      }

      vector<vector<double>> continuous_base(
          count_base[level], vector<double>(1));
      i_level = 0;
      for (int i = 0; i < current_samplesNotNA; i++) {
        if (mixedDiscrete[i] == level) {
          continuous_base[i_level][0] = mixedContinuous[i];
          i_level++;
        }
      }

      vector<vector<double>> continuous_nonNA(
          count_nonNA[level], vector<double>(1));
      int i_level_nonNA = 0;
      for (uint i = 0; i < environment.numSamples; i++) {
        if (samplesToEvaluate[i] == 1 &&
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
    delete[] curr_samplesToEvaluate;
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

}  // namespace utility
}  // namespace miic
