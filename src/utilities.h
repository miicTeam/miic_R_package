#ifndef MIIC_UTILITIES_H_
#define MIIC_UTILITIES_H_

#include "environment.h"
#include "structure.h"

namespace miic {
namespace utility {

void createMemorySpace(structure::Environment&, structure::MemorySpace&);
void deleteMemorySpace(structure::Environment&, structure::MemorySpace&);
std::vector<std::vector<std::string>> getEdgesInfoTable(
    structure::Environment&);
std::string toNameString(
    const structure::Environment&, const std::vector<int>&);
std::vector<std::vector<int>> getAdjMatrix(const structure::Environment&);
int sign(double val);
void sort2arraysConfidence(int len, int a[], int brr[]);
void sort2arrays(int len, int a[], int brr[], int bridge[]);
double get_wall_time();
double ramanujan(int n);
int printProgress(double percentage, double startTime, int prg_numSearchMore);
// KL divergence functions
double compute_kl_divergence(const std::vector<int>& posArray,
    structure::Environment& environment, int samplesNotNA,
    const std::vector<int>& AllLevels_red,
    const std::vector<int>& sample_is_not_NA);
double compute_kl_divergence_continuous(
    std::vector<std::vector<double>>& space1,
    std::vector<std::vector<double>>& space2, int n1, int n2, int ndims, int k,
    bool* flag_break_ties, int* map_samples, double* noise_vec);

double kl(double** freqs1, double** freqs2, int nrows, int ncols);
double kl(int** counts1, double** freqs2, int nrows, int ncols);
double kl(const std::vector<int>& freqs1, const std::vector<int>& freqs2,
    int n1, int n2);

void getJointMixed(const structure::Environment&, int i, int j,
    int* mixedDiscrete, double* mixedContinuous, int* curr_sample_is_not_NA);
double** getJointFreqs(const structure::Environment&, int i, int j,
    const std::vector<int>& sample_is_not_NA = std::vector<int>());
void getJointSpace(const structure::Environment&, int i, int j,
    double** jointSpace, int* curr_sample_is_not_NA);
int getNumSamplesNonNA(const structure::Environment&, int i, int j);

int count_non_NAs(int nbrUi, std::vector<int> &sample_is_not_NA,
    std::vector<int> &NAs_count, const std::vector<int>& posArray,
    structure::Environment& environment, int z=-1);

bool filter_NAs(int nbrUi, std::vector<int> &AllLevels, std::vector<int> &cnt,
    std::vector<int> &posArray_red, const std::vector<int>& posArray,
    std::vector<std::vector<int> > &dataNumeric,
    std::vector<std::vector<int> > &dataNumericIdx,
    std::vector<double> &sample_weights,
    const std::vector<int> &sample_is_not_NA,
    const std::vector<int> &NAs_count,
    structure::Environment& environment, int z=-1);

bool checkInterrupt(bool check = true);

double lookupScore(const std::vector<int> &posArray, int nbrUi, int z,
  structure::Environment& environment);
void lookupScore(const std::vector<int>& posArray, int nbrUi, int z,
    double* score, structure::Environment& environment);
void saveScore(const std::vector<int>& posArray, int nbrUi, int z, double score,
    structure::Environment& environment);
void saveScore(const std::vector<int>& posArray, int nbrUi, int z,
    double* score, structure::Environment& environment);

bool SampleHasNoNA(const structure::Environment& env, int row, int i, int j);

class EdgeSorter {
  const structure::Environment& env;

 public:
  EdgeSorter(const structure::Environment& env) : env(env) {}
  bool operator()(
      const structure::EdgeID& e1, const structure::EdgeID& e2) const {
    const auto info1 = env.edges[e1.i][e1.j].shared_info;
    const auto info2 = env.edges[e2.i][e2.j].shared_info;
    // connected can be 1 or 0
    if (info1->connected != info2->connected)
      return info1->connected > info2->connected;

    if (info1->connected == 0) {
      if (info1->Rxyz_ui == 0 || info2->Rxyz_ui == 0)
        return info2->Rxyz_ui != 0;
      else
        return info1->Rxyz_ui > info2->Rxyz_ui;
    } else {
      return info1->Ixy_ui > info2->Ixy_ui;
    }
  }
};

}  // namespace utility
}  // namespace miic

#endif  // MIIC_UTILITIES_H_
