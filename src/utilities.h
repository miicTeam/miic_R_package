#ifndef MIIC_UTILITIES_H_
#define MIIC_UTILITIES_H_
#include <chrono>

#include "environment.h"
#include "structure.h"

namespace miic {
namespace utility {
using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

std::vector<std::vector<std::string>> getEdgesInfoTable(
    const structure::Environment&);
std::string toNameString(
    const structure::Environment&, const std::vector<int>&);
std::vector<std::vector<int>> getAdjMatrix(const structure::Environment&);
void sort2arraysConfidence(
    int len, const std::vector<int>& a, std::vector<int>& brr);
void sort2arrays(int len, structure::TempVector<int>& a,
    structure::TempVector<int>& brr, structure::TempVector<int>& bridge);
double ramanujan(int n);
TimePoint getLapStartTime();
double getLapInterval(TimePoint);
void printProgress(double percentage, TimePoint, int& n_unsettled);
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

double lookupScore(const std::vector<int>& posArray, int nbrUi, int z,
    structure::Environment& environment);
void lookupScore(const std::vector<int>& posArray, int nbrUi, int z,
    double* score, structure::Environment& environment);
void saveScore(const std::vector<int>& posArray, int nbrUi, int z, double score,
    structure::Environment& environment);
void saveScore(const std::vector<int>& posArray, int nbrUi, int z,
    double* score, structure::Environment& environment);

bool SampleHasNoNA(const structure::Environment& env, int row, int i, int j);

}  // namespace utility
}  // namespace miic

#endif  // MIIC_UTILITIES_H_
