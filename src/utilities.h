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
void sort2arrays(int len, structure::TempVector<int>& a,
    structure::TempVector<int>& brr, structure::TempVector<int>& bridge);
double ramanujan(int n);
TimePoint getLapStartTime();
double getLapInterval(TimePoint);
void printProgress(double percentage, TimePoint, int& n_unsettled);
// KL divergence functions
double compute_kl_divergence(int X, int Y, structure::Environment& environment,
    int samplesNotNA, const structure::TempVector<int>& AllLevels_red,
    const structure::TempVector<int>& sample_is_not_NA);
double compute_kl_divergence_continuous(
    std::vector<std::vector<double>>& space1,
    std::vector<std::vector<double>>& space2, int n1, int n2, int ndims, int k,
    bool* flag_break_ties, int* map_samples, double* noise_vec);

double kl(const structure::TempGrid2d<int>& counts1,
    const structure::TempGrid2d<double>& freqs2);

void getJointMixed(const structure::Environment&, int i, int j,
    int* mixedDiscrete, double* mixedContinuous, int* curr_sample_is_not_NA);
structure::TempGrid2d<double> getJointFreqs(const structure::Environment&,
    int i, int j,
    const structure::TempVector<int>& sample_is_not_NA =
        structure::TempVector<int>());
void getJointSpace(const structure::Environment&, int i, int j,
    double** jointSpace, int* curr_sample_is_not_NA);
int getNumSamplesNonNA(const structure::Environment&, int i, int j);

int count_non_NAs(int X, int Y, const std::vector<int>& ui_list,
    structure::TempVector<int>& sample_is_not_NA,
    structure::TempVector<int>& NAs_count, structure::Environment& environment,
    int z = -1);

bool filter_NAs(int X, int Y, const std::vector<int>& ui_list,
    structure::TempVector<int>& AllLevels, structure::TempVector<int>& cnt,
    structure::TempVector<int>& posArray_red,
    structure::TempGrid2d<int>& dataNumeric,
    structure::TempGrid2d<int>& dataNumericIdx,
    structure::TempVector<double>& sample_weights,
    const structure::TempVector<int>& sample_is_not_NA,
    const structure::TempVector<int>& NAs_count,
    structure::Environment& environment, int z = -1);

bool checkInterrupt(bool check = true);

double lookupScore(int X, int Y, const std::vector<int>& ui_list, int z,
    structure::Environment& environment);
void lookupScore(int X, int Y, const std::vector<int>& ui_list, int z,
    double* score, structure::Environment& environment);
void saveScore(int X, int Y, const std::vector<int>& ui_list, int z,
    double score, structure::Environment& environment);
void saveScore(int X, int Y, const std::vector<int>& ui_list, int z,
    double* score, structure::Environment& environment);

bool SampleHasNoNA(const structure::Environment& env, int row, int i, int j);

}  // namespace utility
}  // namespace miic

#endif  // MIIC_UTILITIES_H_
