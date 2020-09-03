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
TimePoint getLapStartTime();
double getLapInterval(TimePoint);
void printProgress(double percentage, TimePoint, int& n_unsettled);
// KL divergence functions
double compute_kl_divergence(int X, int Y, structure::Environment& environment,
    int samplesNotNA, const structure::TempVector<int>& AllLevels_red,
    const structure::TempVector<int>& sample_is_not_NA);

double kl(const structure::TempGrid2d<int>& counts1,
    const structure::TempGrid2d<double>& freqs2);

structure::TempGrid2d<double> getJointFreqs(const structure::Environment&,
    int i, int j,
    const structure::TempVector<int>& sample_is_not_NA =
        structure::TempVector<int>());
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

bool checkInterrupt();

}  // namespace utility
}  // namespace miic

#endif  // MIIC_UTILITIES_H_
