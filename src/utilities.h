#ifndef MIIC_UTILITIES_H_
#define MIIC_UTILITIES_H_
#include <chrono>

#include "environment.h"
#include "structure.h"

namespace miic {
namespace utility {
using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

inline TimePoint getLapStartTime() { return std::chrono::steady_clock::now(); }
inline double getLapInterval(TimePoint start_time) {
  using second = std::chrono::duration<double>;
  return second(std::chrono::steady_clock::now() - start_time).count();
}

std::vector<std::vector<int>> getAdjMatrix(const structure::Environment&);

std::string toNameString(
    const structure::Environment&, const std::vector<int>&);

std::vector<std::vector<std::string>> getEdgesInfoTable(
    const structure::Environment&);

bool checkInterrupt();

void printProgress(double percentage, TimePoint, int& n_unsettled);

double compute_kl_divergence(int X, int Y, structure::Environment& environment,
    int samplesNotNA, const structure::TempVector<int>& AllLevels_red,
    const structure::TempVector<int>& sample_is_not_NA);

int getNumSamplesNonNA(const structure::Environment&, int i, int j);

int countNonNA(int X, int Y, int Z, const std::vector<int>& ui_list,
    const std::vector<std::vector<int>>& data_numeric,
    structure::TempVector<int>& sample_is_not_NA,
    structure::TempVector<int>& NAs_count);

bool filterNA(int X, int Y, int Z, const std::vector<int>& ui_list,
    const std::vector<std::vector<int>>& data_numeric,
    const std::vector<std::vector<int>>& data_numeric_idx,
    const std::vector<int>& levels,
    const std::vector<int>& is_continuous,
    const std::vector<double>& sample_weights,
    const structure::TempVector<int>& sample_is_not_NA,
    const structure::TempVector<int>& NAs_count,
    structure::TempGrid2d<int>& data_numeric_red,
    structure::TempGrid2d<int>& data_numeric_idx_red,
    structure::TempVector<int>& levels_red,
    structure::TempVector<int>& is_continuous_red,
    structure::TempVector<int>& posArray_red,
    structure::TempVector<double>& sample_weights_red, bool any_na);

size_t getLinearAllocatorSize(int n_samples, int n_nodes, int maxbins,
    int initbins, const std::vector<int>& is_continuous,
    const std::vector<int>& levels);

}  // namespace utility
}  // namespace miic

#endif  // MIIC_UTILITIES_H_
