#ifndef MIIC_UTILITIES_H_
#define MIIC_UTILITIES_H_
#include <chrono>
#include <string>
#include <vector>

#include "structure.h"

namespace miic {
namespace utility {
using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

inline TimePoint getLapStartTime() { return std::chrono::steady_clock::now(); }
inline double getLapInterval(TimePoint start_time) {
  using second = std::chrono::duration<double>;
  return second(std::chrono::steady_clock::now() - start_time).count();
}

std::vector<std::vector<int>> getAdjMatrix(
    const structure::Grid2d<structure::Edge>&);

std::string toNameString(
    const std::vector<structure::Node>&, const std::vector<int>&);

std::vector<std::vector<std::string>> getEdgesInfoTable(
    const structure::Grid2d<structure::Edge>&,
    const std::vector<structure::Node>&);

bool checkInterrupt();

void printProgress(double percentage, TimePoint, int& n_unsettled);

double compute_kl_divergence(const structure::Grid2d<int>& data_numeric,
    const structure::Grid2d<double>& data_double, int X, int Y,
    const std::vector<int>& ui_list, const std::vector<int>& levels,
    const std::vector<int>& is_continuous, int samplesNotNA,
    const structure::TempVector<int>& AllLevels_red,
    const structure::TempVector<int>& sample_is_not_NA,
    const std::vector<double>& noise_vec);

int getNumSamplesNonNA(const structure::Grid2d<int>& data_numeric, int X, int Y,
    const std::vector<int>& ui_list);

int countNonNA(int X, int Y, int Z, const std::vector<int>& ui_list,
    const structure::Grid2d<int>& data_numeric,
    structure::TempVector<int>& sample_is_not_NA,
    structure::TempVector<int>& NAs_count);

bool filterNA(int X, int Y, int Z, const std::vector<int>& ui_list,
    const structure::Grid2d<int>& data_numeric,
    const structure::Grid2d<int>& data_numeric_idx,
    const std::vector<int>& levels, const std::vector<int>& is_continuous,
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
