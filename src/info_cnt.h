// Module to compute the conditional mutual information for mixed continuous
// and discrete variables.
#ifndef MIIC_INFO_CNT_H_
#define MIIC_INFO_CNT_H_

#include <algorithm>
#include <functional>
#include <vector>

#include "structure.h"

namespace miic {
namespace computation {

double* compute_Ixy_alg1(int** data, int** sortidx, int* ptr_cnt,
    int* ptrVarIdx, int* AllLevels, int n, structure::Environment& environment,
    bool saveIterations = false);

double* compute_mi_cond_alg1(int** data, int** sortidx, int* AllLevels,
    int* ptr_cnt, int* ptrVarIdx, int nbrUi, int n,
    structure::Environment& environment, bool saveIterations = false);

double* compute_Rscore_Ixyz_alg5(int** data, int** sortidx, int* AllLevels,
    int* ptr_cnt, int* ptrVarIdx, int nbrUi, int ptrZiIdx, int n,
    structure::Environment& environment, bool saveIterations = false);

double* compute_Rscore_Ixyz_new_alg5(int** data, int** sortidx, int* AllLevels,
    int* ptr_cnt, int* ptrVarIdx, int nbrUi, int ptrZiIdx, int n,
    structure::Environment& environment, bool saveIterations = false);

void optfun_onerun_kmdl_coarse(int* sortidx_var, int* data, int nbrV,
    int** factors, int* r, double sc, int sc_levels1, int previous_levels,
    int n, int nnr, int* cut, int* r_opt, structure::Environment& environment);

void reset_u_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx,
    int init_nbin, int maxbins, int lbin, int* r, int* AllLevels, int n);

void reset_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx,
    int init_nbin, int maxbins, int lbin, int* r, int* AllLevels, int n);

namespace computation_impl {

using std::size_t;
using std::vector;
template <typename T>
struct Grid2d {
 private:
  size_t rows_, cols_;
  vector<T> data_vec_;

 public:
  Grid2d(size_t rows, size_t cols)
      : rows_(rows), cols_(cols), data_vec_(vector<T>(rows * cols)) {}

  Grid2d(size_t rows, size_t cols, T&& init)
      : rows_(rows), cols_(cols), data_vec_(vector<T>(rows * cols, init)) {}

  Grid2d(const Grid2d&) = delete;
  Grid2d& operator=(const Grid2d&) = delete;

  T& operator()(size_t row, size_t col) { return data_vec_[row * cols_ + col]; }
  const T& operator()(size_t row, size_t col) const {
    return data_vec_[row * cols_ + col];
  }

  void add_row(vector<T>& counts, size_t row) {
    std::transform(counts.begin(), counts.end(),
        data_vec_.begin() + row * cols_, counts.begin(), std::plus<T>());
  }

  void subtract_row(vector<T>& counts, size_t row) {
    std::transform(counts.begin(), counts.end(),
        data_vec_.begin() + row * cols_, counts.begin(), std::minus<T>());
  }

  auto begin() { return data_vec_.begin(); }
  auto end() { return data_vec_.end(); }
  auto cbegin() const { return data_vec_.cbegin(); }
  auto cend() const { return data_vec_.cend(); }
};

}  // namespace computation_impl
using computation_impl::Grid2d;
}  // namespace computation
}  // namespace miic

#endif  // MIIC_INFO_CNT_H_
