#include <string>
#include <map>
#include <tuple>
#include <algorithm>

#include "utilities.h"

using namespace std;


//////////////////////////////////////////////////////////////////////////////////
double* compute_Ixy_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels,
						 int n, Environment& environment, bool saveIterations = false);

double* compute_mi_cond_alg1(int** data, int** sortidx, int* AllLevels, int* ptr_cnt, int* ptrVarIdx,
							 int nbrUi, int n, Environment& environment, bool saveIterations = false);

double* compute_Rscore_Ixyz_alg5(int** data, int** sortidx, int* AllLevels, int* ptr_cnt, int* ptrVarIdx,
								 int nbrUi, int ptrZiIdx, int n, Environment& environment, bool saveIterations = false);

double* compute_Rscore_Ixyz_new_alg5(int** data, int** sortidx, int* AllLevels, int* ptr_cnt, int* ptrVarIdx,
									 int nbrUi, int ptrZiIdx, int n, Environment& environment, bool saveIterations = false);


template<typename T>
struct Grid2d {
private:
    std::size_t rows_, cols_;
    std::vector<T> data_vec_;

public:
    Grid2d(std::size_t rows, std::size_t cols):
        rows_(rows), cols_(cols),
        data_vec_(std::vector<T>(rows * cols))
    {}

    Grid2d(std::size_t rows, std::size_t cols, T&& init):
        rows_(rows), cols_(cols),
        data_vec_(std::vector<T>(rows * cols, init))
    {}

    Grid2d(const Grid2d&) = delete;
    Grid2d& operator=(const Grid2d&) = delete;

    T& operator()(size_t row, size_t col) { return data_vec_[row * cols_ + col]; }
    const T& operator()(size_t row, size_t col) const { return data_vec_[row * cols_ + col]; }

    void add_row(vector<T>& counts, size_t row) {
        std::transform(counts.begin(), counts.end(),
                       data_vec_.begin() + row * cols_, counts.begin(),
                       std::plus<T>());
    }

    void subtract_row(vector<T>& counts, size_t row) {
        std::transform(counts.begin(), counts.end(),
                       data_vec_.begin() + row * cols_, counts.begin(),
                       std::minus<T>());
    }

    auto begin() { return data_vec_.begin(); }
    auto end() { return data_vec_.end(); }
    auto cbegin() const { return data_vec_.cbegin(); }
    auto cend() const { return data_vec_.cend(); }
};
