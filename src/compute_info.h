#ifndef MIIC_COMPUTE_INFO_H_
#define MIIC_COMPUTE_INFO_H_

#include "computation_cache.h"
#include "structure.h"

namespace miic {
namespace computation {

double* computeCondMutualInfoDiscrete(const structure::TempGrid2d<int>& data,
    const structure::TempVector<int>& levels,
    const structure::TempVector<int>& var_idx,
    const structure::TempVector<double>& weights, int cplx,
    std::shared_ptr<CtermCache> cache);
double* getAllInfoNEW(const std::vector<std::vector<int>>& data,
    const std::vector<int>& ptrAllLevels, int id_x, int id_y, int id_z,
    const std::vector<int>& ui_list, int sampleSizeEff, int modCplx, int k23,
    const std::vector<double>& weights, std::shared_ptr<CtermCache> cache);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_INFO_H_
