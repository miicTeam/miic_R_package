#ifndef MIIC_COMPUTE_INFO_H_
#define MIIC_COMPUTE_INFO_H_

#include "computation_cache.h"
#include "structure.h"

namespace miic {
namespace computation {

structure::InfoBlock computeCondMutualInfoDiscrete(
    const structure::TempGrid2d<int>& data,
    const structure::TempVector<int>& r_list,
    const structure::TempVector<int>& var_idx,
    const structure::TempVector<double>& weights, int cplx,
    bool negative_info, std::shared_ptr<CtermCache> cache);
structure::Info3PointBlock computeInfo3PointAndScoreDiscrete(
    const structure::TempGrid2d<int>& data,
    const structure::TempVector<int>& r_list,
    const structure::TempVector<int>& var_idx,
    const structure::TempVector<double>& weights, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_INFO_H_
