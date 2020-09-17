// Module to compute the conditional mutual information for mixed continuous
// and discrete variables.
#ifndef MIIC_INFO_CNT_H_
#define MIIC_INFO_CNT_H_

#include "computation_cache.h"
#include "structure.h"

namespace miic {
namespace computation {

structure::InfoBlock computeCondMutualInfo(
    const structure::TempGrid2d<int>& data,
    const structure::TempGrid2d<int>& data_idx,
    const structure::TempVector<int>& levels,
    const structure::TempVector<int>& is_continuous,
    const structure::TempVector<int>& var_idx,
    const structure::TempVector<double>& sample_weights,
    bool flag_sample_weights, int initbins, int maxbins, int cplx,
    std::shared_ptr<CtermCache>,
    std::shared_ptr<structure::CutPointsInfo> = nullptr);

structure::Info3PointBlock computeInfo3PointAndScore(
    const structure::TempGrid2d<int>& data,
    const structure::TempGrid2d<int>& data_idx,
    const structure::TempVector<int>& levels,
    const structure::TempVector<int>& is_continuous,
    const structure::TempVector<int>& var_idx,
    const structure::TempVector<double>& sample_weights,
    bool flag_sample_weights, int initbins, int maxbins, int cplx,
    std::shared_ptr<CtermCache> cache);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_INFO_CNT_H_
