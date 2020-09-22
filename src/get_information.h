#ifndef MIIC_COMPUTE_ENS_INFORMATION_H_
#define MIIC_COMPUTE_ENS_INFORMATION_H_

#include "environment.h"
#include "structure.h"

namespace miic {
namespace computation {

structure::InfoBlock getCondMutualInfo(int X, int Y,
    const std::vector<int>& ui_list, const structure::Grid2d<int>& data_numeric,
    const structure::Grid2d<int>& data_numeric_idx, structure::Environment&);
double getInfo3PointOrScore(structure::Environment&, int X, int Y, int Z,
    const std::vector<int>& ui_list, bool get_info);
void searchForBestContributingNode(
    structure::Environment&, int X, int Y, bool parallelizable);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_ENS_INFORMATION_H_
