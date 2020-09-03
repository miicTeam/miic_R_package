#ifndef MIIC_COMPUTE_ENS_INFORMATION_H_
#define MIIC_COMPUTE_ENS_INFORMATION_H_

#include "environment.h"
#include "structure.h"

namespace miic {
namespace computation {

double* computeEnsInformationNew(
    structure::Environment&, int X, int Y, const std::vector<int>& ui_list);
void SearchForNewContributingNodeAndItsRank(structure::Environment&, int, int);
// Continuous data
double* computeEnsInformationContinuous(structure::Environment&, int X, int Y,
    const std::vector<int>& ui_list);
double getInfo3PointOrScore(structure::Environment&, int X, int Y, int Z,
    const std::vector<int>& ui_list, bool get_info);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_ENS_INFORMATION_H_
