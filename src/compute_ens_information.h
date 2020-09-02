#ifndef MIIC_COMPUTE_ENS_INFORMATION_H_
#define MIIC_COMPUTE_ENS_INFORMATION_H_

#include "environment.h"
#include "structure.h"

namespace miic {
namespace computation {

double* computeEnsInformationNew(structure::Environment&, int X, int Y,
    const std::vector<int>& ui_list, const structure::TempVector<int>& zi_list,
    int cplx);
void SearchForNewContributingNodeAndItsRank(structure::Environment&, int, int);
// Continuous data
double* computeEnsInformationContinuous(structure::Environment&, int X, int Y,
    const std::vector<int>& ui_list, const std::vector<int>& zi_list, int cplx);
double get3PointInfo(structure::Environment&, int X, int Y, int Z,
    const std::vector<int>& ui_list);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_ENS_INFORMATION_H_
