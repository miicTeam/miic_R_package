#ifndef MIIC_COMPUTE_ENS_INFORMATION_H_
#define MIIC_COMPUTE_ENS_INFORMATION_H_

#include "structure.h"

namespace miic {
namespace computation {

void computeContributingScores(structure::Environment&, int* ziContPosIdx,
    int iz, int* myZi, int myNbrUi, unsigned int numSamples_nonNA,
    const std::vector<int>& posArray, double* scoresZ,
    structure::MemorySpace m);
double* computeEnsInformationNew(structure::Environment&, int*, int, int*, int,
    int, int, int, int, structure::MemorySpace&);
void SearchForNewContributingNodeAndItsRank(
    structure::Environment&, int, int, structure::MemorySpace&);
// Continuous data
double* computeEnsInformationContinuous(structure::Environment&, int* myCond,
    int myNbrUi, int* myZi, unsigned int myNbrZi, int myZiPos, int myVarIdxX,
    int myVarIdxY, int cplx, structure::MemorySpace&);
double* computeEnsInformationContinuous_Orientation(structure::Environment&,
    int* myCond, int myNbrUi, int* myZi, int myVarIdxX, int myVarIdxY, int cplx,
    structure::MemorySpace&);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_ENS_INFORMATION_H_
