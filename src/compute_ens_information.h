#ifndef MIIC_COMPUTE_ENS_INFORMATION_H_
#define MIIC_COMPUTE_ENS_INFORMATION_H_

#include "structure.h"

namespace miic {
namespace computation {

void computeContributingScores(structure::Environment&, int* ziContPosIdx,
    int iz, int* myZi, int myNbrUi, unsigned int numSamples_nonNA,
    int* posArray, double* scoresZ, structure::MemorySpace m);
double* computeDifference(double*, double*, int);
double* computeDifferenceByStep(double*, double*, int, int);
double* computeSum(double*, double*, int);
double* computeEnsInformation(
    structure::Environment&, int*, int, int*, int, int, int, int, int);
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
// Gaussian case
double* corrMutInfo(structure::Environment&, double** dataset, int*, int, int*,
    int, int, int, int);
void SearchForNewContributingNodeAndItsRankGaussian(
    structure::Environment& environment, const int posX, const int posY,
    structure::MemorySpace&);
double computeEnsInformationContinuous_Gaussian(
    structure::Environment& environment, const int posX, const int posY,
    const int posZ);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_ENS_INFORMATION_H_
