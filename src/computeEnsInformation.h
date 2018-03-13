#ifndef COMPUTEENSINFORMATION_H
#define COMPUTEENSINFORMATION_H

double* getFromArray(double*, int, int);
double* computeDifference(double*, double*, int);
double* computeDifferenceByStep(double*, double*, int, int);
double* computeDifference(double*, double*, int);
double* computeSum(double*, double*, int);
//double* computeEnsInformation(Environment, int*, int, int*, int, int, const int, const int, int);
double* computeEnsInformationNew(Environment&, int*, int, int*, int, int, const int, const int, int);
double* computeEnsInformationNewThread(Environment&, int*, int, int*, int, int, const int, const int, int, MemorySpace);
bool SearchForNewContributingNodeAndItsRank(Environment&, const int, const int);
void* SearchForNewContributingNodeAndItsRankThread(void*);

#endif
