#ifndef COMPUTEENSINFORMATION_H
#define COMPUTEENSINFORMATION_H

double* getFromArray(double*, int, int);
double* computeDifference(double*, double*, int);
double* computeDifferenceByStep(double*, double*, int, int);
double* computeDifference(double*, double*, int);
double* computeSum(double*, double*, int);
double* computeEnsInformation(Environment&, int*, int, int*, int, int, const int, const int, int);
double* computeEnsInformationNew(Environment&, int*, int, int*, int, int, const int, const int, int, MemorySpace&);
void SearchForNewContributingNodeAndItsRank(Environment&, const int, const int, MemorySpace&);

//Continuous data
double* computeEnsInformationContinuous(Environment&, int*, int, int*, int, int, const int, const int, int, MemorySpace&);
double* computeEnsInformationContinuous_Orientation(Environment& environment, int* myCond, int myNbrUi,  int* myZi, 
													const int myVarIdxX, const int myVarIdxY, const int cplx, MemorySpace&);

//Gaussian case
double* corrMutInfo(Environment&, double** dataset, int*,int,int*,int, int, int, int);
void SearchForNewContributingNodeAndItsRankGaussian(Environment& environment, const int posX, const int posY, MemorySpace&);
double computeEnsInformationContinuous_Gaussian(Environment& environment, const int posX, const int posY,const int posZ);

#endif
