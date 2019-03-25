#ifndef _COMPUTEINFO_H_
#define _COMPUTEINFO_H_

int getDVect( int* myPtrAllLevels, int* myPtrVarIdx, int myNbrAllVar, int* ioPtrDVectVar );
int computeInfoAndCplx(  int *myNxyui, int *myDVect, int *myLevels, int *myVarIdx, int mySSize, int mySSizeEff, int myNbrVar, double* infoAndCplx, int myCplxType, double* );
double computeLogC( int N,  int r, double* looklog, double* c2terms);
double computeLogC( int N,  int r, double* looklog, double** cterms);
double computeLogHDC(int N, int r, int reff, double* looklog, double** cterms);
double compute_rescaledLogC( int sampleSize, int sampleSizeEff, int N, int r, double* );
double compute_LogC_C2( int N, int r, double* c2terms);
unsigned long binomialCoeff(int n, int k);
double logchoose(int n, int k, double* looklog);

#endif /* _COMPUTEINFO_H_ */
