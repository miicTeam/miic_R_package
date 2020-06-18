#ifndef MIIC_COMPUTE_INFO_H_
#define MIIC_COMPUTE_INFO_H_

#include "structure.h"

namespace miic {
namespace computation {

double computeLogC(int N, int r, double* looklog, double* c2terms);
double computeLogC(int N, int r, double* looklog, double** cterms);
double* getAllInfoNEW(int* ptrAllData, unsigned int* ptrAllLevels,
    const std::vector<int>& ptrVarIdx, int nbrUi, int* ptrZiIdx, int nbrZi,
    int ziPos, int sampleSize, int sampleSizeEff, int modCplx, int k23,
    double* looklog, double* c2terms, structure::MemorySpace* memory,
    double* weights, double** freqs1, bool testDistribution);
double lnfactorial(int n, double* looklog);
double logchoose(int n, int k, double* looklog);
double logchoose(int n, int k, double* looklog, double** lookchoose);

}  // namespace computation
}  // namespace miic

#endif  // MIIC_COMPUTE_INFO_H_
