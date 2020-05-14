#ifndef MIIC_UTILITIES_H_
#define MIIC_UTILITIES_H_

#include "structure.h"

namespace miic {
namespace utility {

void createMemorySpace(structure::Environment&, structure::MemorySpace&);
void deleteMemorySpace(structure::Environment&, structure::MemorySpace&);
void deleteStruct(structure::Environment&);
std::string printNodesName(const structure::Environment&);
void printEdges (structure::Environment&, bool filter_status=true, bool half_only=true);
void printNoMoreAdress (structure::Environment&);
void printAdjacencyMatrix (structure::Environment&, std::string status_field="status"); 
void printMatrix(const structure::Environment&, std::string);
void readData(structure::Environment&);
// FRS 6 may 2020 parseCommandLine is not used anymore
// bool parseCommandLine(structure::Environment&, int, char**);
void printEnvironment(const structure::Environment&);
void setEnvironment(structure::Environment&);
int** copyMatrix(int**, int, int);
void setNumberLevels(structure::Environment&);
bool existsTest(const std::string&);
bool checkNA(int**, int, int);
void saveAdjMatrix(const structure::Environment&, const std::string);
void saveAdjMatrixState(const structure::Environment&, const std::string);
std::vector<std::vector<std::string> > saveEdgesListAsTable(
    structure::Environment&);
void saveExecTime(const structure::Environment&, const std::string);
std::string arrayToString(const double*, const int);
std::string vectorToString(const std::vector<int>&);
std::string vectorToStringNodeName(
    const structure::Environment&, const std::vector<int>&);
void readTime(structure::Environment&, std::string);
// FRS 6 may 2020 readFilesAndFillStructures is not used anymore
// void readFilesAndFillStructures(
//    std::vector<std::string> edgesVectorOneLine, structure::Environment&);
bool readBlackbox(std::vector<std::string>, structure::Environment&);
std::vector<std::vector<std::string> > getAdjMatrix(
    const structure::Environment&);
int sign(double val);
void transformToFactors(structure::Environment&, int);
void transformToFactorsContinuous(structure::Environment&, int);
void transformToFactorsContinuousIdx(structure::Environment&, int);
void computeMeansandStandardDeviations(structure::Environment&);
void computeCorrelations(structure::Environment&);
void sort2arraysConfidence(int len, int a[], int brr[]);
bool allVariablesDiscrete(int* array, int* posArray, int num);
void sort2arrays(int len, int a[], int brr[], int bridge[]);
double get_wall_time();
double ramanujan(int n);
int printProgress(double percentage, double startTime, std::string outdir,
    int prg_numSearchMore);
// KL divergence functions
double compute_kl_divergence(int* posArray, structure::Environment&,
    int samplesNotNA, int** dataNumeric_red, int* AllLevels_red,
    int* samplesToEvaluate);

double kl(double** freqs1, double** freqs2, int nrows, int ncols);
double kl(int** counts1, double** freqs2, int nrows, int ncols);
double kl(double* freqs1, double* freqs2, int nlevels);
double kl(int* freqs1, int* freqs2, int n1, int n2, int nlevels);

void getJointMixed(structure::Environment&, int i, int j, int* mixedDiscrete,
    double* mixedContinuous, int* curr_samplesToEvaluate);
double** getJointFreqs(
    structure::Environment&, int i, int j, int numSamples_nonNA);
void getJointSpace(structure::Environment&, int i, int j, double** jointSpace,
    int* curr_samplesToEvaluate);
int getNumSamples_nonNA(structure::Environment&, int i, int j);

bool checkInterrupt(bool check = true);

}  // namespace utility
}  // namespace miic

#endif  // MIIC_UTILITIES_H_
