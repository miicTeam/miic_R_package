#include "structure.h"
void createMemorySpace(Environment&, MemorySpace&);
void createMemorySpaceThreads(Environment&, ContainerMemory&);
void deleteMemorySpace(Environment&, MemorySpace&);
std::string printNodesName(Environment&);
void printMatrix(Environment&, std::string);
void readData(Environment&);
bool parseCommandLine(Environment&, int, char**);
void setEnvironment(Environment&);
int** copyMatrix(int**, int, int);
void setNumberLevels(Environment&);
bool checkNA(int**, int, int);
void saveAdjMatrix(const Environment&, const std::string);
void saveAdjMatrixState(const Environment&, const std::string);
void saveEdgesListAsTable(Environment&, const std::string);
std::vector< std::vector <std::string> > saveEdgesListAsTable1(Environment&);
void saveExecTime(const Environment&, ExecutionTime, const std::string);
double findAvg(const Environment&, double**, int);
std::string arrayToString(Environment&, const int* , const int );
std::string arrayToString1(const double*, const int);
std::string vectorToString(const std::vector<int> vec);
void readTime(std::string, ExecutionTime&);
void readFilesAndFillStructures(std::vector<std::string>, Environment&);
bool readBlackbox(Environment&);
bool removeBlackboxEdges(Environment&);
bool readBlackbox1(std::vector<std::string>, Environment&);
std::vector< std::vector <std::string> > getAdjMatrix(const Environment&);
bool readFileTypeData(std::vector<std::string>, Environment&);
void sort2arrays(int len, int a[], int brr[], int bridge[]);
void deleteStruct(Environment& environment);
void deleteMemorySpaceThreads(Environment& environment, ContainerMemory& m);
double logchoose(int n, int k, double* looklog);
double ramanujan(int n);
bool isContinuousDiscrete(int** dataNumeric_red, int nsamplesNotNA, int pos);
void sort2arraysConfidence(int len, int a[], int brr[]);
double compute_parametric_complexity(int n, int K, double** sc_look);