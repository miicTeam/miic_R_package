#include "structure.h"
#include "log.h"
void createMemorySpace(Environment&, MemorySpace&);
// bool createMemorySpaceThreads(Environment&, ContainerMemory&);
void deleteMemorySpace(Environment&, MemorySpace&);
std::string printNodesName(Environment&);
void printMatrix(Environment&, std::string);
void readData(Environment&);
bool parseCommandLine(Environment&, int, char**);
void printEnvironment(Environment&, Log*);
void setEnvironment(Environment&);
int** copyMatrix(int**, int, int);
void setNumberLevels(Environment&);
bool existsTest(const std::string&);
bool checkNA(int**, int, int);
void saveAdjMatrix(const Environment&, const std::string);
void saveAdjMatrixState(const Environment&, const std::string);
std::vector< std::vector <std::string> > saveEdgesListAsTable(Environment&);
void saveExecTime(const Environment&, const std::string);
double findAvg(const Environment&, double**, int);
std::string arrayToString(Environment&, const int* , const int );
std::string arrayToString1(const double*, const int);
std::string vectorToString(const std::vector<int> vec);
void readTime(Environment&, std::string);
void readFilesAndFillStructures(vector<string> edgesVectorOneLine, Environment& environment);
bool readBlackbox1(std::vector<std::string>, Environment&);
std::vector< std::vector <std::string> > getAdjMatrix(const Environment&);
int sign(double val);
bool readBlackbox(Environment&, std::string);
bool removeBlackboxEdges(Environment& environment);
void transformToFactors(Environment&, int);
void transformToFactorsContinuous(Environment&, int);
void transformToFactorsContinuousIdx(Environment&, int);
void computeMeansandStandardDeviations(Environment& environment);
void computeCorrelations(Environment& environment);
void sort2arraysConfidence(int len, int a[], int brr[]);
void updateNumberofLevelsAndContinuousVariables(int** dataNumeric_red, int* allLevels_red, int* cnt_red, int myNbrUi, int nsamplesNotNA);
bool isNumberOfLevelsLessTwo(int** dataNumeric_red, int myNbrUi, int nsamplesNotNA, int vars);
void updateContinuousVariables(int** dataNumeric_red, int* cnt_red, int*, int myNbrUi, int nsamplesNotNA, int vars);
bool isContinuousDiscrete(int** dataNumeric_red, int nsamplesNotNA, int pos);
bool allVariablesDiscrete(int* array, int* posArray, int num);
bool isContinuousDiscrete(int** dataNumeric_red, int nsamplesNotNA, int pos);
void sort2arrays(int len, int a[], int brr[], int bridge[]);
void deleteStruct(Environment& environment);
double get_wall_time();
double get_cpu_time();
double ramanujan(int n);
int printProgress (double percentage, double startTime, string outdir, int prg_numSearchMore);
//KL divergence functions
//double compute_kl_divergence_continuous(double** space1, double** space2, int n1, int n2, int ndims, int k,
//										bool* flag_break_ties, int* map_samples);
//double compute_kl_divergence_continuous(double* dist1, double* dist2, int n1, int n2, int k,
//										bool flag_break_ties, int* map_samples);
//double compute_k_nearest_distance(double* point, double** space, int ndims, int n, int k);
//double compute_k_nearest_distance(double point, double* dist, int n, int k);
double compute_kl_divergence(int* posArray, Environment& environment, int samplesNotNA, int** dataNumeric_red,
							 int* AllLevels_red, int* samplesToEvaluate);

// static double igf(double S, double Z);
// double approx_gamma(double Z);
// double chisqr(int Dof, double Cv);
// bool testDistribution(double* propObserved, double* propExpected,int length);
//void setJointCounts(Environment& environment);
double kl(double** freqs1, double** freqs2, int nrows, int ncols);
double kl(int** counts1, double** freqs2, int nrows, int ncols);
double kl(double* freqs1, double* freqs2, int nlevels);
double kl(int* freqs1, int* freqs2, int n1, int n2, int nlevels);
//void updateJointCounts(Environment& environment, int i , int j);

void getJointMixed(Environment& environment, int i, int j, int* mixedDiscrete, double* mixedContinuous,
				   int* curr_samplesToEvaluate);
double** getJointFreqs(Environment& environment, int i, int j, int numSamples_nonNA);
void getJointSpace(Environment& environment, int i, int j, double** jointSpace, int* curr_samplesToEvaluate);
int getNumSamples_nonNA(Environment& environment, int i, int j);

bool checkInterrupt(bool check=true);
