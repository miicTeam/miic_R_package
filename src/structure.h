#include "memory.h"
#ifndef STRUCTURE_H
#define STRUCTURE_H
#include <string>
#include <vector>
#include <map>
#include <tuple>

using namespace std;
using uint=unsigned int;

struct StructWithOrtToPropagate{
	int name;
	int idxStruct;
	int vStruct;
	int idxEdge;
	int dir;
};

struct Struct{
	int xi;
	int xk;
	int xj;

	int originalPositionStruct;

	double rv;
	int sv;
	double Ibase;
	double Istruct;

	short int d1;
	short int d2;
	short int isOut;
	short int orgOut;
	short int orgD1;
	short int orgD2;
};

struct EdgeStructure {
	int z_name_idx; // index of the last best contributor
	std::vector<int> ui_vect_idx; // index of ui
	std::vector<int> zi_vect_idx; // index of ui
	double Rxyz_ui; // score of the best contributor

	double Ixy_ui;
	double cplx;
	int Nxy_ui;
	short int status;
	// Keeping track of jointCounts (discrete) and jointSpace(continuous) for KL divergence when
	//adding Us to the conditioning set.
	double mutInfo;  // mutual information without conditioning
	double cplx_noU;  // complexity without conditioning
	int Nxy;  // count of joint factors without NA

	std::vector<int> indexStruct; // memorize the corresponding structure in the structures list in the environment by its index
	std::vector<int> edgesInSpeTpl_list; // memorize the index of open structures

public:
	EdgeStructure():
		z_name_idx(-1),
		ui_vect_idx(std::vector<int>()),
		zi_vect_idx(std::vector<int>()),
		Rxyz_ui(0.0),

		Ixy_ui(0.0),
		cplx(0.0),
		status(-1),
		Nxy_ui(-1),
		mutInfo(0.0),
		cplx_noU(0.0),
		Nxy(-1),

		indexStruct(std::vector<int>()),
		edgesInSpeTpl_list(std::vector<int>())
	{}
};

struct Node{
	std::string name;
	int level;
};

struct ExecutionTime{
	double startTimeInit;
	double startTimeIter;
	long double init;
	long double iter;
	long double initIter;
	long double ort;
	long double cut;
	long double ort_after_cut;
	long double total;
};

struct XJAddress{
	int i;
	int j;
};

struct Edge{
	short int isConnected;
	short int isConnectedAfterInitialization;
	short int areNeighboursAfterIteration;
	EdgeStructure* edgeStructure;
};

//
 /* Structure for all the needed parameters (input plus state variables)
 */
struct Environment {
	ExecutionTime execTime;
	bool consistentPhase;
	bool iterationStepEven;
	// for gaussian case
	double** rho;
	double* standardDeviations;
	double* means;
	double** dataDouble;
	double ** pMatrix;
	int ** nSamples;
	int ** iterative_cuts;
	double* sampleWeights;
	std::vector<double> sampleWeightsVec;


	bool testDistribution;
	int seed;
	uint nThreads;
	MemorySpace m;
	std::vector<int> steps;
	MemorySpace* memoryThreads;
	double** shuffleListNumEdges; // matrix to keep the number of edges for each eta and shuffle
	double* c2terms;
	double** cterms;
	double** lookchoose;
	int* columnAsContinuous;
	int* columnAsGaussian;
	std::vector<int> cntVarVec;
	int* oneLineMatrix;

	std::string myVersion; // -i parameter

	std::string outDir; // -i parameter
	std::string inData; // -o parameter
	std::string cplxType; // -c parameter
	std::string blackbox_name;
	std::string edgeFile;
	std::string dataTypeFile;
	//std::string sampleWeightsFile;// -w parameter

	uint numNodes;
	uint numSamples;
	bool firstIterationDone;

	std::vector<XJAddress*> searchMoreAddress;
	std::vector<XJAddress*> noMoreAddress;
	int numSearchMore;
	int numNoMore;

	int typeOfData;
	int isAllGaussian;
	int atLeastTwoGaussian;
	int atLeastOneContinuous;
	int atLeastTwoDiscrete;

	// keep a trace of the number of edges for every state
	int phantomEdgenNum;

	Node* nodes;
	Edge** edges;
	std::vector<std::string>  vectorData;
	std::vector< std::vector <std::string> > data;
	//std::string** data;
	int** dataNumeric;
	int** dataNumericIdx;
	int* allLevels;

	double** proportions;

	int** cut;

	double logEta;
	double l;

	bool myTest;
	bool isDegeneracy; // -d parameter
	bool isVerbose; // -v parameter
	bool isLatent; // -l parameter
	bool isLatentOnlyOrientation; // -l parameter
	bool isNoInitEta; // -f parameter
	bool propag;
	bool isK23; // -k parameter
	bool isPropagation;

	int isTplReuse; // -r parameter
	int cplx;
	int halfVStructures;

	int numberShuffles; // -s parameter
	double confidenceThreshold; // -e parameter

	uint effN; // -n parameter
	int minN;
	//int nDg;
	int thresPc;

	int maxbins;
	int initbins;
	double* looklog;
	double* lookH;

	double* noiseVec;

	int iCountStruct;
	std::vector<Struct*> globalListOfStruct;
	std::vector<int> vstructWithContradiction;
};

#endif
