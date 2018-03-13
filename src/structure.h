#include "memory.h"
#ifndef STRUCTURE_H
#define STRUCTURE_H
#include <string>
#include <vector>

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
	int Nxyz_ui; // Ixy_uiz of the best z
	double Rxyz_ui; // score of the best contributor
	// string key;
	double Ixy_ui;
	double cplx;
	int Nxy_ui;
	short int status;

	std::vector<int> indexStruct; // memorize the corresponding structure in the structures list in the environment by its index
	std::vector<int> edgesInSpeTpl_list; // memorize the index of open structures
};

struct Node{
	std::string name;
	int level;
};

struct ExecutionTime{
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
	EdgeStructure* edgeStructure;
};

//
 /* Structure for all the needed parameters (input plus state variables)
 */
struct Environment {
	ExecutionTime execTime;
	int seed;
	
	int nThreads;
	MemorySpace m;
	std::vector<int> steps;
	ContainerMemory* memoryThreads;

	double** shuffleListNumEdges; // matrix to keep the number of edges for each eta and shuffle
	double* c2terms;
	int* oneLineMatrix;

	std::string myVersion; // -i parameter

	std::string outDir; // -i parameter
	std::string inData; // -o parameter
	std::string cplxType; // -c parameter
	std::string blackbox_name;
	std::string edgeFile;
		
	int numNodes;
	int numSamples;
	
	int typeOfData;

	std::vector<XJAddress*> searchMoreAddress;
	std::vector<XJAddress*> noMoreAddress;
	int numSearchMore;
	int numNoMore;

	// keep a trace of the number of edges for every state
	int phantomEdgenNum;


	
	Node* nodes;
	Edge** edges;
	std::vector< std::vector <std::string> > data;
	std::vector<std::string>  vectorData;
	int** dataNumeric;

	double** dataDouble;


	int* allLevels;
	
	double logEta;
	double l;

	int numberShuffles;
	double confidenceThreshold;

	bool myTest;
	bool isDegeneracy; // -d parameter
	bool isVerbose; // -v parameter
	bool isLatent; // -l parameter
	bool isNoInitEta; // -f parameter
	bool propag;
	bool isK23; // -k parameter
	bool isPropagation;

	int isTplReuse; // -r parameter
	int cplx;
	int eta; // -e parameter
	int shuffle; // -s parameter
	int halfVStructures;
	int etaIterator; // -e parameter
	int shuffleIterator; // -s parameter

	int effN; // -n parameter
	int minN;
	//int nDg;
	int thresPc;

	int iCountStruct;
	std::vector<Struct*> globalListOfStruct;
	std::vector<int> vstructWithContradiction;
};


struct Container {
	Environment* environment;
	int start;
	int stop; 
	
	MemorySpace m;
};

struct ContainerInit {
	Environment* environment;
	int i;
	int j; 
	int steps;

	MemorySpace m;
};

#endif
