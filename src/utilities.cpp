#include <fstream>
#include <sstream>
#include <string>
#include <getopt.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <unordered_set>

#include "utilities.h"


using namespace std;

class sort_indices
{
   private:
     int* mparr;
   public:
     sort_indices(int* parr) : mparr(parr) {}
     bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
};





// sort the indexed array brr using the array a
void sort2arrays(int len, int a[], int brr[], int bridge[]){
    int i;

    int *pArray = &a[1];
    int *pArray2 = &brr[1];

    std::sort(pArray2, pArray2+len, sort_indices(pArray));

    for(i=1; i < len+1; i++){
    	bridge[i] = pArray2[i-1];
    }

    brr = bridge;
    brr[0] = 0;
}


void createMemorySpace(Environment& environment, MemorySpace& m){

	int maxLevel = 0;
	for(int i =0; i<environment.numNodes; i++)
		if(environment.allLevels[i] > maxLevel)
			maxLevel = environment.allLevels[i];
	m.maxlevel = maxLevel;
	int nrow=environment.numSamples+1;
	int sampleSize = environment.numSamples;
	int ncol=7;
	int bin_max=maxLevel;
	int iii;

	m.sample = (int **)calloc(nrow, sizeof(int*));
	for(iii = 0; iii < nrow; iii++)
		m.sample[iii] = (int *)calloc(ncol, sizeof(int));

	m.sortedSample = (int **)calloc(nrow, sizeof(int*));
	for(iii = 0; iii < nrow; iii++)
	 	m.sortedSample[iii] = (int *)calloc(ncol, sizeof(int));

	m.Opt_sortedSample = (int **)calloc(nrow, sizeof(int*));
	for(iii = 0; iii < nrow; iii++)
	 	m.Opt_sortedSample[iii] = (int *)calloc(ncol, sizeof(int));

	m.Nxuiz = (int **)calloc(bin_max+1, sizeof(int*));
	for(iii = 0; iii < bin_max+1; iii++)
		m.Nxuiz[iii] = (int *)calloc(bin_max+1, sizeof(int));


	m.orderSample = (int *)calloc((sampleSize+2), sizeof(int));
	m.sampleKey = (int *)calloc((sampleSize+2), sizeof(int));

	m.Nxyuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nyuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nz = (int *)calloc((bin_max+1), sizeof(int));

	m.Ny = (int *)calloc((bin_max+1), sizeof(int));
	m.Nxui = (int *)calloc((bin_max+1), sizeof(int));
	m.Nx = (int *)calloc((bin_max+1), sizeof(int));

	m.bridge = (int *)calloc(sampleSize+2, sizeof(int));
}



void createMemorySpaceThreads(Environment& environment, ContainerMemory& m){

	int maxLevel = 0;
	for(int i =0; i<environment.numNodes; i++)
		if(environment.allLevels[i] > maxLevel)
			maxLevel = environment.allLevels[i];

	int nrow=environment.numSamples+1;
	int sampleSize = environment.numSamples;
	int ncol=7;
	int bin_max=maxLevel;
	int iii;


	m.sortedSample = (int **)calloc(nrow, sizeof(int*));
	for(iii = 0; iii < nrow; iii++) 
	 	m.sortedSample[iii] = (int *)calloc(ncol, sizeof(int));


	m.Nxuiz = (int **)calloc(bin_max+1, sizeof(int*));
	for(iii = 0; iii < bin_max+1; iii++) 
		m.Nxuiz[iii] = (int *)calloc(bin_max+1, sizeof(int));


	m.Nxyuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nyuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nuiz = (int *)calloc((bin_max+1), sizeof(int));
	m.Nz = (int *)calloc((bin_max+1), sizeof(int));

	m.Ny = (int *)calloc((bin_max+1), sizeof(int));
	m.Nxui = (int *)calloc((bin_max+1), sizeof(int));
	m.Nx = (int *)calloc((bin_max+1), sizeof(int));	
	m.bridge = (int *)calloc(sampleSize+2, sizeof(int));
}

void deleteMemorySpaceThreads(Environment& environment, ContainerMemory& m){

	int maxLevel = 0;
	for(int i =0; i<environment.numNodes; i++)
		if(environment.allLevels[i] > maxLevel)
			maxLevel = environment.allLevels[i];

	int nrow=environment.numSamples+1;
	int bin_max=maxLevel;
	int iii;


	for(iii = 0; iii < nrow; iii++) 
	 	free(m.sortedSample[iii]);
	free(m.sortedSample);

	
	for(iii = 0; iii < bin_max+1; iii++) 
		free(m.Nxuiz[iii]);
	free(m.Nxuiz);

	free(m.Nxyuiz);
	free(m.Nyuiz);
	free(m.Nuiz);
	free(m.Nz);

	free(m.Ny);
	free(m.Nxui);
	free(m.Nx);
	free(m.bridge);
}



void deleteMemorySpace(Environment& environment, MemorySpace& m){

	int maxLevel = 0;
	for(int i =0; i<environment.numNodes; i++)
		if(environment.allLevels[i] > maxLevel)
			maxLevel = environment.allLevels[i];

	int nrow=environment.numSamples+1;
	int bin_max=maxLevel;
	int i;

	for(i=0; i<nrow;i++)
		free(m.sample[i]);
	free(m.sample);

	for(i=0; i<nrow;i++)
		free(m.sortedSample[i]);
	free(m.sortedSample);

	for(i=0; i<bin_max+1;i++)
		free(m.Nxuiz[i]);
	free(m.Nxuiz);

	for(i=0; i<bin_max+1;i++)
		free(m.Opt_sortedSample[i]);
	free(m.Opt_sortedSample);

	free(m.orderSample);
	free(m.sampleKey);

	free(m.Nxyuiz);
	free(m.Nyuiz);
	free(m.Nuiz);
	free(m.Nz);

	free(m.Ny);
	free(m.Nxui);
	free(m.Nx);
	free(m.bridge);
}




void deleteStruct(Environment& environment){

	for(int i = 0 ; i < environment.numNoMore; i++)
		delete environment.noMoreAddress[i];


	delete [] environment.oneLineMatrix;
	delete [] environment.allLevels;
	for(int i = 0; i < environment.numSamples; i++){
		delete [] environment.dataNumeric[i];
	}
	delete [] environment.dataNumeric;
	delete [] environment.c2terms;
	delete [] environment.nodes;


	for(int i = 0; i < environment.numNodes - 1; i++){
		for(int j = i + 1; j < environment.numNodes; j++){
			delete environment.edges[i][j].edgeStructure;
		}
	}


	for(int i = 0; i < environment.numNodes; i++){
		delete [] environment.edges[i];
	}
	delete [] environment.edges;
}




bool SortFunctionNoMore1(const XJAddress* a, const XJAddress* b, Environment& environment) {
	 return environment.edges[a->i][a->j].edgeStructure->Ixy_ui > environment.edges[b->i][b->j].edgeStructure->Ixy_ui;
}

class sorterNoMore {
	  Environment& environment;
		public:
	  sorterNoMore(Environment& env) : environment(env) {}
	  bool operator()(XJAddress const* o1, XJAddress const* o2) const {
			return SortFunctionNoMore1(o1, o2, environment );
	  }
};

bool comparatorPairs ( const pair<double,int>& l, const pair<double,int>& r)
{ return l.first < r.first; }

bool isOnlyDouble(const char* str) {
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;
    return true;
}

bool SortFunction(const XJAddress* a, const XJAddress* b, Environment& environment) {

	if(environment.edges[a->i][a->j].edgeStructure->status < environment.edges[b->i][b->j].edgeStructure->status)
		return true;
	else if(environment.edges[a->i][a->j].edgeStructure->status > environment.edges[b->i][b->j].edgeStructure->status)
		return false;

	if(environment.edges[a->i][a->j].edgeStructure->status == 1 && environment.edges[b->i][b->j].edgeStructure->status == 1){
		if(environment.edges[a->i][a->j].edgeStructure->Rxyz_ui == 0 && environment.edges[b->i][b->j].edgeStructure->Rxyz_ui != 0)
			return true;
		else if(environment.edges[a->i][a->j].edgeStructure->Rxyz_ui != 0 && environment.edges[b->i][b->j].edgeStructure->Rxyz_ui == 0)
			return false;

		if(environment.edges[a->i][a->j].edgeStructure->Rxyz_ui > environment.edges[b->i][b->j].edgeStructure->Rxyz_ui)
			return true;
		else if(environment.edges[a->i][a->j].edgeStructure->Rxyz_ui < environment.edges[b->i][b->j].edgeStructure->Rxyz_ui)
			return false;
	}

	if(environment.edges[a->i][a->j].edgeStructure->status == 3 && environment.edges[b->i][b->j].edgeStructure->status == 3){
		if(environment.edges[a->i][a->j].edgeStructure->Ixy_ui > environment.edges[b->i][b->j].edgeStructure->Ixy_ui)
			return true;
		else if(environment.edges[a->i][a->j].edgeStructure->Ixy_ui < environment.edges[b->i][b->j].edgeStructure->Ixy_ui)
			return false;
	}
	return false;
}

class sorter {
	  Environment& environment;
		public:
	  sorter(Environment& env) : environment(env) {}
	  bool operator()(XJAddress const* o1, XJAddress const* o2) const {
			return SortFunction(o1, o2, environment );
	  }
};



// bool readTime(Environment& environment, string name){
// 	const char * c = name.c_str();
// 	ifstream input (c);
// 	string lineData;
// 	string s;
// 	int row = 0;
// 	int col = 0;
// 	while(getline(input, lineData))
// 	{
// 		if(row == 1){
// 		istringstream f(lineData);
// 			while (getline(f, s, '\t')) {
// 				if(col == 0)
// 					environment.execTime.init = atof(s.c_str());
// 				else if(col == 1)
// 					environment.execTime.iter = atof(s.c_str());
// 				else if(col == 2)
// 					environment.execTime.initIter = atof(s.c_str());
// 				else if(col == 3)
// 					environment.execTime.ort = atof(s.c_str());
// 				else if(col == 4)
// 					environment.execTime.cut = atof(s.c_str());
// 				else if(col == 5)
// 					environment.execTime.ort_after_cut = atof(s.c_str());
// 				else if(col == 6)
// 					environment.execTime.total = atof(s.c_str());
// 				col++;
// 			}
// 		}
// 		row++;
// 	}
// }


vector< vector <string> > getAdjMatrix(const Environment& environment){
	stringstream ss;
	vector< vector <string> > adjMatrix;
	vector<string> vec;
	for (int i=0; i < environment.numNodes; i++){
		vec.push_back(environment.nodes[i].name);
	}

	adjMatrix.push_back(vec);

	for (int i = 0; i < environment.numNodes; i++)
	{
		vec.clear();
		vec.push_back(environment.nodes[i].name);
		for (int j = 0; j < environment.numNodes; j++){
			ss.str("");
			ss << environment.edges[i][j].isConnected;
			vec.push_back(ss.str());
		}
		adjMatrix.push_back(vec);
	}

	return adjMatrix;
}

/*
 * Transform a vector to a string
 */
string vectorToStringNodeName(const Environment& environment, const vector<int> vec){
	stringstream ss;
	int length = vec.size();
	if(length > 0){
	  	for (int temp = 0; temp < length; temp++){
	  		if(vec[temp] != -1)
				ss << environment.nodes[vec[temp]].name;
			if(temp+1 < length)
				ss << ",";
		}
	} else {
		ss << "NA";
	}
  	return ss.str();
}

string vectorToString(const vector<int> vec){
	stringstream ss;
	int length = vec.size();
	if(length > 0){
	  	for (int temp = 0; temp < length; temp++){
			ss << vec[temp];
			if(temp+1 < length)
				ss << ",";
		}
	}
  	return ss.str();
}

string arrayToString1(const double* int_array, const int length){
	stringstream ss;
	if(length > 0){
	  	for (int temp = 0; temp < length; temp++){
	  		if(int_array[temp] != -1)
				ss << int_array[temp] << ", ";
		}
	} else {
		ss << "NA";
	}
  	return ss.str();
}

string zNameToString(const Environment& environment, vector<int> vec, int pos){
	stringstream ss;
	if(pos != -1)
		ss << environment.nodes[vec[pos]].name;
	else
		ss << "NA";
	return ss.str();
}



vector< vector <string> > saveEdgesListAsTable1(Environment& environment){

	vector< vector <string> > data;

	vector<XJAddress*> allEdges;

	for(int i = 0; i < environment.numNodes -1; i++){
	 	for(int j = i + 1; j < environment.numNodes; j++){
	 		XJAddress* s = new XJAddress();
			s->i=i;
			s->j=j;
			allEdges.push_back(s);
	 	}
	}

	vector<string> row;

	std::sort(allEdges.begin(), allEdges.end(), sorter(environment));

	row.push_back("x");
	row.push_back("y");
	row.push_back("z.name");
	row.push_back("ai.vect");
	row.push_back("zi.vect");
	row.push_back("Ixy_ai");
	row.push_back("cplx");
	row.push_back("Rxyz_ai");
	row.push_back("category");
	row.push_back("Nxy_ai");

	data.push_back(row);


	for(int i = 0; i < allEdges.size();i++){

		stringstream output;
		row.clear();
		for(int j = 0; j < environment.numNodes;j++){
			if(j == allEdges[i]->i || j == allEdges[i]->j)
				output << "1";
			else
				output << "0";
		}
		row.push_back(output.str());
		row.push_back(environment.nodes[allEdges[i]->i].name);
		row.push_back(environment.nodes[allEdges[i]->j].name);
 		row.push_back(zNameToString(environment, environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->zi_vect_idx, environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->z_name_idx));
 		row.push_back(vectorToStringNodeName(environment, environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->ui_vect_idx));
 		row.push_back(vectorToStringNodeName(environment, environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->zi_vect_idx));

 		output.str("");
 		output << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->Ixy_ui;
 		row.push_back(output.str());

 		output.str("");
 		output << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->cplx;
 		row.push_back(output.str());

 		output.str("");
 		output << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->Rxyz_ui;
 		row.push_back(output.str());

 		output.str("");
 		output << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->status;
 		row.push_back(output.str());

 		output.str("");
 		output << environment.edges[allEdges[i]->i][allEdges[i]->j].edgeStructure->Nxy_ui;
 		row.push_back(output.str());

 		data.push_back(row);
	}

	for(int i = 0; i < allEdges.size();i++){
		delete allEdges[i];
	}


	return data;
}


void copyValue(Environment& environment, int i){

 	for(int j = 0; j < environment.numSamples; j++){

			environment.dataNumeric[j][i] = atof(environment.data[j][i].c_str());
	}
}

bool checkNA(int** data, int numRows, int numColumns){
	for(int i = 0; i < numRows; i++)
		for(int j = 0; j < numColumns; j++)
			if(data[i][j] == -1)
				return true;

	return false;
}

/*
 * Initialize all the elements of the array to the given value
 */
bool setArrayValuesInt(int* array, int length, int value){
	for(int i = 0; i < length; i++){
			array[i] = value;
		}
		return true;
}

/*
 * Check if a value is of type integer
 */
bool isInteger(const string &s)
{
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

   char * p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}

/*
 * Find the average of a vector
*/
double findAvg(const Environment& environment, double** shuffleListNumEdges, int i){
	stringstream ss;
	double mean = 0;
	for(int j = 1; j <= environment.shuffle; j++){
		mean += shuffleListNumEdges[i][j-1];
	}

	mean /= environment.shuffle;
	ss >> mean;
	return mean;
}


/*
 * Read data matrix from file
 */
 bool readData(Environment& environment, bool& isNA){

 	vector <string> vec;
	environment.nodes = new Node[environment.numNodes];

	//convert input data
	for(int i = 0; i < environment.vectorData.size(); i++){
		if(i < environment.numNodes){
			environment.nodes[i].name = environment.vectorData[i];
		}
		else{
			if(i % environment.numNodes == 0){
				if(i != environment.numNodes){
					environment.data.push_back(vec);
					vec.clear();
				}
			}
			vec.push_back(environment.vectorData[i]);
		}
	}

	environment.data.push_back(vec);

	environment.numSamples = environment.data.size();

	// set effN if not set before to the number of rows from the input environment.data
	if(environment.effN == -1)
		environment.effN = environment.numSamples;

	return true;
}

/*
 * Remove all the lines that contain only NA
 */
bool removeRowsAllNA(Environment& environment){
	int* indexNA = new int[environment.numSamples];
	setArrayValuesInt(indexNA, environment.numSamples, -1);
	int pos = 0;
	for(int i = 0; i < environment.numSamples; i++){
		bool isNA = true;
		for(int j = 0; j < environment.numNodes && isNA; j++){
			if((environment.data[i][j].compare("NA")  != 0) && (environment.data[i][j].compare("") != 0)){
				isNA = false;
			}
		}
		if(isNA){
			indexNA[pos] = i;
			pos++;
		}
	}

	// if there are rows of NA value
	if(pos != 0){
		//correct variable numSamples

		// save the values
		int pos = 0;
		for(int i = 0; i < environment.numSamples; i++){
			if(i != indexNA[pos]){
				for(int j = 0; j < environment.numNodes; j++){
					environment.data[i-pos][j] = environment.data[i][j];
				}
			} else {
				pos++;
			}
		}
		environment.numSamples -= pos;
	}

	delete [] indexNA;
	return true;
}
/*
 * Transforms the string into factors
 */
void transformToFactors(Environment& environment, int i){
	 // create a dictionary to store the factors of the strings
 	map<string,int> myMap;

	//clean the dictionary since it is used column by column
	myMap.clear();
	myMap["NA"] = -1;
	myMap[""] = -1;
	int factor = 0;

 	for(int j = 0; j < environment.numSamples; j++){

		map<string,int>::iterator it = myMap.find(environment.data[j][i]);
		if ( it != myMap.end() ){
			environment.dataNumeric[j][i] = it->second;
		}
		else {
			myMap[environment.data[j][i]] = factor;
			environment.dataNumeric[j][i] = factor;
			factor++;
		}
	}
}


/*
 * Set the number of levels for each node (the maximum level of each column)
 */
void setNumberLevels(Environment& environment){
	int max;
 	environment.allLevels = new int[environment.numNodes];

	for(int i = 0; i < environment.numNodes;i++){
		max = 0;
	 	for(int j = 0; j < environment.numSamples; j++){
	 		if(environment.dataNumeric[j][i] > max)
	 			max = environment.dataNumeric[j][i];
	 	}


	 	environment.allLevels[i] = max+1;
	 }
}

/*
 * Set the variables in the environment structure
 */
void setEnvironment(Environment& environment){

	// Load the data
	// ----	
	environment.noMoreAddress.clear();
	environment.numNoMore = 0;
	environment.searchMoreAddress.clear();
	environment.numSearchMore = 0;

	environment.globalListOfStruct.clear();
	environment.vstructWithContradiction.clear();

	bool isNA = false;

	readData(environment, isNA);


	if(isNA){
		//// Remove the lines that are all 'NA'
		removeRowsAllNA(environment);
	}


	//create the data matrix for factors
	 environment.dataNumeric = new int*[environment.numSamples];
	 for(int i = 0; i < environment.numSamples; i++){
	 	environment.dataNumeric[i] = new int[environment.numNodes];
	 }

	for(int i = 0; i < environment.numNodes; i++){
		transformToFactors(environment, i);
	}

	//printMatrix(environment, "factors");

	//// Set the effN if not already done
	if(environment.effN == -1 )
		environment.effN = environment.numSamples;

	//// Set a variables with all properties name and levels
	setNumberLevels(environment);

	// create the 1000 entries to store c2 values
	environment.c2terms = new double[environment.numSamples+1];
	for(int i = 0; i < environment.numSamples+1; i++){
		environment.c2terms[i] = -1;
	}

	//// Set the number of digits for the precision while using round( ..., digits = ... )
	//// Make sure the min levels for the data is 0
	environment.minN = 1;

	//// Set the probability threshold for the rank
	environment.thresPc = 0;	// if the contribution probability is the min value
	environment.l = (environment.numNodes*(environment.numNodes-1)/2);

	// Stats test correction
	environment.logEta = std::log(static_cast<double>(environment.eta));//log(environment.eta);

	// create the edge structure and keep track of how many searchMore we have
	environment.edges = new Edge*[environment.numNodes];

	for(int i = 0; i < environment.numNodes; i++)
		environment.edges[i] = new Edge[environment.numNodes];

	for(int i = 0; i < environment.numNodes; i++){
		for(int j = 0; j < environment.numNodes; j++){
			environment.edges[i][j].isConnected = 1;
		}
	}

}

void readFilesAndFillStructures(vector<string> edgesVectorOneLine, Environment& environment){
	//fill nodes
	setEnvironment(environment);

	//create the one line matrix
	environment.oneLineMatrix = new int[environment.numSamples*environment.numNodes];
	for(int i = 0; i < environment.numSamples;i++){
		for(int j = 0; j < environment.numNodes;j++){
			environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
		}
	}

	//create edges
	environment.edges = new Edge*[environment.numNodes];

	for(int i = 0; i < environment.numNodes; i++)
		environment.edges[i] = new Edge[environment.numNodes];

	for(int i = 0; i < environment.numNodes; i++)
		environment.edges[i][i].isConnected = 0;


	for(int i = 0; i < environment.numNodes - 1; i++){
		for(int j = i + 1; j < environment.numNodes; j++){
			// create a structure for the nodes that need to store information about them
			environment.edges[i][j].edgeStructure = new EdgeStructure();
			environment.edges[j][i].edgeStructure = environment.edges[i][j].edgeStructure ;
			// initialize the structure
			environment.edges[j][i].edgeStructure->z_name_idx = -1;
			environment.edges[j][i].edgeStructure->status = -1;
		}
	}

	string lineData;
	string s;
	int posX = -1;
	int posY = -1;
	int numCols = 10;

	vector< vector <string> > vec;
	vector <string> v;

	for(int i = 0; i < edgesVectorOneLine.size(); i++){
		v.push_back(edgesVectorOneLine[i]);
		if((i + 1)% numCols == 0 && i != 0){
			vec.push_back(v);
			v.clear();
		}
	}

	for(int row = 0; row < vec.size(); row++)
	{
		v = vec[row];
		for(int col = 0; col < v.size(); col++) {
			string s = vec[row][col];

			if(col == 0){
				for(int i = 0; i < environment.numNodes;i++)
					if(environment.nodes[i].name.compare(s) == 0)
						posX=i;

			}
			else if(col == 1){
				for(int i = 0; i < environment.numNodes;i++)
					if(environment.nodes[i].name.compare(s) == 0)
						posY=i;
			}
			else if(col == 2){
			}
			else if(col == 3){
				if(s.compare("NA") != 0){
					stringstream ss(s); // Turn the string into a stream.
					string tok;
					char delimiter = ',';
					while(getline(ss, tok, delimiter)) {
						int ival;
						for(int i = 0; i < environment.numNodes;i++){
							if(environment.nodes[i].name.compare(tok) == 0){
								ival=i;
								break;
							}
						}
						environment.edges[posX][posY].edgeStructure->ui_vect_idx.push_back(ival);
					}
				}
			} else if(col == 4){
				if(s.compare("NA") != 0){
					stringstream ss(s); // Turn the string into a stream.
					string tok;
					char delimiter = ',';
					while(getline(ss, tok, delimiter)) {
						int ival;
						for(int i = 0; i < environment.numNodes;i++){
							if(environment.nodes[i].name.compare(tok) == 0){
								ival=i;
								break;
							}
						}
						environment.edges[posX][posY].edgeStructure->zi_vect_idx.push_back(ival);
					}
				}
			} else if(col == 5){
				environment.edges[posX][posY].edgeStructure->Ixy_ui = atof(s.c_str());
			} else if(col == 6){
				environment.edges[posX][posY].edgeStructure->cplx = atof(s.c_str());
			} else if(col == 7){
				environment.edges[posX][posY].edgeStructure->Rxyz_ui = atof(s.c_str());
			} else if(col == 8){
				int state = atoi(s.c_str());
				environment.edges[posX][posY].edgeStructure->status = state;
				if(state == 3){
					environment.edges[posX][posY].isConnected = 1;
					environment.edges[posY][posX].isConnected = 1;
					// add the edge to Nomore
					XJAddress* ij = new XJAddress();
					ij->i = posX;
					ij->j = posY;
					environment.noMoreAddress.push_back(ij);
				}
				else{
					environment.edges[posX][posY].isConnected = 0;
					environment.edges[posY][posX].isConnected = 0;
				}
			} else if(col == 9){
				environment.edges[posX][posY].edgeStructure->Nxy_ui = atof(s.c_str());
			}
		}


	}

	environment.numNoMore = environment.noMoreAddress.size();

	std::sort(environment.noMoreAddress.begin(), environment.noMoreAddress.end(), sorterNoMore(environment));
}

bool readBlackbox1(vector<string> v, Environment& environment){
	string lineData;
	string s;
	string s1;
	string s2;
	int posX;
	int posY;

	for(int pos = 0; pos < v.size(); pos++){
		posX = -1;
		posY = -1;

		s1 = v[pos];
		for(int i = 0; i < environment.numNodes;i++){
			if(environment.nodes[i].name.compare(s1) == 0)
				posX=i;
		}

		pos++;
		s2 = v[pos];
		for(int i = 0; i < environment.numNodes;i++){
			if(environment.nodes[i].name.compare(s2) == 0)
					posY=i;
		}

		if(posX != -1 && posY != -1){
			environment.edges[posX][posY].isConnected = 0;
			environment.edges[posY][posX].isConnected = 0;
		}
	}

	return true;
}

bool isContinuousDiscrete(int** dataNumeric_red, int nsamplesNotNA, int pos){
	
	int j = pos;
 	std::unordered_set<int> s(dataNumeric_red[j], dataNumeric_red[j] + nsamplesNotNA);
 	int size = s.size();
	// cout << size << endl << flush;

	if((size < 0.3 * nsamplesNotNA && size < 40) || size < 2){
		// cout << "true" << flush;
		return true;
	}
	// cout << "false" << flush;
	return false;
}	

void sort2arraysConfidence(int len, int a[], int brr[]){

    std::sort(brr, brr+len, sort_indices(a));
}

double ramanujan(int n){
	// Returns log(fac(n)) with Ramanujan's approximation.
    if(n==0){
      return(1);
    }
    double N = n*log(n) - n + log(1.0*n*(1 + 4*n*(1+2*n)))/6 + log(M_PI)/2L;
    return(N);
}

double logchoose(int n, int k){
	// Returns the log of n choose k with Ramanujan's approximation.
	double res = ramanujan(n) - ramanujan(k) - ramanujan(n-k);
	return(res);
}

double compute_parametric_complexity(int n, int K, double** sc_look){

  if(sc_look[n-1][K-1] != 0){
    return(sc_look[n-1][K-1]);
  }

  double res;
  if(K==1){
    res = 1;
  } else if(K==2){
    if(n<1000){
      res = 0;
      for(int i=0 ; i <= n; i++){
        int h1 = i;
        int h2 = n - h1;
        res = res + exp(ramanujan(n) - ramanujan(h1) - ramanujan(h2)) * pow((1.0*h1/n),h1) * pow((1.0*h2/n),h2);
      }
    } else{
      res = sqrt((n*M_PI)/2) * exp(sqrt(8./(9.*n*M_PI)) + (3*M_PI-16)/(36.*n*M_PI));
    }
  } else {
    res = compute_parametric_complexity(n, K-1, sc_look) + 1.0*n/(K-2) * compute_parametric_complexity(n, K-2, sc_look);
  }

  sc_look[n-1][K-1] = res;
  return(res);
}
