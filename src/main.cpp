//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <string>
//#include <sys/stat.h>
//#include <algorithm>
//#include <ctime>
//#include <unistd.h>
//#include <string.h>
//#include "utilities.h"
//#include "log.h"
//#include "skeleton.h"
//#include "orientation.h"
//#include "orientationProbability.h"
//#include "confidenceCut.h"
//#include <time.h>
//#include <stdlib.h>
//
//
//
//using namespace std;
//#ifdef __unix__		 
//#elif defined(_WIN32) || defined(WIN32) 
//	#define OS_Windows
//#endif
//
//inline bool exists_test (const std::string& name) {
//  struct stat buffer;   
//  return (stat (name.c_str(), &buffer) == 0); 
//}
//
//int main(int argc, char* argv[]){
//
//	string slash;
//	#ifdef OS_Windows
//	 /* Windows code */
//		slash = "\\";
//	#else
//	 /* GNU/Linux code */
//		slash = "/";
//	#endif 
//
//
//	double startTimeWhole = get_wall_time();
//	double startTime;
//
//	// define the environment
//	Environment environment;
//	environment.numNodes=0;
//	environment.numSamples=0;
//	environment.myVersion="V79";
//
//	//set test to FALSE
//	environment.myTest = false;
//
//	if(environment.myTest){
//		cout << "# -- !! -- TEST TEST TEST -- !! -- #\n" ;
//		environment.outDir = "/home/nadir/Desktop/Trial/resSkeletonDelete";
//		environment.inData= "/home/nadir/Desktop/Trial/asia_catData_10000.txt";
//		environment.effN = -1;
//		environment.cplx = -1;
//		environment.isVerbose = true;
//		environment.numberShuffles = 0;
//		environment.isLatent = false;
//		environment.isTplReuse = true;
//		environment.isK23 = true;
//		environment.isDegeneracy = false;
//		environment.isNoInitEta = false;
//	}
//
//	//// Parse command line arguments
//	//// ----
//
//	if(!environment.myTest){
//		parseCommandLine(environment, argc, argv);
//		if(environment.isVerbose) {
//			cout << "# --------\n# Parsing Commmand Line arguments...\n";
//		cout << "# --------> Arguments from any variable inside the environment\n";
//		}
//	} else {
//		if(environment.isVerbose) {
//			cout << "# --------\n# Parsing Commmand Line arguments...\n";
//			cout << "# --------> No Arguments from a variable inside the environment\n";
//		}
//	}
//
//	srand(environment.seed);
//
//	if(find(environment.steps.begin(), environment.steps.end(), 1) != environment.steps.end()){
//		if(environment.isVerbose) cout << "# --------\n# -> START miic Skeleton...\n";
//		//perform the skeleton evaluation
//		skeleton(environment, slash, startTimeWhole, argc, argv);	
//
//
//	} else if(find(environment.steps.begin(), environment.steps.end(), 5) != environment.steps.end() || find(environment.steps.begin(), environment.steps.end(), 2) != environment.steps.end()){
//		readFilesAndFillStructures(environment, slash);
//		createMemorySpace(environment, environment.m);
//	}
//
//	if(environment.isLatentOnlyOrientation == true)
//		environment.isLatent = true;
//
//
//	stringstream ss;
//	ss.str("");
//	//// Save the execTime
//	ofstream output;
//	ss << environment.outDir << slash << "progress.txt";
//	string filename_progress = ss.str();
//
//	// write progress
//	if(!(environment.numberShuffles > 0)){
//		if( (find(environment.steps.begin(), environment.steps.end(), 5) != environment.steps.end()) || 
//			(find(environment.steps.begin(), environment.steps.end(), 2) != environment.steps.end()) ){
//			output.open(filename_progress.c_str());
//			output << "Orientation";
//			output.close();
//		}
//	}
//
//	ss.str("");
//	ss << environment.outDir << slash << "adjacencyMatrix.miic.orientProba.txt";
//	string filename = ss.str();
//
//	environment.testDistribution = false;
//	startTime = get_wall_time();
//	
//	if(!exists_test(filename) && environment.noMoreAddress.size() > 0) {
//		if(find(environment.steps.begin(), environment.steps.end(), 5) != environment.steps.end()){
//			//perform the edge orientation
//			orientation(environment, slash, environment.isVerbose);// environment.isVerbose);
//		} else if(find(environment.steps.begin(), environment.steps.end(), 2) != environment.steps.end()){
//			//perform the edge orientation
//			orientationProbability(environment, slash, environment.isVerbose);// environment.isVerbose);
//		}
//	}
//
//	long double spentTime = (get_wall_time() - startTime);
//	environment.execTime.ort = spentTime;
//
//	if(environment.numberShuffles > 0){
//		startTime = get_wall_time();
//		output.open(filename_progress.c_str());
//		output << "Confidence tool";
//		output.close();
//
//		//change output directory
//		ss.str("");
//		ss << environment.outDir << slash << "shuffle_" << environment.numberShuffles << slash << "filtered_network_" << environment.confidenceThreshold;
//		string od = environment.outDir;
//		environment.outDir = ss.str();
//
//		cout << "# --------\n# -> START Confidence Cut...\n";
//
//		confidenceCut(environment, slash, od, startTime);
//
//		cout << "# -> END Confidence Cut\n";
//
//		ss.str("");
//		ss << environment.outDir << slash << "adjacencyMatrix.miic.txt";
//		saveAdjMatrix(environment, ss.str());
//
//
//		long double spentTime = (get_wall_time() - startTime);
//		environment.execTime.cut = spentTime;
//
//
//		ss.str("");
//		//// Make a table from the list
//		ss << environment.outDir << slash << "edgesList.miic.txt";
//		saveEdgesListAsTable(environment, ss.str());
//
//
//
//		//readFilesAndFillStructures(environment, slash);
//
//		output.open(filename_progress.c_str());
//		output << "Orientation";
//		output.close();
//		startTime = get_wall_time();
//		if(environment.noMoreAddress.size() > 0){
//			if(find(environment.steps.begin(), environment.steps.end(), 5) != environment.steps.end()){
//				//perform the edge orientation
//				orientation(environment, slash, environment.isVerbose);// environment.isVerbose);
//			} else if(find(environment.steps.begin(), environment.steps.end(), 2) != environment.steps.end()){
//				//perform the edge orientation
//				orientationProbability(environment, slash, environment.isVerbose);// environment.isVerbose);
//			}
//		}
//
//		spentTime = (get_wall_time() - startTime);
//		environment.execTime.ort_after_cut = spentTime;
//	}
//
//
//
//	// ss.str("");
//	// //// Save the execTime
//	// ss << environment.outDir << slash << "edges.txt";
//	// filename = ss.str();
//	// output.open(filename.c_str());
//	// for(int i = 0; i < environment.numNodes -1; i++){
//	// 	for (int j = i+1; j < environment.numNodes; j++)
//	// 	{
//	// 		if(environment.edges[i][j].isConnected == -2 || environment.edges[i][j].isConnected == 2 || environment.edges[i][j].isConnected == 1)
//	// 			output << environment.nodes[i].name << "\t" << environment.nodes[j].name << "\n";
//	// 	}
//	// }
//	// output.close();
//
//	spentTime = (get_wall_time() - startTimeWhole);
//	environment.execTime.total = spentTime;
//	ss.str("");
//	ss << environment.outDir << slash << "execTime.miic.txt";
//	saveExecTime(environment, ss.str());
//
//
//	//delete pointers, free memory
//	deleteMemorySpace(environment, environment.m);
//	
//	delete [] environment.allLevels;
//	delete [] environment.oneLineMatrix;
//
//	for(int i = 0; i < environment.numSamples -1 ;i++ ){
//		for (int j = i+1; j < environment.numNodes; j++){
//			delete environment.edges[i][j].edgeStructure;
//		}
//	}
//
//	for(int i = 0; i < environment.numSamples;i++ ){	
//		delete [] environment.data[i];
//		delete [] environment.dataNumeric[i];
//	}
//
//	delete [] environment.columnAsContinuous;
//	delete [] environment.c2terms;
//	delete [] environment.edges;
//	delete [] environment.nodes;
//	delete [] environment.data;
//	delete [] environment.dataNumeric;
//
//	long double elapsedTime = (get_wall_time() - startTimeWhole);
//
//	//usleep(1000000000000000000);
//
//	if(environment.isVerbose == true){ cout << "\n# ----> Total execution time:" << elapsedTime << "millisec\n\n"; }
//}
//