#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include <string.h>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "structure.h"
#include "utilities.h"
#include "log.h"
#include "skeletonInitialization.h"
#include "skeletonIteration.h"
#include "confidenceCut.h"



using namespace std;
using namespace Rcpp;

bool skeletonChanged(::Environment& environment){
	for(uint i = 0; i < environment.numNodes; i++){
		for(uint j = 0; j < environment.numNodes; j++){
			if(environment.edges[i][j].isConnected != environment.edges[i][j].areNeighboursAfterIteration){
				return true;
			}
		}
	}
	return false;
}

List empty_results(){
	List result;
	result = List::create(
		_["interrupted"] = true
	) ;
	return(result);
}

extern "C" SEXP skeleton(SEXP inputDataR, SEXP typeOfDataR, SEXP cntVarR, SEXP numNodesR, SEXP nThreadsR, SEXP blackBoxR, SEXP effNR, SEXP cplxR,
						 SEXP etaR, SEXP isLatentR, SEXP isTplReuseR, SEXP isK23R, SEXP isDegeneracyR, SEXP isNoInitEtaR, SEXP confidenceShuffleR,
						 SEXP confidenceThresholdR, SEXP sampleWeightsR, SEXP consistentR, SEXP testDistR, SEXP verboseR)
{

	vector< vector <double> > outScore;
	vector< vector <string> > edgesMatrix;
	stringstream output;

	// define the environment
	::Environment environment;
	environment.myVersion="V79";

	//set test to FALSE
	environment.myTest = false;

	environment.numNodes = Rcpp::as<int> (numNodesR);
	environment.nThreads = Rcpp::as<int> (nThreadsR);
	#ifdef _OPENMP
		int threads = environment.nThreads;
		if (threads < 0)
			threads = omp_get_num_procs();
		omp_set_num_threads(threads);
	#endif
	environment.consistentPhase = Rcpp::as<bool> (consistentR);
	environment.testDistribution = Rcpp::as<bool> (testDistR);

	environment.vectorData = Rcpp::as< vector <string> > (inputDataR);

	vector<string> v;
	v = Rcpp::as< vector <string> > (blackBoxR);


	string cplx;

	environment.effN = Rcpp::as<int> (effNR);
	environment.typeOfData = Rcpp::as<int> (typeOfDataR);
	cplx = Rcpp::as<string> (cplxR);
	//environment.eta = Rcpp::as<int> (etaR);
	//environment.shuffle = 0;
	environment.isLatent = Rcpp::as<bool> (isLatentR);
	environment.isTplReuse = Rcpp::as<bool> (isTplReuseR);
	environment.isK23 = Rcpp::as<bool> (isK23R);
	environment.isDegeneracy = Rcpp::as<bool> (isDegeneracyR);
	environment.isNoInitEta = Rcpp::as<bool> (isNoInitEtaR);

	environment.confidenceThreshold = Rcpp::as<double> (confidenceThresholdR);
	environment.numberShuffles = Rcpp::as<int> (confidenceShuffleR);

	environment.sampleWeightsVec = Rcpp::as< vector <double> > (sampleWeightsR);
	environment.cntVarVec = Rcpp::as< vector <int> > (cntVarR);

	environment.isVerbose = Rcpp::as<bool> (verboseR);

	if(cplx.compare("nml") == 0)
		environment.cplx = 1;
	else if(cplx.compare("mdl") == 0)
		environment.cplx = 0;

	double startTime;


	vector< vector <string> > retShuffle;
	vector< vector <string> > retShuffleAvg;
	vector <string> shfVec;

	string filePath;

	// set iterator to the bound, i.e only 1 iteration
	//environment.etaIterator = environment.eta;
	//environment.shuffleIterator = environment.shuffle;
	environment.isAllGaussian = false;

	// set the environment
	setEnvironment(environment);
		if(v.size() > 1)
		readBlackbox1(v, environment);

	startTime = get_wall_time();

	// ----
	environment.memoryThreads = new MemorySpace[environment.nThreads];
	for(uint i = 0; i < environment.nThreads; i++){
		createMemorySpace(environment, environment.memoryThreads[i]);
	}
	

	if(!skeletonInitialization(environment))
	{
		// structure the output
	    List result = List::create(
		 	_["error"] = "error during skeleton initialization",
		 	_["interrupted"] = true
	    ) ;
	    return result;
	}
	
	// ----
	long double spentTime = (get_wall_time() - startTime);
	environment.execTime.init = spentTime;
	if( environment.isVerbose == true ){ cout << "\n# ----> First contributing node elapsed time:" << spentTime << "sec\n\n"; }
	

	if(!firstStepIteration(environment)) return(empty_results());
	if( environment.numNoMore == 0 && environment.numSearchMore == 0 ) { 
		if( environment.isVerbose == true ){ cout << "# ------| Only phantom edges found.\n"; }
	} else if( environment.numSearchMore > 0 ) {
	
		//// Search for other Contributing node(s) (possible only for the edges still in 'searchMore', ie. 2)
		if( environment.isVerbose == true ){ cout << "\n# ---- Other Contributing node(s) ----\n\n"; }
		startTime = get_wall_time();
			
		if(!skeletonIteration(environment)) return(empty_results());
		long double spentTime = (get_wall_time() - startTime);
		environment.execTime.iter = spentTime;
		environment.execTime.initIter = environment.execTime.init + environment.execTime.iter;
	}
	// delete memory
	for(uint i = 0; i < environment.nThreads; i++){
		deleteMemorySpace(environment, environment.memoryThreads[i]);
	}
	delete [] environment.memoryThreads;
	cout << endl;

	int NEdgesBeforeIter = environment.numNoMore;
	int maxConsistentIter = 10;
	
	//run of the skeleton iteration phase for the consistent part
	if(environment.consistentPhase){

		cout << "\n======	Consistent phase iterations	=====" << endl;
		do{
			cout << "Consistent phase " << maxConsistentIter-10 << " (maximum 10)" << endl;
			int count = 0;
			//save the neighbours in the areNeighboursAfterIteration structure
			for(uint i = 0; i < environment.numNodes; i++){
				for(uint j = 0; j < environment.numNodes; j++){
					environment.edges[i][j].areNeighboursAfterIteration = environment.edges[i][j].isConnected;
					if(environment.edges[i][j].areNeighboursAfterIteration)
						count++;
				}
			}

            bool graph_is_consistent = true;
			for(uint i = 0; i < environment.numNodes - 1; i++){
				for(uint j = i + 1; j < environment.numNodes; j++){
                    if (!environment.edges[i][j].isConnectedAfterInitialization || environment.edges[i][j].edgeStructure->status != 1)
                        continue;
                    if (!is_consistent(environment, i, j, environment.edges[i][j].edgeStructure->ui_vect_idx))
                        graph_is_consistent = false;
                }
            }
            if (graph_is_consistent && maxConsistentIter == 10)
                break;

			//return back with structures at the moment of initialization
			environment.countSearchMore = 0;
			for(uint i = 0; i < environment.numNodes; i++){
				for(uint j = 0; j < environment.numNodes; j++){
					environment.edges[i][j].isConnected = environment.edges[i][j].isConnectedAfterInitialization; 
					if(environment.edges[i][j].isConnected && i>j)
						environment.countSearchMore++;
				}
			}

			for(uint i = 0; i < environment.numNodes - 1; i++){
				for(uint j = i + 1; j < environment.numNodes; j++){
					environment.edges[i][j].edgeStructure->zi_vect_idx.clear();
					environment.edges[i][j].edgeStructure->ui_vect_idx.clear();
					environment.edges[i][j].edgeStructure->z_name_idx = -1;
					environment.edges[i][j].edgeStructure->Rxyz_ui = 0;
				}
			}

			for(uint i = 0; i < environment.searchMoreAddress.size(); i++){
				delete  environment.searchMoreAddress[i];
			}

			for(uint i = 0; i < environment.noMoreAddress.size(); i++){
				delete  environment.noMoreAddress[i];
			}

			environment.noMoreAddress.clear();
			environment.searchMoreAddress.clear();

			firstStepIteration(environment);

			//// Search for other Contributing node(s) (possible only for the edges still in 'searchMore', ie. 2)
			if( environment.isVerbose == true ){ cout << "\n# ---- Other Contributing node(s) ----\n\n"; }
			startTime = get_wall_time();

			skeletonIteration(environment);

			long double spentTime = (get_wall_time() - startTime);
			environment.execTime.iter = spentTime;
			environment.execTime.initIter = environment.execTime.init + environment.execTime.iter;
			maxConsistentIter --;
		} while(skeletonChanged(environment) && maxConsistentIter>0);
		cout << "Returning consistent graph with " << environment.numNoMore << " edges (" << NEdgesBeforeIter << " before consistency checks)." << endl;
	}


	startTime = get_wall_time();
	vector< vector <string> >confVect;
	if(environment.numberShuffles > 0){
		cout << "Computing confidence cut with permutations..." << flush;
		confVect = confidenceCut(environment);
		long double spentTime = (get_wall_time() - startTime);
		environment.execTime.cut = spentTime;
		cout << " done." << endl;
	}
	else{
		environment.execTime.cut = 0;
	}


	if( environment.numNoMore > 0){
		cout << "Saving results for orientation..." << flush;
		edgesMatrix = saveEdgesListAsTable1(environment);
		cout << " done." << endl;
	}


	vector< vector <string> >  adjMatrix;
	adjMatrix = getAdjMatrix(environment);

	vector<double>time;
	time.push_back(environment.execTime.init);
	time.push_back(environment.execTime.iter);
	time.push_back(environment.execTime.cut);
	time.push_back(environment.execTime.initIter+environment.execTime.cut);


	List result;
	if(environment.numberShuffles > 0){
		// structure the output
	    result = List::create(
	    	_["confData"] = confVect,
		 	_["adjMatrix"] = adjMatrix,
	        _["edges"] = edgesMatrix,
	        _["time"] = time,
			_["interrupted"] = false
	    ) ;
	} else {
		result = List::create(
		 	_["adjMatrix"] = adjMatrix,
	        _["edges"] = edgesMatrix,
	        _["time"] = time,
			_["interrupted"] = false
	    ) ;
	}
	
	deleteMemorySpace(environment, environment.m);
	deleteStruct(environment);

    return result;
}