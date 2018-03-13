#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include <string.h>
#include "confidenceCut.h"

#include "structure.h"
#include "utilities.h"
#include "skeletonInitialization.h"
#include "skeletonIteration.h"

#include <Rcpp.h>

using namespace std;

using namespace Rcpp;

extern "C" SEXP skeleton(SEXP inputDataR, SEXP numNodesR, SEXP nThreadsR, SEXP blackBoxR, SEXP effNR, SEXP cplxR, SEXP etaR,
	SEXP isLatentR, SEXP isTplReuseR, SEXP isK23R, SEXP isDegeneracyR, SEXP isNoInitEtaR, SEXP continuousR,
	SEXP confidenceShuffleR, SEXP confidenceThresholdR, SEXP verboseR){

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

	environment.vectorData = Rcpp::as< vector <string> > (inputDataR);

	vector<string> v;
	v = Rcpp::as< vector <string> > (blackBoxR);


	string cplx;

	environment.effN = Rcpp::as<int> (effNR);
	environment.typeOfData = Rcpp::as<int> (continuousR);
	cplx = Rcpp::as<string> (cplxR);
	environment.eta = Rcpp::as<int> (etaR);
	environment.shuffle = 0;
	environment.isLatent = Rcpp::as<bool> (isLatentR);
	environment.isTplReuse = Rcpp::as<bool> (isTplReuseR);
	environment.isK23 = Rcpp::as<bool> (isK23R);
	environment.isDegeneracy = Rcpp::as<bool> (isDegeneracyR);
	environment.isNoInitEta = Rcpp::as<bool> (isNoInitEtaR);

	environment.confidenceThreshold = Rcpp::as<double> (confidenceThresholdR);
	environment.numberShuffles = Rcpp::as<int> (confidenceShuffleR);

	environment.isVerbose = Rcpp::as<bool> (verboseR);

	if(cplx.compare("nml") == 0)
		environment.cplx = 1;
	else if(cplx.compare("mdl") == 0)
		environment.cplx = 0;

	clock_t startTime;


	vector< vector <string> > retShuffle;
	vector< vector <string> > retShuffleAvg;
	vector <string> shfVec;

	string filePath;

	// set iterator to the bound, i.e only 1 iteration
	environment.etaIterator = environment.eta;
	environment.shuffleIterator = environment.shuffle;

	// set the environment
	setEnvironment(environment);

	if(v.size() > 1)
		readBlackbox1(v, environment);

	startTime = std::clock();

	// ----

	if(!skeletonInitialization(environment))
	{
		// structure the output
	    List result = List::create(
		 	_["error"] = "error during skeleton initialization"

	    ) ;
	    return result;
	}
			// ----
	long double spentTime = (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000) /1000;
	environment.execTime.init = spentTime;


	if( environment.numNoMore == 0 && environment.numSearchMore == 0 ) {
		// if( environment.isVerbose == true ){ cout << "# ------| Only phantom edges found.\n"; }
	} else if( environment.numSearchMore > 0 ) {

		//// Search for other Contributing node(s) (possible only for the edges still in 'searchMore', ie. 2)

		startTime = std::clock();


		if(!skeletonIteration(environment))
		{
			// structure the output
		    List result = List::create(
			 	_["error"] = "error during skeleton iteration"

		    ) ;
		    return result;
		}

		long double spentTime = (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000) /1000;
		environment.execTime.iter = spentTime;
		environment.execTime.initIter = environment.execTime.init + environment.execTime.iter;
	}

	startTime = std::clock();
	vector< vector <string> >confVect;
	if(environment.numberShuffles > 0){
		confVect = confidenceCut(environment);
		spentTime = (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000) /1000;
		environment.execTime.cut = spentTime;
	} else {
		environment.execTime.cut = 0;
	}

	

	if( environment.numNoMore > 0)
		edgesMatrix = saveEdgesListAsTable1(environment);


	vector< vector <string> >  adjMatrix;
	adjMatrix = getAdjMatrix(environment);

	vector<double>time;
	time.push_back(environment.execTime.init);
	time.push_back(environment.execTime.iter);
	time.push_back(environment.execTime.cut);
	time.push_back(environment.execTime.initIter+environment.execTime.cut);


	if(environment.nThreads > 1){
		// delete meory space for thread execution
		for(int i = 0; i < environment.nThreads; i++){
			deleteMemorySpaceThreads(environment, environment.memoryThreads[i]);
		}
		delete [] environment.memoryThreads ;
	}
	
	
	List result;
	if(environment.numberShuffles > 0){
		// structure the output
	    result = List::create(
	    	_["confData"] = confVect,
		 	_["adjMatrix"] = adjMatrix,
	        _["edges"] = edgesMatrix,
	        _["time"] = time
	    ) ;
	} else {
		result = List::create(
		 	_["adjMatrix"] = adjMatrix,
	        _["edges"] = edgesMatrix,
	        _["time"] = time
	    ) ;
	}
	
	deleteMemorySpace(environment, environment.m);
	deleteStruct(environment);


    return result;
}
