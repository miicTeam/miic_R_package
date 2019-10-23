#include "skeleton.h"

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

	int maxConsistentIter = 10;
	//run the skeleton iteration phase if consistency is required
	environment.iterationStepEven = true;
	auto cycle_tracker = CycleTracker(environment);
	do{
		environment.iterationStepEven = !environment.iterationStepEven;
		//save the neighbours in the areNeighboursAfterIteration structure
		//and return back with structures at the moment of initialization
		for(int i = 0; i < environment.numNodes; i++){
			for(int j = 0; j < environment.numNodes; j++){
				environment.edges[i][j].areNeighboursAfterIteration = environment.edges[i][j].isConnected;
				environment.edges[i][j].isConnected = environment.edges[i][j].isConnectedAfterInitialization;
			}

		firstStepIteration(environment);
		if (environment.numNoMore == 0 && environment.numSearchMore == 0) {
			if (environment.isVerbose == true)
				cout << "# ------| Only phantom edges found.\n";
		} else if (environment.numSearchMore > 0) {
			//// Search for other Contributing node(s) (possible only for the edges still in 'searchMore', ie. 2)
			if( environment.isVerbose == true ){ cout << "\n# ---- Other Contributing node(s) ----\n\n"; }
			startTime = get_wall_time();

			skeletonIteration(environment);

			long double spentTime = (get_wall_time() - startTime);
			environment.execTime.iter = spentTime;
			environment.execTime.initIter = environment.execTime.init + environment.execTime.iter;
		}
		cout << "Number of edges: " << environment.numNoMore << endl;
		maxConsistentIter --;
	} while (environment.consistentPhase && !cycle_tracker.hasCycle());

	int union_n_edges = 0;
	for(int i = 0; i < environment.numNodes -1; i++){
		for(int j = i+1; j < environment.numNodes; j++){
			if(environment.edges[i][j].isConnected){
				union_n_edges ++;
			}
		}
	}
	environment.numNoMore = union_n_edges;


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

bool CycleTracker::hasCycle() {
    uint n_edge = env_.numNoMore;
    // before saving the current iteration, search among previous iterations
    // those with the same number of edges
    auto range = edge_index_map_.equal_range(n_edge);
    bool no_cycle_found = range.first == range.second;
    // indices of iteration that is possibly the end point of a cycle
    // example: suppose that iteration #1, #3, #6 have the same
    // number of edges as the current iteration, and each of them is a
    // possible start point of a cycle, therefore #2, #4, #7 are the
    // corresponding possible end points of the cycle.
    std::deque<uint> iter_indices;
    for (auto it = range.first; it != range.second; ++it)
        iter_indices.push_back(it->second + 1);
    // save the current iteration
    saveIteration();
    if (no_cycle_found)
        return false;
    // backtracking requires starting from the largest index first
    std::sort(iter_indices.begin(), iter_indices.end(), std::greater<uint>());
    // set of edges that are to be marked as connected
    std::set<uint> edges_union;
    // check if an edge is changed. vector is chosen over map for
    // quicker access and simpler syntax, at the cost of extra memory trace
    // and (possible) extra time complexity (in practice there are very few
    // changes between each pair of iterations).
    std::vector<uint> changed(env_.numNodes * (env_.numNodes - 1) / 2, 0);
    // vector keeping track of index of the iteration with maximum information
    // initialized to the last iteration (push_front)
    std::vector<EdgeInfo> max_info(env_.numNodes * (env_.numNodes - 1) / 2);
    // backtracking over iteration to get changed_edges
    uint cycle_size = 0;
    for (auto& iter : iterations_) {
        ++cycle_size;
        if (cycle_size > max_cycle_size) {
            // FIXME: decide what to do if cycle_size > max_cycle_size
            std::cout << "cycle size out of limit size: "
                << max_cycle_size << '\n';
            return true;
        }
        for (auto& k : iter.changed_edges) {
            edges_union.insert(k);
            // changed twice == unchanged
            changed[k] = 1 - changed[k];
            // set max_info for each edge
            if (max_info[k].is_empty ||
                    iter.edge_info_list[k].Ixy_ui - iter.edge_info_list[k].cplx
                    > max_info[k].Ixy_ui - max_info[k].cplx) {
                max_info[k] = iter.edge_info_list[k];
            }
        }
        if (iter.index != iter_indices.front())
            continue;
        iter_indices.pop_front();
        // if any edge has been changed
        if (std::any_of(changed.begin(), changed.end(),
                    [](uint j){ return j != 0;})) {
            // no cycle
            if (iter_indices.empty())
                break;
            continue;
        }
        for (auto& k : edges_union)
            setEdgeState(k, max_info[k]);

        std::cout << "cycle found of size " << cycle_size << std::endl;
        return true;
    }
    return false;
}

void CycleTracker::setEdgeState(edge_index_1d index, const EdgeInfo& edge_info) {
    std::pair<uint, uint> p = getEdgeIndex2D(index);
    env_.edges[p.first][p.second].isConnected = true;
    env_.edges[p.second][p.first].isConnected = true;
    EdgeStructure& eij = *(env_.edges[p.first][p.second].edgeStructure);
    EdgeStructure& eji = *(env_.edges[p.second][p.first].edgeStructure);
    // ugly but avoid touching Edge or EdgeStructure class
    eij.Ixy_ui = edge_info.Ixy_ui;
    eij.cplx = edge_info.cplx;
    eij.Nxy_ui = edge_info.Nxy_ui;

    eji.Ixy_ui = edge_info.Ixy_ui;
    eji.cplx = edge_info.cplx;
    eji.Nxy_ui = edge_info.Nxy_ui;
}
