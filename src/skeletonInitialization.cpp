#include "structure.h"
#include "computeEnsInformation.h"
#include "utilities.h"
#include <string>
#include <ctime>
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;


/*
 * Initialize the edges of the network
 */
int initEdgeElt(Environment& environment, int i, int j, MemorySpace& m){

	//// Compute the mutual information and the corresponding CPLX
	double* res = NULL;
	if(environment.columnAsContinuous[i] == 0 && environment.columnAsContinuous[j] == 0){
		res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx, m);
		environment.edges[i][j].shared_info->Ixy_ui = res[1];
		environment.edges[i][j].shared_info->cplx = res[2];
		environment.edges[i][j].shared_info->Nxy_ui = res[0];
		free(res);
	}

	else if(environment.columnAsGaussian[i] == 1 && environment.columnAsGaussian[j] == 1){
		res = corrMutInfo(environment, environment.dataDouble, NULL, 0, NULL, 0, i, j, -2);
		int N = environment.nSamples[i][j];
		environment.edges[i][j].shared_info->Ixy_ui = res[0];
		environment.edges[i][j].shared_info->cplx = 0.5 * (environment.edges[i][j].shared_info->ui_vect_idx.size() + 2) * log(N);
		environment.edges[i][j].shared_info->Nxy_ui = N;
		delete [] res;
		// cout << "GAUSS: " << environment.nodes[i].name << "-" << environment.nodes[j].name << "\t" <<  res[0] << "\t" <<   environment.edges[i][j].shared_info->cplx << "\n" << flush;

	}
	else {
		res = computeEnsInformationContinuous(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx, m);
		environment.edges[i][j].shared_info->Ixy_ui = res[1];
		environment.edges[i][j].shared_info->cplx = res[2];
		environment.edges[i][j].shared_info->Nxy_ui = res[0];
		delete [] res;

		//cout<< "GUIDO: " << environment.nodes[i].name << "-" << environment.nodes[j].name << "\t" <<  res[1] << "\t" <<   res[2] << "\n" << flush;

	}


	if(environment.isVerbose)
	{
		cout << "# --> Ixy_ui = " << environment.edges[i][j].shared_info->Ixy_ui/environment.edges[i][j].shared_info->Nxy_ui << "[Ixy_ui*Nxy_ui =" << environment.edges[i][j].shared_info->Ixy_ui << "]\n"
			 << "# --> Cplx = " << environment.edges[i][j].shared_info->cplx << "\n"
			 << "# --> Nxy_ui = " << environment.edges[i][j].shared_info->Nxy_ui << "\n"
			 << "# --> nbrProp P = " << environment.numNodes << "\n";
	}

	double myTest = 0;
	string category;
	environment.edges[i][j].shared_info->mutInfo = environment.edges[i][j].shared_info->Ixy_ui;
	environment.edges[i][j].shared_info->cplx_noU = environment.edges[i][j].shared_info->cplx;
	environment.edges[i][j].shared_info->Nxy = environment.edges[i][j].shared_info->Nxy_ui;

	if(environment.isNoInitEta)
		myTest = environment.edges[i][j].shared_info->Ixy_ui - environment.edges[i][j].shared_info->cplx;
	 else
	 	myTest = environment.edges[i][j].shared_info->Ixy_ui - environment.edges[i][j].shared_info->cplx - environment.logEta;

	if (myTest <= 0) {
		// Unconditional independence
		environment.edges[i][j].shared_info->connected = 1;
		environment.edges[i][j].status = 0;
		environment.edges[j][i].status = 0;
		category = "phantom";
	} else {
		environment.edges[i][j].shared_info->connected = 3;
		environment.edges[i][j].status = 1;
		environment.edges[j][i].status = 1;
		category= "searchMore";
	}

	if(environment.isVerbose)
        cout << "# --> Category = " << category << "\n";

	return environment.edges[i][j].status;
}


bool skeletonInitialization(Environment& environment){

	environment.oneLineMatrix = new int[environment.numSamples*environment.numNodes];
	for(uint i = 0; i < environment.numSamples;i++){
		for(uint j = 0; j < environment.numNodes;j++){
			// cout << j * environment.numSamples + i << " ";

			environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
		}
	}

	int threadnum = 0;

	cout << "Computing pairwise independencies...";
	fflush (stdout);

	bool interrupt = false;

	#ifdef _OPENMP
	#pragma omp parallel for shared(interrupt) firstprivate(threadnum) schedule(dynamic)
	#endif
	for(uint i = 0; i < environment.numNodes - 1; i++){

		if (interrupt) {
			continue; // will continue until out of for loop
		}
		#ifdef _OPENMP
			threadnum = omp_get_thread_num();
		#endif
			if(checkInterrupt(threadnum == 0)) {
				interrupt = true;
			}
		for(uint j = i + 1; j < environment.numNodes && !interrupt; j++){

			if(environment.isVerbose) { cout << "\n# Edge " << environment.nodes[i].name << "," << environment.nodes[j].name << "\n" ; }

			// create a structure for the nodes that need to store information about them
			environment.edges[i][j].shared_info = std::make_shared<EdgeSharedInfo>();

			environment.edges[j][i].shared_info = environment.edges[i][j].shared_info ;
			// Reserve space for vectors
			if(environment.edges[i][j].status){
				if(initEdgeElt(environment, i, j, environment.memoryThreads[threadnum]) == 1){
					#ifdef _OPENMP
					#pragma omp critical
					#endif
					environment.numSearchMore++;
				}
			}
		}
	}

	if(interrupt) return false;

	cout << " done." << endl;
	for(uint i = 0; i < environment.numNodes; i++){
		for(uint j = 0; j < environment.numNodes; j++){
			environment.edges[i][j].status_init = environment.edges[i][j].status;
		}
	}

	return true;
}
