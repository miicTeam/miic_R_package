#include "structure.h"
#include "computeEnsInformation.h"
#include "utilities.h"
#include <string>
#include <ctime>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unistd.h>

using namespace std;


/*
 * Initialize the edges of the network
 */
int initEdgeElt(Environment& environment, int i, int j){

	//// Compute the mutual information and the corresponding CPLX
	double* res = NULL;
	if(environment.columnAsContinuous[i] == 0 && environment.columnAsContinuous[j] == 0){
		res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx, environment.m);
		environment.edges[i][j].edgeStructure->Ixy_ui = res[1]; 
		environment.edges[i][j].edgeStructure->cplx = res[2];
		environment.edges[i][j].edgeStructure->Nxy_ui = res[0];
		free(res);
	} 

	else if(environment.columnAsGaussian[i] == 1 && environment.columnAsGaussian[j] == 1){
		res = corrMutInfo(environment, environment.dataDouble, NULL, 0, NULL, 0, i, j, -2);
		int N = environment.nSamples[i][j];
		environment.edges[i][j].edgeStructure->Ixy_ui = res[0]; 
		environment.edges[i][j].edgeStructure->cplx = 0.5 * (environment.edges[i][j].edgeStructure->ui_vect_idx.size() + 2) * log(N);
		environment.edges[i][j].edgeStructure->Nxy_ui = N;
		delete [] res;
		// cout << "GAUSS: " << environment.nodes[i].name << "-" << environment.nodes[j].name << "\t" <<  res[0] << "\t" <<   environment.edges[i][j].edgeStructure->cplx << "\n" << flush;

	}
	else {
		//cout<< "GUIDO: " << environment.nodes[i].name << "-" << environment.nodes[j].name << flush;

		res = computeEnsInformationContinuous(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx,environment.m);
		environment.edges[i][j].edgeStructure->Ixy_ui = res[1]; 
		environment.edges[i][j].edgeStructure->cplx = res[2];
		environment.edges[i][j].edgeStructure->Nxy_ui = res[0];
		delete [] res;

		//cout<< "GUIDO: " << environment.nodes[i].name << "-" << environment.nodes[j].name << "\t" <<  res[1] << "\t" <<   res[2] << "\n" << flush;

	} 


	if(environment.isVerbose) 
	{ 
		cout << "# --> Ixy_ui = " << environment.edges[i][j].edgeStructure->Ixy_ui/environment.edges[i][j].edgeStructure->Nxy_ui << "[Ixy_ui*Nxy_ui =" << environment.edges[i][j].edgeStructure->Ixy_ui << "]\n"
			 << "# --> Cplx = " << environment.edges[i][j].edgeStructure->cplx << "\n"
			 << "# --> Nxy_ui = " << environment.edges[i][j].edgeStructure->Nxy_ui << "\n"
			 << "# --> nbrEdges L = " << environment.l << "\n"
			 << "# --> nbrProp P = " << environment.numNodes << "\n";
	}

	double myTest = 0;
	string category;
	environment.edges[i][j].mutInfo = environment.edges[i][j].edgeStructure->Ixy_ui;
	environment.edges[i][j].cplx_noU = environment.edges[i][j].edgeStructure->cplx;
	
	if(environment.isNoInitEta)
		myTest = environment.edges[i][j].edgeStructure->Ixy_ui - environment.edges[i][j].edgeStructure->cplx;	 
	 else 
	 	myTest = environment.edges[i][j].edgeStructure->Ixy_ui - environment.edges[i][j].edgeStructure->cplx - environment.logEta;	 	
	
	//// set the edge status
	if(myTest <= 0){
		// the node is a phantom one
		environment.edges[i][j].edgeStructure->status = 1; 
		environment.edges[i][j].isConnected = 0;
		environment.edges[j][i].isConnected = 0;
		category = "phantom";
	} else {
		// the node is a searchMore one
		environment.edges[i][j].edgeStructure->status = 3; 
		environment.edges[i][j].isConnected = 1;
		environment.edges[j][i].isConnected = 1;
		category= "searchMore";
	}

	// cout << myTest << "\n";

	if(environment.isVerbose) { cout << "# --> Category = " << category << "\n" ; }

	return environment.edges[i][j].isConnected;
}


void* initEdgeEltThread(void* container){

	ContainerInit* cont = static_cast<ContainerInit*>(container);

	Environment* environment1 = cont->environment;
	Environment environment = *environment1;

	MemorySpace m = cont->m;

	int startI = cont->i;
	int startJ= cont->j;
	int steps= cont->steps;

	int i = startI;
	int j = startJ;

	for(int step = 0; step < steps && i < environment.numNodes && j < environment.numNodes; step++)
	{
		if(environment.edges[i][j].isConnected){
			// cout << "i: " << i << " j: " << j << " steps: " <<steps << endl;
			double* res;
			//// Compute the mutual information and the corresponding CPLX
			if(environment.columnAsContinuous[i] == 0 && environment.columnAsContinuous[j] == 0){
				res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx, m);
				environment.edges[i][j].edgeStructure->cplx = res[2];
				environment.edges[i][j].edgeStructure->Nxy_ui = res[0];
				environment.edges[i][j].edgeStructure->Ixy_ui = res[1]; 
				free(res);
			
			} else if(environment.columnAsGaussian[i] == 1 && environment.columnAsGaussian[j] == 1){
				res = corrMutInfo(environment, environment.dataDouble, NULL, 0, NULL, 0, i, j, -2);
				int N = environment.nSamples[i][j];
				environment.edges[i][j].edgeStructure->cplx = 0.5 * (environment.edges[i][j].edgeStructure->ui_vect_idx.size() + 2) * log(N);
				environment.edges[i][j].edgeStructure->Nxy_ui = N;
				environment.edges[i][j].edgeStructure->Ixy_ui = res[0]; 
				delete [] res;

				// cout << "GAUSS: " << environment.nodes[i].name << "-" << environment.nodes[j].name << "\t" <<  res[0] << "\t" <<   environment.edges[i][j].edgeStructure->cplx << "\n" << flush;

			} else {
				res = computeEnsInformationContinuous(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx,m);
				environment.edges[i][j].edgeStructure->Nxy_ui = res[0];
				environment.edges[i][j].edgeStructure->Ixy_ui = res[1]; 
				environment.edges[i][j].edgeStructure->cplx = res[2];
				delete [] res;
			}

			if(environment.isVerbose) 
			{ 
				cout << "# --> Ixy_ui = " << environment.edges[i][j].edgeStructure->Ixy_ui/environment.edges[i][j].edgeStructure->Nxy_ui << "[Ixy_ui*Nxy_ui =" << environment.edges[i][j].edgeStructure->Ixy_ui << "]\n"
					 << "# --> Cplx = " << environment.edges[i][j].edgeStructure->cplx << "\n"
					 << "# --> Nxy_ui = " << environment.edges[i][j].edgeStructure->Nxy_ui << "\n"
					 << "# --> nbrEdges L = " << environment.l << "\n"
					 << "# --> nbrProp P = " << environment.numNodes << "\n";
			}

			double myTest = 0;
			string category;

			environment.edges[i][j].mutInfo = environment.edges[i][j].edgeStructure->Ixy_ui;

			
			if(environment.isNoInitEta)
				myTest = environment.edges[i][j].edgeStructure->Ixy_ui - environment.edges[i][j].edgeStructure->cplx;	 
			 else 
			 	myTest = environment.edges[i][j].edgeStructure->Ixy_ui - environment.edges[i][j].edgeStructure->cplx - environment.logEta;	 	
			
			//// set the edge status
			if(myTest <= 0){
				// the node is a phantom one
				environment.edges[i][j].edgeStructure->status = 1; 
				environment.edges[i][j].isConnected = 0;
				environment.edges[j][i].isConnected = 0;
				category = "phantom";
			} else {
				// the node is a searchMore one
				environment.edges[i][j].edgeStructure->status = 3; 
				environment.edges[i][j].isConnected = 1;
				environment.edges[j][i].isConnected = 1;
				category= "searchMore";
			}

			if(environment.isVerbose)  { cout << "# --> Category = " << category << "\n" ; }
		}

		if(j + 1 < environment.numNodes){
			j++;
		} else {
			i++;
			j = i + 1;
		}
	}
	return NULL;
}


bool skeletonInitialization(Environment& environment){

	createMemorySpace(environment, environment.m);

	environment.oneLineMatrix = new int[environment.numSamples*environment.numNodes];
	for(int i = 0; i < environment.numSamples;i++){
		for(int j = 0; j < environment.numNodes;j++){
			// cout << j * environment.numSamples + i << " ";

			environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
		}
	}

	environment.countSearchMore = 0;

	int nthreadsMax = environment.nThreads;
	// cout << "nthreadsMax: " << nthreadsMax << endl;

	// IF THREADS
	if(nthreadsMax > 1) {
		int step;

		int length = environment.numNodes*(environment.numNodes-1)/2;
		step = (length)/(nthreadsMax);
		if((length)%(nthreadsMax) != 0)
			step ++;

		if((length)%(step) == 0)
			nthreadsMax = (length)/(step);
		else 
			nthreadsMax = ((length)/(step))+1;
		int iterator = 0;


		// if((environment.numNodes*(environment.numNodes-1)/2)%(nthreadsMax-1) == 0)
		// 		nthreadsMax--;

		cout << "Number of threads : " <<  nthreadsMax << endl;

		

		pthread_t* pt;
		pt = (pthread_t*)(malloc(nthreadsMax * sizeof(pthread_t)));
		ContainerInit* cInit = new ContainerInit[nthreadsMax];
		
		

		int counter = 0;
		
		for(int i = 0; i < environment.numNodes - 1; i++){
			for(int j = i + 1; j < environment.numNodes; j++){
				
				// create a structure for the nodes that need to store information about them
				environment.edges[i][j].edgeStructure = new EdgeStructure();
				environment.edges[j][i].edgeStructure = environment.edges[i][j].edgeStructure ;

				// initialize the structure
				environment.edges[j][i].edgeStructure->z_name_idx = -1;
				environment.edges[j][i].edgeStructure->status = -1;


				//reserve space for vectors
				// environment.edges[j][i].edgeStructure->ui_vect_idx.reserve(environment.numNodes);
				// environment.edges[j][i].edgeStructure->zi_vect_idx.reserve(environment.numNodes);

				if(counter%step == 0){
					// cout << "i: " << i << " j: " << j << "iter:"<< iterator << endl;

					cInit[iterator].environment = &environment;
					cInit[iterator].i=i;
					cInit[iterator].j=j;
					cInit[iterator].steps=step;
					createMemorySpace(environment, cInit[iterator].m);

					if (pthread_create(&pt[iterator], NULL, initEdgeEltThread, (void *)&cInit[iterator]) < 0)
					{ 
						return false;
			  		}

					iterator++;
				}

				counter++;
			}
		}

		// cout << "threadsMax" << nthreadsMax;

		for(int pos = 0; pos < nthreadsMax; pos++)
			pthread_join( pt[pos], NULL);

		for(int pos = 0; pos < nthreadsMax; pos++){
			deleteMemorySpace(environment, cInit[pos].m);
		}

		free(pt);
		delete [] cInit;



		for(int i = 0; i < environment.numNodes - 1; i++){
			for(int j = i + 1; j < environment.numNodes; j++){
				if(environment.edges[i][j].isConnected)
					environment.countSearchMore++;
			}
		}
	} else {
		double sum = 0;
		for(int i = 0; i < environment.numNodes - 1; i++){
			for(int j = i + 1; j < environment.numNodes; j++){
				if(environment.isVerbose) { cout << "\n# Edge " << environment.nodes[i].name << "," << environment.nodes[j].name << "\n" ; }
				
				// create a structure for the nodes that need to store information about them
				environment.edges[i][j].edgeStructure = new EdgeStructure();


				environment.edges[j][i].edgeStructure = environment.edges[i][j].edgeStructure ;
				// cout << environment.edges[i][j].edgeStructure << endl;exit(0);
				// cout << environment.edges[j][i].edgeStructure << endl;exit(0);

				// initialize the structure
				environment.edges[j][i].edgeStructure->z_name_idx = -1;
				environment.edges[j][i].edgeStructure->status = -1;

				//reserve space for vectors
				// environment.edges[j][i].edgeStructure->ui_vect_idx.reserve(environment.numNodes);
				// environment.edges[j][i].edgeStructure->zi_vect_idx.reserve(environment.numNodes);
				if(environment.edges[i][j].isConnected){
					if(initEdgeElt(environment, i, j) == 1){
						environment.countSearchMore++;			
					}
				}
				// sum += sizeof(environment.edges[i][j].edgeStructure) +  sizeof(environment.edges[i][j].edgeStructure->z_name_idx) + 
				// sizeof(environment.edges[i][j].edgeStructure->ui_vect_idx) +  sizeof(environment.edges[i][j].edgeStructure->zi_vect_idx) + 
				// sizeof(environment.edges[i][j].edgeStructure->Nxyz_ui) +  sizeof(environment.edges[i][j].edgeStructure->Rxyz_ui) + 
				// sizeof(environment.edges[i][j].edgeStructure->Ixy_ui) +  sizeof(environment.edges[i][j].edgeStructure->cplx) + 
				// sizeof(environment.edges[i][j].edgeStructure->Nxy_ui) +  sizeof(environment.edges[i][j].edgeStructure->status) +
				// sizeof(environment.edges[i][j].edgeStructure->indexStruct) +  sizeof(environment.edges[i][j].edgeStructure->edgesInSpeTpl_list);
			}
		}
		// cout << "Structs edges initedgeEld: " << sum/1024/1024/1024 << "Gb" << endl;
	}

	for(int i = 0; i < environment.numNodes; i++){
		for(int j = 0; j < environment.numNodes; j++){
			environment.edges[i][j].isConnectedAfterInitialization = environment.edges[i][j].isConnected;
		}
	}
	// if(environment.atLeastOneContinuous == 1)
	// 	environment.nThreads = 1;

	return true;
}
