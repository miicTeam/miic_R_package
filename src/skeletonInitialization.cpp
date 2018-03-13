#include "structure.h"
#include "computeEnsInformation.h"
#include "utilities.h"
#include <string>
#include <ctime>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdio.h>

using namespace std;
/*
 * Sort the ranks of the function
 */

bool SortFunction1(const XJAddress* a, const XJAddress* b, const Environment& environment) {
	 return environment.edges[a->i][a->j].edgeStructure->Rxyz_ui > environment.edges[b->i][b->j].edgeStructure->Rxyz_ui;
}

class sorter1 {
	  Environment& environment;
		public:
	  sorter1(Environment& env) : environment(env) {}
	  bool operator()(XJAddress const* o1, XJAddress const* o2) const {
			return SortFunction1(o1, o2, environment );
	  }
};

/*
 * Search for the neighbors of x or y, without duplicates
 */
void searchAndSetZi(Environment& environment, vector<int>& vec, int& numZiPos, const int posX, const int posY){
	//search for neighbors of x
	numZiPos = 0;

	if(environment.isLatent){
		for(int c = 0; c < environment.numNodes; c++){
			//chech if the node is a neighbour of y
			if(c!=posX && c!= posY){
				vec.push_back(c);
				numZiPos++;
			}
		}
	}
	else {

		for(int c = 0; c < environment.numNodes; c++){
		//chech if the node is a neighbour of y
			if(environment.edges[posX][c].isConnected && c!=posX && c!= posY){
				vec.push_back(c);
				numZiPos++;
			}
		}

		for(int c = 0;  c < environment.numNodes; c++){
		//chech if the node is a neighbour of x
			if(environment.edges[posY][c].isConnected && !environment.edges[posX][c].isConnected && c!=posX && c!= posY){
				vec.push_back(c);
				numZiPos++;
			}
		}
	}
}

/*
 * Initialize the edges of the network
 */
int initEdgeElt(Environment& environment, int i, int j){
	//// Compute the mutual information and the corresponding CPLX
	//double* res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx);
	double* res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx);

	environment.edges[i][j].edgeStructure->Ixy_ui = res[1];

	environment.edges[i][j].edgeStructure->cplx = res[2];

	environment.edges[i][j].edgeStructure->Nxy_ui = res[0];

	double myTest = 0;
	string category;

	free(res);

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

	return environment.edges[i][j].isConnected;
}

/*
 * Initialize the edges of the network, for threads
 */
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
			//// Compute the mutual information and the corresponding CPLX
			//double* res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx);
			double* res = computeEnsInformationNewThread(environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx, m);

			environment.edges[i][j].edgeStructure->Ixy_ui = res[1];
			environment.edges[i][j].edgeStructure->cplx = res[2];
			environment.edges[i][j].edgeStructure->Nxy_ui = res[0];

			free(res);

			double myTest = 0;
			string category;

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
		}

		if(j + 1 < environment.numNodes){
			j++;
		} else {
			i++;
			j = i + 1;
		}
	}
	return (void*) 0;
}


bool skeletonInitialization(Environment& environment){

	createMemorySpace(environment, environment.m);

	environment.oneLineMatrix = new int[environment.numSamples*environment.numNodes];
	for(int i = 0; i < environment.numSamples;i++){
		for(int j = 0; j < environment.numNodes;j++){

			environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
		}
	}

	long double startTime = time(0)*1000;
	int countSearchMore = 0;

	int nthreadsMax = environment.nThreads;

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
				environment.edges[j][i].edgeStructure->status = 1;

				if(counter%step == 0){
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
					countSearchMore++;
			}
		}
	} else {
		for(int i = 0; i < environment.numNodes - 1; i++){
			for(int j = i + 1; j < environment.numNodes; j++){

				// create a structure for the nodes that need to store information about them
				environment.edges[i][j].edgeStructure = new EdgeStructure();
				environment.edges[j][i].edgeStructure = environment.edges[i][j].edgeStructure ;

				// initialize the structure
				environment.edges[j][i].edgeStructure->z_name_idx = -1;
				environment.edges[j][i].edgeStructure->status = -1;

				if(environment.edges[i][j].isConnected){
					if(initEdgeElt(environment, i, j) == 1){
						countSearchMore++;
					}
				}
			}
		}
	}

	environment.numSearchMore = countSearchMore;
	environment.numNoMore= 0;

	//set the diagonal of the adj matrix to 0
	for(int i = 0; i < environment.numNodes; i++){
		environment.edges[i][i].isConnected = 0;
	}

	// create and fill the searchMoreAddress struct, that keep track of i and j positions of searchMore Edges
	for(int i = 0; i < environment.numNodes - 1; i++){
		for(int j = i + 1; j < environment.numNodes; j++){
			if(environment.edges[i][j].isConnected){
				XJAddress* s = new XJAddress();
				s->i=i;
				s->j=j;
				environment.searchMoreAddress.push_back(s);
			}
		}
	}

	if(nthreadsMax > 1){
		// create meory space for thread execution
		environment.memoryThreads = new ContainerMemory[environment.nThreads];
		for(int i = 0; i < environment.nThreads; i++){
			createMemorySpaceThreads(environment, environment.memoryThreads[i]);
		}
	}

	// if there exists serachMore edges
	if( countSearchMore > 0 ) {
		startTime = time(0)*1000;

		//// Define the pairs among which the candidate {zi} should be searched
		for(int i = 0; i < environment.numSearchMore; i++){
			int posX = environment.searchMoreAddress[i]->i;
			int posY = environment.searchMoreAddress[i]->j;

			//// Find all zi, such that xzi or yzi is not a phantom
			// search in the row of the edge matrix

			int numZiPos = 0;
			searchAndSetZi(environment, environment.edges[posX][posY].edgeStructure->zi_vect_idx, numZiPos,posX, posY);
		}


		if(nthreadsMax > 1){
			int step1;

			step1 = environment.numSearchMore/(nthreadsMax);
			if(environment.numSearchMore%(nthreadsMax) != 0)
				step1++;

			if(environment.numSearchMore%(step1) == 0)
				nthreadsMax=environment.numSearchMore/(step1);
			else
				nthreadsMax=(environment.numSearchMore/(step1) + 1);
			int iterator1 = 0;

			// if(environment.numSearchMore%(nthreadsMax-1) == 0)
			// 	nthreadsMax--;

			pthread_t* pt1;
			pt1 = (pthread_t*)(malloc(nthreadsMax * sizeof(pthread_t)));
			Container* c1 = new Container[nthreadsMax];



			int i = 0;
			while(i < environment.numSearchMore){
				c1[iterator1].environment = &environment;
				c1[iterator1].start=i;
				c1[iterator1].stop= i + step1;
				createMemorySpace(environment, c1[iterator1].m);


				//// Search for new contributing node and its rank
				if (pthread_create(&pt1[iterator1], NULL, SearchForNewContributingNodeAndItsRankThread, (void *)&c1[iterator1]) < 0)
				{
		  			return false;
		  		}

				iterator1++;
				i += step1;
			}


			for(int pos = 0; pos < nthreadsMax; pos++)
				pthread_join( pt1[pos], NULL);

			for(int pos = 0; pos < nthreadsMax; pos++){
				deleteMemorySpace(environment, c1[pos].m);
			}

			free(pt1);
			delete [] c1;


			for(int i = 0; i < environment.numSearchMore; i++){
				int posX = environment.searchMoreAddress[i]->i;
				int posY = environment.searchMoreAddress[i]->j;
				if(environment.edges[posX][posY].edgeStructure->z_name_idx == -1)
				{
					environment.noMoreAddress.push_back(environment.searchMoreAddress[i]);
					environment.numNoMore++;

					//// Remove the element from searchMore
					environment.searchMoreAddress.erase(environment.searchMoreAddress.begin() + i);
					environment.numSearchMore--;
					i--;

					//update the status
					environment.edges[posX][posY].edgeStructure->status = 3;
				}

			} //// end for non-ph edges
		} else {
			for(int i = 0; i < environment.numSearchMore; i++){
				int posX = environment.searchMoreAddress[i]->i;
				int posY = environment.searchMoreAddress[i]->j;
				if(environment.edges[posX][posY].edgeStructure->zi_vect_idx.size() > 0 ){
					//// Search for new contributing node and its rank
					SearchForNewContributingNodeAndItsRank(environment, posX, posY);
				}


				if(environment.edges[posX][posY].edgeStructure->z_name_idx == -1)
				{
					//// Put this edge element to "noMore" and update the vector of status
					environment.noMoreAddress.push_back(environment.searchMoreAddress[i]);
					environment.numNoMore++;

					//// Remove the element from searchMore
					environment.searchMoreAddress.erase(environment.searchMoreAddress.begin() + i);
					environment.numSearchMore--;
					i--;

					//update the status
					environment.edges[posX][posY].edgeStructure->status = 3;
				}
			}
		}
		
		// sort the ranks
		std::sort(environment.searchMoreAddress.begin(), environment.searchMoreAddress.end(), sorter1(environment));
	}

	return true;
}
