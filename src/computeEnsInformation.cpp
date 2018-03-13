#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "structure.h"
#include <iostream>
#include "utilities.h"
#include "computeInfo_interface.h"
#include <cstdlib>

using namespace std;

/*
/ This functions prepare all data for the call of getAllInfoNEW or getAllInfoNEWThreads and returns its result
*/
double* computeEnsInformationNew(Environment& environment, int* myCond, int myNbrUi,  int* myZi, int myNbrZi, int myZiPos,
	const int myVarIdxX, const int myVaridxY, const int cplx){

	int* posArray = new int[2 + environment.edges[myVarIdxX][myVaridxY].edgeStructure->ui_vect_idx.size()];
	posArray[0] = myVarIdxX;
	posArray[1] = myVaridxY;

	if(myCond != NULL){
		if(myNbrUi  > 0)
		{
			/// The index of the variables in this dataset
			for (int i = 0; i < myNbrUi; ++i)
				posArray[i+2] = myCond[i];
		}
	}

	//// Compute the mutual information
	double* res_new;
	if(myNbrZi>5 && myNbrUi>0 && environment.nThreads > 1){
	 	int nthreads = int(myNbrZi /5) +1;
	 	if(nthreads > environment.nThreads)
	 		nthreads = environment.nThreads;

	 	res_new = getAllInfoNEWThreads(environment.oneLineMatrix, environment.allLevels, posArray,
	 	myNbrUi, myZi, myNbrZi, myZiPos, environment.numSamples, environment.effN, cplx, environment.isK23, environment.c2terms,
	 	&environment.m, nthreads, environment.memoryThreads);
	}
	else {
	 	res_new = getAllInfoNEW(environment.oneLineMatrix, environment.allLevels, posArray,
	 	myNbrUi, myZi, myNbrZi, myZiPos, environment.numSamples, environment.effN, cplx, environment.isK23, environment.c2terms, &environment.m);
	}
	int nbrRetValues = 3;

	// If nbrZi > 0, return {nSample[z1]*I(..|{ui})[z1], NML(..|{ui})[z1], nSample[z1],nSample[z2]*I(..|{ui})[z2], NML(..|{ui})[z2], nSample[z2], ... }
	if( myNbrZi > 0 )
		 nbrRetValues = 9;

 	for(int i = 0; i < nbrRetValues; i++){
		if(res_new[i] > -0.0000000001 && res_new[i] < 0.0000000001){
			res_new[i] = 0.0;
		}
	}

	delete [] posArray;

	return res_new;
}


/*
 * Wrapper for the evaluation of the mutual information, multithreaded. It calls the principal function passing the correct parameters
 */
double* computeEnsInformationNewThread(Environment& environment, int* myCond, int myNbrUi,  int* myZi, int myNbrZi, int myZiPos,
	const int myVarIdxX, const int myVaridxY, const int cplx, MemorySpace m){

	int* posArray = new int[2 + environment.edges[myVarIdxX][myVaridxY].edgeStructure->ui_vect_idx.size()];
	posArray[0] = myVarIdxX;
	posArray[1] = myVaridxY;

	if(myCond != NULL){
		if(myNbrUi  > 0)
		{
			/// The index of the variables in this dataset
			for (int i = 0; i < myNbrUi; ++i)
				posArray[i+2] = myCond[i];
		}
	}

	double *res_new = getAllInfoNEW(environment.oneLineMatrix, environment.allLevels, posArray,
	 	myNbrUi, myZi, myNbrZi, myZiPos, environment.numSamples, environment.effN, cplx, environment.isK23, environment.c2terms, &m);

	int nbrRetValues = 3;

	// If nbrZi > 0, return {nSample[z1]*I(..|{ui})[z1], NML(..|{ui})[z1], nSample[z1],nSample[z2]*I(..|{ui})[z2], NML(..|{ui})[z2], nSample[z2], ... }
	if( myNbrZi > 0 )
		 nbrRetValues = 9;

 	for(int i = 0; i < nbrRetValues; i++){
		if(res_new[i] > -0.0000000001 && res_new[i] < 0.0000000001){
			res_new[i] = 0.0;
		}
	}

	delete [] posArray;

	return res_new;
}


/*
 * If zx and zy are phantom, remove the z from the zi
 */
void removeifNA(Environment& environment, vector<int>& vec, const int posX, const int posY){
	vector<int>::iterator it = vec.begin();

	while(it != vec.end()) {
		if(*it == -1){
			vec.erase(it);
		}
		else
			++it;
	}
}
/*
 * If zx and zy are phantom, remove the z from the zi
 */
void removeifBothPhantomAndNA(Environment& environment, vector<int>& vec, const int posX, const int posY){
	vector<int>::iterator it = vec.begin();

	while(it != vec.end()) {
		if(!environment.edges[posX][*it].isConnected && !environment.edges[posY][*it].isConnected){
			vec.erase(it);
		}
		else
			++it;
	}
}

/*
 * search for the best z and find the rank
 */
bool SearchForNewContributingNodeAndItsRank(Environment& environment, const int posX, const int posY) {

	if(environment.edges[posX][posY].edgeStructure->zi_vect_idx.size() == 0)
	 	return true;

	//// If needed, remove the NA (-1) elements
	removeifNA(environment, environment.edges[posX][posY].edgeStructure->zi_vect_idx,
	  	posX, posY);

	if(!environment.isLatent)
		removeifBothPhantomAndNA(environment, environment.edges[posX][posY].edgeStructure->zi_vect_idx,
	  	posX, posY);

	int nbrZi = environment.edges[posX][posY].edgeStructure->zi_vect_idx.size();

	if(nbrZi == 0)
	 	return true;

	int* ui;
	int* zi;

	if(environment.edges[posX][posY].edgeStructure->ui_vect_idx.empty())
		ui = NULL;
	else
		ui = &environment.edges[posX][posY].edgeStructure->ui_vect_idx[0];

	if(environment.edges[posX][posY].edgeStructure->zi_vect_idx.empty())
		zi = NULL;
	else
		zi = &environment.edges[posX][posY].edgeStructure->zi_vect_idx[0];


	int argEnsInfo = -1;
	if(environment.isK23 == true)
		argEnsInfo = environment.cplx;

	double* vect = computeEnsInformationNew(environment, ui, environment.edges[posX][posY].edgeStructure->ui_vect_idx.size(), zi, environment.edges[posX][posY].edgeStructure->zi_vect_idx.size(), environment.edges[posX][posY].edgeStructure->ui_vect_idx.size()+2,  posX, posY, argEnsInfo);
	//// There can be more than one zi with the same rank... so, arbitrarly take the first one
	if(vect[6] - environment.edges[posX][posY].edgeStructure->Rxyz_ui > 0 ){
		//// The order matters: set first the z.name.idx, than get the corresponding zi from the original vect
		//// Doing this way, we make sure that the z.name has the right bin xyzi key
		environment.edges[posX][posY].edgeStructure->z_name_idx = vect[3];
		//environment.edges[posX][posY].edgeStructure->z_name = environment.nodes[z_name_idx];

		environment.edges[posX][posY].edgeStructure->Rxyz_ui = vect[6];
		environment.edges[posX][posY].edgeStructure->Nxyz_ui = vect[6];
	}

	free(vect);
	
	return true;
}


 /*
 * search for the best z and find the rank, threaded
 */
void* SearchForNewContributingNodeAndItsRankThread(void* container){

	Container* cont = static_cast<Container*>(container);

	Environment* environment1 = cont->environment;
	Environment environment = *environment1;

	int start = cont->start;
	int stop= cont->stop;

	MemorySpace m = cont->m;


	// // Environment& environment, const int posX, const int posY) {
	for(int pos = start; pos < stop && pos < environment.numSearchMore; pos++){
		int posX = environment.searchMoreAddress[pos]->i;
		int posY = environment.searchMoreAddress[pos]->j;

		if(environment.edges[posX][posY].edgeStructure->zi_vect_idx.size() != 0)
		{
			//// If needed, remove the NA (-1) elements
			removeifNA(environment, environment.edges[posX][posY].edgeStructure->zi_vect_idx,
			  	posX, posY);

			if(!environment.isLatent)
				removeifBothPhantomAndNA(environment, environment.edges[posX][posY].edgeStructure->zi_vect_idx,
			  	posX, posY);

			int nbrZi = environment.edges[posX][posY].edgeStructure->zi_vect_idx.size();

			if(nbrZi != 0)
			{

				int* ui;
				int* zi;

				if(environment.edges[posX][posY].edgeStructure->ui_vect_idx.empty())
					ui = NULL;
				else
					ui = &environment.edges[posX][posY].edgeStructure->ui_vect_idx[0];

				if(environment.edges[posX][posY].edgeStructure->zi_vect_idx.empty())
					zi = NULL;
				else
					zi = &environment.edges[posX][posY].edgeStructure->zi_vect_idx[0];

				int argEnsInfo = -1;
				if(environment.isK23 == true)
					argEnsInfo = environment.cplx;
				double* vect = computeEnsInformationNewThread(environment, ui, environment.edges[posX][posY].edgeStructure->ui_vect_idx.size(), zi,
					environment.edges[posX][posY].edgeStructure->zi_vect_idx.size(), environment.edges[posX][posY].edgeStructure->ui_vect_idx.size()+2,
					posX, posY, argEnsInfo, m);

				//// There can be more than one zi with the same rank... so, arbitrarly take the first one
				if(vect[6] - environment.edges[posX][posY].edgeStructure->Rxyz_ui > 0 ){

					//// The order matters: set first the z.name.idx, than get the corresponding zi from the original vect
					//// Doing this way, we make sure that the z.name has the right bin xyzi key
					environment.edges[posX][posY].edgeStructure->z_name_idx = vect[3];
					environment.edges[posX][posY].edgeStructure->Rxyz_ui = vect[6];
					environment.edges[posX][posY].edgeStructure->Nxyz_ui = vect[3];
				}

				free(vect);
			}
		}
	}

	return (void*) 0;
}
