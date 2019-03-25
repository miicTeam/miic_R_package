#include <math.h>
#include <algorithm>
#include <iostream>
#include <sstream>

#include "structure.h"
#include "computeEnsInformation.h"
#include "utilities.h"
using namespace std;

void shuffle_lookup(int *array, int *array2, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
          array2[t] = i;
        }
    }
}

//void saveConfidenceVector(Environment& environment, int** inferredEdges_tab, double* confVect,  string outDir, string slash){
//	
//
//	stringstream ss;
//	ss.str("");
//	ss << outDir << slash << "confRatios.txt";
//	string filename = ss.str();
//
//	if(environment.isVerbose)
//		cout << "Saving confidence matrix\n";
//	ofstream output;
//	output.open(filename.c_str());
//	output << "x" << "\t" << "y" << "\t" << "confidence_ratio" << endl;
//	for (int i=0; i < environment.numNoMore; i++){
//		int X = inferredEdges_tab[i][0] ;
// 		int Y = inferredEdges_tab[i][1] ;
//
//		output << environment.nodes[X].name << "\t" << environment.nodes[Y].name << "\t" << confVect[i] << endl;
//	}
//	
//	output.close();
//}

bool SortFunctionNoMore2(const XJAddress* a, const XJAddress* b, const Environment& environment) {
	 return (environment.edges[a->i][a->j].edgeStructure->Ixy_ui) > (environment.edges[b->i][b->j].edgeStructure->Ixy_ui);
}

class sorterNoMore2 {
	  Environment environment;
		public:
	  sorterNoMore2(Environment& env) : environment(env) {}
	  bool operator()(XJAddress const* o1, XJAddress const* o2) const {
			return SortFunctionNoMore2(o1, o2, environment );
	  }
};

// class sort_indicesLookup
// {
//    private:
//      int* parr;
//    public:
//      sort_indicesLookup(int* parr) : parr(parr) {}
//      bool operator()(int i, int j) const { return parr[i]<parr[j]; }
// };

// void sort2arraysLookup(int len, int a[], int brr[]){

//     std::sort(a, brr+len, sort_indicesLookup(a));

// }

vector< vector <string> > confidenceCut(Environment& environment){

	int** safe_state;
	int** safe_stateIdx;
	double** safe_stateDouble;

	int* lookup = new int[environment.numSamples];
	int* lookup2 = new int[environment.numSamples];
	for(int i=0; i<environment.numSamples; i++)
		lookup[i] = i;

	double noMore = environment.numNoMore;

	// Allocate the true edges table
	int** inferredEdges_tab;
	inferredEdges_tab = new int*[environment.numNoMore];
	for(int i=0; i < environment.numNoMore; i++)
		inferredEdges_tab[i] = new int[2];

	int pos = 0;

	for(int i = 0; i < environment.numNodes -1; i++){
		for(int j = i+1; j < environment.numNodes; j++){
			if(environment.edges[i][j].isConnected){
				inferredEdges_tab[pos][0] = i;
				inferredEdges_tab[pos][1] = j;
				pos++;
			}
		}
	}

	if(environment.atLeastTwoGaussian == 1){
		safe_stateDouble = new double*[environment.numSamples];
		for(int i = 0; i < environment.numSamples; i++)
			safe_stateDouble[i] = new double[environment.numNodes];
	}

	// Create a back up of the data, for later randomization
	safe_state = new int*[environment.numSamples];
	for(int i=0; i<environment.numSamples; i++)
		safe_state[i] = new int[environment.numNodes];

	if(environment.atLeastOneContinuous == 1){
		safe_stateIdx = new int*[environment.numNodes];
		for(int i = 0; i < environment.numNodes; i++)
			safe_stateIdx[i] = new int[environment.numSamples];
	}

	//copy to safe state
	float p;
	double* safe_weights;
	if(environment.sampleWeightsVec[0]==-1){
		if(environment.effN != environment.numSamples){
			p=environment.effN*1.0/environment.numSamples;
			safe_weights=new double[environment.numSamples];
		}
	}

	for(int i=0; i<environment.numSamples; i++)
	{
		for(int j=0; j<environment.numNodes;j++){
			safe_state[i][j]=environment.dataNumeric[i][j];
			if(environment.columnAsContinuous[j] == 1)
				safe_stateIdx[j][i]=environment.dataNumericIdx[j][i];
			if(environment.atLeastTwoGaussian)
				safe_stateDouble[i][j]=environment.dataDouble[i][j];
		}
		if(environment.sampleWeightsVec[0]==-1){
			//re-init weigths
			if(environment.effN != environment.numSamples){
				safe_weights[i]=environment.sampleWeights[i];
				environment.sampleWeights[i]=p;
			}
		}
	}

	int* nodes_toShf = new int[environment.numNodes];
  	for(int i=0; i<environment.numNodes; i++)
  		nodes_toShf[i]=0;
  	
  	//indexes of nodes to shuffle
  	for(int i=0; i<environment.numNoMore; i++)
	{
		nodes_toShf[inferredEdges_tab[i][0]] = 1;
	}

	double* confVect = new double[environment.numNoMore];
	for(int nb=0;nb<environment.numNoMore;nb++)
		confVect[nb] = 0;


	int* ptrVarIdx = new int[2];

	int prg_progress = -1;

	// loop on the number of shuffling
	for(int nb=1; nb<=environment.numberShuffles; nb++)
	{
		//prg_progress = printProgress(nb/(1.0*environment.numberShuffles), startTime, outdir, prg_progress);

		//Shuffle the dataset only for the variables present in nodes_toShf
		for(int col=0; col<environment.numNodes; col++)
		{
			if( nodes_toShf[col] == 1)
			{
				int row2=0;
				shuffle_lookup(lookup, lookup2, environment.numSamples);
				if(environment.columnAsContinuous[col] != 0){
					for(int i = 0; i < environment.numSamples; i++){
						lookup2[i] = i;
					}

					sort2arraysConfidence(environment.numSamples, lookup, lookup2);
				}
				for(int row = 0; row < environment.numSamples; row++){
					environment.dataNumeric[row][col] = safe_state[lookup[row]][col];

					if(environment.columnAsGaussian[col] == 1 && environment.atLeastTwoGaussian == 1){
						environment.dataDouble[row][col] = safe_stateDouble[lookup[row]][col];
					}
				}

				for(int row = 0; row < environment.numSamples; row++){
					if(environment.columnAsContinuous[col] != 0){

						// cout << environment.dataNumericIdx[col][row] << endl;
						// environment.dataNumericIdx[col][row] = lookup2[environment.dataNumericIdx[col][row]];
						if( environment.dataNumeric[lookup2[row]][col] != -1 ){
							// cout << lookup2[safe_stateIdx[col][row]] << " ";
							environment.dataNumericIdx[col][row2] = lookup2[row];
							
							row2++;
						}
					}					
				}
				if(environment.columnAsContinuous[col] != 0){
					// cout << "\n" << col << " " << row2  << flush;

					while(row2<environment.numSamples){

						environment.dataNumericIdx[col][row2] = -1;
						row2++;
					}
				}
			}
			// if( nodes_toShf[col] == 1)
			// {
			// 	int row2=0;
			// 	shuffle_lookup(lookup, lookup2, environment.numSamples);
			// 	if(environment.columnAsContinuous[col] != 0){
			// 		for(int i = 0; i < environment.numSamples; i++){
			// 			lookup2[i] = i;
			// 		}

			// 		sort2arraysConfidence(environment.numSamples, lookup, lookup2);
			// 	}
			// 	for(int row = 0; row < environment.numSamples; row++){
			// 		environment.dataNumeric[row][col] = safe_state[lookup[row]][col];

			// 		if(environment.columnAsGaussian[col] == 1){
			// 			environment.dataDouble[row][col] = safe_stateDouble[lookup[row]][col];
			// 		}

			// 		if(environment.columnAsContinuous[col] != 0){

			// 			// cout << environment.dataNumericIdx[col][row] << endl;
			// 			// environment.dataNumericIdx[col][row] = lookup2[environment.dataNumericIdx[col][row]];
			// 			if( lookup2[safe_stateIdx[col][row]]!=-1 ){
			// 				// cout << lookup2[safe_stateIdx[col][row]] << " ";
			// 				environment.dataNumericIdx[col][row2] = lookup2[safe_stateIdx[col][row]];
							
			// 				row2++;
			// 			}
			// 		}					
			// 	}

			// 	if(environment.columnAsContinuous[col] != 0){
			// 		while(row2<environment.numSamples){
			// 			environment.dataNumericIdx[col][row2] = -1;
			// 			row2++;
			// 		}
			// 	}
			// }
		}

		for(int i = 0; i < environment.numSamples;i++){
			for(int j = 0; j < environment.numNodes;j++){
				environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
			}
		}

		if(environment.atLeastTwoGaussian == 1){
			for (int i = 0; i < environment.numNodes; i++) {
				environment.means[i] = 0.0;
				environment.standardDeviations[i] = 0.0;
			}
			//compute the means and standard deviations
			computeMeansandStandardDeviations(environment);

			
			//compute the correlations coefficients
			computeCorrelations(environment);
		}

		int X,Y;
		int N_xy_ui;
		double NIxy_ui,k_xy_ui;
		// evaluate the mutual information for every edge
		for(int i=0; i<environment.numNoMore; i++)
		{
		 	X = inferredEdges_tab[i][0] ;
	 		Y = inferredEdges_tab[i][1] ;

			double* res; 

			// cout << environment.nodes[X].name << "( " << X << " " <<  nodes_toShf[X] << ")" << environment.nodes[Y].name << "( " << Y << " " <<  
					// nodes_toShf[Y] << ")" << endl;
			//discrete case
			if(environment.columnAsContinuous[X] == 0 && environment.columnAsContinuous[Y] == 0){
				res = computeEnsInformationNew(environment, NULL, 0, NULL, 0, -1, X, Y, environment.cplx, environment.m);
				N_xy_ui = res[0]; 
				NIxy_ui = res[1];
				k_xy_ui = res[2];
				free(res);

				// cout << "DISC: " << NIxy_ui << " " << k_xy_ui << endl;

			}  else if(environment.columnAsGaussian[X] == 1 && environment.columnAsGaussian[Y] == 1){
				int N = environment.numSamples;
				N_xy_ui = N; 
				double NIxy_ui1 = environment.rho[X][Y];
				k_xy_ui = log(N);

				NIxy_ui = (-log(1 - pow(NIxy_ui1,2))/2) * environment.numSamples;
			}

			// continuous case non gaussian, mixed
			else {
				// cout << "X: " << X  << " " << environment.columnAsContinuous[X] << "		Y: " << Y << " " << environment.columnAsContinuous[Y]<< endl;
				res = computeEnsInformationContinuous(environment, NULL, 0, NULL, 0, -1, X, Y, environment.cplx, environment.m);
				N_xy_ui = res[0]; 
				NIxy_ui = res[1];
				k_xy_ui = res[2];
				free(res);

				// cout << "GUIDO: " << N_xy_ui << " " << NIxy_ui << " " << k_xy_ui << endl;
				

			} 

			double ni = NIxy_ui-k_xy_ui;
			if(ni <= 0){ ni = 0 ;}
			confVect[i] += exp(-ni);
	    }
	}

	//evaluate the average confidence
	for(int nb=0;nb<environment.numNoMore;nb++){
		confVect[nb] /= environment.numberShuffles;	
	}

	//put values > 1 to 1
	for(int nb=0;nb<environment.numNoMore;nb++)
		if(confVect[nb] > 1) 
			confVect[nb] = 1;	

	//remove edges based on confidence cut
	double confidence;
	vector<int> toDelete ;

	for(int i = 0; i < environment.numNoMore; i++){
		int X = inferredEdges_tab[i][0] ;
 		int Y = inferredEdges_tab[i][1] ;
 		confidence = exp (- (environment.edges[X][Y].edgeStructure->Ixy_ui - environment.edges[X][Y].edgeStructure->cplx));
 		confVect[i] = confidence / confVect[i];
 		if(confVect[i] > environment.confidenceThreshold){
 			if(X > Y)
 				environment.edges[X][Y].edgeStructure->status = 1; 
 			else
 				environment.edges[Y][X].edgeStructure->status = 1; 
			environment.edges[X][Y].isConnected = 0;
			environment.edges[Y][X].isConnected = 0;
			// environment.edges[Y][X].edgeStructure->Ixy_ui = 0;

			// cout << environment.nodes[X].name << "\t" <<environment.nodes[Y].name << "\t" << confidence << "\t" << 
			// exp (- (environment.edges[X][Y].edgeStructure->Ixy_ui - environment.edges[X][Y].edgeStructure->cplx)) << "\t" <<
			// confVect[i] << endl;

			toDelete.push_back(i);
 		}
	}

	cout << "# -- number of edges cut: " << toDelete.size() << "\n";

	//delete from vector
	environment.noMoreAddress.clear();
	for(int i = 0; i < environment.numNoMore; i++){
		if(!(std::find(toDelete.begin(), toDelete.end(), i) != toDelete.end())) {
			int X = inferredEdges_tab[i][0] ;
 			int Y = inferredEdges_tab[i][1] ;
			XJAddress* s = new XJAddress();
			s->i=X;
			s->j=Y;
			environment.noMoreAddress.push_back(s);
		}
	}



		
	for(int X = 0; X < environment.numNodes -1; X++){
		for(int Y = X+1; Y < environment.numNodes; Y++){
			if(environment.edges[X][Y].isConnected == -2 || environment.edges[X][Y].isConnected == 2 || environment.edges[X][Y].isConnected == 6){
	 			environment.edges[X][Y].isConnected = 1;
	 			environment.edges[Y][X].isConnected = 1;
	 		}
	 	}
	 }



	//////////////////////////////////////////////// copy data back
	for(int i=0; i<environment.numSamples; i++)
	{
		if(environment.sampleWeightsVec[0]==-1){
			if(environment.effN != environment.numSamples)
				environment.sampleWeights[i]=safe_weights[i];
		}


		for(int j=0; j<environment.numNodes;j++){
			environment.dataNumeric[i][j]=safe_state[i][j];
			if(environment.atLeastTwoGaussian)
				environment.dataDouble[i][j]=safe_stateDouble[i][j];
		}
	}

	for(int i = 0; i < environment.numSamples;i++){
		for(int j = 0; j < environment.numNodes;j++){
			environment.oneLineMatrix[j * environment.numSamples + i] = environment.dataNumeric[i][j];
		}
	}

	if(environment.atLeastOneContinuous){
		//create the data matrix for factors indexes
	
		for(int i = 0; i < environment.numNodes; i++){
			for(int j = 0; j < environment.numSamples; j++) 
				environment.dataNumericIdx[i][j]=-1;
		}


		for(int j = 0; j < environment.numNodes; j++){
			if(environment.columnAsContinuous[j] != 0){
				transformToFactorsContinuous(environment, j);
				transformToFactorsContinuousIdx(environment, j);
				transformToFactors(environment, j);
			}
		}
	} else if(environment.atLeastTwoGaussian){
		for (int i = 0; i < environment.numNodes; i++) {
			environment.means[i] = 0.0;
			environment.standardDeviations[i] = 0.0;
		}
		//compute the means and standard deviations
		computeMeansandStandardDeviations(environment);

		
		//compute the correlations coefficients
		computeCorrelations(environment);
	}
	 //////////////////////////////////////////////// end copy data back

	//saveConfidenceVector(environment, inferredEdges_tab, confVect, environment.outDir, slash);


	std::sort(environment.noMoreAddress.begin(), environment.noMoreAddress.end(), sorterNoMore2(environment));
	environment.numNoMore = environment.noMoreAddress.size();

	delete[] ptrVarIdx;

	for(int i=0; i<environment.numSamples; i++)
		delete safe_state[i];
	delete[] safe_state;

	if(environment.atLeastTwoGaussian){
		for(int i=0; i<environment.numSamples; i++)
			delete safe_stateDouble[i];
		delete[] safe_stateDouble;
	}
	delete[] lookup;
	delete[] lookup2;

	if(environment.sampleWeightsVec[0]==-1){
		if(environment.effN != environment.numSamples){
			delete[] safe_weights;
		}
	}
	delete[] nodes_toShf;


	vector< vector <string> > confVect1;
	vector <string> v;
	v.push_back("x");
	v.push_back("y");
	v.push_back("confidence_ratio");
	confVect1.push_back(v);
	for(int i =0; i < noMore; i++){
		vector <string> v;
		v.push_back(environment.nodes[inferredEdges_tab[i][0]].name);
		v.push_back(environment.nodes[inferredEdges_tab[i][1]].name);
		v.push_back(to_string(confVect[i]));
		confVect1.push_back(v);
	}

	for(int i=0; i < noMore; i++)
		delete inferredEdges_tab[i];
	delete[] inferredEdges_tab;
	delete[] confVect; 
	
	return confVect1;
}