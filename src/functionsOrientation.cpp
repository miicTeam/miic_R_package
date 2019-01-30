#include <iostream>
#include <vector>
#include <cmath>
#include "utilities.h"
#include "computeEnsInformation.h"
#include <sys/stat.h>


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include <string.h>

using namespace std;


void saveTableOfOrientations(Environment environment, string filename){
	ofstream output;
	output.open(filename.c_str());

	output << "b1" << "\t" << "d1" << "\t" << "v" << "\t" << "d2" << "\t" << "b2" << "\t" 
				<< "RV" << "\t" << "SV" << "\t" << "isOut" << "\t" << "orgOut" << "\t" << 
				"orgD1" << "\t" << "orgD2" << "\t" << "Ibase" << "\t" << "Istruct\n";

	for(int i = 0; i < environment.globalListOfStruct.size(); i++){
		output << i + 1 << "\t" << 
		environment.nodes[environment.globalListOfStruct[i]->xi].name << "\t" <<
		environment.globalListOfStruct[i]->d1 << "\t" << 
		environment.nodes[environment.globalListOfStruct[i]->xk].name << "\t" << 
		environment.globalListOfStruct[i]->d2 << "\t" << 
		environment.nodes[environment.globalListOfStruct[i]->xj].name << "\t" << 
		environment.globalListOfStruct[i]->rv << "\t" <<
		environment.globalListOfStruct[i]->sv << "\t" <<
		environment.globalListOfStruct[i]->isOut << "\t" <<
		environment.globalListOfStruct[i]->orgOut << "\t" <<
		environment.globalListOfStruct[i]->orgD1 << "\t" <<
		environment.globalListOfStruct[i]->orgD2 << "\t" <<
		environment.globalListOfStruct[i]->Ibase << "\t" <<
		environment.globalListOfStruct[i]->Istruct << "\n";
	}

	output.close();
}

int sign(double val){
	if(val < 0) 
		return -1;
	else if(val > 0) 
		return 1;
	else 
		return 0;
}

bool getSingleStructure(Environment& environment, const int posX, const int posY, const int posZ, bool isVerbose){
	//// Value to return
	bool structFound = false;	
	//// To check if xk belongs to the {ui} of the base
	bool isXkInUi = false;
	//// Get the {ui} of the base edge

	vector<int> u(environment.edges[posX][posZ].edgeStructure->ui_vect_idx);
	if(environment.edges[posX][posZ].edgeStructure->ui_vect_idx.size() > 0){
		//// Remove xk from the {ui} if needed
		for(int i = 0; i < environment.edges[posX][posZ].edgeStructure->ui_vect_idx.size(); i++){
			if(u[i] == posY){
				u.erase(u.begin() + i);
				isXkInUi = true;
				break;
			}
		}
	}

	int* ui;
	int* zi = NULL;

	if(u.empty())
		ui = NULL;
	else
		ui = &u[0];

	vector<int> z;
	z.clear();
	z.push_back(posY);

	zi = &z[0];

	//// Compute the conditional mutual interaction of the tpl
	//// ----
	//// Compute I(xi xj|ui)[xi,xj,xk,ui]


	double* res = computeEnsInformationNew(environment, ui, u.size(), zi, z.size(), -1, posX, posZ, environment.cplx, environment.m);

	double Is = res[7];
	double Cs = res[8];

	// double* Ixixj_ui_xk_list = computeEnsInformation(environment, ui, u.size(), zi, z.size(), -1, posX, posZ, environment.cplx);

	// //// Compute I(xi xj|ui,xk)[xi,xj,xk,ui]
	// double* Ixixj_uixk_list = computeEnsInformation(environment, ui, u.size(), zi, z.size(), u.size()+2, posX, posZ, environment.cplx);

	//// Is = N.I(xi; xj; xk) = N.I(xi; xj|{ui})[xk] - N.I(xi; xj|{xk, ui})
	// int nbrRetValues = 3 * environment.edges[posX][posZ].edgeStructure->zi_vect_idx.size();

	////  Compute I(xyz|ui)[xyuiz] = I(xy|ui)[xyuiz] - I(xz|ui,z)[xyuiz]

	if(environment.isK23){

		if(environment.isDegeneracy){
				Cs += log(3);
		}

		Is = Is + Cs;
	} 

	double myRV_tmp = Is;

	//// Check for the sign of Is
	//// If Is > 0: we have a non v-structure
	//// If Is < 0 AND xk in not in {ui}: we have a v-structure
	if(myRV_tmp < 0 && isXkInUi){ cout << "# !!!!!! xk (" << posY << "belongs to {ui}\n" ; }

	if(myRV_tmp > 0 || (myRV_tmp < 0 && !isXkInUi)){
		structFound = true;

		// Get the sign and the absolute scoreCompute the score
		int mySV_tmp = sign(myRV_tmp);
		myRV_tmp = abs( myRV_tmp );



		Struct* s = new Struct();
		s->xi = posX;
		s->xk = posY;
		s->xj = posZ;
		s->rv = myRV_tmp;
		s->sv = mySV_tmp;
		// cout <<"SV: " << s->sv << " myRV_tmp :" << myRV_tmp << " sign(myRV_tmp):" << sign(myRV_tmp) << endl;
		s->Ibase = res[1];
		s->Istruct = Is;
		s->originalPositionStruct = environment.globalListOfStruct.size();

		s->d1 = 1;
		s->d2 = 1;

		environment.globalListOfStruct.push_back(s);

		//// Add also the tpl to the list of edges involved in structure
		environment.edges[posX][posY].edgeStructure->indexStruct.push_back(environment.iCountStruct);
		environment.edges[posY][posZ].edgeStructure->indexStruct.push_back(environment.iCountStruct);

		environment.edges[posX][posY].edgeStructure->edgesInSpeTpl_list.push_back(environment.iCountStruct);
		environment.edges[posY][posZ].edgeStructure->edgesInSpeTpl_list.push_back(environment.iCountStruct);

		if(isVerbose)
		{ 
			string mySign;
			if(mySV_tmp > 0)
				mySign = "positive";
			else
				mySign = "negative";

		   	cout << "\n# --[" <<
		   		environment.nodes[posX].name << ", " << environment.nodes[posY].name << ", " << environment.nodes[posZ].name <<
	   			"]-- is a " << mySign << " structure \n# ( RV = " << myRV_tmp << "; score_xyz = " << Is << ")\n";
		}  

		//// Add the structure to the global list
		environment.iCountStruct++;

	}
	delete res;
	return structFound;
}


void getAllStructures(Environment& environment, bool isVerbose){
	
	environment.iCountStruct = 0;

	//// For all non phantom "xy"
	for(int pos = 0; pos < environment.numNoMore; pos++){
		//// Get the current non phantom edge
		int posX = environment.noMoreAddress[pos]->i;
		int posY = environment.noMoreAddress[pos]->j;

		if(isVerbose){ cout << "\n# (" << pos << ") ---- Edge " << environment.nodes[posX].name << "--" << environment.nodes[posY].name << "\n"; }

		//// Prepare a list that will contain the neighbors of "x" and the neighbors of "y"
		vector<int> neighboursX;
		vector<int> neighboursY;

		// vector<int> neighboursCommon;
		//// Get other links (.z), with:
		//// - x_z or y_z (or both) not phantom
		//// - x_y < x_z and x_y < y_z (in the given order of edges)
		//// (below, we have the index from which we search starting at iNotPhEdge+1)
		//// ----
		//// Get the neighbours of x and y

		for(int i = pos + 1; i < environment.noMoreAddress.size(); i++){
			int posX1 = environment.noMoreAddress[i]->i;
			int posY1 = environment.noMoreAddress[i]->j;

			if(posY1 == posX && !environment.edges[posY][posX1].isConnected)
				neighboursX.push_back(posX1);
			else if(posY1 == posY && !environment.edges[posX][posX1].isConnected)
				neighboursY.push_back(posX1);
		}

		for(int i = pos + 1; i < environment.noMoreAddress.size(); i++){
			int posX1 = environment.noMoreAddress[i]->i;
			int posY1 = environment.noMoreAddress[i]->j;
			
			if(posX1 == posX && !environment.edges[posY][posY1].isConnected)
				neighboursX.push_back(posY1);
			else if(posX1 == posY && !environment.edges[posX][posY1].isConnected)
				neighboursY.push_back(posY1);
		}

		int sizeX = neighboursX.size();
		int sizeY = neighboursY.size();

		if(sizeX == 0 && sizeY == 0)
			continue;

		//// Is the tpl a v-/non v-structure
		//// ----
		// iterate on neighbours of x
		for(int i = 0; i < sizeX; i++){
			if(isVerbose){ cout << "\n# (" << pos << ") -----| Test tpl (" << environment.nodes[posY].name << ", " << environment.nodes[posX].name << ", " << environment.nodes[neighboursX[i]].name << ")\n" ; } 
			//// Get the structure if any
		   	getSingleStructure(environment, posY, posX, neighboursX[i], isVerbose);
		}

		// iterate on neighbours of y
		for(int i = 0; i < sizeY; i++){
			if(isVerbose){ cout << "\n# (" << pos << ") -----| Test tpl (" << environment.nodes[posX].name << ", " << environment.nodes[posY].name << ", " << environment.nodes[neighboursY[i]].name << ")\n" ; } 
			//// Get the structure if any
			getSingleStructure(environment, posX, posY, neighboursY[i], isVerbose);
		}

		// cout << "\n Size of struct: " << environment.edges[posX][posY].edgeStructure->edgesInSpeTpl_list.size() << endl;
	}
}


void prepareDFforPropagation(Environment& environment, string slash){
	stringstream ss;
	ss.str("");
	//// Save the execTime
	ss << environment.outDir << slash << "structnotOrt.txt";
	saveTableOfOrientations(environment, ss.str());
}

bool simpleContradictionORrestartAll(Environment& environment, StructWithOrtToPropagate* myCurrStruct, int myInvTpl, 
									 bool isVerbose)
{
	//// If myInvTpl is less likely than myCurrStruct that induced the contradiction
	//// then, myInvTpl is contradicted
	//// else
	//// all the orientations are deleted, the myCurrStruct is ignore and we restart the propagation		
	bool doBreak = false;

	if(isVerbose) {
		cout << "# ------! The new orientation comes from the vStruct" << myCurrStruct->vStruct << "\n";
	}

	if( myCurrStruct->vStruct < myInvTpl )
	{

		if(isVerbose) 
		{ cout << "# ------! Wrong tpl: " << myInvTpl << "\n";  }

		//// Set the involved tpl as contradicted
		environment.globalListOfStruct[myInvTpl]->isOut = 3;
		environment.globalListOfStruct[myInvTpl]->orgOut = myCurrStruct->vStruct;

	} else {
		if(isVerbose) 
		{ cout << "# ------! Wrong tpl is the vStruct: " << myCurrStruct->vStruct <<  " => restart without      \n";}

		//// Set the v-struct corresponding to the current tpl as contradicted
		environment.vstructWithContradiction.push_back(myCurrStruct->vStruct);


		//// Remove all the orientations
		//// Set the vStruct as contradicted in the reinitialized table
		for(int i = 0; i < environment.globalListOfStruct.size(); i++){
			environment.globalListOfStruct[i]->d1 = 1;
			environment.globalListOfStruct[i]->d2 = 1;
			environment.globalListOfStruct[i]->isOut = 0;
		}

		for(int i = 0; i < environment.vstructWithContradiction.size(); i++){
			environment.globalListOfStruct[environment.vstructWithContradiction[i]]->isOut = 3;
		}

		//// And finish here the loop on invTpl
		doBreak = true;



	}
	return doBreak;
}

bool setNewOrt(Environment& environment, int myInvTpl, StructWithOrtToPropagate*& myCurrStructToPpg, 
			   int type, int myOrtToBeSet, bool isVerbose)
{
	//
	//
	bool isNewOrt = false;
	int val;
	// If the edge to orient is not oriented
	if(type == 1){
		val = environment.globalListOfStruct[myInvTpl]->d1;
	} else if(type == 2){
		val = environment.globalListOfStruct[myInvTpl]->d2;
	}
	if(val == 1 ){
		// Set the new orientation
		if(type == 1){
			environment.globalListOfStruct[myInvTpl]->d1 = myOrtToBeSet;
			environment.globalListOfStruct[myInvTpl]->orgD1 = myCurrStructToPpg->vStruct;
		}
		else if(type == 2){
			environment.globalListOfStruct[myInvTpl]->d2 = myOrtToBeSet;
			environment.globalListOfStruct[myInvTpl]->orgD2 = myCurrStructToPpg->vStruct;
		}	
		isNewOrt = true;
	} else { 
		// If the edge to orient is already oriented
		if(environment.isLatent){
			// If the existing orientation is 2 (-2) and the new one is -2 (2), then it becomes 6!
			if( abs(val) == 2 && abs(myOrtToBeSet) == 2 && sign(val) != sign(myOrtToBeSet) ){
				// Set the new orientation (--> here means bi-directed)
				if(type == 1){
					environment.globalListOfStruct[myInvTpl]->d1 = 6;
					environment.globalListOfStruct[myInvTpl]->orgD1 = myCurrStructToPpg->vStruct;
				}
				else if(type == 2){
					environment.globalListOfStruct[myInvTpl]->d2 = 6;
					environment.globalListOfStruct[myInvTpl]->orgD2 = myCurrStructToPpg->vStruct;
				}
				isNewOrt = true;
			} else if( abs(val) == 2 & abs(myOrtToBeSet) == 4 && sign(val) == sign(myOrtToBeSet)){
				// Set the new orientation

				if(type == 1){
					environment.globalListOfStruct[myInvTpl]->d1 = myOrtToBeSet;
					environment.globalListOfStruct[myInvTpl]->orgD1 = myCurrStructToPpg->vStruct;
				}
				else if(type == 2){
					environment.globalListOfStruct[myInvTpl]->d2 = myOrtToBeSet;
					environment.globalListOfStruct[myInvTpl]->orgD2 = myCurrStructToPpg->vStruct;
				}

				isNewOrt = true;
			} else {
				if(isVerbose)
				{
					cout << "# --WAR--> Third case in 'setNewOrt()'\n" <<
						  "# -------> Existing Ort.(" << val << ")\n" <<
						  "# -------> New Ort.(" << myOrtToBeSet << ")\n" <<
						  "# --------\n";
				}
			}
		}
	}

	//cout << "d1: " << environment.globalListOfStruct[myInvTpl]->d1 << " d2: " << environment.globalListOfStruct[myInvTpl]->d2 << endl;

	return isNewOrt;
}

void learnNewOrtLocally(Environment& environment, int myInvTpl, StructWithOrtToPropagate*& myCurrStructToPpg, int type, 
						int myOrtToBeSet, vector<StructWithOrtToPropagate*>& structWithOrtToPropagate, bool isVerbose)
{
	// If the new orientation is convergent, there may be a local propagation
	// if( myOrtToBeSet %in% c(2,4) )
	// cout << "learnNewOrtLocally myOrt: " << myOrtToBeSet << endl;
	if( myOrtToBeSet >= 2 ){
		// By default, consider the case without allowing for latent variables
		int myLearnedDir = -1;
		if(environment.isLatent){
			myLearnedDir = environment.globalListOfStruct[myInvTpl]->sv > 0 ? -4 : 2;
		} else { 
			myLearnedDir = environment.globalListOfStruct[myInvTpl]->sv > 0 ? -2 : 2;
		}

		// Get the orientation of the edge targeted by the local propagation
		int myTargetEdgeDir;
		if(type == 1){
			myTargetEdgeDir = environment.globalListOfStruct[myInvTpl]->d2;
		} else if(type == 2){
			myTargetEdgeDir = environment.globalListOfStruct[myInvTpl]->d1;
		}

		// If the learned orientation is different from the existing one, keep for future propagation

		// cout << "myTargetEdgeDir: " << myTargetEdgeDir << " myLearnedDir: " << myLearnedDir << endl;
		if( myTargetEdgeDir != myLearnedDir ){
			//// Keep the index of the current structure to propagate the new direction later to the second edge
			// environment.structWithOrtToPropagate

			// cout << "SIZE: " << structWithOrtToPropagate.size() << endl;
			int idx = -1;

			for(int i = 0; i < structWithOrtToPropagate.size();i++){
				if(structWithOrtToPropagate[i]->name == myInvTpl){
					idx = i;
					break;
				}
			}

			if(idx == -1){
				idx = structWithOrtToPropagate.size();
				StructWithOrtToPropagate* s = new StructWithOrtToPropagate();
				structWithOrtToPropagate.push_back(s);
			}

			structWithOrtToPropagate[idx]->name = myInvTpl;
			structWithOrtToPropagate[idx]->idxStruct = myInvTpl;
			structWithOrtToPropagate[idx]->vStruct = myCurrStructToPpg->vStruct;
			structWithOrtToPropagate[idx]->dir =  myLearnedDir;
			if(type == 1)
				structWithOrtToPropagate[idx]->idxEdge = 2;
			else
				structWithOrtToPropagate[idx]->idxEdge = 1;
		
			if(isVerbose){
			cout << "# --------> Keep edge " << environment.nodes[environment.globalListOfStruct[myInvTpl]->xi].name << "-" << environment.nodes[environment.globalListOfStruct[myInvTpl]->xk].name <<
			" for future propagation (" << myLearnedDir << "), Out:" << environment.globalListOfStruct[myInvTpl]->isOut << ")\n"; 
			
				
			}
		}
	}
}


bool hasLowerIndex(vector<StructWithOrtToPropagate*> vec, int index){
	for(int i = 0; i < vec.size(); i++){
		if(vec[i]->name <= index)
			return true;
	}
	return false;
}
