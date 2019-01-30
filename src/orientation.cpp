#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include <string.h>

#include "structure.h"
#include "utilities.h"

#include "functionsOrientation.h"

using namespace std;

bool sortFunctionListStruct(const Struct* a, const Struct* b) {
     return a->rv > b->rv;
}

void orientation(Environment environment, string slash, bool isVerbose){
	cout << "# -> START miic orientation...\n";
	std::stringstream ss;
	ss.str("");
	ss << environment.outDir << slash << "log.miic.orient.txt";

	char *addr = new char[ss.str().length()+1];
	strcpy(addr,ss.str().c_str());
	Log* pLog = new Log(addr);

	delete[] addr;

	// print the environment
	printEnvironment(environment, pLog);

	if(environment.isVerbose){ cout << "\n# ---- Get all the structures (v-structures or non-v-structures) ----\n" ; }
	
	// ALGO START TIME
	double startTime_algo = get_wall_time();
	double startTime = get_wall_time();
	getAllStructures(environment, isVerbose);

	//// ORDER THEM
	//// Order the RV_xyz index of the vstruct by decreasing value
    std::sort(environment.globalListOfStruct.begin(), environment.globalListOfStruct.end(), sortFunctionListStruct);

    int* orig = new int[environment.globalListOfStruct.size()];
    for(int i = 0; i < environment.globalListOfStruct.size(); i++){
    	// cout << i << " " << environment.globalListOfStruct[i]->originalPositionStruct << "\n";
    	orig[environment.globalListOfStruct[i]->originalPositionStruct] = i;
    }

	long double spentTime = (get_wall_time() - startTime);
	if(environment.isVerbose){ cout << "\n# ----> Get all the structures elapsed time:" << spentTime << "sec.\n\n" ; }

	////// PREPARE DATA FRAME for PROPAGATION  //////
	// ------------
	prepareDFforPropagation(environment, slash);
	// ------------

	////// PROPAGATE ORIENTATIONS #######
	if(isVerbose){ cout << "\n# ---- Propagate Orientations ----\n" ; }
	startTime = get_wall_time();

	vector<StructWithOrtToPropagate*> structWithOrtToPropagate;

    // cout << "STRUCTURES #: " << environment.iCountStruct << "\n";

	int iCurrTpl = 0;
	while( ( iCurrTpl < environment.iCountStruct ) || ( !structWithOrtToPropagate.empty() ) ){

		if(isVerbose){
			cout << "\n# Current tpl " << iCurrTpl +1 << ": " << environment.nodes[environment.globalListOfStruct[iCurrTpl]->xi].name << "." << 
			environment.globalListOfStruct[iCurrTpl]->d1 <<  "." <<
			environment.nodes[environment.globalListOfStruct[iCurrTpl]->xk].name << "." << environment.globalListOfStruct[iCurrTpl]->d2 << "." <<
			environment.nodes[environment.globalListOfStruct[iCurrTpl]->xj].name << " (SV:" << environment.globalListOfStruct[iCurrTpl]->sv << 
				", Out:" << environment.globalListOfStruct[iCurrTpl]->isOut << ")\n"; 
 		}

 		//// Is it a v-structure without a collider?
		//// ----------------------------------------
	    if(	( environment.globalListOfStruct[iCurrTpl]->sv < 0 ) && ( environment.globalListOfStruct[iCurrTpl]->d1 < 2 ) &&
	     ( environment.globalListOfStruct[iCurrTpl]->d2 < 2 ) ){
	    	if( environment.globalListOfStruct[iCurrTpl]->isOut == 3 )
       		{
				if(isVerbose) { cout << "# --> is a contradicted v-structure!\n" ; }

			} else {
		
				if(isVerbose) { cout << "# --> is a v-structure without convergent orientations => keep b1-v for propagation\n" ; }

    	    	// Keep the index of the current structure to globally propagate the direction of b1-->v
    	    	StructWithOrtToPropagate* s = new StructWithOrtToPropagate();
    	    	s->name = iCurrTpl;
    	    	s->idxStruct = iCurrTpl;
    	    	s->vStruct = iCurrTpl;
    	    	s->idxEdge = 1;
    	    	s->dir = 2;
				structWithOrtToPropagate.push_back(s);
			}
	    }

	    // ----------------------------------------
    
    	// Is there any possible 'global' propagation?
    	// ----------------------------------------
	    // While there exists a upper structure from where a new orientation can be propagated

		StructWithOrtToPropagate* currStructToPpg;
	    while(hasLowerIndex(structWithOrtToPropagate, iCurrTpl))
    	{
			//if(isVerbose) 
			//{ cout << "\n# ----> ( currTpl:" << iCurrTpl +1 << ") all pending tpl idx :" << vectorToString(structWithOrtToPropagateIdx) << "\n" ; }
		
			//// Order the index by increasing order and take the first struct
			int min = 0;
			for(int i = 0; i < structWithOrtToPropagate.size(); i++)
				if(structWithOrtToPropagate[i]->name < structWithOrtToPropagate[min]->name)
					min = i;

			//// Read the first structure having an orientation that can be propagated, and remove it from the list
			currStructToPpg = structWithOrtToPropagate[min];
			structWithOrtToPropagate.erase(structWithOrtToPropagate.begin() + min);

			if(isVerbose) 
			{
				int type = currStructToPpg->idxEdge;
				string ort;
				if(type == 1)
					ort = "d1";
				else
					ort = "d2"; 

				cout << "\n# --> Ort. " << ort << " to propagate from tpl " <<
					currStructToPpg->idxStruct << ": " << environment.nodes[environment.globalListOfStruct[currStructToPpg->idxStruct]->xi].name << "." <<
					environment.globalListOfStruct[currStructToPpg->idxStruct]->d1 << "." << environment.nodes[environment.globalListOfStruct[currStructToPpg->idxStruct]->xk].name << "." <<
					environment.globalListOfStruct[currStructToPpg->idxStruct]->d2 << "." << environment.nodes[environment.globalListOfStruct[currStructToPpg->idxStruct]->xj].name <<
					" (SV:" << environment.globalListOfStruct[currStructToPpg->idxStruct]->sv << ")\n"; 
			}

			//// Check if this structure bringing a new orientation has not been contradicted
			if(environment.globalListOfStruct[currStructToPpg->idxStruct]->isOut == 3 ) {
				if(isVerbose) { cout << "\n# ----> Previously contradicted tpl!\n" ; }

				//// If yes, go to the next element of the loop
				//// But before, reinit the vector that contains all the 'pending' structure
				continue;
			}

			//// Get the newly oriented edge + direction from the structure
			int newOrtEdgeX;
        	int newOrtEdgeY;
        	int type;
			if(currStructToPpg->idxEdge == 1){
		        newOrtEdgeX = environment.globalListOfStruct[currStructToPpg->idxStruct]->xi;
	        	newOrtEdgeY = environment.globalListOfStruct[currStructToPpg->idxStruct]->xk;
	        } else if(currStructToPpg->idxEdge == 2){
	        	newOrtEdgeX = environment.globalListOfStruct[currStructToPpg->idxStruct]->xj;
	        	newOrtEdgeY = environment.globalListOfStruct[currStructToPpg->idxStruct]->xk;
	        }

	        // cout << "newOrtEdgeX " << environment.nodes[newOrtEdgeX].name << " newOrtEdgeY" << environment.nodes[newOrtEdgeY].name << endl;

	        int newOrtDir = currStructToPpg->dir;


			//// Process the open tpl list
			// cout << "iCurrTpl: " << iCurrTpl << endl;
			// cout << "X : " << environment.nodes[newOrtEdgeX].name << "\nY: " << environment.nodes[newOrtEdgeY].name << endl;


			vector<int> orgInvTpl = environment.edges[newOrtEdgeX][newOrtEdgeY].edgeStructure->edgesInSpeTpl_list;

			// sort
			vector<int> invTpl;
			for(int i = 0; i < orgInvTpl.size(); i++){
				invTpl.push_back(orig[orgInvTpl[i]]);
			}	

			//// Finally, sort the index of the structures concerned by this new orientation
			sort(invTpl.begin(), invTpl.end());

			for(int pos = 0; pos < invTpl.size(); pos++) {
				// cout << "POS: " << pos << endl;
        		int iInvTpl = invTpl[pos];
    			if(isVerbose) {
    				cout << "\n# --> iInvTpl " << iInvTpl << ": " <<
    					environment.nodes[environment.globalListOfStruct[iInvTpl]->xi].name << "." << environment.globalListOfStruct[iInvTpl]->d1 << "." <<
    					environment.nodes[environment.globalListOfStruct[iInvTpl]->xk].name << "." << environment.globalListOfStruct[iInvTpl]->d2 << "." <<
    					environment.nodes[environment.globalListOfStruct[iInvTpl]->xj].name << " (SV:" << environment.globalListOfStruct[iInvTpl]->sv << ", Out: " << 
    								environment.globalListOfStruct[iInvTpl]->isOut << ")\n";
				}	

				//// Make sure first that this structure has not been contradicted
            	if( environment.globalListOfStruct[iInvTpl]->isOut == 3 ){ continue; }

            	//// Find the edge in the iInvTpl, that should receive the new orientation
            	int edgeToBeSetX;
            	int edgeToBeSetY;

            	if(environment.globalListOfStruct[iInvTpl]->xi == newOrtEdgeX){
					edgeToBeSetX = newOrtEdgeX;
					edgeToBeSetY = newOrtEdgeY;
				} else if(environment.globalListOfStruct[iInvTpl]->xk == newOrtEdgeX){
					edgeToBeSetX = newOrtEdgeY;
					edgeToBeSetY = newOrtEdgeX;
				} else if(environment.globalListOfStruct[iInvTpl]->xj == newOrtEdgeX){
					edgeToBeSetX = newOrtEdgeX;
					edgeToBeSetY = newOrtEdgeY;
				}

				type = 1;
				if(environment.globalListOfStruct[iInvTpl]->xj == edgeToBeSetX || environment.globalListOfStruct[iInvTpl]->xj == edgeToBeSetY)
					type=2;

            	int ortToBeSet;

            	// cout << "edgeToBeSetX: " << environment.nodes[edgeToBeSetX].name << " edgeToBeSetY" << environment.nodes[edgeToBeSetY].name << endl;

            	if(newOrtEdgeY == edgeToBeSetY)
            		ortToBeSet = newOrtDir;
            	else
            		ortToBeSet = -1 * newOrtDir;

            	//// Display what should be done
				if(isVerbose) {
					cout << "# ----> ?Orient " << environment.nodes[edgeToBeSetX].name << "-" <<
					environment.nodes[edgeToBeSetY].name << " as " << ortToBeSet << "\n";
				}

				// If the invTpl is a non v-structure
				// cout << " environment.globalListOfStruct[iInvTpl]->sv " <<  environment.globalListOfStruct[iInvTpl]->sv << endl;
				if( environment.globalListOfStruct[iInvTpl]->sv > 0 ){
					// If the other edge is convergent, ie in (2,4,6)
					bool test = false;
					if(type == 1){
						if(environment.globalListOfStruct[iInvTpl]->d2 >= 2 )
							test=true;
						}
					else if(type == 2){
						if(environment.globalListOfStruct[iInvTpl]->d1 >= 2 )
							test=true;
					}
					// cout << "TEST COLL: " << test << " ortToBeSet: " << ortToBeSet << endl;

					if(test){
						// If the new orientation is convergent, ie in (2,4)
			        	//if( ortToBeSet %in% c(2,4) )
			        	if( ortToBeSet >= 2 ){
	                        // Display message
	           				if(isVerbose) { cout << "# -- COLL > 0 --! Trying to set a collider in a non vStruct!\n"; }

	           				//// Contradict the iInvTpl or contradict the vStruct and restart all
	    	                bool reinit_ret = simpleContradictionORrestartAll(environment, currStructToPpg, iInvTpl, isVerbose);

			                if(reinit_ret) 
			                { 
			                	if(isVerbose)
			                		cout << "#### RESTART!\n";
			                	//// Empty the list of orientation to propagate, set the current tpl at 0 and abort the loop on invTpl
			                	structWithOrtToPropagate.clear(); 
			                	iCurrTpl = 0; 
			                	break;
			                }
				        } else {
				        	// The new orientation is divergent, ie in -2 or -4
				        	// Set the new orientation
			            	setNewOrt(environment, iInvTpl, currStructToPpg, type, ortToBeSet, isVerbose);
				        }
					} else { 
						// The other edge is non oriented or divergent, ie in (1,-2,-4)
				    
				        // Set the new orientation
			            bool retValue = setNewOrt(environment, iInvTpl, currStructToPpg, type, ortToBeSet, isVerbose );
				        
				        // If the orientation has been set, check for local propagation
				        if(retValue){
				        	learnNewOrtLocally(environment, iInvTpl, currStructToPpg, type, ortToBeSet, structWithOrtToPropagate, isVerbose);
				        }
				    }
				}  else {
					// The invTpl is a v-structure
					// If the new orientation is uniquely divergent
			    	// if( ortToBeSet == (-4) )

			    	if( ( environment.isLatent && ortToBeSet == (-4) ) || ( !environment.isLatent && ortToBeSet == (-2) ) ){
	    			    // Display message
	       				if(isVerbose) { cout << "# -- MED < 0 --! Trying to set a mediation in a vStruct!\n"; }

	       				//// Contradict the iInvTpl or contradict the vStruct and restart all
	  	                bool reinit_ret = simpleContradictionORrestartAll( environment, currStructToPpg, iInvTpl, isVerbose);

		                if(reinit_ret) 
		                { 
		                	cout << "#### RESTART!\n";
		                	//// Empty the list of orientation to propagate, set the current tpl at 0 and abort the loop on invTpl
		                	structWithOrtToPropagate.clear(); 
		                	iCurrTpl = 0; 
		                	break;
		                }

				    } else {
					    // The new orientation can be (2,-2,4) if latent, or (2) if no latent
				    
	    			    // Set the new orientation
			            bool retValue = setNewOrt( environment, iInvTpl, currStructToPpg, type, ortToBeSet, isVerbose );

				        // If the orientation has been set, check for local propagation
				        if(retValue)
				        {
				            learnNewOrtLocally(environment, iInvTpl, currStructToPpg, type, ortToBeSet, structWithOrtToPropagate, isVerbose );
				        }
				    }
				} // END The invTpl is a v-structure
			} // END for each of the structure in which the edge is involved (iInvTpl)

		} // end while( length( structWithOrtToPropagateIdx[ structWithOrtToPropagateIdx <= iCurrTpl ] ) > 0 )

		iCurrTpl++;
	} // end while iCurrTpl
	spentTime = get_wall_time() - startTime;
	if(isVerbose){ cout << "\n# ----> Propagate Orientations elapsed time:" << spentTime << "sec.\n\n" ; }

	// ALGO STOP TIME
	long double spentTime_algo = (get_wall_time() - startTime_algo); 

	// Save the completed table of orientations
	ss.str("");
	//// Save the execTime
	ss << environment.outDir << slash << "struct_ort_tplOpen_ppgClosed.txt";
	saveTableOfOrientations(environment, ss.str());
	//write.table( gV$allOrderedStruct, file = "struct_ort_tplOpen_ppgClosed.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE )

	// ------

	////// UPDATE ADJ MATRIX #######
	// ------
	if(isVerbose){ cout << "\n# ---- Update Adjacency Matrix----\n"; }
	startTime = get_wall_time();

	for(int i = 0; i < environment.globalListOfStruct.size();i++ ){
		if(environment.globalListOfStruct[i]->isOut == 3)
			continue;
		int myOrt;
		int myOrt_opposite;

		myOrt = environment.globalListOfStruct[i]->d1;

		if(myOrt == 6 || myOrt == 1)
			myOrt_opposite = myOrt;
		else
		 	myOrt_opposite = -1*myOrt;

		environment.edges[environment.globalListOfStruct[i]->xi][environment.globalListOfStruct[i]->xk].isConnected = myOrt;
		environment.edges[environment.globalListOfStruct[i]->xk][environment.globalListOfStruct[i]->xi].isConnected = myOrt_opposite;

		myOrt = environment.globalListOfStruct[i]->d2;
		if(myOrt == 6 || myOrt == 1)
			myOrt_opposite = myOrt;
		else 
			myOrt_opposite = -1*myOrt;

		environment.edges[environment.globalListOfStruct[i]->xj][environment.globalListOfStruct[i]->xk].isConnected = myOrt;
		environment.edges[environment.globalListOfStruct[i]->xk][environment.globalListOfStruct[i]->xj].isConnected = myOrt_opposite;

		//// If the tpl is closed, copy the third orientation
	    if( environment.globalListOfStruct[i]->sv == 0 )
	    {
	        myOrt = environment.globalListOfStruct[i]->xj;

	        environment.edges[environment.globalListOfStruct[i]->xi][environment.globalListOfStruct[i]->xj].isConnected = myOrt;

	        if(myOrt == 1)
				environment.edges[environment.globalListOfStruct[i]->xj][environment.globalListOfStruct[i]->xi].isConnected = 1;
			else
				environment.edges[environment.globalListOfStruct[i]->xj][environment.globalListOfStruct[i]->xi].isConnected = -1*myOrt;
	    }
	}

	//// Save the adjacency matrix
	ss.str("");
	    
	//// Make an adjacency matrix from the list
	ss << environment.outDir << slash << "adjacencyMatrix.miic.orient.txt";
	saveAdjMatrix(environment, ss.str());

	spentTime = (get_wall_time() - startTime) / (double)(CLOCKS_PER_SEC / 1000) /1000;
	if(isVerbose)
		cout << "\n# ----> Update Adjacency Matrix elapsed time:" << spentTime << "sec\n\n";


	delete [] orig;
}
