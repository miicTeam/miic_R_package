#include <iostream>
#include <string>
#include <cmath>
#include "structure.h"
#include "computeEnsInformation.h"
#include "skeletonInitialization.h"
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;



// sort the noMore edges according to the three point mutual information
bool SortFunctionNoMore(const XJAddress* a, const XJAddress* b, const Environment& environment) {
	 return environment.edges[a->i][a->j].edgeStructure->Ixy_ui > environment.edges[b->i][b->j].edgeStructure->Ixy_ui;
}

class sorterNoMore {
	  Environment& environment;
		public:
	  sorterNoMore(Environment& env) : environment(env) {}
	  bool operator()(XJAddress const* o1, XJAddress const* o2) const {
			return SortFunctionNoMore(o1, o2, environment );
	  }
};


bool skeletonIteration(Environment& environment){
	int iIteration_count = 0;
	int max = 0;

	while( environment.numSearchMore > 0 )
	{
		iIteration_count++;
		//// Get the first edge
		int posX = environment.searchMoreAddress[max]->i;
		int posY = environment.searchMoreAddress[max]->j;

		
		EdgeStructure* topEdgeElt = environment.edges[posX][posY].edgeStructure;

		//// Keep the previous z.name for this edge
		string accepted_z_name = environment.nodes[topEdgeElt->z_name_idx].name;

		//// Ui
		topEdgeElt->ui_vect_idx.push_back(topEdgeElt->zi_vect_idx[topEdgeElt->z_name_idx]);

		//// Zi			
		// ttopEdgeElt->z.name = NA
	   	topEdgeElt->zi_vect_idx[topEdgeElt->z_name_idx] = -1;
	   	topEdgeElt->z_name_idx = -1;
		//// Nxyui			
		topEdgeElt->Nxy_ui = topEdgeElt->Nxyz_ui;
		topEdgeElt->Nxyz_ui = -1;

		//// Set Ixy_ui
		double* v = computeEnsInformationNew(environment, &environment.edges[posX][posY].edgeStructure->ui_vect_idx[0], environment.edges[posX][posY].edgeStructure->ui_vect_idx.size(), 
					NULL, 0, -1,  posX, posY, environment.cplx);


		topEdgeElt->Ixy_ui = v[1];
		topEdgeElt->Nxy_ui = v[0];
		topEdgeElt->cplx = v[2];
		
		free(v);
		double topEdgeElt_kxy_ui = topEdgeElt->cplx;

	 	if(environment.isDegeneracy)
			topEdgeElt_kxy_ui = topEdgeElt->cplx + (topEdgeElt->ui_vect_idx.size()*std::log(static_cast<double>(3)));//log(3));


		if( topEdgeElt->Ixy_ui - topEdgeElt_kxy_ui - environment.logEta <= 0 )
		{   

			//// Move this edge from the list searchMore to phantom
			delete environment.searchMoreAddress[max] ;
			environment.searchMoreAddress.erase(environment.searchMoreAddress.begin() + max);
			environment.numSearchMore--;

			// set the connection to 0 on the adj matrix
			environment.edges[posX][posY].isConnected = 0;
			environment.edges[posY][posX].isConnected = 0;
			//// Save the phantom status
			topEdgeElt->status = 1;	   
		} else {
			//// Reinit Rxyz_ui
			topEdgeElt->Rxyz_ui = environment.thresPc;

			if( topEdgeElt->zi_vect_idx.size() > 0 )
			{
				//// Search for a new contributing node and its rank
				SearchForNewContributingNodeAndItsRank(environment, posX, posY);
			}

			//// Update the information about the edge
			if( topEdgeElt->z_name_idx == -1)
			{			
				//// Move this edge from the list searchMore to noMore
				environment.noMoreAddress.push_back(environment.searchMoreAddress[max]);
				environment.numNoMore++;
			 	environment.searchMoreAddress.erase(environment.searchMoreAddress.begin() + max);
				environment.numSearchMore--;
				// environment.edges[posX][posY].isConnected = 1;
				// environment.edges[posY][posX].isConnected = 1;
				//// Update the status of the edge
				topEdgeElt->status = 3;
			}				
		}

		//// Sort all pairs xy with a contributing node z in decreasing order of their ranks, R(xy;z| )

		max = 0;
		for(int i = 0; i < environment.numSearchMore; i++){
			if(environment.edges[environment.searchMoreAddress[i]->i][environment.searchMoreAddress[i]->j].edgeStructure->Rxyz_ui > 
				environment.edges[environment.searchMoreAddress[max]->i][environment.searchMoreAddress[max]->j].edgeStructure->Rxyz_ui)
			max = i;
		}
	}
	std::sort(environment.noMoreAddress.begin(), environment.noMoreAddress.end(), sorterNoMore(environment));
	return true;
}
