//*****************************************************************************
// Filename   : tmiic.cpp                   Namespace: tmiic                                   
//
// Author     : Franck SIMON                Creation date: 07 may 2020 
//
// Description: Store functions for temporal mode of miic (tmiic)
//
// Changes history:
// - 07 may 2020 : initial version
//*****************************************************************************

//=============================================================================
// INCLUDE SECTION
//=============================================================================
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "environment.h"
#include "structure.h"
#include "tmiic.h"

//=============================================================================
// CONSTANTS
//=============================================================================

//=============================================================================
// NAMESPACES
//=============================================================================
namespace tmiic {

using std::string;
using namespace miic::structure;
using namespace miic::utility;

//-----------------------------------------------------------------------------
// repeatEdgesOverHistory
//-----------------------------------------------------------------------------
// Description: Duplicates edges over the temporal graph
//
// Detail: For the orientation step when latent variable discovery is enabled, 
// we duplicate the edges over the history, assuming stationarity, to improve 
// the orientation step
//
// Params: 
// - Environment&: the environment structure
//
// Returns: none
//-----------------------------------------------------------------------------
void repeatEdgesOverHistory (Environment& environment) {
  //
  // Now, we iterate over computed edges to duplicate (if needed)
  // the edges over the history
  //
  int n_nodes_not_lagged = environment.n_nodes / (environment.tau + 1);
  auto& edge_list = environment.connected_list;
  auto edge_end = end(edge_list);
  for (auto iter0 = begin(edge_list); iter0 != edge_end; ++iter0) {
    //
    // Find the Edge structure that will be duplicated 
    //
    int node1_pos = iter0->X;
    int node2_pos = iter0->Y;
    const Edge& edge_orig = environment.edges(node1_pos, node2_pos);
    //
    // To create the same edge over history, we switch the nodes pos by n_nodes_not_lagged
    // until node pos > nb total nodes
    //
    while (true) {
      node1_pos += n_nodes_not_lagged;
      node2_pos += n_nodes_not_lagged;
      if ( (node1_pos >= environment.n_nodes) || (node2_pos >= environment.n_nodes) )
        break;
     //
      // Duplicate the edge info into environment.edges array
      //
      Edge& edge_to_modif = environment.edges(node1_pos,node2_pos);
      edge_to_modif.status = edge_orig.status;
      edge_to_modif.status_init = edge_orig.status_init;
      edge_to_modif.status_prev = edge_orig.status_prev;
      //
      // Add the nodes into connected_list
      //
      environment.connected_list.emplace_back (node1_pos, node2_pos, environment.edges(node1_pos,node2_pos) );

      Edge& edge_to_modif_inverse = environment.edges(node2_pos,node1_pos);
      edge_to_modif_inverse.status = edge_orig.status;
      edge_to_modif_inverse.status_init = edge_orig.status_init;
      edge_to_modif_inverse.status_prev = edge_orig.status_prev;
    }
  }
}

//-----------------------------------------------------------------------------
// dropPastEdges
//-----------------------------------------------------------------------------
// Description: 
// Drop past edges (the edges having no node of the final timestep)
//
// Detail: In orientation step, when latent variable discovery is enabled, 
// we duplicated stationary edges over history. Here, we ensure that the edges
// are restored to their state prior to orientation. 
//
// Params: 
// - Environment&: the environment structure
//
// Returns: none
//--------------------------------------------------------------------------------
void dropPastEdges (Environment& environment) {
  int n_nodes_not_lagged = environment.n_nodes / (environment.tau + 1);
  //
  // We iterate over computed edges to find edges previously duplicated using stationary.
  // All the edges duplicated are disconnected and removed from environment.connected_list
  //
  auto it = begin(environment.connected_list); 
  while ( it != end(environment.connected_list) ) {
    // 
    // For each edge, we compute the minimum lag of the 2 nodes
    //    
    int node1_pos = it->X;
    int node2_pos = it->Y;
    int node1_layer = node1_pos / n_nodes_not_lagged;
    int node2_layer = node2_pos / n_nodes_not_lagged;
    // 
    // When the two lag are greater than 0, the edge is not contemporanous
    // and is removed
    //    
    if ((node1_layer > 0) && (node2_layer > 0))
      it = environment.connected_list.erase(it);
    else
      it++;
  }  
  //
  // We remove from the edges all those having pos > n_nodes_not_lagged
  //
  for (int node1_pos = n_nodes_not_lagged; node1_pos < environment.n_nodes; node1_pos++)
    for (int node2_pos = n_nodes_not_lagged; node2_pos < environment.n_nodes; node2_pos++) {
      environment.edges(node1_pos,node2_pos).status = 0;
      environment.edges(node1_pos,node2_pos).status_init = 0;
      environment.edges(node1_pos,node2_pos).status_prev = 0;
    }
}

}  // namespace tmiic
