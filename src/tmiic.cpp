//*****************************************************************************
// Filename   : tmiic.cpp                   Namespace: tmiic                                   
//
// Author     : Franck SIMON                Creation date: 07 may 2020 
//
// Description: Store functions for temporal mode of miic (tmiic)
//
// Changes history:
// - 07 may 2020 : initial version
// - 13 oct 2020 : add of completeOrientationUsingTime
//*****************************************************************************

//=============================================================================
// INCLUDE SECTION
//=============================================================================
#include "orientation.h"
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

      Edge& edge_to_modif_inverse = environment.edges(node2_pos,node1_pos);
      edge_to_modif_inverse.status = edge_orig.status;
      edge_to_modif_inverse.status_init = edge_orig.status_init;
      edge_to_modif_inverse.status_prev = edge_orig.status_prev;
      //
      // Add the nodes into connected_list
      //
      environment.connected_list.emplace_back (node1_pos, node2_pos, environment.edges(node1_pos,node2_pos) );
    }
  }
}

//-----------------------------------------------------------------------------
// completeOrientationUsingTime
//-----------------------------------------------------------------------------
// Description: Complete the orientations with the orientation of temporal  
// edges that were not previously oriented 
//
// Detail: completeOrientationUsingTime will look in the list of connected
// edges the ones that have not been oriented using the unshielded triples, 
// are lagged and have a lag0 node (edges not having a lag0 are only a 
// consequence of duplicating edges over history to maximize the number of 
// unshielded triples and we are not interested to orient them as they will
// be removed at the end of the orientation step). Edges matching these
// criteria will be oriented using time from oldest node to newest node.
//
// Params:
// - Environment&: the environment structure
// - std::vector<Triple>& : list of unshielded triples
//--------------------------------------------------------------------------------
void completeOrientationUsingTime (Environment& environment, 
                                   const std::vector<Triple>& triples) {
  int n_nodes_not_lagged = environment.n_nodes / (environment.tau + 1);
  const auto& edge_list = environment.connected_list;
  for (auto iter0 = begin(edge_list); iter0 != end(edge_list); ++iter0) {
    int posX = iter0->X, posY = iter0->Y;
    //
    // If edges has no lag0 node, it is a duplicated one => skip it
    //
    if ( ! (   (posX < n_nodes_not_lagged)
            || (posY < n_nodes_not_lagged) ) )
      continue;
    //
    // If the edge is not lagged, no information from time => skip it
    //
    int layer_x = int(posX / n_nodes_not_lagged);
    int layer_y = int(posY / n_nodes_not_lagged);
    if (layer_x == layer_y)
      continue;
    //
    // If edge is in triple, head/tail probas have already been computed
    // and applied in adjacency matrix
    //
    bool is_in_triple = false;
    for (unsigned int i = 0; i < triples.size() ; i++)
      if (   ( (triples[i][0] == posX) && (triples[i][1] == posY) )
          || ( (triples[i][0] == posY) && (triples[i][1] == posX) )
          || ( (triples[i][1] == posX) && (triples[i][2] == posY) )
          || ( (triples[i][1] == posY) && (triples[i][2] == posX) ) ) {
        is_in_triple = true;
        break;
      }
    if (is_in_triple)
      continue;
    //
    // The edge was not in open triples, has a lag0 node and is lagged
    // => we can and need to orient it using time 
    // As time goes from past to present: edge orientation is max layer -> min layer
    //
    if (layer_x > layer_y)
      miic::reconstruction::updateAdj(environment, posX, posY, 0, 1);
    else
      miic::reconstruction::updateAdj(environment, posX, posY, 1, 0);
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
    // For each edge, we compute the layer
    //    
    int node1_pos = it->X;
    int node2_pos = it->Y;
    int node1_layer = node1_pos / n_nodes_not_lagged;
    int node2_layer = node2_pos / n_nodes_not_lagged;
    // 
    // When the two lags are greater than 0, the edge is not contemporaneous
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
