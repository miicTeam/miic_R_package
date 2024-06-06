//*****************************************************************************
// Filename   : tmiic.cpp                           Creation date: 07 may 2020
//
// Author     : Franck SIMON
//
// Description: Store functions for temporal mode of miic (tmiic)
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
// getListLaggedEdges
//-----------------------------------------------------------------------------
// Description: in temporal mode, find the past lagged counterparts of an edge
// assuming stationarity
//
// Params:
// - Environment&: the environment structure
// - int: edge's first node index
// - int: edge's second node index
// Return:
// - std::vector< std::pair<int>, <int> > : list of lagged edges
//--------------------------------------------------------------------------------
std::vector< std::pair<int, int> > getListLaggedEdges
  (Environment& environment, int node1_pos, int node2_pos) {
  std::vector< std::pair<int, int> > list_ret;
  if (  (node1_pos >= environment.n_nodes_not_lagged)
     && (node2_pos >= environment.n_nodes_not_lagged) )
    //
    // The edge is duplicated, this function deals only with the original ones
    //
    return (list_ret);
  //
  // If one of the variables is not lagged, the lag is not constant
  //
  int sav_lag = environment.nodes_lags[node1_pos] - environment.nodes_lags[node2_pos];
  bool same_lag_needed = true;
  if (   (node1_pos < environment.n_nodes_not_lagged)
      && (environment.list_n_layers[node1_pos] <= 1) )
    same_lag_needed = false;
  else if (  (node2_pos < environment.n_nodes_not_lagged)
          && (environment.list_n_layers[node2_pos] <= 1) )
    same_lag_needed = false;
  //
  // Look for the same edge lagged over all layers of history
  //
  while (true) {
    //
    // We shift the nodes positions using pre-computed vector nodes_shifts
    //
    int node1_shift = environment.nodes_shifts[node1_pos];
    int node2_shift = environment.nodes_shifts[node2_pos];
    if ( (node1_shift <= 0) && (node2_shift <= 0) )
      break;
    node1_pos += node1_shift;
    node2_pos += node2_shift;
    //
    // Ensure if both variable are lagged than we keep the same lag when duplicating
    //
    bool same_lag_impossible = false;
    if (same_lag_needed) {
      int new_lag = environment.nodes_lags[node1_pos] - environment.nodes_lags[node2_pos];
      while (sav_lag != new_lag) {
        if (sav_lag < new_lag) {
          int node2_shift = environment.nodes_shifts[node2_pos];
          if (node2_shift <= 0) {
            same_lag_impossible = true;
            break;
          }
          node2_pos += node2_shift;
        }
        else { // sav_lag > new_lag
          int node1_shift = environment.nodes_shifts[node1_pos];
          if (node1_shift <= 0) {
            same_lag_impossible = true;
            break;
          }
          node1_pos += node1_shift;
        }
        new_lag = environment.nodes_lags[node1_pos] - environment.nodes_lags[node2_pos];
      }
    }
    if (same_lag_impossible)
      break;
    //
    // A lagged edge has been found
    //
    list_ret.push_back (std::make_pair (node1_pos, node2_pos) );
  }
  return (list_ret);
}

//-----------------------------------------------------------------------------
// repeatEdgesOverHistory
//-----------------------------------------------------------------------------
// Description: Duplicates edges over the temporal graph
//
// Detail: For consistency and the orientation step when latent variable
// discovery is enabled, we duplicate the edges over the history, assuming
// stationarity, to improve the consistency assessment and the orientations
//
// Params:
// - Environment&: the environment structure
//-----------------------------------------------------------------------------
void repeatEdgesOverHistory (Environment& environment) {
  //
  // We iterate over computed edges to duplicate (if needed)
  // the edges over the history
  //
  auto& edge_list = environment.connected_list;
  std::vector<EdgeID>::size_type size = edge_list.size();
  for (std::vector<EdgeID>::size_type i = 0; i < size; ++i) {
    const Edge& edge_orig = environment.edges (edge_list[i].X, edge_list[i].Y);
    std::vector< std::pair<int, int> > list_lagged = getListLaggedEdges
      (environment, edge_list[i].X, edge_list[i].Y);
    for (auto const& it_lagged : list_lagged) {
      //
      // Duplicate the edge info into environment.edges array
      //
      Edge& edge_to_modif = environment.edges (it_lagged.first, it_lagged.second);
      edge_to_modif.status = edge_orig.status;
      edge_to_modif.status_init = edge_orig.status_init;
      edge_to_modif.status_prev = edge_orig.status_prev;
      edge_to_modif.proba_head = edge_orig.proba_head;

      Edge& edge_to_modif_inverse = environment.edges (it_lagged.second, it_lagged.first);
      edge_to_modif_inverse.status = edge_orig.status;
      edge_to_modif_inverse.status_init = edge_orig.status_init;
      edge_to_modif_inverse.status_prev = edge_orig.status_prev;
      edge_to_modif_inverse.proba_head = edge_orig.proba_head;
      //
      // Add the nodes into connected_list
      //
      environment.connected_list.emplace_back (it_lagged.first, it_lagged.second,
                                environment.edges (it_lagged.first, it_lagged.second) );
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
// edges for the ones that have not been oriented using the unshielded triples,
// that are lagged and have a node on the layer 0.
// Edges matching these criteria will be oriented using time from the oldest
// node to the newest.
//
// N.B.: Edges not having a node on the layer 0 are edges "past only" and are a
// consequence of duplicating edges over history to maximize the number of
// unshielded triples, There is no interest to orient "past only" edges
// using time as they will be removed at the end of the orientation step.
//
// Params:
// - Environment&: the environment structure
// - std::vector<Triple>& : list of unshielded triples
//--------------------------------------------------------------------------------
void completeOrientationUsingTime (Environment& environment,
                                   const std::vector<Triple>& triples)
  {
  // Tail probability to use for lagged edges differs if latent variables are
  // authorized or not:
  // - 0 if no latent var, we are sure that the oldest node is the cause
  // - 0.5 with latent var as we can not be sure that the oldest node is the cause
  //
  double tail_proba = 0;
  if (environment.latent_orientation)
    tail_proba = 0.5;
  //
  // Loop over edges to find edges that were not considered when orienting with
  // open triples but can be oriented using time
  //
  const auto& edge_list = environment.connected_list;
  for (auto iter0 = begin(edge_list); iter0 != end(edge_list); ++iter0)
    {
    int posX = iter0->X, posY = iter0->Y;
    //
    // If edges has no node on the layer 0, it is a duplicated one => skip it
    //
    if ( ! (  (posX < environment.n_nodes_not_lagged)
           || (posY < environment.n_nodes_not_lagged) ) )
      continue;
    //
    // If the edge is contemporaneous, no information from time => skip it
    //
    if (environment.nodes_lags[posX] == environment.nodes_lags[posY])
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
          || ( (triples[i][1] == posY) && (triples[i][2] == posX) ) )
        {
        is_in_triple = true;
        break;
        }
    if (is_in_triple)
      continue;
    //
    // The edge was not in open triples, has a node on layer 0, and is not
    // contemporaneous => we can and need to orient it using time
    // As time goes from past to present: edge orientation is max lag -> min lag
    //
    if (environment.nodes_lags[posX] > environment.nodes_lags[posY])
      {
      if (environment.is_contextual[posX])
        miic::reconstruction::updateAdj(environment, posX, posY, 0, 1);
      else
        miic::reconstruction::updateAdj(environment, posX, posY, tail_proba, 1);
      }
    else
      {
      if (environment.is_contextual[posY])
        miic::reconstruction::updateAdj(environment, posX, posY, 1, 0);
      else
        miic::reconstruction::updateAdj(environment, posX, posY, 1, tail_proba);
      }
    }
  }

//-----------------------------------------------------------------------------
// dropPastEdges
//-----------------------------------------------------------------------------
// Description:
// Drop past edges (the edges having no node of the final time step)
//
// Detail: for consistency assessment or orientation step with latent variable
// discovery is enabled, we duplicated edges over history. Here, we ensure that
// the edges are restored to their state prior before duplication.
//
// Params:
// - Environment&: the environment structure
//--------------------------------------------------------------------------------
void dropPastEdges (Environment& environment) {
  //
  // We iterate over computed edges to find edges previously duplicated using stationary.
  // All the edges duplicated are disconnected and removed from environment.connected_list
  //
  auto it = begin(environment.connected_list);
  while ( it != end(environment.connected_list) ) {
    //
    // When the two nodes are not on the layer 0, the edge is removed
    //
    if (  (it->X >= environment.n_nodes_not_lagged)
       && (it->Y >= environment.n_nodes_not_lagged) )
      it = environment.connected_list.erase (it);
    //
    // When one of the nodes is not lagged (i.e. contextual)
    // and the other not on the layer 0, the edge is removed
    //
    else if (  (it->X < environment.n_nodes_not_lagged)
            && (environment.list_n_layers[it->X] <= 1)
            && (it->Y >= environment.n_nodes_not_lagged) )
      it = environment.connected_list.erase (it);
    else if (  (it->Y < environment.n_nodes_not_lagged)
            && (environment.list_n_layers[it->Y] <= 1)
            && (it->X >= environment.n_nodes_not_lagged) )
      it = environment.connected_list.erase (it);
    else
      it++;
  }
  //
  // We remove from the edges all those having pos > n_nodes_not_lagged
  //
  for (int node1_pos = environment.n_nodes_not_lagged; node1_pos < environment.n_nodes; ++node1_pos)
    for (int node2_pos = environment.n_nodes_not_lagged; node2_pos < environment.n_nodes; ++node2_pos) {
      environment.edges(node1_pos,node2_pos).status = 0;
      environment.edges(node1_pos,node2_pos).status_init = 0;
      environment.edges(node1_pos,node2_pos).status_prev = 0;
      environment.edges(node1_pos,node2_pos).proba_head = -1;
    }
  //
  // In addition, remove lagged edges added on non lagged variable (i.e.:contextual)
  //
  for (int i = 0; i < environment.n_nodes_not_lagged; i++)
    if (environment.is_contextual[i])
      for (int j = environment.n_nodes_not_lagged; j < environment.n_nodes; j++) {
        environment.edges(i, j).status = 0;
        environment.edges(i, j).status_init = 0;
        environment.edges(i, j).status_prev = 0;
        environment.edges(i, j).proba_head = -1;
        environment.edges(j, i).status = 0;
        environment.edges(j, i).status_init = 0;
        environment.edges(j, i).status_prev = 0;
        environment.edges(j, i).proba_head = -1;
        }
}

}  // namespace tmiic
