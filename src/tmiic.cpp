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

#include "debug.h"
//=============================================================================
// CONSTANTS
//=============================================================================
#define _DEBUG 1

//=============================================================================
// NAMESPACES
//=============================================================================
namespace tmiic {

using std::string;
using namespace miic::structure;
using namespace miic::utility;

//=============================================================================
// FUNCTIONS
//=============================================================================
// countNbNodesNotLagged
//-----------------------------------------------------------------------------
// Description: Counts the number of nodes not lagged
//
// Detail: Counts the number of nodes ending with lag0
//
// Params: 
// - Environment&: the environment structure
//
// Returns: int, the number of nodes not lagged
//-----------------------------------------------------------------------------
int countNbNodesNotLagged (Environment& environment) 
  {
  std::regex lag_expr(".*lag");
  int n_nodes_not_lagged = 0;
  while (n_nodes_not_lagged < environment.n_nodes)
    {
    string node_name = environment.nodes[n_nodes_not_lagged].name;
    int node_lag = std::stoi (std::regex_replace(node_name, lag_expr, "") );
    if (node_lag > 0)
      break;
    n_nodes_not_lagged++;
    }
  return (n_nodes_not_lagged);
  }

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
void repeatEdgesOverHistory (Environment& environment) 
  {
#if _DEBUG
  Rcpp::Rcout << "\ntmiic::RepeatEdgesOverHistory:\n";
  Rcpp::Rcout << "Connected list:\n";
  miic::debug::printConnectedList  (environment);
  Rcpp::Rcout << "Edges:\n";
  miic::debug::printEdges (environment);
  Rcpp::Rcout << "\nEdge node tot=" << environment.n_nodes 
              << " not lagged nodes count=" << environment.n_nodes_not_lagged << "\n";
#endif
  //
  // Now, we iterate over computed edges to duplicate (if needed)
  // the edges over the history
  //
  int n_layers_tot = environment.n_nodes / environment.n_nodes_not_lagged;
  auto& edge_list = environment.connected_list;
  Rcpp::Rcout << edge_list.size() << "\n";
  auto edge_end = end(edge_list);
  for (auto iter0 = begin(edge_list); iter0 != edge_end; ++iter0) 
    {
    // 
    // For each edge, we compute the layer of each node
    //    
    int node1_pos = iter0->X;
    int node2_pos = iter0->Y;
    int node1_layer = node1_pos / environment.n_nodes_not_lagged;
    int node2_layer = node2_pos / environment.n_nodes_not_lagged;
    // 
    // When the layer is smaller than the total layers in the temporal graph,
    // We can add lagged edges in the history until total layers is reached
    //    
    int min_layer = std::min (node1_layer, node2_layer); // Normally 0 as one node should be lag0
    int n_layers_between_nodes = std::abs (node1_layer - node2_layer);
#if _DEBUG
    Rcpp::Rcout << "Edge " << node1_pos << " (" << environment.nodes[node1_pos].name << ") - " 
                    << node2_pos << " (" << environment.nodes[node2_pos].name << ")" 
                    << " min_layer=" << min_layer  
                    << " layers_between_nodes=" << n_layers_between_nodes 
                    << " => will add " << n_layers_tot - 1 - min_layer - n_layers_between_nodes 
                    << " into history\n";
#endif
    //
    // Find the Edge structure that will be duplicated 
    //
    const Edge& edge_orig = environment.edges(node1_pos, node2_pos);
    //
    // To create the same edge over history, we switch the nodes pos by n_nodes_not_lagged
    // until node pos > nb total nodes
    //
    while (true)
      {
      node1_pos += environment.n_nodes_not_lagged;
      node2_pos += environment.n_nodes_not_lagged;
      if ( (node1_pos >= environment.n_nodes) || (node2_pos >= environment.n_nodes) )
        break;
#if _DEBUG
      Rcpp::Rcout << "Edge " << node1_pos << " (" << environment.nodes[node1_pos].name << ") - "
                      << node2_pos << " (" << environment.nodes[node2_pos].name << ") to add\n";
#endif
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
#if _DEBUG
      Rcpp::Rcout << "Edge " << node1_pos << " (" << environment.nodes[node1_pos].name << ") - "
           << node2_pos << " (" << environment.nodes[node2_pos].name << ") modified:"
           << edge_to_modif.status << "\n";
#endif

      Edge& edge_to_modif_inverse = environment.edges(node2_pos,node1_pos);
      edge_to_modif_inverse.status = edge_orig.status;
      edge_to_modif_inverse.status_init = edge_orig.status_init;
      edge_to_modif_inverse.status_prev = edge_orig.status_prev;
#if _DEBUG
      Rcpp::Rcout << "Edge " << node2_pos << " (" << environment.nodes[node2_pos].name << ") - "
           << node1_pos << " (" << environment.nodes[node1_pos].name << ") modified:"
           << edge_to_modif_inverse.status << "\n";
#endif
      }
    }

#if _DEBUG
  Rcpp::Rcout << "\n";
  Rcpp::Rcout << "\nEnd of tmiic::RepeatEdgesOverHistory:\n";
  Rcpp::Rcout << "Connected list:\n";
  miic::debug::printConnectedList  (environment);
  Rcpp::Rcout << "Edges:\n";
  miic::debug::printEdges (environment);
#endif
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
void dropPastEdges (Environment& environment)
  {
#if _DEBUG
  Rcpp::Rcout << "\ntmiic::dropPastEdges:\n";
  miic::debug::printEdges (environment);
#endif
  //
  // We iterate over computed edges to find edges previously duplicated using stationary.
  // All the edges duplicated are disconnected and removed from environment.connected_list
  //
  std::regex lag_expr(".*lag");
  auto it = begin(environment.connected_list); 
  while ( it != end(environment.connected_list) ) 
	  {
    // 
    // For each edge, we compute the minimum lag of the 2 nodes
    //    
    int node1_pos = it->X;
    int node2_pos = it->Y;
    int node1_layer = node1_pos / environment.n_nodes_not_lagged;
    int node2_layer = node2_pos / environment.n_nodes_not_lagged;
    // 
    // When the two lag are greater than 0, the edge is not contemporanous
    // and is removed
    //    
    if ((node1_layer > 0) && (node2_layer > 0))
      {
      it = environment.connected_list.erase(it);
#if _DEBUG
      Rcpp::Rcout << "Edge " << node1_pos << " (" << environment.nodes[node1_pos].name << ") - " 
                      << node2_pos << " (" << environment.nodes[node2_pos].name << ")" 
                      << " => edge removed from environment.connected_list (new size="
                      << environment.connected_list.size() << ")\n";
#endif
      }
    else
			it++;
	  }  
  //
  // We remove from the edges all those having pos > n_nodes_not_lagged
  //
  for (int node1_pos = environment.n_nodes_not_lagged; node1_pos < environment.n_nodes; node1_pos++)
    for (int node2_pos = environment.n_nodes_not_lagged; node2_pos < environment.n_nodes; node2_pos++)
  	  {
#if _DEBUG
      if (environment.edges(node1_pos,node2_pos).status != 0)
        Rcpp::Rcout << "Edge " << node1_pos << " (" << environment.nodes[node1_pos].name << ") - " 
                        << node2_pos << " (" << environment.nodes[node2_pos].name << ")" 
                        << " => status set to 0 in environment.edges\n";
#endif
      environment.edges(node1_pos,node2_pos).status = 0;
      environment.edges(node1_pos,node2_pos).status_init = 0;
      environment.edges(node1_pos,node2_pos).status_prev = 0;
      }
    
#if _DEBUG
  Rcpp::Rcout << "\ntmiic::dropPastEdges end:\n";
  Rcpp::Rcout << "Connected list:\n";
  miic::debug::printConnectedList  (environment);
  Rcpp::Rcout << "Edges:\n";
  miic::debug::printEdges (environment);
#endif
  }

}  // namespace tmiic
