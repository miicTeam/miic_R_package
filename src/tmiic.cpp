//*****************************************************************************
// Filename   : tmiic.cpp                   Namespace: tmiic                                   
//
// Author     : Franck SIMON                Creation date: 07 may 2020 
//
// Description: Store functions for temporal version of miic (tmiic)
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

#include "structure.h"
#include "utilities.h"
#include "tmiic.h"

//=============================================================================
// CONSTANTS
//=============================================================================
#define _DEBUG 0

//=============================================================================
// NAMESPACES
//=============================================================================
namespace tmiic {

using uint = unsigned int;
using std::string;
using std::cout;
using namespace miic::structure;
using namespace miic::utility;

//=============================================================================
// FUNCTIONS
//=============================================================================
// repeatEdgesOverHistory
//-----------------------------------------------------------------------------
// Description: Repeats edges over the lagged graph
//
// Detail: For orientation step when latent variable discovery is enabled, 
// we need to include in the graph the past edges. So we duplicate over the 
// history all the edges found on nodes of the final timestep (assuming 
// stationarity).
//
// Params: 
// - Environment&: the environment structure
//
// Returns: none
//-----------------------------------------------------------------------------
void repeatEdgesOverHistory (Environment& environment) 
  {
#if _DEBUG
    cout << "\ntmiic::RepeatEdgesOverHistory:\n";
    // printNoMoreAdress  (environment);
    cout << "\n";
    // printEdges (environment);
#endif
  //
  // First, we count how much non lagged nodes we have 
  // => we look for the first lag1 node
  //
  std::regex lag_expr(".*lag");
  uint nodes_cnt = 0;
  while (nodes_cnt < environment.n_nodes)
    {
    string node_name = environment.nodes[nodes_cnt].name;
    int node_lag = std::stoi (std::regex_replace(node_name, lag_expr, "") );
    if (node_lag > 0)
      break;
    nodes_cnt++;
    }
#if _DEBUG
  cout << "\nEdge node tot=" << environment.n_nodes << " not lagged nodes count=" << nodes_cnt << "\n";
#endif
  //
  // Now, we iterate over computed edges to duplicate (if needed)
  // the edges over the history
  //
//TODO
//   for (auto edge_iter = begin(edge_list); edge_iter != end(edge_list); ++edge_iter) 
//     {
//     // For each edge, we compute the lag between the 2 nodes
//     //    
//     int node1_pos = edge_iter->i;
//     int node2_pos = edge_iter->j;
//   
//     string node1_name = environment.nodes[node1_pos].name;
//     string node2_name = environment.nodes[node2_pos].name;
//       
//     std::regex lag_expr(".*lag");
//     int node1_lag = std::stoi (std::regex_replace(node1_name, lag_expr, "") );
//     int node2_lag = std::stoi (std::regex_replace(node2_name, lag_expr, "") );
//     // 
//     // When the lag is smaller than the global lag (tau) in the temporal graph,
//     // We can add lagged edges in the history until tau is reached
//     //    
//     int min_lag = std::min (node1_lag, node2_lag); // Normally 0 as one edge should be lag0
//     int lag_between_nodes = std::abs (node1_lag - node2_lag);
// #if _DEBUG
//     cout << "Edge " << node1_pos << " (" << node1_name << ") - " 
//                     << node2_pos << " (" << node2_name << ")" 
//                     << " min_lag=" << min_lag << " lag_between_nodes=" << lag_between_nodes 
//                     << " => will add " << environment.tau - lag_between_nodes << " into history\n";
// #endif
//     //
//     // Find the Edge structure that will be duplicated 
//     //
//     const Edge& edge_orig = environment.edges[node1_pos][node2_pos];
//     //
//     // To create the same edge over history, we switch the nodes pos by nodes_cnt
//     // until node pos > nb total nodes
//     //
//     while (true)
//       {
//       node1_pos += nodes_cnt;
//       node2_pos += nodes_cnt;
//       if ( (node1_pos >= environment.n_nodes) || (node2_pos >= environment.n_nodes) )
//         break;
// #if _DEBUG
//       cout << "Edge " << node1_pos << " (" << environment.nodes[node1_pos].name << ") - " 
//                       << node2_pos << " (" << environment.nodes[node2_pos].name << ") to add\n";
// #endif
//       //
//       // Add the nodes into noMoreAddress
//       //
// //TODO      environment.noMoreAddress.emplace_back ( new EdgeID (node1_pos, node2_pos) );
//       //
//       // Duplicate the edge info into environment.edges array
//       //
//       Edge& edge_to_modif = environment.edges[node1_pos][node2_pos];
//       
//       edge_to_modif.status = edge_orig.status;
//       edge_to_modif.status_init = edge_orig.status_init;
//       edge_to_modif.status_prev = edge_orig.status_prev;
// #if _DEBUG
//       cout << "Edge " << node1_pos << " (" << environment.nodes[node1_pos].name << ") - " 
//            << node2_pos << " (" << environment.nodes[node2_pos].name << ") modified:"
//            << edge_to_modif.status << "\n";
// #endif
//       
//       Edge& edge_to_modif_inverse = environment.edges[node2_pos][node1_pos];
//       edge_to_modif_inverse.status = edge_orig.status;
//       edge_to_modif_inverse.status_init = edge_orig.status_init;
//       edge_to_modif_inverse.status_prev = edge_orig.status_prev;
// #if _DEBUG
//       cout << "Edge " << node2_pos << " (" << environment.nodes[node2_pos].name << ") - " 
//            << node1_pos << " (" << environment.nodes[node1_pos].name << ") modified:"
//            << edge_to_modif_inverse.status << "\n";
// #endif
//       }
//     }
//TODO

#if _DEBUG
  cout << "\n";
  cout << "\nEnd of tmiic::RepeatEdgesOverHistory:\n";
  printNoMoreAdress  (environment);
  cout << "\n";
  printEdges (environment);
#endif
  }

//-----------------------------------------------------------------------------
// dropPastEdges
//-----------------------------------------------------------------------------
// Description: 
// Drop past edges (the edges having no node of the final timestep)
//
// Detail: For skeleton construction, we use in temporal miic only edges having
// at least one node on the final timestep.
// As for orientation step, when latent variable discovery is enabled, 
// we need to include in the graph the past edges, the orientation step adds
// extra historical edges. So, this function ensures that we restore the edges
// to their state prior to orientation for the next skeleton iteration. 
//
// Params: 
// - Environment&: the environment structure
//
// Returns: none
//--------------------------------------------------------------------------------
void dropPastEdges (Environment& environment)
  {
#if _DEBUG
  cout << "\ntmiic::dropPastEdges:\n";
  printNoMoreAdress  (environment);
  cout << "\n";
  printEdges (environment);
#endif
  //
  // We iterate over computed edges to find edges previously duplicated using stationary.
  // All the edges duplicated are removed from environment.noMoreAddress
  //
  std::regex lag_expr(".*lag");
//TODO	std::vector<EdgeID*>::iterator it = environment.noMoreAddress.begin();
//TODO	while ( it != environment.noMoreAddress.end() ) 
// 	  {
//     // 
//     // For each edge, we compute the minimum lag of the 2 nodes
//     //    
//     int node1_pos = (*it)->i;
//     int node2_pos = (*it)->j;
//   
//     string node1_name = environment.nodes[node1_pos].name;
//     string node2_name = environment.nodes[node2_pos].name;
//       
//     std::regex lag_expr(".*lag");
//     int node1_lag = std::stoi (std::regex_replace(node1_name, lag_expr, "") );
//     int node2_lag = std::stoi (std::regex_replace(node2_name, lag_expr, "") );
//     int min_lag = std::min (node1_lag, node2_lag);
//     // 
//     // When the two lag are greater than 0, the edge is not comptemporanous
//     // and is removed
//     //    
//     if (min_lag > 0)
//       {
//       delete (*it);
// 			it = environment.noMoreAddress.erase(it);
// #if _DEBUG
//       cout << "Edge " << node1_pos << " (" << node1_name << ") - " 
//                       << node2_pos << " (" << node2_name << ")" 
//                       << " min_lag=" << min_lag << " => edge removed from environment.noMoreAddress (new size="
//                       << environment.noMoreAddress.size() << ")\n";
// #endif
//       }
//     else
// 			it++;
// 	  }  
//TODO
  //
  // We iterate over the list of edges to find edges previously duplicated using stationary.
  // For all the edges duplicated, we set the status to 0 (no connection) in environment.edges
  //
  for (int node1_pos = 0; node1_pos < environment.n_nodes; node1_pos++)
    {
    for (int node2_pos = 0; node2_pos < environment.n_nodes; node2_pos++)
  	  {
      // 
      // For each edge, we compute the minimum lag of the 2 nodes
      //    
      string node1_name = environment.nodes[node1_pos].name;
      string node2_name = environment.nodes[node2_pos].name;
        
      std::regex lag_expr(".*lag");
      int node1_lag = std::stoi (std::regex_replace(node1_name, lag_expr, "") );
      int node2_lag = std::stoi (std::regex_replace(node2_name, lag_expr, "") );
      int min_lag = std::min (node1_lag, node2_lag);
      // 
      // When the min lag is greater than 0, the edge is not contemporanous
      // and status is set to 0
      //    
      if (min_lag > 0)
        {
#if _DEBUG
        if (environment.edges[node1_pos][node2_pos].status != 0)
          cout << "Edge " << node1_pos << " (" << node1_name << ") - " 
                          << node2_pos << " (" << node2_name << ")" 
                          << " min_lag=" << min_lag << " => status set to 0 in environment.edges\n";
#endif
        environment.edges[node1_pos][node2_pos].status = 0;
        environment.edges[node1_pos][node2_pos].status_init = 0;
        environment.edges[node1_pos][node2_pos].status_prev = 0;
        }
      }
	  }
    
#if _DEBUG
  cout << "\ntmiic::dropPastEdges end:\n";
  //printNoMoreAdress  (environment);
  cout << "\n";
  //printEdges (environment);
#endif
  }

}  // namespace tmiic
