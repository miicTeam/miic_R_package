#include "utilities.h"

#include <Rcpp.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "debug.h"
#include "environment.h"

using namespace miic::structure;

namespace miic {
namespace debug {

using uint = unsigned int;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;

// string printNodesName(const Environment& environment) {
//   string s = "";
//   for (uint i = 0; i < environment.numNodes; i++) {
//     Rcpp::Rcout << environment.nodes[i].name;
//     if (i + 1 < environment.numNodes) Rcpp::Rcout << " ";
//   }
//   Rcpp::Rcout << "\n";
//   return s;
// }
// 
//-----------------------------------------------------------------------------
// debugEdges
//-----------------------------------------------------------------------------
// Description: print the list of edges in environment.edges
// (only half of them as the matrix is symetrical)
//
// Params:
// - Environment&: the environment structure
// - filter_status: boolean. Optional, true by default.
//   When true, displays only edges with status != 0.
//   When false, displays all edges.
// - half_only: boolean. Optional, true by default.
//   When true, displays only edges with col > row (matrix is symetrical).
//   When false, displays all edges.
//
// Returns: None
//-----------------------------------------------------------------------------
void debugEdges (Environment &environment, std::string title, 
                 bool filter_status, bool half_only)
  {
  string half_or_full_str = "";
  if (half_only)
    half_or_full_str = "half list";
  else
    half_or_full_str = "full list";
  if (title != "")
    Rcpp::Rcout << title << ": ";
  if (filter_status)
    Rcpp::Rcout << "List of edges in environment.edges (" << half_or_full_str << " filtered on status != 0):\n";
  else
    Rcpp::Rcout << "List of edges in environment.edges (" << half_or_full_str << "):\n";

  Rcpp::Rcout << "node\tnode\tstat\tstat\tstat\tconnect.  Nxy      mutInfo      cplx_noU   Nxy_ui cplx Ixy_ui Rxyz_ui z_name_idx zi_vect_idx  ui_vect_idx\n";
  Rcpp::Rcout << "  1 \t  2 \t    \tinit\tprev\t        nb joint mutual info   Complexity                      Score  Index last candid.nodes  Indice of\n";
  Rcpp::Rcout << "    \t    \t    \t   \t    \t         factors    without       without                       best     best     contributing  separating\n";
  Rcpp::Rcout << "    \t    \t    \t   \t    \t          not NA  conditioning conditioning                    contrib  contrib   cond. indep.   nodes\n";

  for (uint i = 0; i < environment.n_nodes; i++)
    {
    uint j;
    if (half_only)
      j = i + 1;
    else
      j = 0;
    for (; j < environment.n_nodes; j++)
      {
      const Edge& one_edge = environment.edges(i,j);
      std::shared_ptr<EdgeSharedInfo> shared_info_ptr = one_edge.shared_info;

      if ( (!filter_status) || (environment.edges(i,j).status) )
        {
        Rcpp::Rcout << environment.nodes[i].name << "\t"
             << environment.nodes[j].name << "\t"
             << environment.edges(i,j).status << "\t"
             << environment.edges(i,j).status_init << "\t"
             << environment.edges(i,j).status_prev << "\t";

        if (shared_info_ptr != NULL)
          Rcpp::Rcout << environment.edges(i,j).shared_info->connected << "\t"
               << environment.edges(i,j).shared_info->Nxy << "\t"
               << environment.edges(i,j).shared_info->Ixy << "\t"
               << environment.edges(i,j).shared_info->exp_shuffle << "\t"
               << environment.edges(i,j).shared_info->Nxy_ui << "\t"
//               << environment.edges(i,j).shared_info->cplx << "\t"
               << environment.edges(i,j).shared_info->Ixy_ui << "\t"
               << environment.edges(i,j).shared_info->Rxyz_ui << "\t"
               << environment.edges(i,j).shared_info->top_z << "\t";
//               << environment.edges(i,j).shared_info->zi_list << "\t"
//               << environment.edges(i,j).shared_info->ui_list << "\t";
        Rcpp::Rcout << "\n";
        }
      }
    }
  }

//-----------------------------------------------------------------------------
// debugConnectedList
//-----------------------------------------------------------------------------
// Desciption: print the list of edges in environment.noMoreAdress
//
// Params:
// - Environment&: the environment structure
//
// Returns: None
//-----------------------------------------------------------------------------
void debugConnectedList (Environment& environment, std::string title)
  {
  if (title != "")
    Rcpp::Rcout << title << ": ";
  Rcpp::Rcout << "List of edges in connected_list:\n";
  Rcpp::Rcout << "Node 1 idx + (name) - Node 2 idx + (name)\n";

  for (auto iter0 = begin(environment.connected_list); iter0 != end(environment.connected_list); ++iter0) 
    {
    int posX = iter0->X;
    int posY = iter0->Y;
    Rcpp::Rcout << posX << " (" << environment.nodes[posX].name << ") - "
         << posY << " (" << environment.nodes[posY].name << ")" << "\n";
    }
  Rcpp::Rcout << "\n";
  }

//-----------------------------------------------------------------------------
// debugAdjacencyMatrix
//-----------------------------------------------------------------------------
// Desciption: 
// print an adjacency matrix from the list of edges environment.edges
//
// Params: 
// - Environment&: the environment structure
// - status_field: string. Optional, "status" by default.
//   Control what information is used to print the adjacency matrix.
//   Possible values are "status", "status_init" and "status_prev"  
//
// Returns: None
//-----------------------------------------------------------------------------
void debugAdjacencyMatrix (Environment& environment, std::string title, 
                           std::string status_field) 
  {
  if (title != "")
    Rcpp::Rcout << title << ": ";
  Rcpp::Rcout << "Adjacency matrix using environment.edges on col " << status_field << "\n";
  Rcpp::Rcout << std::setw (20);
  for (uint i = 0; i < environment.n_nodes; i++) 
    Rcpp::Rcout  <<  "\t" << environment.nodes[i].name;
  Rcpp::Rcout << "\n";
  for (uint i = 0; i < environment.n_nodes; i++) 
    {
    Rcpp::Rcout << std::setw (20) << environment.nodes[i].name;
    for (uint j = 0; j < environment.n_nodes; j++) 
      {
      if (status_field == "status_init")      
        Rcpp::Rcout << "\t" << environment.edges(i, j).status_init;
      else if (status_field == "status_prev")      
        Rcpp::Rcout << "\t" << environment.edges(i, j).status_prev;
      else
        Rcpp::Rcout << "\t" << environment.edges(i, j).status;
      }
    Rcpp::Rcout << "\n";
    }
  }


// void printMatrix(const Environment& environment, string type) {
//   if (type.compare("string") == 0) {
//     Rcpp::Rcout << "Data matrix of strings\n";
//     printNodesName(environment);
//     for (uint i = 0; i < environment.numSamples; i++) {
//       for (uint j = 0; j < environment.numNodes; j++) {
//         Rcpp::Rcout << environment.data(i,j) << " ";
//       }
//       Rcpp::Rcout << "\n";
//     }
//   } else if (type.compare("factors") == 0) {
//     Rcpp::Rcout << "Data matrix of factors\n";
//     printNodesName(environment);
//     for (uint i = 0; i < environment.numSamples; i++) {
//       for (uint j = 0; j < environment.numNodes; j++) {
//         Rcpp::Rcout << environment.dataNumeric(i,j) << " ";
//       }
//       Rcpp::Rcout << "\n";
//     }
//   }
// }

//-----------------------------------------------------------------------------
// debugTriples
//-----------------------------------------------------------------------------
// Description: print the list of triples
//
// Params:
// - const std::vector<Triple>&: the triples list
// Returns: None
//-----------------------------------------------------------------------------
void debugTriples (const structure::Environment& environment, 
                   const std::vector<Triple>& triples, 
                   std::string title)
  {
  if (title != "")
    Rcpp::Rcout << title << ": ";
  Rcpp::Rcout << "List of triples=\n";
  Rcpp::Rcout << "X\tZ\tY\n";
  for (int i = 0; i < triples.size(); i++)
    {
    if (triples[i][0] != 0)
      continue;
    if (triples[i][1] != 1)
      continue;
    if (triples[i][2] != 2)
      continue;
    std::string node1 = environment.nodes[triples[i][0]].name;
    std::string node2 = environment.nodes[triples[i][1]].name;
    std::string node3 = environment.nodes[triples[i][2]].name;
    Rcpp::Rcout << node1 << "\t" << node2 << "\t" << node3 << "\n";
    }
  }

//-----------------------------------------------------------------------------
// debugProbaArray
//-----------------------------------------------------------------------------
// Description: print the list of triples
//
// Params:
// - const std::vector<Triple>&: the triples list
// Returns: None
//-----------------------------------------------------------------------------
void debugProbaArray (const structure::Environment& environment, 
                   const std::vector<Triple>& triples, 
                   const std::vector<ProbaArray>& prob_array, 
                   std::string title)
  {
  if (title != "")
    Rcpp::Rcout << title << ": ";
  Rcpp::Rcout << "List of triples=\n";
  Rcpp::Rcout << "X\tZ\tY\n";
  for (int i = 0; i < triples.size(); i++)
    {
    if (triples[i][0] != 0)
      continue;
    if (triples[i][1] != 1)
      continue;
    if (triples[i][2] != 2)
      continue;
    std::string node1 = environment.nodes[triples[i][0]].name;
    std::string node2 = environment.nodes[triples[i][1]].name;
    std::string node3 = environment.nodes[triples[i][2]].name;
    Rcpp::Rcout << "triple index= " << i << ": " << node1 << "\t" << prob_array[i][0] 
                << " - " << prob_array[i][1] << "\t" << node2 << "\t" << prob_array[i][2] 
                << " - " << prob_array[i][3] << "\t" << node3 << "\n";
    }
  }

//-----------------------------------------------------------------------------
// debugProbaArray
//-----------------------------------------------------------------------------
// Description: print the list of triples
//
// Params:
// - const std::vector<Triple>&: the triples list
// Returns: None
//-----------------------------------------------------------------------------
void debugProbaMap (const structure::Environment& environment, 
                   const std::map<std::pair<int, int>, double>& proba_map, 
                   std::string title)
  {
  if (title != "")
    Rcpp::Rcout << title << ": ";
  Rcpp::Rcout << "Proba map\n";
  Rcpp::Rcout << "Node1\tNode2\tProba\n";
  for (auto it = proba_map.begin(); it != proba_map.end(); it++) 
    {
    int node1_pos = (it->first).first;
    int node2_pos = (it->first).second;
    if (node1_pos > 2)
      continue;
    if (node2_pos > 2)
      continue;
    Rcpp::Rcout << environment.nodes[node1_pos].name 
                << "\t" << environment.nodes[node2_pos].name 
                << "\t" << it->second << "\n";      
    }
  }
}  // namespace debug
}  // namespace miic
