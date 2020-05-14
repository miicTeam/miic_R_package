#include "orientation_probability.h"

#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <string>

#include "compute_ens_information.h"
#include "proba_orientation.h"
#include "structure.h"
#include "utilities.h"
#include "tmiic.h"

#define _DEBUG 0

namespace miic {
namespace reconstruction {

using uint = unsigned int;
using std::string;
using std::vector;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;
using namespace tmiic;

void getSStructure(Environment& environment, const int posX, const int posY,
    const int posZ, bool isVerbose, vector<vector<int> >& allTpl,
    vector<double>& allI3) {
  // Check if xk belongs to the {ui} of the base
  vector<int> u(environment.edges[posX][posZ].shared_info->ui_vect_idx);
  if (environment.edges[posX][posZ].shared_info->ui_vect_idx.size() > 0) {
    // Remove xk from the {ui} if needed
    for (uint i = 0;
         i < environment.edges[posX][posZ].shared_info->ui_vect_idx.size();
         i++) {
      if (u[i] == posY) {
        u.erase(u.begin() + i);
        break;
      }
    }
  }
  int* ui;
  int* zi = NULL;
  if (u.empty())
    ui = NULL;
  else
    ui = &u[0];

  vector<int> z;
  z.clear();
  z.push_back(posY);

  zi = &z[0];

  double* res = NULL;
  double Is = -1;
  double Cs = -1;
  if (environment.columnAsContinuous[posX] == 0 &&
      environment.columnAsContinuous[posZ] == 0 &&
      environment.columnAsContinuous[posY] == 0) {
    res = computeEnsInformationNew(environment, ui, u.size(), zi, z.size(), -1,
        posX, posZ, environment.cplx, environment.m);

    Is = res[7];
    Cs = res[8];
    if (environment.isK23) {
      if (environment.isDegeneracy) {
        Cs += log(3);
      }
      // to fit eq(20) and eq(22) in BMC Bioinfo 2016
      Is = Is + Cs;
    }
  } else if (environment.typeOfData == 2 ||
             (environment.typeOfData == 1 && environment.isAllGaussian == 0)) {
    res = computeEnsInformationContinuous_Orientation(environment, ui, u.size(),
        zi, posX, posZ, environment.cplx, environment.m);

    Is = res[1];
    Cs = res[2];
    if (environment.isK23) {
      if (environment.isDegeneracy) {
        Cs += log(3);
      }
      // I(x;y;z|ui) - cplx(x;y;z|ui)
      Is = Is - Cs;
    }
  } else {
    // THIS PART IS TO DO ?
    Is =
        computeEnsInformationContinuous_Gaussian(environment, posX, posY, posZ);
    Cs = 0.5 *
         (environment.edges[posX][posY].shared_info->ui_vect_idx.size() + 2) *
         log(environment.numSamples);
    if (environment.isK23) {
      if (environment.isDegeneracy) {
        Cs += log(3);
      }
      Is = Is - Cs;
    }
  }
  delete (res);

  allTpl.emplace_back(std::initializer_list<int>{posX + 1, posY + 1, posZ + 1});
  allI3.push_back(Is);
}

//----------------------------------------------------------------------------
vector<vector<string> > orientationProbability(Environment& environment) {
  vector<vector<string> > orientations;  // output of orientation table
  vector<vector<int> > allTpl;
  vector<double> allI3;
  double* ptrRetProbValues;
  bool isVerbose = environment.isVerbose;
  uint nodes_cnt_not_lagged = 0; // for tmiic, number of nodes (not lagged)

  if (environment.tau > 0) 
    {
    //
    // For temporal graph, repeat edges found over history if latent variable 
    // discovery is activated (TODO add why ...)
    //
    if ( (environment.isLatent) || (environment.isLatentOnlyOrientation) )
      tmiic::repeatEdgesOverHistory (environment);
    //
    // We count how much non lagged nodes we have 
    // => we look for the first lag1 node
    //
    std::regex lag_expr(".*lag");
    while (nodes_cnt_not_lagged < environment.numNodes)
      {
      string node_name = environment.nodes[nodes_cnt_not_lagged].name;
      int node_lag = std::stoi (std::regex_replace(node_name, lag_expr, "") );
      if (node_lag > 0)
        break;
      nodes_cnt_not_lagged++;
      }
    }
  //  
  // GET ALL TPL that could be V/NON-V-STRUCTURES #######
  //
  // In classic (non temporal) mode of miic, only open triplets are considered. 
  // Closed triplets (forming a triangle in the graph) can not be used to 
  // determine the orientation of edges. 
  //
  // In temporal mode, the oriention of temporal edges can be determined
  // using the time. So next to open triplets, we should prepate a list
  // of temporal edges that are, by nature, oriented. 
  // As all the orientation code has been designed for triplets, it is not
  // suited to deal with a temporal edges. As a quick fix, we add closed 
  // triplets as long as they include a temporal information.
  //
  for (uint pos = 0; pos < environment.noMoreAddress.size(); pos++) 
    {
    int posX = environment.noMoreAddress[pos]->i;
    int posY = environment.noMoreAddress[pos]->j;
    //
    // Prepare a list that will contain the neighbors of "x" and the neighbors
    // of "y"
    //
    vector<int> neighboursX;
    vector<int> neighboursY;
    //
    // First loop to identify possible triplets on the Y1 side 
    // of the edges in the loop
    //
    for (uint i = pos + 1; i < environment.noMoreAddress.size(); i++) 
      {
      int posX1 = environment.noMoreAddress[i]->i;
      int posY1 = environment.noMoreAddress[i]->j;
      
      if (environment.tau <= 0) 
        {
        // Classic (non temporal) mode
        //
        if (posY1 == posX && !environment.edges[posY][posX1].status)
          neighboursX.push_back(posX1);
        else if (posY1 == posY && !environment.edges[posX][posX1].status)
          neighboursY.push_back(posX1);
        continue;
        }
      //
      // For temporal mode of miic : we are interested in orienting only
      // the triplets having at least one node contemporaneous (lag = 0)
      //
      int lag_nodeX = posX / nodes_cnt_not_lagged;
      int lag_nodeY = posY / nodes_cnt_not_lagged;
      int lag_nodeX1 = posX1 / nodes_cnt_not_lagged;
      bool has_0_level = ( (lag_nodeX == 0) || (lag_nodeY == 0) || (lag_nodeX1 == 0) );
      
      if (!has_0_level)
        continue;
      // 
      // Triplet has at least one node contemporanous (lag = 0)
      // First, consider open triplets as in normal mode of miic 
      //
      if (posY1 == posX && !environment.edges[posY][posX1].status)
        neighboursX.push_back(posX1);
      else if (posY1 == posY && !environment.edges[posX][posX1].status)
        neighboursY.push_back(posX1);
      else
        {
        //
        // Then, quick fix to orient temporal edges not part of open triplets: 
        // compute the lag level of all the nodes in the triplet and 
        // add the closed triplet if lagged
        //
        bool is_lagged = ( (lag_nodeX != lag_nodeY) || (lag_nodeX != lag_nodeX1) || (lag_nodeY != lag_nodeX1) );
        if ( (posY1 == posX) && is_lagged )
          {
          neighboursX.push_back(posX1);
#if _DEBUG
          std::cout << "Side1: Add closed triplet X=" << posX << " Y=" << posY << " X1=" << posX1 << " in neighboursX\n";         
#endif
          }
        else if ( (posY1 == posY) && is_lagged )
          {
          neighboursY.push_back(posX1);
#if _DEBUG
          std::cout << "Side1: Add closed triplet X=" << posX << " Y=" << posY << " X1=" << posX1 << " in neighboursY\n";         
#endif
          }
        }
      }
    //
    // Second loop to identify possible triplets on the other side 
    // of the edges in the loop, the X1 side
    //
    for (uint i = pos + 1; i < environment.noMoreAddress.size(); i++) 
      {
      int posX1 = environment.noMoreAddress[i]->i;
      int posY1 = environment.noMoreAddress[i]->j;
      
      if (environment.tau <= 0) 
        {
        // Classic (non temporal) case
        //
        if (posX1 == posX && !environment.edges[posY][posY1].status)
          neighboursX.push_back(posY1);
        else if (posX1 == posY && !environment.edges[posX][posY1].status)
          neighboursY.push_back(posY1);
        continue;
        }
      //
      // For temporal mode of miic : we are interested in orienting only
      // the triplets having at least one node contemporanous (lag = 0)
      //
      int lag_nodeX = posX / nodes_cnt_not_lagged;
      int lag_nodeY = posY / nodes_cnt_not_lagged;
      int lag_nodeY1 = posY1 / nodes_cnt_not_lagged;
      bool has_0_level = ( (lag_nodeX == 0) || (lag_nodeY == 0) || (lag_nodeY1 == 0) );
      
      if (!has_0_level)
        continue;
      // 
      // Triplet has at least one node contemporanous (lag = 0)
      // First, consider open triplets as normal mode of miic 
      //
      if (posX1 == posX && !environment.edges[posY][posY1].status)
        neighboursX.push_back(posY1);
      else if (posX1 == posY && !environment.edges[posX][posY1].status)
        neighboursY.push_back(posY1);
      else
        {
        //
        // Quick fix to orient temporal edges not in open temporal triplets: 
        // compute the lag level of all the nodes in the triplet and 
        // add the closed triplet if lagged
        //
        bool is_lagged = ( (lag_nodeX != lag_nodeY) || (lag_nodeX != lag_nodeY1) || (lag_nodeY != lag_nodeY1) );
          
        if ( (posX1 == posX) && is_lagged )
          {
          neighboursX.push_back(posY1);
#if _DEBUG
          std::cout << "Side2: Add closed triplet X=" << posX << " Y=" << posY << " Y1=" << posY1 << " in neighboursX\n";         
#endif
          }
        else if ( (posX1 == posY) && is_lagged )
          {
          neighboursY.push_back(posY1);
#if _DEBUG
          std::cout << "Side2: Add closed triplet X=" << posX <<  " Y=" << posY << " Y1=" << posY1 << " in neighboursY\n";         
#endif
          }
        }
      }
    //
    // The lists of neighbours on both sides is ready for the edge
    //
    int sizeX = neighboursX.size();
    int sizeY = neighboursY.size();
    if (sizeX == 0 && sizeY == 0) continue;

    for (int i = 0; i < sizeX; i++) {
      // Get the structure if any
      getSStructure(
          environment, posY, posX, neighboursX[i], isVerbose, allTpl, allI3);
    }
    // iterate on neighbours of y
    for (int i = 0; i < sizeY; i++) {
      //// Get the structure if any
      getSStructure(
          environment, posX, posY, neighboursY[i], isVerbose, allTpl, allI3);
    }
  }
  
  int* oneLineMatrixallTpl = new int[allTpl.size() * 3];
  for (uint i = 0; i < allTpl.size(); i++) {
    for (int j = 0; j < 3; j++) {
      oneLineMatrixallTpl[j * allTpl.size() + i] = allTpl[i][j];
      allTpl[i][j]--;
    }
  }
  // Compute the arrowhead probability of each edge endpoint
  int myNbrTpl = allTpl.size();
  
#if _DEBUG
    std::cout << "\noneLineMatrixallTpl Init:\n";
    for (int i = 0; i < myNbrTpl; i++) 
      {
      // CAUTION : in oneLineMatrixallTpl nodes indices start from 1 !!!!!!!!
      // The indice shift occurs line 164-165
      int node0 = oneLineMatrixallTpl[0 * myNbrTpl + i] - 1;
      int node1 = oneLineMatrixallTpl[1 * myNbrTpl + i] - 1;
      int node2 = oneLineMatrixallTpl[2 * myNbrTpl + i] - 1;
      std::cout << "!!! Tpl " << node0 << "=" << environment.nodes[node0].name;
      std::cout << " " << node1 << "=" << environment.nodes[node1].name;
      std::cout << " " << node2 << "=" << environment.nodes[node2].name;
      std::cout << " -- I3=" << allI3[i] << "\n";
     }
#endif
  
  if (myNbrTpl > 0) {
    int propag = 0;
    if (environment.isPropagation) propag = 1;
    int degeneracy = 0;
    if (environment.isDegeneracy) degeneracy = 1;
    int latent = 0;
    if (environment.isLatent || environment.isLatentOnlyOrientation) latent = 1;

    ptrRetProbValues = getOrientTplLVDegPropag(environment, myNbrTpl, oneLineMatrixallTpl,
        &allI3[0], latent, degeneracy, propag, environment.halfVStructures);
    // update ptrRetProbValues for possible inconsistencies
    std::map<string, double> probabsMap;
    string s;
    for (uint i = 0; i < allTpl.size(); i++) {
      // 0 -> 1
      s = environment.nodes[allTpl[i][0]].name +
          environment.nodes[allTpl[i][1]].name;
      double proba = ptrRetProbValues[i + (1 * allTpl.size())];
      if (probabsMap.find(s) == probabsMap.end()) {
        // not found
        probabsMap[s] = proba;
      } else {
        // found
        if (abs(probabsMap.find(s)->second - 0.5) < abs(proba - 0.5)) {
          probabsMap[s] = proba;
        }
      }
      // 1 -> 0
      s = environment.nodes[allTpl[i][1]].name +
          environment.nodes[allTpl[i][0]].name;
      proba = ptrRetProbValues[i + (0 * allTpl.size())];
      if (probabsMap.find(s) == probabsMap.end()) {
        // not found
        probabsMap[s] = proba;
      } else {
        // found
        if (abs(probabsMap.find(s)->second - 0.5) < abs(proba - 0.5)) {
          probabsMap[s] = proba;
        }
      }
      // 1 -> 2
      s = environment.nodes[allTpl[i][1]].name +
          environment.nodes[allTpl[i][2]].name;
      proba = ptrRetProbValues[i + (3 * allTpl.size())];
      if (probabsMap.find(s) == probabsMap.end()) {
        // not found
        probabsMap[s] = proba;
      } else {
        // found
        if (abs(probabsMap.find(s)->second - 0.5) < abs(proba - 0.5)) {
          probabsMap[s] = proba;
        }
      }
      // 2 -> 1
      s = environment.nodes[allTpl[i][2]].name +
          environment.nodes[allTpl[i][1]].name;
      proba = ptrRetProbValues[i + (2 * allTpl.size())];
      if (probabsMap.find(s) == probabsMap.end()) {
        // not found
        probabsMap[s] = proba;
      } else {
        // found
        if (abs(probabsMap.find(s)->second - 0.5) < abs(proba - 0.5)) {
          probabsMap[s] = proba;
        }
      }
    }
    for (uint i = 0; i < allTpl.size(); i++) {
      s = environment.nodes[allTpl[i][0]].name +
          environment.nodes[allTpl[i][1]].name;
      ptrRetProbValues[i + (1 * allTpl.size())] = probabsMap.find(s)->second;

      s = environment.nodes[allTpl[i][1]].name +
          environment.nodes[allTpl[i][0]].name;
      ptrRetProbValues[i + (0 * allTpl.size())] = probabsMap.find(s)->second;

      s = environment.nodes[allTpl[i][1]].name +
          environment.nodes[allTpl[i][2]].name;
      ptrRetProbValues[i + (3 * allTpl.size())] = probabsMap.find(s)->second;

      s = environment.nodes[allTpl[i][2]].name +
          environment.nodes[allTpl[i][1]].name;
      ptrRetProbValues[i + (2 * allTpl.size())] = probabsMap.find(s)->second;
    }
    // set results to file
    vector<string> vec;
    vec.push_back("source1");
    vec.push_back("p1");
    vec.push_back("p2");
    vec.push_back("target");
    vec.push_back("p3");
    vec.push_back("p4");
    vec.push_back("source2");
    vec.push_back("NI3");
    vec.push_back("Error");

    orientations.push_back(vec);

    for (uint i = 0; i < allTpl.size(); i++) {
      vec.clear();
      int error = 0;
      int info_sign = round(allI3[i] - 0.5);
      int proba_sign = 0;

      if (ptrRetProbValues[i + (1 * allTpl.size())] > 0.5 &&
          ptrRetProbValues[i + (2 * allTpl.size())] > 0.5)
        proba_sign = -1;
      else
        proba_sign = 1;

      if ((sign(info_sign) != proba_sign && info_sign != 0) ||
          (info_sign == 0 && proba_sign == -1))
        error = 1;

      vec.push_back(environment.nodes[allTpl[i][0]].name);

      std::stringstream output;

      output.str("");
      output << ptrRetProbValues[i + (0 * allTpl.size())];
      vec.push_back(output.str());

      output.str("");
      output << ptrRetProbValues[i + (1 * allTpl.size())];
      vec.push_back(output.str());

      vec.push_back(environment.nodes[allTpl[i][1]].name);

      output.str("");
      output << ptrRetProbValues[i + (2 * allTpl.size())];
      vec.push_back(output.str());

      output.str("");
      output << ptrRetProbValues[i + (3 * allTpl.size())];
      vec.push_back(output.str());

      vec.push_back(environment.nodes[allTpl[i][2]].name);

      output.str("");
      output << allI3[i];
      vec.push_back(output.str());

      output.str("");
      output << error;
      vec.push_back(output.str());

      orientations.push_back(vec);
    }
  }

  vector<double> myProba;
  //// UPDATE ADJ MATRIX
  if (myNbrTpl > 0) {
    for (uint i = 0; i < allTpl.size(); i++) {
      myProba.clear();
      myProba.push_back(ptrRetProbValues[i + (0 * allTpl.size())]);
      myProba.push_back(ptrRetProbValues[i + (1 * allTpl.size())]);
      // If there is AT LEAST ONE arrowhead
      if ((*max_element(myProba.begin(), myProba.end()) - 0.5) > 0) {
        // If there is ONLY ONE arrowhead
        if ((*min_element(myProba.begin(), myProba.end()) - 0.5) <= 0) {
          double a = myProba[0] - *max_element(myProba.begin(), myProba.end());
          // If p1 is max: n1 <-- n2
          if (a == 0) {
            environment.edges[allTpl[i][0]][allTpl[i][1]].status = -2;
            environment.edges[allTpl[i][1]][allTpl[i][0]].status = 2;
          } else {
            environment.edges[allTpl[i][0]][allTpl[i][1]].status = 2;
            environment.edges[allTpl[i][1]][allTpl[i][0]].status = -2;
          }
        } else {
          environment.edges[allTpl[i][0]][allTpl[i][1]].status = 6;
          environment.edges[allTpl[i][1]][allTpl[i][0]].status = 6;
        }
      }
      myProba.clear();
      myProba.push_back(ptrRetProbValues[i + (2 * allTpl.size())]);
      myProba.push_back(ptrRetProbValues[i + (3 * allTpl.size())]);
      // If there is AT LEAST ONE arrowhead
      if ((*max_element(myProba.begin(), myProba.end()) - 0.5) > 0) {
        // If there is ONLY ONE arrowhead
        if ((*min_element(myProba.begin(), myProba.end()) - 0.5) <= 0) {
          double a = myProba[0] - *max_element(myProba.begin(), myProba.end());
          // If p1 is max: n1 <-- n2
          if (a == 0) {
            environment.edges[allTpl[i][1]][allTpl[i][2]].status = -2;
            environment.edges[allTpl[i][2]][allTpl[i][1]].status = 2;
          } else {
            environment.edges[allTpl[i][1]][allTpl[i][2]].status = 2;
            environment.edges[allTpl[i][2]][allTpl[i][1]].status = -2;
          }
        } else {
          environment.edges[allTpl[i][1]][allTpl[i][2]].status = 6;
          environment.edges[allTpl[i][2]][allTpl[i][1]].status = 6;
        }
      }
    }
  }
  delete[] oneLineMatrixallTpl;
  delete[] ptrRetProbValues;

  if (environment.tau > 0) 
    {
    //
    // For temporal graph, when latent variable discovery is enabled, 
    // we duplicated the edges over the history at the beginning 
    // of the function => we need now to drop these extra edges
    //
    if ( (environment.isLatent) || (environment.isLatentOnlyOrientation) )
      tmiic::dropPastEdges (environment);
    //
    // Quick fix for orientation of temporal edges not part of a triplet. 
    // Isolated edges are not oriented by the orientation step as it is based 
    // only on triplets => In temporal mode, orient isolated edges using time 
    //
    for (uint edge_idx = 0; edge_idx < environment.noMoreAddress.size(); edge_idx++) 
      {
      int egde_node0 = environment.noMoreAddress[edge_idx]->i;
      int egde_node1 = environment.noMoreAddress[edge_idx]->j;
      //
      // If the edge is not temporal, nothing that we can do
      //
      int lag_node0 = egde_node0 / nodes_cnt_not_lagged;
      int lag_node1 = egde_node1 / nodes_cnt_not_lagged;
      if (lag_node0 == lag_node1)
        continue;
      //
      // The edge is temporal, check if edge was in the triplet list
      //
      for (uint tpl_idx = 0; tpl_idx < allTpl.size(); tpl_idx++) 
        {
        int tpl_node0 = allTpl[tpl_idx][0];
        int tpl_node1 = allTpl[tpl_idx][1];
        int tpl_node2 = allTpl[tpl_idx][2];
        if (  (egde_node0 == tpl_node0) && (egde_node1 == tpl_node1)
           || (egde_node0 == tpl_node1) && (egde_node1 == tpl_node0)
           || (egde_node0 == tpl_node1) && (egde_node1 == tpl_node2)
           || (egde_node0 == tpl_node2) && (egde_node1 == tpl_node1) )
           continue;
        }
      //
      // Edge is temporal and was not in the triplet list, look if not oriented
      //
      int orient = environment.edges[egde_node0][egde_node1].status;
      if (orient != 1)
        continue;
      //
      // Edge is temporal and is not oriented => orient it
      //
      if (lag_node0 < lag_node1)
        orient = -2;
      else
        orient = 2;
      environment.edges[egde_node0][egde_node1].status = orient;
      environment.edges[egde_node1][egde_node0].status = -orient;
#if _DEBUG
      std::cout << "Edge " << environment.nodes[egde_node0].name
                       << "-" << environment.nodes[egde_node1].name
                       << " oriented using time, orientation=" << orient << "\n";
#endif
      }
    //
    // Check if orientation found by miic is aligned with time
    //
    for (uint edge_idx = 0; edge_idx < environment.noMoreAddress.size(); edge_idx++) 
      {
      int egde_node0 = environment.noMoreAddress[edge_idx]->i;
      int egde_node1 = environment.noMoreAddress[edge_idx]->j;

      int lag_node0 = egde_node0 / nodes_cnt_not_lagged;
      int lag_node1 = egde_node1 / nodes_cnt_not_lagged;
      if (lag_node0 == lag_node1)
        continue;
      
      int orient_from_time;
      if (lag_node0 < lag_node1)
        orient_from_time = -2;
      else
        orient_from_time = 2;
      
      int orient_from_miic = environment.edges[egde_node0][egde_node1].status;
      
      if (orient_from_miic != orient_from_time)
        std::cout << "Warning: Edge " << environment.nodes[egde_node0].name
                       << "-" << environment.nodes[egde_node1].name
                       << ", miic orientation=" << orient_from_miic
                       << " differs from time orientation=" << orient_from_time << "\n";
      }
    }
  
#if _DEBUG
  std::cout << "\norientationProbability end:\n";
  printEdges (environment);
#endif

  return orientations;
}

}  // namespace reconstruction
}  // namespace miic
