#include "orientation.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define _USE_MATH_DEFINES
#include <cmath>
#include <tuple>  // std::tie
#include <map>

#include "get_information.h"
#include "linear_allocator.h"
#include "proba_orientation.h"
#include "tmiic.h"

#include "debug.h"
#define _DEBUG 1

namespace miic {
namespace reconstruction {

using std::string;
using std::vector;
using std::fabs;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;

namespace {

bool acceptProba(double proba, double ori_proba_ratio) {
  if (proba <= 0.5) return false;
  return (1 - proba) / proba < ori_proba_ratio;
}

// y2x: probability that there is an arrow from node y to x
// x2y: probability that there is an arrow from node x to y
void updateAdj(Environment& env, int x, int y, double y2x, double x2y) {
  double lower, higher;
  std::tie(lower, higher) = std::minmax(y2x, x2y);
  // No arrowhead
  if (higher <= 0.5) return;
  // Only one arrowhead
  if (lower <= 0.5) {
    if (y2x == higher && acceptProba(y2x, env.ori_proba_ratio)) {
      env.edges(x, y).status = -2;
      env.edges(y, x).status = 2;
    } else if (acceptProba(x2y, env.ori_proba_ratio)) {
      env.edges(x, y).status = 2;
      env.edges(y, x).status = -2;
    }
  } else if (acceptProba(x2y, env.ori_proba_ratio) &&
             acceptProba(y2x, env.ori_proba_ratio)) {
    env.edges(x, y).status = 6;
    env.edges(y, x).status = 6;
  }
}

}  // anonymous namespace

vector<vector<string>> orientationProbability(Environment& environment) {
  vector<Triple> triples;
  if (environment.tau <= 0)
    {
    // In regular mode , get all unshielded triples X -- Z -- Y
    //
    const auto& edge_list = environment.connected_list;
    for (auto iter0 = begin(edge_list); iter0 != end(edge_list); ++iter0) 
      {
      int posX = iter0->X, posY = iter0->Y;
      for (auto iter1 = iter0 + 1; iter1 != end(edge_list); ++iter1) 
        {
        int posX1 = iter1->X, posY1 = iter1->Y;
        if (posY1 == posX && !environment.edges(posY, posX1).status)
          triples.emplace_back(Triple{posY, posX, posX1});
        else if (posY1 == posY && !environment.edges(posX, posX1).status)
          triples.emplace_back(Triple{posX, posY, posX1});
        if (posX1 == posX && !environment.edges(posY, posY1).status)
          triples.emplace_back(Triple{posY, posX, posY1});
        else if (posX1 == posY && !environment.edges(posX, posY1).status)
          triples.emplace_back(Triple{posX, posY, posY1});
        }
      }
    }
  else
    {
#if _DEBUG
  Rcpp::Rcout << "\norientationProbability start:\n";
  Rcpp::Rcout << "Connected list\n";
  miic::debug::printConnectedList (environment);
  Rcpp::Rcout << "Edges\n";
  miic::debug::printEdges (environment);
#endif  
    // In temporal mode: 
    //
    // When latent variable discovry is activated, duplicate edges over history 
    // assuming stationarity to increasethe number of possible unshielded triples.
    // Get only unshielded triples X -- Z -- Y having a node at lag0
    // (the past only triples are not interesting as edges forming these past
    // only triples are the result of the duplication and these edges will
    // be removed at the end of this function)
    // In addition of unshielded triples, get also edges not already included
    // in the triples having a lag0 node and which can be oriented using time
    //
    if (environment.latent || environment.latent_orientation)
      tmiic::repeatEdgesOverHistory (environment);
    //
    // Get unshielded triples with at least one lag0 node
    //
    const auto& edge_list = environment.connected_list;
    for (auto iter0 = begin(edge_list); iter0 != end(edge_list); ++iter0)
      {
      int posX = iter0->X, posY = iter0->Y;
      bool edge_has_lag0 =    (posX < environment.n_nodes_not_lagged)
                           || (posY < environment.n_nodes_not_lagged);

      for (auto iter1 = iter0 + 1; iter1 != end(edge_list); ++iter1)
        {
        int posX1 = iter1->X, posY1 = iter1->Y;
        bool edge1_has_lag0 =    (posX1 < environment.n_nodes_not_lagged)
                              || (posY1 < environment.n_nodes_not_lagged);
        if (! (edge_has_lag0 || edge1_has_lag0) )
          {
            
#if _DEBUG
          Rcpp::Rcout << "Couple of edges not lag0 discarded "
                      << environment.nodes[posX].name  << "-" << environment.nodes[posY].name << " "
                      << environment.nodes[posX1].name << "-" << environment.nodes[posY1].name <<"\n";
#endif  
          continue;
          }

        if (posY1 == posX && !environment.edges(posY, posX1).status)
          triples.emplace_back(Triple{posY, posX, posX1});
        else if (posY1 == posY && !environment.edges(posX, posX1).status)
          triples.emplace_back(Triple{posX, posY, posX1});
        if (posX1 == posX && !environment.edges(posY, posY1).status)
          triples.emplace_back(Triple{posY, posX, posY1});
        else if (posX1 == posY && !environment.edges(posX, posY1).status)
          triples.emplace_back(Triple{posX, posY, posY1});
        }
      }
    //
    // Add as fake triple the lagged edges with at least one lag0 node
    // not already in triples
    //
#if _DEBUG
    Rcpp::Rcout << "\norientationProbability: triples with a lag0 node list=\n";
    for (int i = 0; i < triples.size() ; i++)
       Rcpp::Rcout << "Triple " << environment.nodes[triples[i][0]].name
                << "-" << environment.nodes[triples[i][1]].name 
                << "-" << environment.nodes[triples[i][2]].name << "\n";
    Rcpp::Rcout << "\nEvaluate edges not in triples\n";
#endif  
    for (auto iter0 = begin(edge_list); iter0 != end(edge_list); ++iter0)
      {
      int posX = iter0->X, posY = iter0->Y;
      
#if _DEBUG
      Rcpp::Rcout << "Evaluate if edge " << environment.nodes[posX].name 
                  << "-" << environment.nodes[posY].name << " is in triples\n";
#endif  
  
      if ( ! (   (posX < environment.n_nodes_not_lagged)
              || (posY < environment.n_nodes_not_lagged) ) )
        {
#if _DEBUG
      Rcpp::Rcout << "Edge not lag0 => not considered\n";
#endif  
        continue;
        }

      bool is_in_triple = false;
      for (int i = 0; i < triples.size() ; i++)
        if (   ( (triples[i][0] == posX) && (triples[i][1] == posY) )
            || ( (triples[i][0] == posY) && (triples[i][1] == posX) )
            || ( (triples[i][1] == posX) && (triples[i][2] == posY) )
            || ( (triples[i][1] == posY) && (triples[i][2] == posX) ) )
          {
#if _DEBUG
      Rcpp::Rcout << "Edge found in triple " << i << "\n";
#endif  
          is_in_triple = true;
          break;
          }
      if (!is_in_triple)
        {
        triples.emplace_back(Triple{posX, posY, posX});
#if _DEBUG
        Rcpp::Rcout << "Edge " << environment.nodes[posX].name 
                    << "-" << environment.nodes[posY].name << " added as fake triplet\n";
#endif  
        }
     }
#if _DEBUG
  Rcpp::Rcout << "\norientationProbability: triples after tmiic\n";
  for (int i = 0; i < triples.size() ; i++)
    Rcpp::Rcout << "Triple " << environment.nodes[triples[i][0]].name
                << "-" << environment.nodes[triples[i][1]].name 
                << "-" << environment.nodes[triples[i][2]].name << "\n";
#endif  
    }
  
  if (triples.empty())
    return vector<vector<string>>();

  // Compute the 3-point mutual info (N * I'(X;Y;Z|{ui})) for each triple
  vector<double> I3_list(triples.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (size_t i = 0; i < triples.size(); ++i) {
    int X{triples[i][0]}, Z{triples[i][1]}, Y{triples[i][2]};
    if (X == Y) // fake triple, put highest rank in the list
      I3_list[i] = 0;
    else
      {
      const auto& ui_list = environment.edges(X, Y).shared_info->ui_list;
      vector<int> ui_no_z(ui_list);
      ui_no_z.erase(remove(begin(ui_no_z), end(ui_no_z), Z), end(ui_no_z));
      I3_list[i] = getInfo3PointOrScore(
          environment, X, Y, Z, ui_no_z, /* get_info = */ true);
      }
  }

  // Compute the arrowhead probability of each edge endpoint
  vector<ProbaArray> probas_list = getOriProbasList(triples, I3_list,
      environment.tau, environment.nodes, environment.latent,                                            
      environment.latent_orientation, environment.degenerate,
      environment.propagation, environment.half_v_structure);
  // update probas_list for possible inconsistencies
  class ProbaArrayMap : public std::map<std::pair<int, int>, double> {
   public:
    auto insert_or_update(std::pair<int, int> key, double proba) {
      auto iter_key = this->find(key);
      if (iter_key == this->end()) {
        return this->insert({std::move(key), proba});
      } else {
        // Update if proba is more probable
        if (fabs(iter_key->second - 0.5) < fabs(proba - 0.5))
          this->at(key) = proba;
        return std::make_pair(iter_key, false);
      }
    }
  } proba_map;

  for (size_t i = 0; i < triples.size(); i++) {
    const auto& triple = triples[i];
    const auto& probas = probas_list[i];
    proba_map.insert_or_update({triple[1], triple[0]}, probas[0]);  // 1 -> 0
    proba_map.insert_or_update({triple[0], triple[1]}, probas[1]);  // 0 -> 1
    proba_map.insert_or_update({triple[2], triple[1]}, probas[2]);  // 2 -> 1
    proba_map.insert_or_update({triple[1], triple[2]}, probas[3]);  // 1 -> 2
  }
  // Update probabilities
  std::transform(begin(triples), end(triples), begin(probas_list),
      [&proba_map](const auto& triple) {
        return ProbaArray{proba_map.at({triple[1], triple[0]}),
            proba_map.at({triple[0], triple[1]}),
            proba_map.at({triple[2], triple[1]}),
            proba_map.at({triple[1], triple[2]})};
      });
  // Update adj matrix
  for (size_t i = 0; i < triples.size(); i++) {
    const auto& triple = triples[i];
    const auto& probas = probas_list[i];
    updateAdj(environment, triple[0], triple[1], probas[0], probas[1]);
    updateAdj(environment, triple[1], triple[2], probas[2], probas[3]);
  }
  // Write output
   vector<vector<string>> orientations{
       {"source1", "p1", "p2", "target", "p3", "p4", "source2", "NI3", "Error"}};
  for (size_t i = 0; i < triples.size(); i++) {
    const auto& triple = triples[i];
    const auto& probas = probas_list[i];

    int info_int = round(I3_list[i] - 0.5);
    int info_sign = (info_int > 0) - (info_int < 0);

    int proba_sign = 1;
    if (probas[1] > 0.5 && probas[2] > 0.5) proba_sign = -1;

    string error = "0";
    if ((info_int != 0 && info_sign != proba_sign) ||
        (info_int == 0 && proba_sign == -1))
      error = "1";

    using std::to_string;
    orientations.emplace_back(std::initializer_list<string>{
        environment.nodes[triple[0]].name,
        to_string(probas[0]), to_string(probas[1]),
        environment.nodes[triple[1]].name,
        to_string(probas[2]), to_string(probas[3]),
        environment.nodes[triple[2]].name,
        to_string(I3_list[i]), error});
  }
  //
  // In temporal mode, when latent variable discovery is enabled, 
  // we duplicated the edges over the history at the beginning 
  // of the function => we need now to drop these extra edges
  //
  if (   (environment.tau >= 1) 
      && (environment.latent || environment.latent_orientation) )
      tmiic::dropPastEdges (environment);
  
#if _DEBUG
  Rcpp::Rcout << "\norientation end:\n";
  Rcpp::Rcout << "Orientations returned\n";    
  for (int i = 0; i < orientations.size(); i++)
    {
    vector<string> lig = orientations[i];
    for (int j = 0; j < lig.size(); j++) 
      Rcpp::Rcout << lig[j] << " ";    
    Rcpp::Rcout << "\n";    
    }
  Rcpp::Rcout << "Connected list\n";
  miic::debug::printConnectedList (environment);
  Rcpp::Rcout << "Edges\n";
  miic::debug::printEdges (environment);
#endif  
  return orientations;
}

}  // namespace reconstruction
}  // namespace miic
