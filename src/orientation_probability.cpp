#include "orientation_probability.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <map>
#include <regex>

#include "compute_ens_information.h"
#include "proba_orientation.h"
#include "tmiic.h"

namespace miic {
namespace reconstruction {

using std::string;
using std::vector;
using std::fabs;
using namespace miic::computation;
using namespace miic::structure;

namespace {

bool acceptProba(double proba, double ori_proba_ratio) {
  if (proba <= 0.5) return false;
  return (1 - proba) / proba < ori_proba_ratio;
}

double getI3(Environment& environment, const Triple& t) {
  int posX{t[0]}, posZ{t[1]}, posY{t[2]};

  vector<int> ui_no_z(environment.edges[posX][posY].shared_info->ui_list);
  ui_no_z.erase(remove(begin(ui_no_z), end(ui_no_z), posZ), end(ui_no_z));
  int* ui = ui_no_z.empty() ? NULL : &ui_no_z[0];

  vector<int> z{posZ};
  int* zi = &z[0];

  double* res = NULL;
  double Ixyz_ui = -1;
  double cplx = -1;
  if (!environment.is_continuous[posX] && !environment.is_continuous[posZ] &&
      !environment.is_continuous[posY]) {
    res = computeEnsInformationNew(environment, ui, ui_no_z.size(), zi,
        z.size(), -1, posX, posY, environment.cplx, environment.m);
    Ixyz_ui = res[7];
    cplx = res[8];
    if (environment.degenerate) cplx += log(3.0);
    // To fit eq(20) and eq(22) in BMC Bioinfo 2016
    if (environment.is_k23) Ixyz_ui += cplx;
  } else {
    res = computeEnsInformationContinuous_Orientation(environment, ui,
        ui_no_z.size(), zi, posX, posY, environment.cplx, environment.m);
    Ixyz_ui = res[1];
    cplx = res[2];
    if (environment.degenerate) cplx += log(3.0);
    if (environment.is_k23) Ixyz_ui -= cplx;
  }
  delete[] res;
  return Ixyz_ui;
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
      env.edges[x][y].status = -2;
      env.edges[y][x].status = 2;
    } else if (acceptProba(x2y, env.ori_proba_ratio)) {
      env.edges[x][y].status = 2;
      env.edges[y][x].status = -2;
    }
  } else if (acceptProba(x2y, env.ori_proba_ratio) &&
             acceptProba(y2x, env.ori_proba_ratio)) {
    env.edges[x][y].status = 6;
    env.edges[y][x].status = 6;
  }
}

}  // anonymous namespace

vector<vector<string>> orientationProbability(Environment& environment) {
  // Get all unshielded triples X -- Z -- Y
  vector<Triple> triples;
  const auto& edge_list = environment.connected_list;
  
  int nodes_cnt_not_lagged = 0; // for tmiic, number of nodes (not lagged)

  if (environment.tau > 0) 
    {
    // For temporal graph, repeat edges found over history if latent variable 
    // discovery is activated (TODO add why ...)
    //
    if (environment.latent || environment.latent_orientation)
      tmiic::repeatEdgesOverHistory (environment);
    //
    // We count how much non lagged nodes we have 
    // => we look for the first lag1 node
    //
    std::regex lag_expr(".*lag");
    while (nodes_cnt_not_lagged < environment.n_nodes)
      {
      string node_name = environment.nodes[nodes_cnt_not_lagged].name;
      int node_lag = std::stoi (std::regex_replace(node_name, lag_expr, "") );
      if (node_lag > 0)
        break;
      nodes_cnt_not_lagged++;
      }
    }
  
  for (auto iter0 = begin(edge_list); iter0 != end(edge_list); ++iter0) {
    int posX = iter0->i, posY = iter0->j;

    for (auto iter1 = iter0 + 1; iter1 != end(edge_list); ++iter1) {
      int posX1 = iter1->i, posY1 = iter1->j;
      //
      // In temporal mode, we are only interested to orient
      // triplets having at least one node contemporaneaous (lag0).
      // (Triplets having only past nodes are induced by the call to 
      // repeatEdgesOverHistory when latent discovery is on)
      //
      bool add_Y_X_X1 = true;
      bool add_X_Y_X1 = true;
      bool add_Y_X_Y1 = true;
      bool add_X_Y_Y1 = true;
      if (  (environment.tau > 0) 
         && (environment.latent || environment.latent_orientation) )
        {
        int lag_nodeX = posX / nodes_cnt_not_lagged;
        int lag_nodeY = posY / nodes_cnt_not_lagged;
        int lag_nodeX1 = posX1 / nodes_cnt_not_lagged;
        int lag_nodeY1 = posY1 / nodes_cnt_not_lagged;
        add_Y_X_X1 = (lag_nodeY == 0) || (lag_nodeX == 0) || (lag_nodeX1 == 0);
        add_X_Y_X1 = (lag_nodeX == 0) || (lag_nodeY == 0) || (lag_nodeX1 == 0);
        add_Y_X_Y1 = (lag_nodeY == 0) || (lag_nodeX == 0) || (lag_nodeY1 == 0);
        add_X_Y_Y1 = (lag_nodeX == 0) || (lag_nodeY == 0) || (lag_nodeY1 == 0);
        }

      if (posY1 == posX && !environment.edges[posY][posX1].status && add_Y_X_X1)
        triples.emplace_back(Triple{posY, posX, posX1});
      else if (posY1 == posY && !environment.edges[posX][posX1].status && add_X_Y_X1)
        triples.emplace_back(Triple{posX, posY, posX1});
      if (posX1 == posX && !environment.edges[posY][posY1].status && add_Y_X_Y1)
        triples.emplace_back(Triple{posY, posX, posY1});
      else if (posX1 == posY && !environment.edges[posX][posY1].status && add_X_Y_Y1)
        triples.emplace_back(Triple{posX, posY, posY1});
    }
  }
  //
  // In temporal mode, the oriention of temporal edges can be determined
  // using the time. So even if a temporal edge is not part of an open
  // triple, we will add it for orientation.
  // 
  // As all the orientation code has been designed for triplets, 
  // we use a trick here by creating a fake triplet with:
  // node1 of edge - node2 of edge - node1 of edge 
  //
  if (environment.tau > 0) 
    {
    //
    // Go over all edges computed
    //
    for (auto edge_iter = begin(edge_list); edge_iter != end(edge_list); ++edge_iter) 
      {
      int egde_node0 = edge_iter->i;
      int egde_node1 = edge_iter->j;
      //
      // If edge is not temporal, we don't add it
      //
      int lag_node0 = egde_node0 / nodes_cnt_not_lagged;
      int lag_node1 = egde_node1 / nodes_cnt_not_lagged;
      if (lag_node0 == lag_node1)
        continue;
      //
      // If edge has no contemporaneous node, we don't add it
      // (this can only happen whe latent discovery is activated
      // because we duplicated edges over history and we are not 
      // interested to orient duplicated past edges)
      //
      if (lag_node0 > 0 && lag_node1 > 0)
        continue;
      //
      // The computed edge is temporal and we want to have it oriented.
      // We now check if the edge is already in the open triplets list
      //
      bool is_edge_in_list = false;
      for (size_t i = 0; i < triples.size(); i++) 
        {
        const auto& triple = triples[i];
    
        if (  (egde_node0 == triple[0] && egde_node1 == triple[1])
           || (egde_node0 == triple[1] && egde_node1 == triple[0])
           || (egde_node0 == triple[1] && egde_node1 == triple[2])
           || (egde_node0 == triple[2] && egde_node1 == triple[1]) )
          {
          is_edge_in_list = true;
          break;
          }
        }
      //
      // If the edge is already in the triplets list, nothing to do
      //
      if (is_edge_in_list)
        continue;
      //
      // If the edge is not in the triplets list, add a fake triplet
      // with node1 of edge - node2 of edge - node1 of edge
      // (and we add 1 to the node indexes like getSStructure does)
      //
      triples.emplace_back (Triple{egde_node0, egde_node1, egde_node0});
      // TODO See why we were setting I3=0 allI3.push_back (0);
      // and how not compute I3 in the next lines below
#if _DEBUG
      Rcout << "Fake triplet: Edges " << environment.nodes[egde_node0].name
            << "-" << environment.nodes[egde_node1].name 
            << "-" << environment.nodes[egde_node0].name << " with I3=0 added\n";
#endif
      }
    }
  
  if (triples.empty())
    return vector<vector<string>>();

  // Compute the 3-point mutual info (N * I'(X;Y;Z|{ui})) for each triple
  vector<double> I3_list(triples.size());
  std::transform(begin(triples), end(triples), begin(I3_list),
      [&environment](const auto& t) { return getI3(environment, t); });

  // Compute the arrowhead probability of each edge endpoint
  vector<ProbaArray> probas_list = getOriProbasList(environment, triples, I3_list,
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
  //
  // In temporal mode, when latent variable discovery is enabled, 
  // we duplicated the edges over the history at the beginning 
  // of the function => we need now to drop these extra edges
  //
  if (   environment.tau > 0 
      && (environment.latent || environment.latent_orientation) )
      tmiic::dropPastEdges (environment);
  
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

  if (environment.tau > 0)
    {
    // Check if orientation found by miic is aligned with time
    //
    for (auto edge_iter = begin(edge_list); edge_iter != end(edge_list); ++edge_iter) 
      {
      int egde_node0 = edge_iter->i;
      int egde_node1 = edge_iter->j;
      int orient_from_miic = environment.edges[egde_node0][egde_node1].status;

      int lag_node0 = egde_node0 / nodes_cnt_not_lagged;
      int lag_node1 = egde_node1 / nodes_cnt_not_lagged;
      if (lag_node0 == lag_node1)
        continue;
      
      int orient_from_time;
      if (lag_node0 < lag_node1)
        orient_from_time = -2;
      else
        orient_from_time = 2;
      
      if (orient_from_miic != orient_from_time)
        Rcpp::Rcout << "Warning: Edge " << environment.nodes[egde_node0].name
                    << "-" << environment.nodes[egde_node1].name
                    << ", miic orientation=" << orient_from_miic
                    << " differs from time orientation=" << orient_from_time << "\n";
      }
    }
  
#if _DEBUG
  Rcout << "\norientationProbability end:\n";
  printEdges (environment);
#endif
    }
  return orientations;
}

}  // namespace reconstruction
}  // namespace miic
