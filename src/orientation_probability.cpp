#include "orientation_probability.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <map>

#include "compute_ens_information.h"
#include "proba_orientation.h"
#include "structure.h"
#include "utilities.h"

namespace miic {
namespace reconstruction {

using std::string;
using std::vector;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;

double getI3(Environment& environment, Triple t) {
  int posX{t[0]}, posZ{t[1]}, posY{t[2]};

  vector<int> ui_no_z(environment.edges[posX][posY].shared_info->ui_list);
  ui_no_z.erase(remove(begin(ui_no_z), end(ui_no_z), posZ), end(ui_no_z));
  int* ui = ui_no_z.empty() ? NULL : &ui_no_z[0];

  vector<int> z{posZ};
  int* zi = &z[0];

  double* res = NULL;
  double Is = -1;
  double Cs = -1;
  if (!environment.is_continuous[posX] &&
      !environment.is_continuous[posZ] &&
      !environment.is_continuous[posY]) {
    res = computeEnsInformationNew(environment, ui, ui_no_z.size(), zi,
        z.size(), -1, posX, posY, environment.cplx, environment.m);
    Is = res[7];
    Cs = res[8];
    if (environment.degenerate)
      Cs += log(3.0);
    if (environment.is_k23)
      Is += Cs;  // To fit eq(20) and eq(22) in BMC Bioinfo 2016
  } else {
    res = computeEnsInformationContinuous_Orientation(environment, ui,
        ui_no_z.size(), zi, posX, posY, environment.cplx, environment.m);
    Is = res[1];
    Cs = res[2];
    if (environment.degenerate)
      Cs += log(3.0);
    if (environment.is_k23)
      Is -= Cs;  // I(x;y;z|ui) - cplx(x;y;z|ui)
  }
  delete[] res;
  return Is;
}

vector<vector<string>> orientationProbability(Environment& environment) {
  vector<vector<string>> orientations;  // output of orientation table
  // Get all unshielded triples
  vector<Triple> triples;
  const auto& edge_list = environment.connected_list;
  for (auto iter0 = begin(edge_list); iter0 != end(edge_list); ++iter0) {
    int posX = iter0->i, posY = iter0->j;

    for (auto iter1 = iter0 + 1; iter1 != end(edge_list); ++iter1) {
      int posX1 = iter1->i, posY1 = iter1->j;
      if (posY1 == posX && !environment.edges[posY][posX1].status)
        triples.emplace_back(Triple{posY, posX, posX1});
      else if (posY1 == posY && !environment.edges[posX][posX1].status)
        triples.emplace_back(Triple{posX, posY, posX1});
      if (posX1 == posX && !environment.edges[posY][posY1].status)
        triples.emplace_back(Triple{posY, posX, posY1});
      else if (posX1 == posY && !environment.edges[posX][posY1].status)
        triples.emplace_back(Triple{posX, posY, posY1});
    }
  }
  // Compute I3 for each triple
  vector<double> I3_list(triples.size());
  std::transform(begin(triples), end(triples), begin(I3_list),
      [&environment](Triple t) { return getI3(environment, t); });

  int* oneLineMatrixallTpl = new int[triples.size() * 3];
  for (size_t i = 0; i < triples.size(); i++) {
    for (int j = 0; j < 3; j++) {
      oneLineMatrixallTpl[j * triples.size() + i] = triples[i][j];
    }
  }
  double* ptrRetProbValues;
  // Compute the arrowhead probability of each edge endpoint
  int myNbrTpl = triples.size();
  if (myNbrTpl > 0) {
    ptrRetProbValues = getOrientTplLVDegPropag(myNbrTpl, oneLineMatrixallTpl,
        &I3_list[0], environment.latent_orientation, environment.degenerate,
        environment.propagation, environment.half_v_structure);
    // update ptrRetProbValues for possible inconsistencies
    std::map<string, double> probabsMap;
    string s;
    for (size_t i = 0; i < triples.size(); i++) {
      // 0 -> 1
      s = environment.nodes[triples[i][0]].name +
          environment.nodes[triples[i][1]].name;
      double proba = ptrRetProbValues[i + (1 * triples.size())];
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
      s = environment.nodes[triples[i][1]].name +
          environment.nodes[triples[i][0]].name;
      proba = ptrRetProbValues[i + (0 * triples.size())];
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
      s = environment.nodes[triples[i][1]].name +
          environment.nodes[triples[i][2]].name;
      proba = ptrRetProbValues[i + (3 * triples.size())];
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
      s = environment.nodes[triples[i][2]].name +
          environment.nodes[triples[i][1]].name;
      proba = ptrRetProbValues[i + (2 * triples.size())];
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
    for (size_t i = 0; i < triples.size(); i++) {
      s = environment.nodes[triples[i][0]].name +
          environment.nodes[triples[i][1]].name;
      ptrRetProbValues[i + (1 * triples.size())] = probabsMap.find(s)->second;

      s = environment.nodes[triples[i][1]].name +
          environment.nodes[triples[i][0]].name;
      ptrRetProbValues[i + (0 * triples.size())] = probabsMap.find(s)->second;

      s = environment.nodes[triples[i][1]].name +
          environment.nodes[triples[i][2]].name;
      ptrRetProbValues[i + (3 * triples.size())] = probabsMap.find(s)->second;

      s = environment.nodes[triples[i][2]].name +
          environment.nodes[triples[i][1]].name;
      ptrRetProbValues[i + (2 * triples.size())] = probabsMap.find(s)->second;
    }

    orientations.emplace_back(std::initializer_list<string>{"source1", "p1",
        "p2", "target", "p3", "p4", "source2", "NI3", "Error"});
    for (size_t i = 0; i < triples.size(); i++) {
      string error = "0";
      int info_sign = round(I3_list[i] - 0.5);
      int proba_sign = 0;

      if (ptrRetProbValues[i + (1 * triples.size())] > 0.5 &&
          ptrRetProbValues[i + (2 * triples.size())] > 0.5)
        proba_sign = -1;
      else
        proba_sign = 1;

      if ((sign(info_sign) != proba_sign && info_sign != 0) ||
          (info_sign == 0 && proba_sign == -1))
        error = "1";

      using std::to_string;
      orientations.emplace_back(std::initializer_list<string>{
          environment.nodes[triples[i][0]].name,
          to_string(ptrRetProbValues[i + (0 * triples.size())]),
          to_string(ptrRetProbValues[i + (1 * triples.size())]),
          environment.nodes[triples[i][1]].name,
          to_string(ptrRetProbValues[i + (2 * triples.size())]),
          to_string(ptrRetProbValues[i + (3 * triples.size())]),
          environment.nodes[triples[i][2]].name,
          to_string(I3_list[i]),
          error
      });
    }

    // Update adj matrix
    vector<double> myProba;
    for (size_t i = 0; i < triples.size(); i++) {
      myProba.clear();
      myProba.push_back(ptrRetProbValues[i + (0 * triples.size())]);
      myProba.push_back(ptrRetProbValues[i + (1 * triples.size())]);
      // If there is AT LEAST ONE arrowhead
      if ((*max_element(myProba.begin(), myProba.end()) - 0.5) > 0) {
        // If there is ONLY ONE arrowhead
        if ((*min_element(myProba.begin(), myProba.end()) - 0.5) <= 0) {
          double a = myProba[0] - *max_element(myProba.begin(), myProba.end());
          // If p1 is max: n1 <-- n2
          if (a == 0) {
            environment.edges[triples[i][0]][triples[i][1]].status = -2;
            environment.edges[triples[i][1]][triples[i][0]].status = 2;
          } else {
            environment.edges[triples[i][0]][triples[i][1]].status = 2;
            environment.edges[triples[i][1]][triples[i][0]].status = -2;
          }
        } else {
          environment.edges[triples[i][0]][triples[i][1]].status = 6;
          environment.edges[triples[i][1]][triples[i][0]].status = 6;
        }
      }
      myProba.clear();
      myProba.push_back(ptrRetProbValues[i + (2 * triples.size())]);
      myProba.push_back(ptrRetProbValues[i + (3 * triples.size())]);
      // If there is AT LEAST ONE arrowhead
      if ((*max_element(myProba.begin(), myProba.end()) - 0.5) > 0) {
        // If there is ONLY ONE arrowhead
        if ((*min_element(myProba.begin(), myProba.end()) - 0.5) <= 0) {
          double a = myProba[0] - *max_element(myProba.begin(), myProba.end());
          // If p1 is max: n1 <-- n2
          if (a == 0) {
            environment.edges[triples[i][1]][triples[i][2]].status = -2;
            environment.edges[triples[i][2]][triples[i][1]].status = 2;
          } else {
            environment.edges[triples[i][1]][triples[i][2]].status = 2;
            environment.edges[triples[i][2]][triples[i][1]].status = -2;
          }
        } else {
          environment.edges[triples[i][1]][triples[i][2]].status = 6;
          environment.edges[triples[i][2]][triples[i][1]].status = 6;
        }
      }
    }
  }
  delete[] oneLineMatrixallTpl;
  if (myNbrTpl > 0) delete[] ptrRetProbValues;

  return orientations;
}

}  // namespace reconstruction
}  // namespace miic
