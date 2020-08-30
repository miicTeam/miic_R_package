#include "skeleton.h"

#include <Rcpp.h>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "compute_ens_information.h"
#include "structure.h"
#include "utilities.h"

using Rcpp::Rcout;
using Rcpp::Rcerr;
using std::endl;
using std::vector;
using namespace miic::computation;
using namespace miic::reconstruction;
using namespace miic::structure;
using namespace miic::utility;

template <bool latent = false>
void searchAndSetZi(
    Environment& environment, const int posX, const int posY) {
  // Search candidate nodes for the separation set of X and Y
  int numZiPos = 0;
  for (int c = 0; c < environment.n_nodes; c++) {
    if (c == posX || c == posY) continue;
    if (!latent && !environment.edges[posX][c].status_prev &&
        !environment.edges[posY][c].status_prev)
      continue;
    environment.edges[posX][posY].shared_info->zi_list.push_back(c);
    numZiPos++;
  }

  if (environment.verbose)
    Rcout << "The number of neighbours is: " << numZiPos << endl;
}

void searchAndSetZi(
    Environment& environment, const int posX, const int posY) {
  if (environment.latent)
    return searchAndSetZi<true>(environment, posX, posY);
  else
    return searchAndSetZi<false>(environment, posX, posY);
}

namespace miic {
namespace reconstruction {

// Initialize the edges of the network
int initializeEdge(Environment& environment, int i, int j) {
  // Compute the mutual information and the corresponding CPLX
  double* res = NULL;
  if (!environment.is_continuous[i] &&
      !environment.is_continuous[j]) {
    res = computeEnsInformationNew(environment, i, j, std::vector<int>(),
        std::vector<int>(), environment.cplx);
    environment.edges[i][j].shared_info->Ixy_ui = res[1];
    environment.edges[i][j].shared_info->cplx = res[2];
    environment.edges[i][j].shared_info->Nxy_ui = res[0];
  } else {
    res = computeEnsInformationContinuous(environment, i, j, std::vector<int>(),
        std::vector<int>(), environment.cplx);
    environment.edges[i][j].shared_info->Ixy_ui = res[1];
    environment.edges[i][j].shared_info->cplx = res[2];
    environment.edges[i][j].shared_info->Nxy_ui = res[0];
  }
  delete[] res;

  double myTest = 0;
  environment.edges[i][j].shared_info->Ixy =
      environment.edges[i][j].shared_info->Ixy_ui;
  environment.edges[i][j].shared_info->cplx_no_u =
      environment.edges[i][j].shared_info->cplx;
  environment.edges[i][j].shared_info->Nxy =
      environment.edges[i][j].shared_info->Nxy_ui;

  if (environment.no_init_eta)
    myTest = environment.edges[i][j].shared_info->Ixy_ui -
             environment.edges[i][j].shared_info->cplx;
  else
    myTest = environment.edges[i][j].shared_info->Ixy_ui -
             environment.edges[i][j].shared_info->cplx - environment.log_eta;

  if (myTest <= 0) {
    // Unconditional independence
    environment.edges[i][j].shared_info->connected = 0;
    environment.edges[i][j].status = 0;
    environment.edges[j][i].status = 0;
  } else {
    environment.edges[i][j].shared_info->connected = 1;
    environment.edges[i][j].status = 1;
    environment.edges[j][i].status = 1;
  }

  return environment.edges[i][j].status;
}

bool skeletonInitialization(Environment& environment) {
  for (int i = 0; i < environment.n_samples; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      environment.oneLineMatrix[j * environment.n_samples + i] =
          environment.data_numeric[i][j];
    }
  }

  bool interrupt = false;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < environment.n_nodes - 1; i++) {
    if (interrupt) {
      continue;  // will continue until out of for loop
    }
    int threadnum = 0;
#ifdef _OPENMP
    threadnum = omp_get_thread_num();
#endif
    if (threadnum == 0) {
      if (checkInterrupt()) interrupt = true;
    }
    for (int j = i + 1; j < environment.n_nodes && !interrupt; j++) {
      environment.edges[i][j].shared_info = std::make_shared<EdgeSharedInfo>();
      environment.edges[j][i].shared_info = environment.edges[i][j].shared_info;
      if (environment.edges[i][j].status) {
        if (initializeEdge(environment, i, j) == 1) {
#ifdef _OPENMP
#pragma omp critical
#endif
          environment.numSearchMore++;
        }
      }
    }
  }
  if (interrupt) return false;
  for (int i = 0; i < environment.n_nodes; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      environment.edges[i][j].status_init = environment.edges[i][j].status;
    }
  }
  return true;
}

bool firstStepIteration(Environment& environment, BiconnectedComponent& bcc) {
  // During first step iteration, search for U contributors is not
  // parallelizable see flag "parallelizable" in
  // computeEnsInformationContinuous() in computeEnsInformation.cpp
  environment.first_iter_done = false;

  environment.connected_list.clear();
  environment.unsettled_list.clear();

  // create and fill the unsettled_list struct, that keep track of i and j
  // positions of searchMore Edges
  environment.numSearchMore = 0;
  environment.numNoMore = 0;
  for (int i = 0; i < environment.n_nodes - 1; i++) {
    for (int j = i + 1; j < environment.n_nodes; j++) {
      // Do dot consider edges removed with unconditional independence
      if (!environment.edges[i][j].status) continue;
      environment.edges[i][j].shared_info->reset();
      environment.unsettled_list.emplace_back(i, j, environment.edges[i][j]);
      environment.numSearchMore++;
    }
  }

  if (environment.numSearchMore > 0) {
    if (environment.verbose)
      Rcout << "\n# -> searchMore edges, to get zi and noMore...\n";

    for (int i = 0; i < environment.numSearchMore; i++) {
      int posX = environment.unsettled_list[i].i;
      int posY = environment.unsettled_list[i].j;
      if (environment.verbose) {
        Rcout << "\n# --------------------\n# ----> EDGE: "
             << environment.nodes[posX].name << "--"
             << environment.nodes[posY].name << "\n# --------------------";
      }
      if (environment.consistent > 0)
        bcc.setCandidateZ(posX, posY);
      else
        searchAndSetZi(environment, posX, posY);
    }

    if (environment.verbose) {
      Rcout << "SEARCH OF BEST Z: ";
    }

    bool interrupt = false;
    int progress_percentile = -1;
    int n_jobs_done{0};
    auto loop_start_time = getLapStartTime();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < environment.numSearchMore; i++) {
      if (interrupt) {
        continue;  // will continue until out of for loop
      }
      int threadnum = 0;
#ifdef _OPENMP
      threadnum = omp_get_thread_num();
#endif
      if (threadnum == 0) {
        if (checkInterrupt()) interrupt = true;
      }
      const auto& edgeid = environment.unsettled_list[i];
      if (edgeid.getEdge().shared_info->zi_list.size() > 0)
        SearchForNewContributingNodeAndItsRank(environment, edgeid.i, edgeid.j);

#ifdef _OPENMP
#pragma omp atomic
      ++n_jobs_done;
#endif
      if (threadnum == 0)
        printProgress(
            static_cast<double>(n_jobs_done) / environment.numSearchMore,
            loop_start_time, progress_percentile);
    }
    // Print finished progress bar
    printProgress(1.0, loop_start_time, progress_percentile);
    Rcerr << '\n';

    if (interrupt) return false;

    for (int i = 0; i < environment.numSearchMore; i++) {
      int posX = environment.unsettled_list[i].i;
      int posY = environment.unsettled_list[i].j;
      if (environment.edges[posX][posY].shared_info->z_name_idx != -1) {
        if (environment.verbose) {
          Rcout << "## ------!!--> Update the edge element in 'searchMore': "
               << environment
                      .nodes[environment.edges[posX][posY].shared_info->zi_list
                                 [environment.edges[posX][posY]
                                         .shared_info->z_name_idx]]
                      .name
               << " is a good zi candidate\n";
        }
      } else {
        if (environment.verbose) {
          Rcout << "## ------!!--> Remove the edge element from searchMore.\n## "
                  "------!!--> Add edge to 'noMore' (no good zi candidate)\n";
        }
        // Move this edge element to "noMore"
        environment.connected_list.push_back(environment.unsettled_list[i]);
        environment.numNoMore++;
        environment.unsettled_list.erase(
            environment.unsettled_list.begin() + i);
        environment.numSearchMore--;
        i--;
        // Update the status
        environment.edges[posX][posY].shared_info->connected = 1;
      }
      if (environment.verbose) Rcout << "\n";
    }

    if (checkInterrupt()) {
      return false;
    }
    // sort the ranks
    std::sort(
        environment.unsettled_list.begin(), environment.unsettled_list.end());
  }
  environment.first_iter_done = true;
  return true;
}

bool skeletonIteration(Environment& environment) {
  int iIteration_count = 0;
  int max = 0;

  if (environment.verbose)
    Rcout << "Number of numSearchMore: " << environment.numSearchMore << endl;

  auto loop_start_time = getLapStartTime();
  int start_numSearchMore = environment.numSearchMore;

  int progress_percentile = -1;

  while (environment.numSearchMore > 0) {
    if (checkInterrupt()) {
      return false;
    }
    iIteration_count++;
    if (environment.verbose) {
      Rcout << "\n# Iteration " << iIteration_count << "\n";
    }
    // Get the first edge
    int posX = environment.unsettled_list[max].i;
    int posY = environment.unsettled_list[max].j;

    if (environment.verbose)
      Rcout << "Pos x : " << posX << " , pos y: " << posY << endl;

    auto topEdgeElt = environment.edges[posX][posY].shared_info;

    if (environment.verbose) Rcout << "# Before adding new zi to {ui}: ";

    // Keep the previous z.name for this edge
    std::string accepted_z_name =
        environment.nodes[topEdgeElt->z_name_idx].name;

    // Reinit ui.vect, z.name, zi.vect, z.name.idx
    if (environment.verbose) {
      Rcout << "# DO: Add new zi to {ui}: " << topEdgeElt->z_name_idx << endl;
    }
    // move top z_name_idx from zi_vect to ui_vect
    topEdgeElt->ui_list.push_back(topEdgeElt->zi_list[topEdgeElt->z_name_idx]);
    topEdgeElt->zi_list.erase(
        topEdgeElt->zi_list.begin() + topEdgeElt->z_name_idx);
    topEdgeElt->z_name_idx = -1;

    double* v = NULL;
    if (!environment.is_continuous[posX] && !environment.is_continuous[posY] &&
        std::all_of(topEdgeElt->ui_list.cbegin(), topEdgeElt->ui_list.cend(),
            [&environment](int i) { return !environment.is_continuous[i]; })) {
      v = computeEnsInformationNew(environment, posX, posY,
          environment.edges[posX][posY].shared_info->ui_list, vector<int>(),
          environment.cplx);

      topEdgeElt->Ixy_ui = v[1];
      topEdgeElt->Nxy_ui = v[0];
      topEdgeElt->cplx = v[2];
    } else {
      v = computeEnsInformationContinuous(environment, posX, posY,
          environment.edges[posX][posY].shared_info->ui_list, vector<int>(),
          environment.cplx);
      topEdgeElt->Nxy_ui = v[0];
      topEdgeElt->Ixy_ui = v[1];
      topEdgeElt->cplx = v[2];
    }
    delete[] v;
    double topEdgeElt_kxy_ui = topEdgeElt->cplx;

    if (environment.degenerate)
      topEdgeElt_kxy_ui =
          topEdgeElt->cplx + (topEdgeElt->ui_list.size() * log(3.0));

    int nRemainingEdges = environment.numSearchMore + environment.numNoMore;

    if (environment.verbose) {
      Rcout << "# --> nbrEdges L = " << nRemainingEdges << "\n";
      Rcout << "# --> nbrProp P = " << environment.n_nodes << "\n\n";
      Rcout << "topEdgeElt->Ixy_ui " << topEdgeElt->Ixy_ui << "\n";
      Rcout << "topEdgeElt_kxy_ui " << topEdgeElt_kxy_ui << "\n";
      Rcout << "environment.log_eta " << environment.log_eta << "\n";
      Rcout << "IsPhantom? "
           << (topEdgeElt->Ixy_ui - topEdgeElt_kxy_ui - environment.log_eta <= 0)
           << endl;
    }
    if (topEdgeElt->Ixy_ui - topEdgeElt_kxy_ui - environment.log_eta <= 0) {
      // Conditional independence found, remove edge
      if (environment.verbose) {
        Rcout << "# PHANTOM" << environment.nodes[posX].name << ","
             << environment.nodes[posY].name << "\n";
      }
      environment.unsettled_list.erase(
          environment.unsettled_list.begin() + max);
      environment.numSearchMore--;
      // Set the connection to 0 on the adj matrix
      environment.edges[posX][posY].status = 0;
      environment.edges[posY][posX].status = 0;
      // Save the phantom status
      topEdgeElt->connected = 0;
    } else {
      // Reinit Rxyz_ui
      topEdgeElt->Rxyz_ui = environment.thresPc;

      if (environment.verbose) {
        Rcout << "# Do SearchForNewContributingNodeAndItsRank\n";
      }

      if (topEdgeElt->zi_list.size() > 0) {
        SearchForNewContributingNodeAndItsRank(environment, posX, posY);
      }

      if (environment.verbose) {
        if (environment.edges[posX][posY].shared_info->z_name_idx == -1)
          Rcout << "# See topEdgeElt[['z.name']]: NA\n";
        else
          Rcout << "# See topEdgeElt[['z.name']]: "
               << environment.nodes[topEdgeElt->zi_list[topEdgeElt->z_name_idx]]
                      .name
               << "\n";
      }
      //// Update the information about the edge
      if (topEdgeElt->z_name_idx != -1) {
        if (environment.verbose) {
          Rcout << "# Do update myAllEdges$searchMore\n";
        }
        // myGv$allEdges[["searchMore"]][[topEdgeElt[["key"]]]] = topEdgeElt

      } else {
        if (environment.verbose) {
          Rcout << "# Do update myAllEdges$noMore\n";
        }
        // Move this edge from the list searchMore to noMore
        environment.connected_list.push_back(environment.unsettled_list[max]);
        environment.numNoMore++;
        environment.unsettled_list.erase(
            environment.unsettled_list.begin() + max);
        environment.numSearchMore--;
        // Update the status of the edge
        topEdgeElt->connected = 1;
      }
    }

    // Sort all pairs xy with a contributing node z in decreasing order of
    // their ranks, R(xy;z| )
    if (environment.verbose) {
      Rcout << "# Do Sort all pairs by Rxyz_ui\n";
    }

    max = 0;
    for (int i = 0; i < environment.numSearchMore; i++) {
      if (environment
              .edges[environment.unsettled_list[i].i]
                    [environment.unsettled_list[i].j]
              .shared_info->Rxyz_ui >
          environment
              .edges[environment.unsettled_list[max].i]
                    [environment.unsettled_list[max].j]
              .shared_info->Rxyz_ui)
        max = i;
    }
    printProgress(1.0 * (start_numSearchMore - environment.numSearchMore) /
                      (start_numSearchMore),
        loop_start_time, progress_percentile);
  }
  Rcerr << "\n";
  std::sort(
      environment.connected_list.begin(), environment.connected_list.end());
  return (true);
}

}  // namespace reconstruction
}  // namespace miic
