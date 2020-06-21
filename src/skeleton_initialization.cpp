#include <math.h>
#include <unistd.h>

#include <algorithm>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include "reconstruct.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "compute_ens_information.h"
#include "structure.h"
#include "utilities.h"

using uint = unsigned int;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;

// Initialize the edges of the network
int initEdgeElt(Environment& environment, int i, int j, MemorySpace& m) {
  // Compute the mutual information and the corresponding CPLX
  double* res = NULL;
  if (!environment.is_continuous[i] &&
      !environment.is_continuous[j]) {
    res = computeEnsInformationNew(
        environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx, m);
    environment.edges[i][j].shared_info->Ixy_ui = res[1];
    environment.edges[i][j].shared_info->cplx = res[2];
    environment.edges[i][j].shared_info->Nxy_ui = res[0];
    free(res);
  } else {
    res = computeEnsInformationContinuous(
        environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx, m);
    environment.edges[i][j].shared_info->Ixy_ui = res[1];
    environment.edges[i][j].shared_info->cplx = res[2];
    environment.edges[i][j].shared_info->Nxy_ui = res[0];
    delete[] res;
  }

  if (environment.verbose) {
    std::cout << "# --> Ixy_ui = "
              << environment.edges[i][j].shared_info->Ixy_ui /
                     environment.edges[i][j].shared_info->Nxy_ui
              << "[Ixy_ui*Nxy_ui ="
              << environment.edges[i][j].shared_info->Ixy_ui << "]\n"
              << "# --> Cplx = " << environment.edges[i][j].shared_info->cplx
              << "\n"
              << "# --> Nxy_ui = "
              << environment.edges[i][j].shared_info->Nxy_ui << "\n"
              << "# --> nbrProp P = " << environment.n_nodes << "\n";
  }

  double myTest = 0;
  std::string category;
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
    category = "phantom";
  } else {
    environment.edges[i][j].shared_info->connected = 1;
    environment.edges[i][j].status = 1;
    environment.edges[j][i].status = 1;
    category = "searchMore";
  }

  if (environment.verbose)
    std::cout << "# --> Category = " << category << "\n";

  return environment.edges[i][j].status;
}

bool miic::reconstruction::skeletonInitialization(Environment& environment) {
  environment.oneLineMatrix =
      new int[environment.n_samples * environment.n_nodes];
  for (uint i = 0; i < environment.n_samples; i++) {
    for (uint j = 0; j < environment.n_nodes; j++) {
      environment.oneLineMatrix[j * environment.n_samples + i] =
          environment.data_numeric[i][j];
    }
  }

  int threadnum = 0;
  std::cout << "Computing pairwise independencies...";
  fflush(stdout);

  bool interrupt = false;

#ifdef _OPENMP
#pragma omp parallel for shared(interrupt) firstprivate(threadnum) \
    schedule(dynamic)
#endif
  for (uint i = 0; i < environment.n_nodes - 1; i++) {
    if (interrupt) {
      continue;  // will continue until out of for loop
    }
#ifdef _OPENMP
    threadnum = omp_get_thread_num();
#endif
    if (checkInterrupt(threadnum == 0)) {
      interrupt = true;
    }
    for (uint j = i + 1; j < environment.n_nodes && !interrupt; j++) {
      if (environment.verbose) {
        std::cout << "\n# Edge " << environment.nodes[i].name << ","
                  << environment.nodes[j].name << "\n";
      }
      environment.edges[i][j].shared_info = std::make_shared<EdgeSharedInfo>();
      environment.edges[j][i].shared_info = environment.edges[i][j].shared_info;
      if (environment.edges[i][j].status) {
        if (initEdgeElt(
                environment, i, j, environment.memoryThreads[threadnum]) == 1) {
#ifdef _OPENMP
#pragma omp critical
#endif
          environment.numSearchMore++;
        }
      }
    }
  }
  if (interrupt) return false;
  std::cout << " done.\n";
  for (uint i = 0; i < environment.n_nodes; i++) {
    for (uint j = 0; j < environment.n_nodes; j++) {
      environment.edges[i][j].status_init = environment.edges[i][j].status;
    }
  }
  return true;
}
