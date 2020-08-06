#include "reconstruct.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "compute_ens_information.h"
#include "structure.h"
#include "utilities.h"

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
  } else {
    res = computeEnsInformationContinuous(
        environment, NULL, 0, NULL, 0, -1, i, j, environment.cplx, m);
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

bool miic::reconstruction::skeletonInitialization(Environment& environment) {
  for (int i = 0; i < environment.n_samples; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      environment.oneLineMatrix[j * environment.n_samples + i] =
          environment.data_numeric[i][j];
    }
  }

  int threadnum = 0;
  bool interrupt = false;
#ifdef _OPENMP
#pragma omp parallel for shared(interrupt) firstprivate(threadnum) \
    schedule(dynamic)
#endif
  for (int i = 0; i < environment.n_nodes - 1; i++) {
    if (interrupt) {
      continue;  // will continue until out of for loop
    }
#ifdef _OPENMP
    threadnum = omp_get_thread_num();
#endif
    if (checkInterrupt(threadnum == 0)) {
      interrupt = true;
    }
    for (int j = i + 1; j < environment.n_nodes && !interrupt; j++) {
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
  for (int i = 0; i < environment.n_nodes; i++) {
    for (int j = 0; j < environment.n_nodes; j++) {
      environment.edges[i][j].status_init = environment.edges[i][j].status;
    }
  }
  return true;
}
