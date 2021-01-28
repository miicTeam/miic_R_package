#include "skeleton.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>

#include <algorithm>  // std::sort, std::max_element
#include <vector>

#include "get_information.h"
#include "structure.h"
#include "utilities.h"

namespace miic {
namespace reconstruction {

using Rcpp::Rcerr;
using Rcpp::Rcout;
using namespace miic::computation;
using namespace miic::structure;
using namespace miic::utility;

constexpr double kMinCheckInterval = 1;

// Initialize the edges of the network
int initializeEdge(Environment& environment, int X, int Y) {
  // Compute the mutual information and the corresponding CPLX
  auto info = environment.edges(X, Y).shared_info;
  auto xy = getCondMutualInfo(X, Y, std::vector<int>(),
      environment.data_numeric, environment.data_numeric_idx, environment);
  info->Nxy = xy.n_samples;
  info->Ixy = xy.I;
  info->kxy = xy.k;

  info->Nxy_ui = info->Nxy;
  info->Ixy_ui = info->Ixy;
  info->kxy_ui = info->kxy;

  double shifted_mi = info->Ixy - info->kxy;
  if (!environment.no_init_eta) shifted_mi -= environment.log_eta;

  if (shifted_mi <= 0) {
    // Unconditional independence
    info->connected = 0;
    environment.edges(X, Y).status = 0;
    environment.edges(Y, X).status = 0;
    environment.edges(X, Y).status_init = 0;
    environment.edges(Y, X).status_init = 0;
    environment.edges(X, Y).proba_head = -1;
    environment.edges(Y, X).proba_head = -1;
  } else {
    info->connected = 1;
    environment.edges(X, Y).status = 1;
    environment.edges(Y, X).status = 1;
    environment.edges(X, Y).status_init = 1;
    environment.edges(Y, X).status_init = 1;
  }

  return environment.edges(X, Y).status;
}

bool initializeSkeleton(Environment& environment) {
  auto& edges = environment.edges;
  bool interrupt{false};
  auto t_last_check = getLapStartTime();
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int threadnum{0};
#ifdef _OPENMP
    threadnum = omp_get_thread_num();
#pragma omp for schedule(dynamic)
#endif
    for (int i = 0; i < environment.n_nodes - 1; ++i) {
      if (interrupt) continue;  // Will continue until out of for loop
      // Check interrupt on main thread, if enough time passed since last check
      if (threadnum == 0 && getLapInterval(t_last_check) > kMinCheckInterval) {
        t_last_check = getLapStartTime();
        if (checkInterrupt()) {
          interrupt = true;
          continue;
        }
      }

      for (int j = i + 1; j < environment.n_nodes; ++j) {
        edges(i, j).shared_info = std::make_shared<EdgeSharedInfo>();
        edges(j, i).shared_info = edges(i, j).shared_info;

        if (edges(i, j).status) initializeEdge(environment, i, j);
      }
    }
  }  // omp parallel
  return interrupt ? false : true;
}

bool setBestContributingNode(
    Environment& environment, BiconnectedComponent& bcc) {
  auto& connected_list = environment.connected_list;
  auto& unsettled_list = environment.unsettled_list;
  auto& edges = environment.edges;
  connected_list.clear();
  unsettled_list.clear();
  for (int i = 0; i < environment.n_nodes - 1; i++) {
    for (int j = i + 1; j < environment.n_nodes; j++) {
      // Do dot consider edges removed with unconditional independence
      if (!edges(i, j).status) continue;
      edges(i, j).shared_info->reset();
      unsettled_list.emplace_back(i, j, edges(i, j));
    }
  }
  int n_jobs_total = unsettled_list.size();
  if (unsettled_list.empty()) return true;

  int progress_percentile{-1}, n_jobs_done{0};
  bool interrupt = false;
  auto loop_start_time = getLapStartTime();
  auto t_last_check = loop_start_time;
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int threadnum{0};
#ifdef _OPENMP
    threadnum = omp_get_thread_num();
#pragma omp for schedule(dynamic)
#endif
    for (size_t i = 0; i < environment.unsettled_list.size(); ++i) {
      if (interrupt) continue;  // Will continue until out of for loop
      // Check interrupt on main thread, if enough time passed since last check
      if (threadnum == 0 && getLapInterval(t_last_check) > kMinCheckInterval) {
        t_last_check = getLapStartTime();
        if (checkInterrupt()) {
          interrupt = true;
          continue;
        }
      }

      const auto& edgeid = unsettled_list[i];
      auto info = edgeid.getEdge().shared_info;
      int X = edgeid.X, Y = edgeid.Y;

      bcc.setCandidateZ(X, Y, info->zi_list);
      if (!info->zi_list.empty())
        searchForBestContributingNode(environment, X, Y, /* parallel */ false);

#ifdef _OPENMP
#pragma omp atomic
#endif
      ++n_jobs_done;
      if (threadnum == 0 && !environment.verbose)
        printProgress(static_cast<double>(n_jobs_done) / n_jobs_total,
            loop_start_time, progress_percentile);
    }
  }  // omp parallel
  if (interrupt) return false;
  if (!environment.verbose) {
    printProgress(1.0, loop_start_time, progress_percentile);  // finish
    Rcerr << '\n';
  }
  auto no_z_found = [&connected_list](auto& id) {
    if (id.getEdge().shared_info->top_z != -1) return false;
    // Move edges without top contributing node to connected_list
    connected_list.push_back(id);
    id.getEdge().shared_info->connected = 1;
    return true;
  };
  unsettled_list.erase(
      remove_if(begin(unsettled_list), end(unsettled_list), no_z_found),
      end(unsettled_list));

  return true;
}

bool searchForConditionalIndependence(Environment& environment) {
  auto& unsettled_list = environment.unsettled_list;

  if (environment.verbose)
    Rcout << "Number of unsettled edges: " << unsettled_list.size() << "\n";

  int iter_count = 0;
  int n_jobs_total = environment.unsettled_list.size();
  int progress_percentile = -1;
  auto loop_start_time = getLapStartTime();
  auto t_last_check = loop_start_time;
  while (!environment.unsettled_list.empty()) {
    if (getLapInterval(t_last_check) > kMinCheckInterval) {
      t_last_check = getLapStartTime();
      if (checkInterrupt()) return false;
    }
    ++iter_count;

    auto it_max = std::max_element(begin(unsettled_list), end(unsettled_list),
        [](const EdgeID& a, const EdgeID& b) {
          return a.getEdge().shared_info->Rxyz_ui <
                 b.getEdge().shared_info->Rxyz_ui;
        });
    int X = it_max->X, Y = it_max->Y;
    auto top_info = it_max->getEdge().shared_info;

    // move top z from zi_vect to ui_vect
    top_info->ui_list.push_back(top_info->top_z);
    top_info->zi_list.erase(remove(begin(top_info->zi_list),
                                end(top_info->zi_list), top_info->top_z),
        end(top_info->zi_list));
    top_info->top_z = -1;

    auto res = getCondMutualInfo(X, Y, top_info->ui_list,
        environment.data_numeric, environment.data_numeric_idx, environment);
    top_info->Nxy_ui = res.n_samples;
    top_info->Ixy_ui = res.I;
    top_info->kxy_ui = res.k;

    if (environment.verbose) {
      Rcout << "Edge " << iter_count << ": " << environment.nodes[X].name
            << " -- " << environment.nodes[Y].name << ":\n";
      Rcout << "ui: {";
      for (auto u : top_info->ui_list)
        Rcout << environment.nodes[u].name << ", ";
      Rcout << "}\n";
      Rcout << "Ixy_ui " << top_info->Ixy_ui << "\n";
      Rcout << "kxy_ui " << top_info->kxy_ui << "\n";
    }
    if (top_info->Ixy_ui - top_info->kxy_ui - environment.log_eta <= 0) {
      // Conditional independence found, remove edge
      unsettled_list.erase(it_max);
      environment.edges(X, Y).status = 0;
      environment.edges(Y, X).status = 0;
      environment.edges(X, Y).proba_head = -1;
      environment.edges(Y, X).proba_head = -1;
      top_info->connected = 0;
    } else {
      // Search for next candidate separating node
      searchForBestContributingNode(environment, X, Y, /* parallel */ true);

      if (environment.verbose) {
        int top_z = top_info->top_z;
        Rcout << "Edge " << environment.nodes[X].name << " -- "
              << environment.nodes[Y].name << ", best contributing node: ";
        Rcout << (top_z == -1 ? "NA" : environment.nodes[top_z].name) << "\n";
      }
      // Update the information about the edge
      if (top_info->top_z == -1) {
        // Move edge to connected list
        environment.connected_list.push_back(*it_max);
        unsettled_list.erase(it_max);
        top_info->connected = 1;
      }
    }
    if (!environment.verbose)
      printProgress(1.0 * (n_jobs_total - unsettled_list.size()) / n_jobs_total,
          loop_start_time, progress_percentile);
  }
  if (!environment.verbose) Rcerr << "\n";

  std::sort(begin(environment.connected_list), end(environment.connected_list));
  return true;
}

}  // namespace reconstruction
}  // namespace miic
