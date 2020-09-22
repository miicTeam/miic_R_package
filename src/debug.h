#ifndef MIIC_DEBUG_H_
#define MIIC_DEBUG_H_

#include <iostream>
#include <array>
#include <Rcpp.h>

#include "environment.h"
#include "structure.h"
#include "linear_allocator.h"

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    out << "{";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];
        if (i != last) 
            out << ", ";
    }
    out << "}";
    return out;
}

namespace miic {
namespace debug {
using namespace miic::structure;
using namespace miic::utility;
using Triple = std::array<int, 3>;
using ProbaArray = std::array<double, 4>;

// std::string printNodesName(const structure::Environment&);
void debugEdges (structure::Environment&, std::string title="", 
                 bool filter_status=true, bool half_only=true);
void debugConnectedList (structure::Environment&, std::string title="");
void debugAdjacencyMatrix (structure::Environment&, std::string title="", 
                           std::string status_field="status"); 
// void printMatrix(const structure::Environment&, std::string);
// void printEnvironment(const structure::Environment&);
void debugTriples (const structure::Environment& environment, 
                   const std::vector<Triple>& triples, 
                   std::string title);
void debugProbaArray (const structure::Environment& environment, 
                   const std::vector<Triple>& triples, 
                   const std::vector<ProbaArray>& prob_array, 
                   std::string title);
void debugProbaMap (const structure::Environment& environment, 
                   const std::map<std::pair<int, int>, double>& proba_map, 
                   std::string title);  
template <typename T>
void debugVector (const std::vector<double> &vect,
                  const std::string &title, 
                  const std::size_t max=15);
void debugTempVectorDouble (const miic::structure::TempVector<double> &vect, 
                  const std::string &title, 
                  const std::size_t max=15);
void debugTempVectorInt (const miic::structure::TempVector<int> &vect, 
                  const std::string &title, 
                  const std::size_t max=15);
template <typename T, typename = IsIntContainer<T>>
void debugIntContainer (const T &vect, const std::string &title, 
                        const std::size_t max=15)
  {
  Rcpp::Rcout << title
              << " (vect size=" << vect.size() << ")";
  for(std::size_t i = 0; i < vect.size(); i++)
    {
    if (i>max)
      {
      Rcpp::Rcout << " ...";
      break;
      }
    Rcpp::Rcout << " " << vect[i];
    }
  Rcpp::Rcout << "\n";
  }


template <typename T>
void debugGrid (const structure::Grid2d<T> &grid,
                const std::string &title,  bool transpose=false,
                const std::size_t max_rows=15, const std::size_t max_cols=15);
void debugTempGridInt (const structure::TempGrid2d<int> &grid,
                const std::string &title,  bool transpose=false,
                const std::size_t max_rows=15, const std::size_t max_cols=15);
void debugTempGridDouble (const structure::TempGrid2d<double> &grid,
                const std::string &title,  bool transpose=false,
                const std::size_t max_rows=15, const std::size_t max_cols=15);
void debugCutPointsInfo (std::shared_ptr<CutPointsInfo> &cuts_info,
                const std::string &title,  bool transpose=false,
                const std::size_t max_rows=15, const std::size_t max_cols=15);
}  // namespace debug
}  // namespace miic

#endif  // MIIC_UTILITIES_H_
