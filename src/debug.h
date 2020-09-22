#ifndef MIIC_DEBUG_H_
#define MIIC_DEBUG_H_

#include <iostream>
#include <array>
#include <Rcpp.h>

#include "structure.h"
#include "environment.h"

template<typename T>
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
}  // namespace debug
}  // namespace miic

#endif  // MIIC_UTILITIES_H_
