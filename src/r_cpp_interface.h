#ifndef MIIC_R_CPP_INTERFACE_H_
#define MIIC_R_CPP_INTERFACE_H_

#include <Rcpp.h>

#include "environment.h"

namespace miic {
namespace utility {
void setEnvironmentFromR(
    const Rcpp::List&, const Rcpp::List&, structure::Environment&);
}  // namespace utility
}  // namespace miic
#endif  // MIIC_R_CPP_INTERFACE_H_
