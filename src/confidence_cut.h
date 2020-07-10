#ifndef MIIC_CONFIDENCE_CUT_H_
#define MIIC_CONFIDENCE_CUT_H_

#include <vector>

#include "environment.h"

namespace miic {
namespace reconstruction {

std::vector<std::vector<std::string> > confidenceCut(structure::Environment&);
}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_CONFIDENCE_CUT_H_
