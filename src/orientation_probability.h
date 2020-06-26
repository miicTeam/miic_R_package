#ifndef MIIC_ORIENTATION_PROBABILITY_H_
#define MIIC_ORIENTATION_PROBABILITY_H_

#include <string>
#include <vector>

#include "environment.h"
#include "structure.h"

namespace miic {
namespace reconstruction {

std::vector<std::vector<std::string> > orientationProbability(
    structure::Environment&);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_ORIENTATION_PROBABILITY_H_
