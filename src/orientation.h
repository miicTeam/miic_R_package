#ifndef MIIC_ORIENTATION_PROBABILITY_H_
#define MIIC_ORIENTATION_PROBABILITY_H_

#include <string>
#include <vector>

#include "environment.h"

namespace miic {
namespace reconstruction {

void updateAdj(structure::Environment&, int, int, double, double);
std::vector<std::vector<std::string>> orientationProbability(
    structure::Environment&);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_ORIENTATION_PROBABILITY_H_
