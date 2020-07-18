#ifndef MIIC_ORIENTATION_PROBABILITY_H_
#define MIIC_ORIENTATION_PROBABILITY_H_

#include <string>
#include <vector>

#include "environment.h"
#include "structure.h"

namespace miic {
namespace reconstruction {

namespace reconstruction_impl {
double getI3(structure::Environment&, const structure::Triple&);
void updateAdj(structure::Environment&, int x, int y, double y2x, double x2y);
}  // namespace reconstruction_impl

std::vector<std::vector<std::string>> orientationProbability(
    structure::Environment&);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_ORIENTATION_PROBABILITY_H_
