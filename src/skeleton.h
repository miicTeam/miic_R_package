#ifndef MIIC_RECONSTRUCT_H_
#define MIIC_RECONSTRUCT_H_

#include "biconnected_component.h"
#include "environment.h"

namespace miic {
namespace reconstruction {

bool skeletonInitialization(structure::Environment&);
bool firstStepIteration(structure::Environment&, BiconnectedComponent&);
bool skeletonIteration(structure::Environment&);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_RECONSTRUCT_H_
