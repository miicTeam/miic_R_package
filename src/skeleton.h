#ifndef MIIC_RECONSTRUCT_H_
#define MIIC_RECONSTRUCT_H_

#include "biconnected_component.h"
#include "environment.h"

namespace miic {
namespace reconstruction {

bool initializeSkeleton(structure::Environment&);
bool setBestContributingNode(structure::Environment&, BiconnectedComponent&);
bool searchForConditionalIndependence(structure::Environment&);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_RECONSTRUCT_H_
