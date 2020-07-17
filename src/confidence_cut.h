#ifndef MIIC_CONFIDENCE_CUT_H_
#define MIIC_CONFIDENCE_CUT_H_

#include "environment.h"

namespace miic {
namespace reconstruction {

void setConfidence(structure::Environment&);
void confidenceCut(structure::Environment&);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_CONFIDENCE_CUT_H_
