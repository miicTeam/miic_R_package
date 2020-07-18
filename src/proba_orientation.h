#ifndef MIIC_PROBA_ORIENTATION_H_
#define MIIC_PROBA_ORIENTATION_H_

#include <vector>

#include "structure.h"

namespace miic {
namespace reconstruction {
using ProbaArray = std::array<double, 4>;

std::vector<ProbaArray> getOriProbasList(const std::vector<structure::Triple>&,
    const std::vector<double>& I3_list, bool latent, bool degenerate,
    bool propagation, bool half_v_structure);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_PROBA_ORIENTATION_H_
