#ifndef MIIC_PROBA_ORIENTATION_H_
#define MIIC_PROBA_ORIENTATION_H_

#include <array>
#include <vector>

namespace miic {
namespace reconstruction {
using ProbaArray = std::array<double, 4>;
using Triple = std::array<int, 3>;

std::vector<ProbaArray> getOriProbasList(const std::vector<Triple>&,
    const std::vector<double>& I3_list, bool latent, bool degenerate,
    bool propagation, bool half_v_structure);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_PROBA_ORIENTATION_H_
