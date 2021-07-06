#ifndef MIIC_PROBA_ORIENTATION_H_
#define MIIC_PROBA_ORIENTATION_H_

#include <array>
#include <vector>

namespace miic {
namespace reconstruction {
// An unshielded Triple (X, Z, Y):
using Triple = std::array<int, 3>;
// An array of orientation probabilities of arrowhead
// ProbaArray[i] > 0.5 means likely to be an arrowhead, ProbaArray[i] < 0.5
// means unlikely to be an arrowhead (thus likely to be a tail)
// For an unshielded Triple (X *-* Z *-* Y):
// ProbaArray[0]: probability of the arrowhead X <-* Z
// ProbaArray[1]: probability of the arrowhead X *-> Z
// ProbaArray[2]: probability of the arrowhead Z <-* Y
// ProbaArray[3]: probability of the arrowhead Z *-> Y
// where * is either arrowhead (< or >) or tail (-)
// X [0]-[1] Z [2]-[3] Y
using ProbaArray = std::array<double, 4>;

struct ProbaScore {
  double value = 0;
  // An unsettled score can either be updated by another score, or be used to
  // update other scores.
  bool settled = false;
};
using ScoreArray = std::array<ProbaScore, 4>;

std::vector<ProbaArray> getOriProbasList(const std::vector<Triple>&,
    const std::vector<double>& I3_list, const std::vector<int>& is_contextual,
    bool latent, bool degenerate, bool propagation, bool half_v_structure);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_PROBA_ORIENTATION_H_
