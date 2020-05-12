#ifndef MIIC_PROBA_ORIENTATION_H_
#define MIIC_PROBA_ORIENTATION_H_

#include "structure.h"

namespace miic {
namespace reconstruction {

double *getOrientTplLVDegPropag(structure::Environment&, int, int *, double *, int, int, int, int);
int OrientTpl_LV_Deg_Propag(int NbTpl, int *Tpl, double *I3,
    double *ProbArrowhead, int LV, int deg, int Propag, int);

}  // namespace reconstruction
}  // namespace miic

#endif  // MIIC_PROBA_ORIENTATION_H_
