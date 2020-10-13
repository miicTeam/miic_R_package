//*****************************************************************************
// Filename   : tmiic.h                     Namespace: tmiic                                   
//
// Author     : Franck SIMON                Creation date: 07 may 2020 
//
// Description: header file of tmiic (temporal miic)
//
// Changes history:
// - 07 may 2020 : initial version
//*****************************************************************************
#ifndef TMIIC_
#define TMIIC_

#include "environment.h"


namespace tmiic {
// An unshielded Triple (X, Z, Y):
using Triple = std::array<int, 3>;

void repeatEdgesOverHistory (miic::structure::Environment&);
void completeOrientationUsingTime (miic::structure::Environment&,
                                   const std::vector<Triple>&);
void dropPastEdges (miic::structure::Environment&);

}  // namespace tmiic

#endif  // TMIIC_
