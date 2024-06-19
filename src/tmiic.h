//******************************************************************************
// Filename   : tmiic.h                              Creation date: 07 may 2020
//
// Description: header file of tmiic (temporal miic)
//
// Author     : Franck SIMON
//******************************************************************************
#ifndef TMIIC_
#define TMIIC_

#include <array>

#include "environment.h"


namespace tmiic {
// An unshielded Triple (X, Z, Y):
using Triple = std::array<int, 3>;

std::vector< std::pair<int, int> > getListLaggedEdges
  (miic::structure::Environment&, int, int);
void repeatEdgesOverHistory (miic::structure::Environment&);
void completeOrientationUsingTime (miic::structure::Environment&,
                                   const std::vector<Triple>&);
void dropPastEdges (miic::structure::Environment&);

}  // namespace tmiic

#endif  // TMIIC_
