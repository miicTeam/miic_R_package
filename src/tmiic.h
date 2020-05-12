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

#include "structure.h"

namespace tmiic {

using miic::structure::Environment;

void repeatEdgesOverHistory (miic::structure::Environment&);
void dropPastEdges (miic::structure::Environment&);

}  // namespace tmiic

#endif  // TMIIC_
