#ifndef _FUNCTIONSORIENTATION_H_
#define _FUNCTIONSORIENTATION_H_
#include "structure.h"

void getAllStructures(Environment&, bool);
void prepareDFforPropagation(Environment&, std::string);
bool hasLowerIndex(std::vector<StructWithOrtToPropagate*>, int);
bool simpleContradictionORrestartAll(Environment&, StructWithOrtToPropagate*, int,bool);
bool setNewOrt(Environment&, int, StructWithOrtToPropagate*&, int, int, bool);
void learnNewOrtLocally(Environment&, int, StructWithOrtToPropagate*&, int, int, std::vector<StructWithOrtToPropagate*>&, bool);
bool saveTableOfOrientations(Environment, std::string);

#endif
