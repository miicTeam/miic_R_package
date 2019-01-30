#ifndef _SKELETONITERATION_H_
#define _SKELETONITERATION_H_
#include "structure.h"

void skeletonIteration(Environment&);
void firstStepIteration(Environment&);
vector<int> bfs(const Environment& environment, int start, int end, const vector<int>& excludes=vector<int>());
bool is_consistent(const Environment& environment, int x, int y, int z);
bool is_consistent(const Environment& environment, int x, int y, const vector<int>& vect_z);

#endif
