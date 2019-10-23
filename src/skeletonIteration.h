#ifndef _SKELETONITERATION_H_
#define _SKELETONITERATION_H_
#include <set>
#include <stack>

#include "structure.h"

bool skeletonIteration(Environment&);
bool firstStepIteration(Environment&);
vector<int> bfs(const Environment& environment, int start, int end, const vector<int>& excludes=vector<int>());
bool is_consistent(const Environment& environment, int x, int y, int z);
bool is_consistent(const Environment& environment, int x, int y, const vector<int>& vect_z);

class BCC {
    Environment& environment;
    vector<int> is_cp;
    vector<int> degree_of;
    vector<int> bc_tree_rep;
    vector<int> bc_tree_inverse_index;
    vector<int> bc_tree_node_is_cp;
    vector<set<int> > bcc_list;
    vector<set<int> > bcc_set_indices;
    vector<set<int> > bc_tree_adj_list;

    void bcc_aux(int u, int& time, vector<int>& parent, vector<int>& lowest,
            vector<int>& depth, stack<pair<int, int> >& st);
    void bcc();
    vector<int> bc_tree_bfs(int start, int end);

    public:

    BCC(Environment& env):
        environment(env),
        is_cp(env.numNodes, 0),
        degree_of(env.numNodes, 0),
        bc_tree_rep(env.numNodes, -1),
        bcc_set_indices(env.numNodes, set<int>())
    {
        if (env.consistentPhase)
            bcc();
    }

    void set_candidate_z(int x, int y);
};

#endif
