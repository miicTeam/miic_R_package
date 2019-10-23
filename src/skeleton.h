#ifndef _SKELETON_H_
#define _SKELETON_H_

#include <cmath>
#include <deque>
#include <set>
#include <vector>
#include <assert.h>

#include "structure.h"

bool skeleton(Environment&, std::string, double, int, char*[]);

class CycleTracker {
    using uint = unsigned int;
    using edge_index_1d = uint;
private:
    // max number of iterations for which we keep the edgeStructure of all
    // changed edges
    static constexpr uint max_cycle_size = 10;
    // internal struct keeping record of edge information
    struct EdgeInfo {
        double Ixy_ui;  // mutual information with conditioning
        double cplx;  // complexity with conditioning
        int Nxy_ui;  // count of joint factors with conditioning and without NA

        bool is_empty = true;

        explicit EdgeInfo(const Edge& e) :
            Ixy_ui(e.edgeStructure->Ixy_ui),
            cplx(e.edgeStructure->cplx), Nxy_ui(e.edgeStructure->Nxy_ui),
            is_empty(false)
        {}
        EdgeInfo() = default;
    };

    struct Iteration {
        uint index;
        std::set<uint> changed_edges;
        std::vector<EdgeInfo> edge_info_list;

        Iteration(const Environment& env, uint i) : index(i),
            // edges undirected, keep track of the lower triangular part
            edge_info_list(env.numNodes * (env.numNodes - 1) / 2)
        {
            for (int i = 1; i < env.numNodes; ++i) {
                for (int j = 0; j < i; ++j) {
                    int index_1d = getEdgeIndex1D(i, j);
                    if (env.edges[i][j].areNeighboursAfterIteration ==
                            env.edges[i][j].isConnected) {
                        continue;
                    }
                    changed_edges.insert(index_1d);
                    edge_info_list[index_1d] = EdgeInfo(env.edges[i][j]);
                }
            }
        }
    };

    // internal class keeping track of iterations, max size is given
    // by max_cycle_size, new iteration is inserted at the front
    class IterationList {
        std::deque<Iteration> iteration_list_;
    public:
        void add(Iteration&& i) {
            iteration_list_.push_front(std::move(i));
            if (iteration_list_.size() > max_cycle_size)
                iteration_list_.pop_back();
        }
        Iteration& get(uint i) {
            return iteration_list_[i];
        }

        size_t size() { return iteration_list_.size(); }

        auto begin() { return iteration_list_.begin(); }
        auto cbegin() const { return iteration_list_.cbegin(); }
        auto end() { return iteration_list_.end(); }
        auto cend() const { return iteration_list_.cend(); }
    };

    Environment& env_;
    IterationList iterations_;
    uint n_saved = 0;  // number of saving operations performed
    // key: number of edges in the graph whose equal relation is a necessary
    // condition for a cycle
    // value: index of iteration, used to access set in changed_edges_list_
    std::multimap<uint, uint> edge_index_map_;
    // save changed edges of an iteration
    void saveIteration() {
        uint n_edge = env_.numNoMore;
        // absolute index of the iteration (since IterationList has fixed size)
        uint index = n_saved++;
        edge_index_map_.insert(std::make_pair(n_edge, index));
        // skip the first iteration as its previous step is the initial graph
        if (index != 0)
            iterations_.add(Iteration(env_, index));
    }
    // set info of edges of interest to its max value;
    void setEdgeState(edge_index_1d index, const EdgeInfo& edge_info);

public:
    CycleTracker(Environment& env) : env_(env) {}

    // convert lower triangular indices to 1d index
    static edge_index_1d getEdgeIndex1D(uint i, uint j) {
        assert(i != j);
        return (j < i ? j + i * (i - 1) / 2 : i + j * (j - 1) / 2);
    }

    // convert 1d index to lower triangular indices
    std::pair<uint, uint> getEdgeIndex2D(uint k) {
        // floor is equivalent to a int cast for positive number
        uint i = (uint) (0.5 + std::sqrt(0.25 + 2 * k));
        uint j = k - i * (i - 1) / 2;
        return std::make_pair(i, j);
    }
    // check if a cycle exists between the current and the past iterations
    bool hasCycle();
};

#endif
