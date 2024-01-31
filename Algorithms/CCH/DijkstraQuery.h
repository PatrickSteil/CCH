/*
 * File: ./Algorithms/CCH/DijkstraQuery.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include "../../DataStructures/Graph/Graphs.h"
#include "../../DataStructures/TypeDefs.h"
#include "../../Helpers/StatisticsCollecter.h"
#include "../Dijkstra/StaticDijkstra.h"

namespace CCH {

// The next two definitions are used in the StatisticsCollecter
enum CCH_DIJKSTRA_METRICS {};

enum CCH_DIJKSTRA_PHASES {
    DIJKSTRA_INIT_DS,
    DIJKSTRA_CLEAR,
    DIJKSTRA_QUERY,
};

static std::vector<std::string> CCH_DIJKSTRA_METRIC_NAMES = {};

static std::vector<std::string> CCH_DIJKSTRA_PHASE_NAMES = {
    "\"Init datastructures\"", "\"Clear datastructures\"", "\"Run Query\""
};

template <typename GRAPH>
class DijkstraQuery {
public:
    DijkstraQuery(Data<GRAPH>& data)
        : forwardDijkstra(getUpGraph(data))
        , backwardDijkstra(getDownGraph(data))
        , tentativeDistance(INFTY)
        , meetingVertex(-1)
        , statsCounter(CCH_DIJKSTRA_METRIC_NAMES, CCH_DIJKSTRA_PHASE_NAMES)
    {
    }

    inline Graph::StaticGraph& getUpGraph(Data<GRAPH>& data)
    {
        // Reserve space for static graph
        upGraph.reserve(data.getChordalGraph().numberOfNodes(),
            data.getChordalGraph().numberOfEdges());
        int upPrefix = 0;
        for (VertexID vertex = 0; vertex < data.getChordalGraph().numberOfNodes();
             ++vertex) {
            upGraph.getToAdj().push_back(upPrefix);

            data.getChordalGraph().doForAllOutgoingEdgeIDs(
                vertex, [&](auto /* vertex */, auto to, auto out_weight, auto /* in_weight */, auto edgeId) {
                    // get edges for upgraph
                    // check if edge is removed or not
                    if (!(data.getRemovalFlags()[edgeId] & UPWARD)) {
                        upPrefix++;
                        upGraph.getToVertex().push_back(to);
                        upGraph.getWeight().push_back(out_weight);
                    }
                });
        }
        upGraph.getToAdj().push_back(upPrefix);
        return upGraph;
    }

    inline Graph::StaticGraph& getDownGraph(Data<GRAPH>& data)
    {
        // Reserve space for static graph
        downGraph.reserve(data.getChordalGraph().numberOfNodes(),
            data.getChordalGraph().numberOfEdges());

        int downPrefix = 0;
        for (VertexID vertex = 0; vertex < data.getChordalGraph().numberOfNodes();
             ++vertex) {
            downGraph.getToAdj().push_back(downPrefix);

            data.getChordalGraph().doForAllOutgoingEdgeIDs(
                vertex, [&](auto /* vertex */, auto to, auto /* out_weight */, auto in_weight, auto edgeId) {
                    // get edges for downgraph
                    // check if edge is removed or not
                    if (!(data.getRemovalFlags()[edgeId] & DOWNWARD)) {
                        downPrefix++;
                        downGraph.getToVertex().push_back(to);
                        downGraph.getWeight().push_back(in_weight);
                    }
                });
        }
        downGraph.getToAdj().push_back(downPrefix);
        return downGraph;
    }

    inline void clear()
    {
        statsCounter.startPhase(DIJKSTRA_CLEAR);
        forwardDijkstra.clear();
        backwardDijkstra.clear();
        tentativeDistance = INFTY;
        meetingVertex = -1;
        statsCounter.stopPhase(DIJKSTRA_CLEAR);
    }

    inline void run(const VertexID source, const VertexID target)
    {
        statsCounter.startPhase(DIJKSTRA_QUERY);
        statsCounter.newRound();

        clear();

        statsCounter.startPhase(DIJKSTRA_INIT_DS);
        forwardDijkstra.initSource(source);
        backwardDijkstra.initSource(target);
        statsCounter.stopPhase(DIJKSTRA_INIT_DS);

        bool fwdSearch = false;

        while (!stopForwardDijkstra() || !stopBackwardDijkstra()) {
            fwdSearch = !fwdSearch;
            VertexID u;

            if ((fwdSearch && !stopForwardDijkstra()) || stopBackwardDijkstra()) {
                u = forwardDijkstra.settleVertex();
            } else {
                u = backwardDijkstra.settleVertex();
            }

            processVertex(u);
        }

        statsCounter.stopPhase(DIJKSTRA_QUERY);

        // std::cout << "source: " << source + 1 << "\t"
        //           << "target: " << target + 1 << "\t";
        // std::cout << "upDownPathRoot: " << meetingVertex + 1 << "\t";
        // std::cout << "minDistance: " << tentativeDistance << std::endl;
    }

    inline void processVertex(const VertexID v)
    {
        if (forwardDijkstra.getDistance(v) == INFTY || backwardDijkstra.getDistance(v) == INFTY)
            return;
        const size_t distance = forwardDijkstra.getDistance(v) + backwardDijkstra.getDistance(v);

        if (distance < tentativeDistance) {
            meetingVertex = v;
            tentativeDistance = distance;
        }
    }

    // stop criteria
    inline bool stopForwardDijkstra()
    {
        // return true <=> stop fwd dijkstra
        return forwardDijkstra.QisEmpty() || tentativeDistance <= forwardDijkstra.getMinKey();
    }

    inline bool stopBackwardDijkstra()
    {
        // return true <=> stop bwd dijkstra
        return backwardDijkstra.QisEmpty() || tentativeDistance <= backwardDijkstra.getMinKey();
    }

    inline size_t getDistance() { return tentativeDistance; }

    inline Graph::StaticGraph& getForwardGraph()
    {
        return forwardDijkstra.getGraph();
    }

    inline Graph::StaticGraph& getBackwardGraph()
    {
        return backwardDijkstra.getGraph();
    }

    // Print all the collected stats
    inline void showStats() const noexcept { statsCounter.printStats(); }

    inline FullStatisticsCollecter& getStatsCollecter() noexcept
    {
        return statsCounter;
    }

private:
    Graph::StaticGraph upGraph;
    Graph::StaticGraph downGraph;
    Dijkstra::StaticDijkstra forwardDijkstra;
    Dijkstra::StaticDijkstra backwardDijkstra;

    size_t tentativeDistance;
    VertexID meetingVertex;

    // StatisticsCollecter, which keeps track of all sort of metrics.
    FullStatisticsCollecter statsCounter;
};
} // namespace CCH
