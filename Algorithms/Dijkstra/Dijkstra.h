/*
 * File: ./Algorithms/Dijkstra/Dijkstra.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <cassert>
#include <vector>

#include "../../DataStructures/Container/ExternalKHeap.h"
#include "../../DataStructures/Graph/Graphs.h"
#include "../../DataStructures/TypeDefs.h"
#include "../../Helpers/StatisticsCollecter.h"

namespace Dijkstra {

// The next two definitions are used in the StatisticsCollecter
enum DIJKSTRA_METRICS {
    NUM_OF_RELAXED_EDGES,
    NUM_EXPLORED_VERTICES,
    NUM_UPDATED_LABELS,
};

enum DIJKSTRA_PHASES {
    INIT_DS,
    CLEAR,
    QUERY,
    RELAX,
};

static std::vector<std::string> DIJKSTRA_METRIC_NAMES = {
    "\"# of relaxed edges:\"",
    "\"# of explored vertices:\"",
    "\"# of updated labels:\"",

};

static std::vector<std::string> DIJKSTRA_PHASE_NAMES = {
    "\"Init datastructures\"", "\"Clear datastructures\"", "\"Run Query\"",
    "\"Relax \""
};

template <typename GRAPH>
class Dijkstra {
public:
    struct VertexLabel : public ExternalKHeapElement {
        VertexLabel()
            : ExternalKHeapElement()
            , distance(INFTY)
            , parent(-1)
            , timeStamp(-1)
        {
        }
        inline void reset(int time)
        {
            distance = INFTY;
            parent = -1;
            timeStamp = time;
        }
        inline bool hasSmallerKey(const VertexLabel* other) const
        {
            return distance < other->distance;
        }

        size_t distance;
        VertexID parent;
        int timeStamp;
    };

    Dijkstra(GRAPH& graph)
        : graph(graph)
        , Q(graph.numberOfNodes())
        , labels(graph.numberOfNodes())
        , lastSource(0)
        , timeStamp(0)
        , statsCounter(DIJKSTRA_METRIC_NAMES, DIJKSTRA_PHASE_NAMES)
    {
    }

    // Query methods

    inline void clear()
    {
        statsCounter.startPhase(CLEAR);
        ++timeStamp;
        if (timeStamp == 0) {
            // reset distance arrays
            VertexLabel emptyLabel = VertexLabel();
            std::fill(labels.begin(), labels.end(), emptyLabel);
        }
        Q.clear();
        statsCounter.stopPhase(CLEAR);
    }

    inline void initSource(const VertexID source)
    {
        statsCounter.startPhase(INIT_DS);
        statsCounter.count(NUM_UPDATED_LABELS);
        lastSource = source;
        VertexLabel& label = getLabel(lastSource);
        label.distance = 0;
        label.parent = lastSource;
        Q.update(&label);
        statsCounter.stopPhase(INIT_DS);
    }

    template <bool OUTGOING = true, typename STOP_CRIT = NO_OPERATION>
    inline VertexID settleVertex(const STOP_CRIT& stopCrit = NoOperation)
    {
        statsCounter.startPhase(RELAX);
        statsCounter.count(NUM_EXPLORED_VERTICES);
        VertexLabel* uLabel = Q.extractFront();
        const VertexID u = uLabel - &(labels[0]);

        if (stopCrit(u, uLabel))
            return u;

        if (OUTGOING) {
            graph.doForAllOutgoingEdges(
                u, [&](auto /* source */, auto target, auto weight) {
                    statsCounter.count(NUM_OF_RELAXED_EDGES);
                    VertexID v = target;
                    VertexLabel& vLabel = getLabel(v);

                    if (uLabel->distance + weight < vLabel.distance) {
                        statsCounter.count(NUM_UPDATED_LABELS);
                        vLabel.distance = uLabel->distance + weight;
                        vLabel.parent = u;

                        Q.update(&vLabel);
                    }
                });
        } else {
            graph.doForAllIncomingEdges(
                u, [&](auto source, auto /* target */, auto weight) {
                    VertexID v = source;
                    VertexLabel& vLabel = getLabel(v);

                    if (uLabel->distance + weight < vLabel.distance) {
                        vLabel.distance = uLabel->distance + weight;
                        vLabel.parent = u;

                        Q.update(&vLabel);
                    }
                });
        }
        statsCounter.stopPhase(RELAX);
        return u;
    }

    inline void run(const VertexID source, const VertexID target)
    {
        statsCounter.startPhase(QUERY);
        statsCounter.newRound();

        clear();
        initSource(source);

        while (!QisEmpty()) {
            settleVertex([&](auto u, auto /* uLabel */) { return (u == target); });
        }
        statsCounter.stopPhase(QUERY);
    }

    // Access

    inline GRAPH& getGraph() { return graph; }

    inline ExternalKHeap<2, VertexLabel>& getQueue() { return Q; }

    inline VertexLabel& getLabel(const VertexID vertex) noexcept
    {
        assert((size_t)vertex < labels.size());
        VertexLabel& result = labels[vertex];
        if (result.timeStamp != timeStamp)
            result.reset(timeStamp);
        return result;
    }

    inline size_t getDistance(const VertexID vertex)
    {
        VertexLabel& label = getLabel(vertex);
        return label.distance;
    }

    inline bool QisEmpty() { return Q.empty(); }

    inline size_t getMinKey() { return Q.front()->distance; }

    // Print all the collected stats
    inline void showStats() const noexcept { statsCounter.printStats(); }

    inline FullStatisticsCollecter& getStatsCollecter() noexcept
    {
        return statsCounter;
    }

private:
    GRAPH& graph;
    ExternalKHeap<2, VertexLabel> Q;
    std::vector<VertexLabel> labels;

    VertexID lastSource;
    int timeStamp;

    // StatisticsCollecter, which keeps track of all sort of metrics.
    FullStatisticsCollecter statsCounter;
};
} // namespace Dijkstra
