/*
 * File: ./Algorithms/Dijkstra/StaticDijkstra.h
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

namespace Dijkstra {

class StaticDijkstra {
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

    StaticDijkstra(Graph::StaticGraph& graph)
        : graph(graph)
        , Q(graph.numberOfNodes())
        , labels(graph.numberOfNodes())
        , lastSource(0)
        , timeStamp(0)
    {
    }

    // Query methods

    inline void clear()
    {
        ++timeStamp;
        if (timeStamp == 0) {
            // reset distance arrays
            VertexLabel emptyLabel = VertexLabel();
            std::fill(labels.begin(), labels.end(), emptyLabel);
        }
        Q.clear();
    }

    inline void initSource(const VertexID source)
    {
        lastSource = source;
        VertexLabel& label = getLabel(lastSource);
        label.distance = 0;
        label.parent = lastSource;
        Q.update(&label);
    }

    template <typename STOP_CRIT = NO_OPERATION>
    inline VertexID settleVertex(const STOP_CRIT& stopCrit = NoOperation)
    {
        VertexLabel* uLabel = Q.extractFront();
        const VertexID u = uLabel - &(labels[0]);
        if (stopCrit(u, uLabel))
            return u;

        graph.doForAllOutgoingEdges(
            u, [&](auto /* source */, auto target, auto weight, auto /* edgeId */) {
                VertexID v = target;
                VertexLabel& vLabel = getLabel(v);

                if (uLabel->distance + weight < vLabel.distance) {
                    vLabel.distance = uLabel->distance + weight;
                    vLabel.parent = u;

                    Q.update(&vLabel);
                }
            });
        return u;
    }

    inline void run(const VertexID source, const VertexID target)
    {
        clear();
        initSource(source);

        while (!QisEmpty()) {
            settleVertex([&](auto u, auto /* uLabel */) { return (u == target); });
        }
    }

    // Access

    inline Graph::StaticGraph& getGraph() { return graph; }

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

private:
    Graph::StaticGraph& graph;
    ExternalKHeap<2, VertexLabel> Q;
    std::vector<VertexLabel> labels;

    VertexID lastSource;
    int timeStamp;
};
} // namespace Dijkstra
