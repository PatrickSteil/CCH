/*
 * File: ./Algorithms/CCH/Preprocessing.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <algorithm>
#include <cassert>

#include "../../DataStructures/CCH/Data.h"
#include "../../DataStructures/Graph/Graphs.h"
#include "../../DataStructures/TypeDefs.h"
#include "../../Helpers/Colors.h"
#include "../../Helpers/ProgressBar.h"
#include "../../Helpers/StatisticsCollecter.h"

// The class is used to compute the chordal graph given a CCH::Data object. It
// sorts the (outgoing / incoming) edges of each vertex, in order to sweep over
// the adjacency vectors to find triangles
namespace CCH {

// ******************************
// The next two definitions are used in the StatisticsCollecter
enum METRICS { NUM_OF_SHORTCUTS_ADDED };

static std::vector<std::string> METRIC_NAMES = { "\"# of shortcuts added:\"" };

enum PHASES {
    COPY_GRAPH_TO_BUILDER,
    CONTRACTTION,
    SORT_EDGES,
    CONVERT_TO_STATIC
};

static std::vector<std::string> PHASE_NAMES = {
    "\"Copy given graph to builder graph\"",
    "\"Contraction\"",
    "\"Sort Edges\"",
    "\"Convert to Static\"",
};
// ******************************

template <typename GRAPH, bool VERBOSE = true, bool DEBUG = false>
class Preprocessing {
public:
    Preprocessing(Data<GRAPH>& data)
        : data(data)
        , bobTheBuilder(data.getCore().numberOfNodes())
        , numberOfEdges(0)
        , statsCounter(METRIC_NAMES, PHASE_NAMES)
    {
        for (VertexID node = 0; node < data.getCore().numberOfNodes(); ++node) {
            // too optimistic, but so we kind of assume no realloc
            bobTheBuilder[node].reserve(data.getCore().getDegree(node));
        }
        copyGraphToBobTheBuilder();
    };

    inline void copyGraphToBobTheBuilder() noexcept
    {
        statsCounter.startPhase(COPY_GRAPH_TO_BUILDER);
        // copies the graph, but flip all edges upwards
        for (VertexID node = 0; node < data.getCore().numberOfNodes(); ++node) {
            for (auto& edge : data.getCore().getEdgesFrom(node)) {
                // which one is smaller?
                VertexID from = std::min(node, edge.otherVertex);
                VertexID to = std::max(node, edge.otherVertex);
                if (!hasEdge(from, to)) {
                    bobTheBuilder[from].push_back(to);
                    ++numberOfEdges;
                }
            }
        }
        statsCounter.stopPhase(COPY_GRAPH_TO_BUILDER);
    }

    inline bool hasEdge(const VertexID from, const VertexID to) noexcept
    {
        assert(from < bobTheBuilder.size());
        assert(to < bobTheBuilder.size());

        return std::find(bobTheBuilder[from].begin(), bobTheBuilder[from].end(),
                   to)
            != bobTheBuilder[from].end();
    }

    // call this method to start the overall contraction
    void runContraction() noexcept
    {
        assert(data.getCore().numberOfNodes() > 0);

        // progressbar progress(data.getCore().numberOfNodes());

        if (VERBOSE)
            std::cout << CLR_BLUE << "[START]" << CLR_RESET << " contraction"
                      << std::endl;

        statsCounter.startPhase(CONTRACTTION);

        data.doForAllVerticesBottomToTop(
            [&](VertexID vertex) { contractVertex(vertex); });

        statsCounter.stopPhase(CONTRACTTION);

        // sort the toVertices

        statsCounter.startPhase(SORT_EDGES);

        for (auto& adj : bobTheBuilder) {
            std::sort(adj.begin(), adj.end());
        }

        statsCounter.stopPhase(SORT_EDGES);

        if (VERBOSE)
            std::cout << CLR_RED << "[DONE]" << CLR_RESET << " contraction"
                      << std::endl;

        statsCounter.startPhase(CONVERT_TO_STATIC);

        // convert the vector<vector<>> to a directed static graph
        data.getChordalGraph().copy(bobTheBuilder, numberOfEdges);
        // set removal flags for all egdes to false
        data.getRemovalFlags().assign(data.getChordalGraph().numberOfEdges(),
            NOTHING);

        statsCounter.stopPhase(CONVERT_TO_STATIC);

        if (DEBUG)
            areAllTrianglesThere();
    }

    // method to contract the given vertex
    inline void contractVertex(VertexID vertex) noexcept
    {
        assert(data.getCore().isVertex(vertex));

        // insert all ToVertices to the smallest neighbour
        statsCounter.newRound();

        if (bobTheBuilder[vertex].size() == 0) [[unlikely]]
            return;

        // find parent
        VertexID parent = *std::min_element(bobTheBuilder[vertex].begin(),
            bobTheBuilder[vertex].end());

        for (auto toVertex : bobTheBuilder[vertex]) {
            if (toVertex != parent && !hasEdge(parent, toVertex)) {
                bobTheBuilder[parent].push_back(toVertex);
                ++numberOfEdges;
                statsCounter.count(NUM_OF_SHORTCUTS_ADDED);
            }
        }
    }

    // Print all the collected stats
    inline void showStats() const noexcept { statsCounter.printStats(); }

    // to use for asserts
    // checks that all edges are in triangles

    inline void areAllTrianglesThere() noexcept
    {
        for (VertexID node(0); node < data.getChordalGraph().numberOfNodes();
             ++node) {
            data.getChordalGraph().doForAllOutgoingEdgeIDs(
                node, [&](VertexID x, VertexID y, auto /* upWeight */, auto /* downWeight */, auto /* id */) {
                    if (node == y)
                        return;
                    // check that for all edges starting at to, from -> (new to vertex)
                    data.getChordalGraph().doForAllOutgoingEdgeIDs(
                        node, [&](VertexID /* x */, VertexID z, auto /* upWeight */, auto /* downWeight */, auto /* id */) {
                            if (y == z || x == z)
                                return;

                            assert(data.getChordalGraph().hasEdge(std::min(y, z),
                                std::max(y, z)));
                        });
                });
        }
    }

    inline FullStatisticsCollecter& getStatsCollecter() noexcept
    {
        return statsCounter;
    }

private:
    Data<GRAPH>& data;
    std::vector<std::vector<VertexID>> bobTheBuilder;
    size_t numberOfEdges;

    // StatisticsCollecter, which keeps track of all sort of metrics.
    FullStatisticsCollecter statsCounter;
};

} // namespace CCH
