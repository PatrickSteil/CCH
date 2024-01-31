/*
 * File: ./DataStructures/CCH/Data.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <omp.h>

#include <cassert>
#include <vector>

#include "../../Extern/vector_io.h"
#include "../../Helpers/TypeChecker.h"
#include "../Graph/Graphs.h"
#include "../TypeDefs.h"
#include "DirectedTree.h"

// This class stores the graph, for which we want to compute the CCH.
// Chordal Graph is the graph that contains the shortcuts, and on which we
// perform the triangle-stuff and so on
namespace CCH {

// use to mark which direction an edge should be removed
enum RemovalFlag { NOTHING = 0,
    UPWARD = 1,
    DOWNWARD = 2 };

template <typename GRAPH>
class Data {
public:
    Data(GRAPH& graph, const int numberOfThreads = 1)
        : core()
        , chordalGraph(graph.numberOfNodes())
        , removalFlags()
        , edgeMapping(graph.numberOfEdges(), std::make_pair(noIndex, false))
        , parents(graph.numberOfNodes(), noVertexID)
        , toAdjOfLevel(0)
        , vertexOfLevel(0)
    {
        omp_set_num_threads(numberOfThreads);
        copyGraphToCore(graph);
    };

    // the 'core' graph needs to be dynamic, hence we copy the entire given graph
    // into this dynamicgraph
    inline void copyGraphToCore(GRAPH& graph) noexcept
    {
        // note: no edge weight is copied, edge weights are set to INFTY
        core.addVertices(graph.numberOfNodes());
        core.clear();

#pragma omp parallel for schedule(dynamic)
        for (VertexID node = 0; node < graph.numberOfNodes(); ++node) {
            core.reserveOutEdges(node, graph.getOutDegree(node));
            core.reserveInEdges(node, graph.getInDegree(node));
        }

        for (VertexID node(0); node < graph.numberOfNodes(); ++node) {
            graph.doForAllOutgoingEdgeIDs(
                node, [&](VertexID from, VertexID to, size_t /* weight */, auto /* edgeId */) { core.addEdge(from, to, INFTY); });
        }

        core.sortEdges();
    }

    // some getter, note we want references, since we change some graph attributes
    // like isolateVertex, ...
    inline Graph::DynamicGraph& getCore() noexcept { return core; }

    // access the chordal graph directly (as reference, since we want to change
    // weights later)
    inline Graph::DiStaticGraph& getChordalGraph() noexcept
    {
        return chordalGraph;
    }

    // returns the parent array
    inline std::vector<VertexID>& getParents() noexcept { return parents; }

    // simple getter + setter
    // inline std::vector<Index>& getRank() noexcept { return rank; }
    // inline std::vector<VertexID>& getOrder() noexcept { return order; }
    inline std::vector<uint8_t>& getRemovalFlags() noexcept
    {
        return removalFlags;
    }

    // sets the edge for removal
    inline void setForRemoval(Index edge, RemovalFlag flag) noexcept
    {
        assert(edge < removalFlags.size());

        removalFlags[edge] |= flag;
    }

    inline void computeEdgeMapping(Graph::DynamicGraph& graph) noexcept
    {
        // read the weights and so on from the graph
        edgeMapping.assign(graph.numberOfEdges(), std::make_pair(noIndex, false));

        doForAllVerticesTopToBottomInParallel([&](const VertexID node) {
            graph.doForAllOutgoingEdgeIDs(
                node,
                [&](auto /* from */, auto target, auto /* weight */, auto edgeIndex) {
                    VertexID from = std::min(node, target);
                    VertexID to = std::max(node, target);
                    Index edge = chordalGraph.findEdge(from, to);
                    assert(edge != noIndex);

                    // (from == node) <=> true iff upwards edge
                    edgeMapping[edgeIndex] = std::make_pair(edge, (from == node));
                });
        });
    }

    // get entire edge mapping
    inline std::vector<std::pair<Index, bool>>& getEdgeMapping() noexcept
    {
        return edgeMapping;
    }

    // or get only the new edge (in the chordal graph) of the original edge id
    inline std::pair<Index, bool>& getEdgeMapping(Index edgeId) noexcept
    {
        assert(edgeId < edgeMapping.size());
        return edgeMapping[edgeId];
    }

    // ******************************
    // The next methods make it easier to loop over lower, middle and upper
    // triangles
    // Note: you need the edges to be sorted
    // ******************************

    // function takes the parameters (x, y, z, edgeindex of xz, edgeindex of yz)
    template <typename FUNC>
    inline void doForAllLowerTriangles(const VertexID x, const VertexID y,
        FUNC&& function) noexcept
    {
        // lower triangle of {x, y} <=> z < x < y
        assert(core.isVertex(x) && core.isVertex(y));
        assert(x < y); // the rank check

        assert(areEdgesSorted(x));
        assert(areEdgesSorted(y));

        // the idea is to sweep over the two adjacency vectors in one "sweep", like
        // e.g., use in hub labeling

        std::vector<Index>& xEdges = chordalGraph.getIncomingEdgeIDs(x);
        std::vector<Index>& yEdges = chordalGraph.getIncomingEdgeIDs(y);

        Index xIndex(0);
        Index yIndex(0);

        while (xIndex < xEdges.size() && yIndex < yEdges.size()) {
            VertexID xNeighbour = chordalGraph.getFromVertex(xEdges[xIndex]);
            VertexID yNeighbour = chordalGraph.getFromVertex(yEdges[yIndex]);

            assert(core.isVertex(xNeighbour));
            assert(core.isVertex(yNeighbour));

            if (xNeighbour == yNeighbour) {
                // we found a common vertex => now process the edges
                // Note: xNeighbour == yNeighbour == "z"

                function(x, y, xNeighbour, xEdges[xIndex], yEdges[yIndex]);
            }

            xIndex += !(xNeighbour > yNeighbour);
            yIndex += !(xNeighbour < yNeighbour);
        }
    }

    // function takes the parameters (x, y, z, edgeindex of xz, edgeindex of yz)
    // enumerates until first lower triangle is found where the function returns
    // true then breaks and returns true itself. Otherwise returns false.
    template <typename FUNC>
    inline bool doForAllLowerTrianglesAndBreak(const VertexID x, const VertexID y,
        FUNC&& function) noexcept
    {
        // lower triangle of {x, y} <=> z < x < y
        assert(core.isVertex(x) && core.isVertex(y));
        assert(x < y); // the rank check

        assert(areEdgesSorted(x));
        assert(areEdgesSorted(y));

        // the idea is to sweep over the two adjacency vectors in one "sweep", like
        // e.g., use in hub labeling

        std::vector<Index>& xEdges = chordalGraph.getIncomingEdgeIDs(x);
        std::vector<Index>& yEdges = chordalGraph.getIncomingEdgeIDs(y);

        Index xIndex(0);
        Index yIndex(0);

        while (xIndex < xEdges.size() && yIndex < yEdges.size()) {
            VertexID xNeighbour = chordalGraph.getFromVertex(xEdges[xIndex]);
            VertexID yNeighbour = chordalGraph.getFromVertex(yEdges[yIndex]);

            assert(core.isVertex(xNeighbour));
            assert(core.isVertex(yNeighbour));

            if (xNeighbour < yNeighbour) {
                ++xIndex;
            } else if (xNeighbour > yNeighbour) {
                ++yIndex;
            } else {
                // we found a common vertex => now process the edges
                // Note: xNeighbour == yNeighbour == "z"
                if (function(x, y, xNeighbour, xEdges[xIndex], yEdges[yIndex])) {
                    return true;
                }
                ++xIndex;
                ++yIndex;
            }
        }
        return false;
    }

    // function takes the parameters (x, y, z, edgeindex of xz, edgeindex of yz)
    template <typename FUNC>
    inline void doForAllMiddleTriangles(const VertexID x, const VertexID y,
        FUNC&& function) noexcept
    {
        // lower triangle of {x, y} <=> x < z < y
        assert(core.isVertex(x) && core.isVertex(y));
        assert(x < y); // the rank check

        assert(areEdgesSorted(x));
        assert(areEdgesSorted(y));

        // the idea is to sweep over the two adjacency vectors in one "sweep", like
        // e.g., use in hub labeling

        // for a middle triangle, we need to scan the outgoing edges of x and the
        // incoming of y

        std::vector<Index>& yEdges = chordalGraph.getIncomingEdgeIDs(y);

        Index xIndex(chordalGraph.getToAdj(x));
        Index yIndex(0);

        while (xIndex < chordalGraph.getToAdj(x + 1) && yIndex < yEdges.size()) {
            VertexID xNeighbour = chordalGraph.getToVertex(xIndex);
            VertexID yNeighbour = chordalGraph.getFromVertex(yEdges[yIndex]);

            assert(core.isVertex(xNeighbour));
            assert(core.isVertex(yNeighbour));

            if (xNeighbour == yNeighbour) {
                // we found a common vertex => now process the edges
                // Note: xNeighbour == yNeighbour == "z"
                function(x, y, xNeighbour, xIndex, yEdges[yIndex]);
            }
            xIndex += !(xNeighbour > yNeighbour);
            yIndex += !(xNeighbour < yNeighbour);
        }
    }

    // function takes the parameters (x, y, z, edgeindex of xz, edgeindex of yz)
    template <typename FUNC>
    inline void doForAllUpperTriangles(const VertexID x, const VertexID y,
        FUNC&& function) noexcept
    {
        // upper triangle of {x, y} <=> x < y < z
        assert(core.isVertex(x) && core.isVertex(y));
        assert(x < y); // the rank check

        assert(areEdgesSorted(x));
        assert(areEdgesSorted(y));

        // the idea is to sweep over the two adjacency vectors in one "sweep", like
        // e.g., use in hub labeling

        Index xIndex(chordalGraph.getToAdj(x));
        Index yIndex(chordalGraph.getToAdj(y));

        while (xIndex < chordalGraph.getToAdj(x + 1) && yIndex < chordalGraph.getToAdj(y + 1)) {
            VertexID xNeighbour = chordalGraph.getToVertex(xIndex);
            VertexID yNeighbour = chordalGraph.getToVertex(yIndex);

            assert(core.isVertex(xNeighbour));
            assert(core.isVertex(yNeighbour));

            if (xNeighbour == yNeighbour) {
                // we found a common vertex => now process the edges
                // Note: xNeighbour == yNeighbour == "z"
                function(x, y, xNeighbour, xIndex, yIndex);
            }
            xIndex += !(xNeighbour > yNeighbour);
            yIndex += !(xNeighbour < yNeighbour);
        }
    }

    // used in asserts to assert that the adjacency vector of vertex is sorted
    inline bool areEdgesSorted(const VertexID vertex) noexcept
    {
        return std::is_sorted(
            chordalGraph.getToVertex().begin() + chordalGraph.getToAdj()[vertex],
            chordalGraph.getToVertex().begin() + chordalGraph.getToAdj()[vertex + 1]);
    }

    // loop over the vertices in order (and also in reversed order)
    // function takes the parameters (index, vertex)
    template <typename FUNC>
    inline void doForAllVerticesBottomToTop(FUNC&& function) noexcept
    {
        for (VertexID node(0); node < core.numberOfNodes(); ++node) {
            function(node);
        }
    }

    // function takes the parameters (index, vertex)
    template <typename FUNC>
    inline void doForAllVerticesTopToBottom(FUNC&& function) noexcept
    {
        for (Index i(0); i < core.numberOfNodes(); ++i) {
            function(core.numberOfNodes() - 1 - i);
        }
    }

    // loop over the vertices in order (and also in reversed order)
    // function takes the parameters (index, vertex)
    // Note: this method runs in parallel, hence pay attention to what your
    // function does
    template <typename FUNC>
    inline void doForAllVerticesBottomToTopInParallel(FUNC&& function) noexcept
    {
#pragma omp parallel for schedule(dynamic)
        for (VertexID node = 0; node < core.numberOfNodes(); ++node) {
            function(node);
        }
    }

    // function takes the parameters (index, vertex)
    // Note: this method runs in parallel, hence pay attention to what your
    // function does
    template <typename FUNC>
    inline void doForAllVerticesTopToBottomInParallel(FUNC&& function) noexcept
    {
#pragma omp parallel for schedule(dynamic)
        for (Index i = 0; i < core.numberOfNodes(); ++i) {
            function(core.numberOfNodes() - 1 - i);
        }
    }

    // ****************************************
    // parallel stuff
    inline void buildLevelDataStructure() noexcept
    {
        buildParentsVector();

        DirectedTree tree(parents);
        // levelAssignment[vertex] == level of vertex
        // note: root is level 0
        std::vector<int> levelAssignment = tree.computeLevels();

        // now build the adj datastructure
        int maxLevel = *std::max_element(levelAssignment.begin(), levelAssignment.end());

        std::vector<size_t> numVerticesPerLevel(maxLevel + 1, 0);

        for (auto& level : levelAssignment) {
            if (level == -1)
                level = 0;
            else
                level = maxLevel - level;
            assert(0 <= level && level <= maxLevel);
            ++numVerticesPerLevel[level];
        }

        toAdjOfLevel.resize(maxLevel + 2);
        vertexOfLevel.resize(levelAssignment.size());

        size_t runningSum(0);

        for (int level(0); level <= maxLevel; ++level) {
            toAdjOfLevel[level] = runningSum;
            runningSum += numVerticesPerLevel[level];
            // trick to use bucket sort later
            numVerticesPerLevel[level] = 0;
        }
        toAdjOfLevel.back() = runningSum;

        for (Index i = 0; i < levelAssignment.size(); ++i) {
            auto level = levelAssignment[i];
            // where to write it?
            Index positionToWrite = toAdjOfLevel[level] + numVerticesPerLevel[level];
            assert(positionToWrite < toAdjOfLevel[level + 1]);
            assert(positionToWrite < vertexOfLevel.size());
            ++numVerticesPerLevel[level];

            vertexOfLevel[positionToWrite] = i;
        }
    }

    inline void buildParentsVector() noexcept
    {
        doForAllVerticesBottomToTopInParallel([&](auto vertex) {
            VertexID newParent = chordalGraph.numberOfNodes() - 1; // the root of the chordal graph

            chordalGraph.doForAllOutgoingEdgeIDs(
                vertex,
                [&newParent](auto /* from */, auto toVertex, auto /* upWeight */,
                    auto /* downWeight */, auto /* i */) {
                    newParent = std::min(toVertex, newParent);
                });

            assert(chordalGraph.isVertex(newParent));
            parents[vertex] = newParent;
        });
    }

    // Getter for the parallel customization
    inline Index getToAdjOfLevel(Index level) noexcept
    {
        assert(level < toAdjOfLevel.size());

        return toAdjOfLevel[level];
    }

    // returns the vertex stored at the index
    inline VertexID getVertexOfLevel(Index index) noexcept
    {
        assert(index < vertexOfLevel.size());

        return vertexOfLevel[index];
    }

    inline size_t numberOfLevels() noexcept { return toAdjOfLevel.size() - 1; }
    // ****************************************

private:
    // stores the original graph
    Graph::DynamicGraph core;

    // the chordal graph is a Directed Static Graph
    Graph::DiStaticGraph chordalGraph;

    // vector keeps track of which edges should be removed by the perfect
    // customization
    std::vector<uint8_t> removalFlags;

    // maps original edge ids to a tuple of <chordal edge ids, UPWARDS?>
    std::vector<std::pair<Index, bool>> edgeMapping;

    // the parent vector, used in the Elimination Tree Query, as well as in the
    // parallel section
    std::vector<VertexID> parents;

    // ****************************************
    // the idea of maybe parallel?
    // maps a level to the first vertex of this level
    std::vector<Index> toAdjOfLevel;
    // the index points on this vector
    std::vector<VertexID> vertexOfLevel;
};

} // namespace CCH
