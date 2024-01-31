/*
 * File: ./DataStructures/Graph/DiStaticGraph.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "../../Extern/vector_io.h"
#include "../../Helpers/Colors.h"
#include "../../Helpers/TypeChecker.h"
#include "../TypeDefs.h"
#include "DynamicGraph.h"

// this class is used for the chordal graph in the CCH::Data class
// it stores edge only once (namely the edge {x, y} is stored at the lower
// ranked vertex, but the edge index is stored at y as an incoming edge)
// ** Attributes: **
// toAdj:             toAdj[x] holds the starting outgoing edge index of x
// toVertex:          toVertex[edgeIndex] stores the toVertex of the edge
// edgeIndex
// fromVertex:        fromVertex[edgeIndex] stores the fromVertex of the edge
// edgeIndex
// upWeight:          upWeight[edgeIndex] stores the upWeight of the edge
// edgeIndex
// downWeight:        downWeight[edgeIndex] stores the downWeight of the edge
// edgeIndex
// incomingEdgeIDs:  incomingEdgeIDs[x] stores the incoming edge idsof vertex x
// (note: this is a vector<vector<>>)
namespace Graph {

class DiStaticGraph {
public:
    explicit DiStaticGraph(size_t vertices = 0)
    {
        toAdj.resize(vertices + 1);
        incomingEdgeIDs.assign(vertices, {});
    }

    inline void reserve(size_t newNumberOfVertices, size_t newNumberOfEdges)
    {
        toAdj.reserve(newNumberOfVertices + 1);
        incomingEdgeIDs.reserve(newNumberOfVertices);
        toVertex.reserve(newNumberOfEdges);
        fromVertex.reserve(newNumberOfEdges);
        upWeight.reserve(newNumberOfEdges);
        downWeight.reserve(newNumberOfEdges);
    }

    inline void copy(DynamicGraph& from)
    {
        clear();
        reserve(from.numberOfNodes(), from.numberOfEdges());
        incomingEdgeIDs.assign(from.numberOfNodes(), {});

        int prefix(0);
        for (int i(0); i < (int)from.numberOfNodes(); ++i) {
            int deg = from.getOutDegree(i);

            toAdj.push_back(prefix);
            prefix += deg;

            from.doForAllOutgoingEdges(
                i, [&](auto /* source */, auto target, auto /* weight */) {
                    incomingEdgeIDs[target].push_back(Index(toVertex.size()));

                    toVertex.push_back(target);
                    fromVertex.push_back(i);
                    // by default, ignore all the weights
                    upWeight.push_back(INFTY);
                    downWeight.push_back(INFTY);
                });
        }

        toAdj.push_back(prefix);
    }

    // from is a vector<vector<>>, and numberOfEdges is the number of edges
    // represented in from
    inline void copy(const std::vector<std::vector<VertexID>>& from,
        size_t numberOfEdges)
    {
        clear();
        reserve(from.size(), numberOfEdges);
        incomingEdgeIDs.assign(from.size(), {});

        int prefix(0);
        for (int i(0); i < (int)from.size(); ++i) {
            int deg = from[i].size();

            toAdj.push_back(prefix);
            prefix += deg;

            for (auto& target : from[i]) {
                incomingEdgeIDs[target].push_back(Index(toVertex.size()));

                toVertex.push_back(target);
                fromVertex.push_back(i);
                // by default, ignore all the weights
                upWeight.push_back(INFTY);
                downWeight.push_back(INFTY);
            };
        }

        toAdj.push_back(prefix);
    }

    // Helpers
    inline bool isVertex(VertexID vertex) const
    {
        return ((size_t)vertex < numberOfNodes());
    }

    inline bool isEdge(Index edge) const
    {
        return ((size_t)edge < numberOfEdges());
    }

    inline bool hasEdge(VertexID from, VertexID to) const
    {
        assert(isVertex(from));
        assert(isVertex(to));

        for (size_t i(toAdj[from]); i < (size_t)toAdj[from + 1]; ++i) {
            if (toVertex[i] == to)
                return true;
        }

        return false;
    }

    inline Index findEdge(VertexID from, VertexID to) const noexcept
    {
        for (size_t i(toAdj[from]); i < (size_t)toAdj[from + 1]; ++i) {
            if (toVertex[i] == to)
                return i;
        }
        return noIndex;
    }

    // getter and hidden setter
    // the following methods can be used without argument => it returns a
    // reference to the complete attribute or by passing either VertexID or Index
    // to return the specified attribute for the given vertex / edge

    inline std::vector<Index>& getToAdj() noexcept { return toAdj; }

    inline Index getToAdj(VertexID vertex) noexcept
    {
        assert(isVertex(vertex) || vertex == numberOfNodes());
        return toAdj[vertex];
    }

    inline std::vector<VertexID>& getToVertex() noexcept { return toVertex; }

    inline VertexID& getToVertex(Index edge) noexcept
    {
        assert(isEdge(edge));
        return toVertex[edge];
    }

    inline std::vector<VertexID>& getFromVertex() noexcept { return fromVertex; }

    inline VertexID& getFromVertex(Index edge) noexcept
    {
        assert(isEdge(edge));
        return fromVertex[edge];
    }

    inline std::vector<size_t>& getUpWeight() noexcept { return upWeight; }

    inline size_t& getUpWeight(Index edge) noexcept
    {
        assert(isEdge(edge));
        return upWeight[edge];
    }

    inline std::vector<size_t>& getDownWeight() noexcept { return downWeight; }

    inline size_t& getDownWeight(Index edge) noexcept
    {
        assert(isEdge(edge));
        return downWeight[edge];
    }

    std::vector<std::vector<Index>>& getIncomingEdgeIDs() noexcept
    {
        return incomingEdgeIDs;
    }

    std::vector<Index>& getIncomingEdgeIDs(VertexID vertex) noexcept
    {
        assert(isVertex(vertex));
        return incomingEdgeIDs[vertex];
    }

    // degrees, outgoing, incoming and sum of both
    inline size_t getOutDegree(VertexID vertex) const
    {
        assert(isVertex(vertex));
        return toAdj[vertex + 1] - toAdj[vertex];
    }

    inline size_t getInDegree(VertexID vertex) const
    {
        assert(isVertex(vertex));
        return incomingEdgeIDs[vertex].size();
    }

    inline size_t getDegree(VertexID vertex) const
    {
        assert(isVertex(vertex));
        return getOutDegree(vertex) + getInDegree(vertex);
    }

    // loop over all outgoing edges of vertex and apply function func(vertex,
    // toVertex, upWeight, downWeight, edge id)
    template <typename FUNC>
    inline void doForAllOutgoingEdgeIDs(VertexID vertex, FUNC func)
    {
        assert(isVertex(vertex));

        for (size_t i(toAdj[vertex]); i < (size_t)toAdj[vertex + 1]; ++i) {
            func(vertex, toVertex[i], upWeight[i], downWeight[i], i);
        }
    }

    // loop over all incoming edges of vertex and apply function func(fromVertex,
    // vertex, upWeight, downWeight, edge id)
    template <typename FUNC>
    inline void doForAllIncomingEdgeIDs(VertexID vertex, FUNC func)
    {
        assert(isVertex(vertex));

        for (size_t i(0); i < (size_t)getInDegree(vertex); ++i) {
            size_t edgeId = incomingEdgeIDs[vertex][i];
            func(fromVertex[edgeId], vertex, upWeight[edgeId], downWeight[edgeId],
                edgeId);
        }
    }

    // permutate nodes
    // order is a vector of VertexID, and the idea is that order[i] is the new
    // i'th node
    inline void reorderNodes(const std::vector<VertexID>& order) noexcept
    {
        assert(order.size() == numberOfNodes());
        assert(std::all_of(std::begin(order), std::end(order),
                   [this](int x) { return isVertex(x); })
            && "Not all elements fulfill the condition isVertex(x)");

        // compute reverse mapping, so that mapping[vertex] = index in order (== new
        // Vertex id)
        std::vector<VertexID> mapping(numberOfNodes());

        for (size_t i = 0; i < numberOfNodes(); ++i) {
            mapping[order[i]] = i;
        }

        // create all the new fields
        std::vector<Index> newToAdj(numberOfNodes() + 1, 0);
        std::vector<std::vector<Index>> newIncomingEdgeIDs(numberOfNodes());
        std::vector<VertexID> newToVertex(numberOfEdges(), 0);
        std::vector<VertexID> newFromVertex(numberOfEdges(), 0);
        std::vector<size_t> newUpWeight(numberOfEdges(), 0);
        std::vector<size_t> newDownWeight(numberOfEdges(), 0);

        // to maintain the prefix for toAdj
        size_t runningEdgeIndex(0);

        for (size_t i = 0; i < order.size(); ++i) {
            auto& vertex = order[i];

            newToAdj[i] = runningEdgeIndex;

            // now for all edges, create new ones
            doForAllOutgoingEdgeIDs(
                vertex, [&](auto /* from */, auto toVertexOfEdge, auto upWeightOfEdge, auto downWeightOfEdge, auto /* edgeIndex */) {
                    auto mappedToVertex = mapping[toVertexOfEdge];

                    newToVertex[runningEdgeIndex] = mappedToVertex;
                    newFromVertex[runningEdgeIndex] = i;
                    newUpWeight[runningEdgeIndex] = upWeightOfEdge;
                    newDownWeight[runningEdgeIndex] = downWeightOfEdge;
                    newIncomingEdgeIDs[mappedToVertex].push_back(runningEdgeIndex);

                    ++runningEdgeIndex;
                });
        }
        newToAdj.back() = runningEdgeIndex;

        clear();
        // assign the newly created fields to attributes
        toAdj = newToAdj;
        incomingEdgeIDs = newIncomingEdgeIDs;
        fromVertex = newFromVertex;
        toVertex = newToVertex;
        upWeight = newUpWeight;
        downWeight = newDownWeight;
    }

    // now some simple methods
    inline size_t numberOfNodes() const { return toAdj.size() - 1; }

    inline size_t numberOfEdges() const { return toVertex.size(); }

    inline void clear()
    {
        toAdj.clear();
        incomingEdgeIDs.clear();
        toVertex.clear();
        fromVertex.clear();
        upWeight.clear();
        downWeight.clear();
    }

    inline void showGraphStatistics() const
    {
        int numVertices = numberOfNodes();
        int numEdges = numberOfEdges();
        int degree0Vertices = 0;
        int degree1Vertices = 0;
        int degree2Vertices = 0;
        int minDegree = numVertices;
        int maxDegree = 0;

        for (int i = 0; i < numVertices; ++i) {
            int degree = getDegree(i);
            minDegree = std::min(minDegree, degree);
            maxDegree = std::max(maxDegree, degree);

            if (degree == 0) {
                ++degree0Vertices;
            } else if (degree == 1) {
                ++degree1Vertices;
            } else if (degree == 2) {
                ++degree2Vertices;
            }
        }

        double averageDegree = static_cast<double>(numberOfEdges()) / numVertices;

        std::cout << "**** DiStaticGraph ****" << std::endl;
        std::cout << "   Number of vertices:           " << numVertices
                  << std::endl;
        std::cout << "   Number of edges:              " << numEdges << std::endl;
        std::cout << "   Average degree:               " << averageDegree
                  << std::endl;
        std::cout << "   Minimum degree:               " << minDegree << std::endl;
        std::cout << "   Maximum degree:               " << maxDegree << std::endl;
        std::cout << "   Vertices with degree 0:    " << degree0Vertices
                  << std::endl;
        std::cout << "   Vertices with degree 1:    " << degree1Vertices
                  << std::endl;
        std::cout << "   Vertices with degree 2:    " << degree2Vertices
                  << std::endl;

        size_t totalMemoryConsumption = sizeof(int) * (toAdj.size());
        totalMemoryConsumption += sizeof(VertexID) * (toVertex.size() + fromVertex.size());
        totalMemoryConsumption += sizeof(size_t) * (upWeight.size() + downWeight.size());
        totalMemoryConsumption += sizeof(Index) * toVertex.size(); // the memory consumption for the incoming edgeIDs

        std::cout << "   ~ total memory consumption:   " << totalMemoryConsumption
                  << std::endl;
    }

    // read from binary (the RoutingKIT format)
    inline void readBinary(const std::string& fileName,
        const std::string& weightName) noexcept
    {
        std::cout << CLR_RED << "Reading DiStatic Graph" << CLR_RESET << " from "
                  << fileName << std::endl;
        clear();

        std::vector<unsigned> first_out = load_vector<unsigned>(fileName + "first_out");
        std::vector<unsigned> head = load_vector<unsigned>(fileName + "head");
        std::vector<unsigned> weight = load_vector<unsigned>(fileName + weightName);

        size_t numEdges = head.size();
        size_t numVertices = first_out.size() - 1;

        toAdj.resize(numVertices + 1);
        toVertex.resize(numEdges);
        fromVertex.resize(numEdges);
        upWeight.resize(numEdges);
        downWeight.resize(numEdges);
        incomingEdgeIDs.assign(numVertices, {});

        for (size_t i(0); i < numVertices; ++i) {
            toAdj[i] = Index(first_out[i]);

            for (Index j(first_out[i]); j < first_out[i + 1]; ++j) {
                assert(j < numEdges);
                toVertex[j] = VertexID(head[j]);
                fromVertex[j] = VertexID(i);
                upWeight[j] = VertexID(weight[j]);
                downWeight[j] = VertexID(weight[j]);

                assert(head[j] < incomingEdgeIDs.size());
                incomingEdgeIDs[toVertex[j]].push_back(j);
            }
        }

        toAdj.back() = first_out.back();
    }

private:
    // maps vertex to first outgoing edge
    std::vector<Index> toAdj;

    // maps edge to it's head
    std::vector<VertexID> toVertex;

    // maps edge to it's tail
    std::vector<VertexID> fromVertex;

    // maps edge to it's weights (up && down)
    std::vector<size_t> upWeight;
    std::vector<size_t> downWeight;

    // maps vertex to a vector of incoming edge ids
    std::vector<std::vector<Index>> incomingEdgeIDs;
};

} // namespace Graph
