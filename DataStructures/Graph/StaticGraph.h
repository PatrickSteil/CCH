/*
 * File: ./DataStructures/Graph/StaticGraph.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../Helpers/TypeChecker.h"
#include "../TypeDefs.h"

namespace Graph {

class StaticGraph {
public:
    StaticGraph() { }
    StaticGraph(size_t vertices) { toAdj.assign(vertices + 1, 0); }

    void reserve(size_t newNumberOfVertices, size_t newNumberOfEdges)
    {
        toAdj.reserve(newNumberOfVertices + 1);
        toVertex.reserve(newNumberOfEdges);
        weight.reserve(newNumberOfEdges);
    }

    // getter && setter
    std::vector<Index>& getToAdj() { return toAdj; }
    std::vector<VertexID>& getToVertex() { return toVertex; }
    std::vector<size_t>& getWeight() { return weight; }

    // Helpers
    bool isVertex(VertexID vertex) const
    {
        return ((size_t)vertex < numberOfNodes());
    }

    bool hasEdge(VertexID from, VertexID to) const
    {
        assert(isVertex(from));
        assert(isVertex(to));

        for (size_t i(toAdj[from]); i < (size_t)toAdj[from + 1]; ++i) {
            if (toVertex[i] == to)
                return true;
        }

        return false;
    }

    size_t getOutDegree(VertexID vertex) const
    {
        assert(isVertex(vertex));
        return toAdj[vertex + 1] - toAdj[vertex];
    }

    size_t getDegree(VertexID vertex) const
    {
        assert(isVertex(vertex));
        return getOutDegree(vertex);
    }

    int getWeight(int edge) const
    {
        assert(0 <= edge && (size_t)edge < toVertex.size());
        return weight[edge];
    }

    // func takes 4 parameters (vertex, toVertex[i], weight[i], i)
    template <typename FUNC>
    void doForAllOutgoingEdges(VertexID vertex, FUNC func)
    {
        assert(isVertex(vertex));

        for (size_t i(toAdj[vertex]); i < (size_t)toAdj[vertex + 1]; ++i) {
            func(vertex, toVertex[i], weight[i], i);
        }
    }

    std::vector<int> getOutgoingNeighbours(VertexID vertex)
    {
        assert(isVertex(vertex));

        std::vector<int> result;
        result.reserve(getOutDegree(vertex));

        doForAllOutgoingEdges(vertex,
            [&](auto /* from */, auto to, auto /* weight */,
                auto /* edgeId */) { result.push_back(to); });
        return result;
    }

    size_t numberOfNodes() const { return toAdj.size() - 1; }

    size_t numberOfEdges() const { return toVertex.size() - 1; }

    void clear()
    {
        toAdj.clear();
        toVertex.clear();
        weight.clear();
    }

    void showGraphStatistics() const
    {
        int numVertices = numberOfNodes();
        int numEdges = toVertex.size();
        int degree0Vertices = 0;
        int degree1Vertices = 0;
        int degree2Vertices = 0;
        int minDegree = numVertices;
        int maxDegree = 0;

        for (int i = 0; i < numVertices; ++i) {
            int degree = getOutDegree(i);
            minDegree = std::min(minDegree, degree);
            maxDegree = std::max(maxDegree, degree);

            if (degree == 0) {
                degree0Vertices++;
            } else if (degree == 1) {
                degree1Vertices++;
            } else if (degree == 2) {
                degree2Vertices++;
            }
        }

        double averageDegree = static_cast<double>(numberOfEdges()) / numVertices;

        std::cout << "**** StaticGraph ****" << std::endl;
        std::cout << "   Number of vertices:           " << numVertices
                  << std::endl;
        std::cout << "   Number of edges:              " << numEdges << std::endl;
        std::cout << "   Average degree:               " << averageDegree
                  << std::endl;
        std::cout << "   Minimum degree:               " << minDegree << std::endl;
        std::cout << "   Maximum degree:               " << maxDegree << std::endl;
        std::cout << "   Vertices with outdegree 0:    " << degree0Vertices
                  << std::endl;
        std::cout << "   Vertices with outdegree 1:    " << degree1Vertices
                  << std::endl;
        std::cout << "   Vertices with outdegree 2:    " << degree2Vertices
                  << std::endl;

        size_t totalMemoryConsumption = sizeof(int) * (toAdj.size() + toVertex.size() + weight.size());

        std::cout << "   ~ total memory consumption:   " << totalMemoryConsumption
                  << std::endl;
    }

private:
    std::vector<Index> toAdj;
    std::vector<VertexID> toVertex;
    std::vector<size_t> weight;
};

} // namespace Graph