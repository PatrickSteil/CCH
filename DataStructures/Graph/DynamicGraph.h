/*
 * File: ./DataStructures/Graph/DynamicGraph.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../../Extern/vector_io.h"
#include "../TypeDefs.h"
#include "Graphs.h"

namespace Graph {
// Dynamic Graph, used for fast insertion of edges.
// uses vector<vector<Edge>> as underlying datatstructure
// this graph also contains the reverse edges in order to iterate over incoming
// edges
class DynamicGraph {
public:
    // init with the given number of vertices
    DynamicGraph(size_t numVertices = 0)
        : outgoing(numVertices)
        , incoming(numVertices)
        , numEdges(0) {};

    // clear the whole graph
    // Note that the capacity does not decrease
    void clear()
    {
        for (VertexID vertex(0); vertex < numberOfNodes(); ++vertex) {
            outgoing[vertex].clear();
            incoming[vertex].clear();
        }
        numEdges = 0;
    }

    // after this method, the graph holds numVertices many vertices
    // Note: we erase previous stuff
    inline void addVertices(size_t numVertices = 0) noexcept
    {
        outgoing.assign(numVertices, {});
        incoming.assign(numVertices, {});
        numEdges = 0;
    }

    // method reserves *size* outgoing edges for the given vertex
    inline void reserveOutEdges(VertexID vertex, size_t size)
    {
        assert(isVertex(vertex));

        outgoing[vertex].reserve(size);
    }

    // method reserves *size* incoming edges for the given vertex
    inline void reserveInEdges(VertexID vertex, size_t size)
    {
        assert(isVertex(vertex));

        incoming[vertex].reserve(size);
    }

    // Adds the edge from -> to and from <- to to the graph even if the edge is
    // already present
    inline void addEdge(VertexID from, VertexID to, size_t weight = INFTY,
        Index edgeId = noIndex)
    {
        assert(isVertex(from) && isVertex(to));

        outgoing[from].emplace_back(to, weight, edgeId);
        incoming[to].emplace_back(from, weight, edgeId);

        ++numEdges;
    }

    // returns true iff from -> to is an edge
    inline bool hasEdge(VertexID from, VertexID to) const
    {
        assert(isVertex(from) && isVertex(to));

        return std::any_of(
            outgoing[from].begin(), outgoing[from].end(),
            [to](const Edge& edge) { return edge.otherVertex == to; });
    }

    // removes the edges from -> to and *not* from <- to
    inline void removeEdge(VertexID from, VertexID to)
    {
        assert(isVertex(from) && isVertex(to));

        auto it = std::find_if(outgoing[from].begin(), outgoing[from].end(),
            [to](const Edge& edge) { return edge.otherVertex == to; });

        if (it != outgoing[from].end()) {
            std::iter_swap(it, outgoing[from].end() - 1);
            outgoing[from].pop_back();
        }
    }

    // returns the edge from -> to
    inline Edge getEdge(VertexID from, VertexID to)
    {
        assert(isVertex(from) && isVertex(to));

        for (auto& edge : outgoing[from]) {
            if (edge.otherVertex == to)
                return edge;
        }

        return Edge(-1, -1);
    }

    // returns the handle to the incoming edges of this vertex
    inline std::vector<Edge>& getEdgesTo(VertexID vertex)
    {
        assert(isVertex(vertex));

        return incoming[vertex];
    }

    // returns the handle to the outgoing edges of this vertex
    inline std::vector<Edge>& getEdgesFrom(VertexID vertex)
    {
        assert(isVertex(vertex));

        return outgoing[vertex];
    }

    // if the edge from -> to already exists, it sets the edgeweight to new_weight
    // (and also modifies the backward edge) if the edge was not present, it will
    // be added with the new weight
    inline void setOrAddEdge(VertexID from, VertexID to, size_t new_weight,
        Index edgeIndex = noIndex)
    {
        assert(isVertex(from) && isVertex(to));

        for (size_t i(0); i < outgoing[from].size(); ++i) {
            Edge& edge = outgoing[from][i];
            if (edge.otherVertex == to) {
                if (edge.weight == new_weight)
                    return;

                edge.weight = new_weight;
                edge.edgeIndex = edgeIndex;

                for (auto& incomEdge : incoming[to]) {
                    if (incomEdge.otherVertex == from) {
                        incomEdge.weight = new_weight;
                        incomEdge.edgeIndex = edgeIndex;
                        break;
                    }
                }
                return;
            }
        }

        outgoing[from].emplace_back(to, new_weight, edgeIndex);

        ++numEdges;

        // also add reverse edge
        incoming[to].emplace_back(from, new_weight, edgeIndex);
    }

    // if the edge is already there, do nothing
    inline void findOrAddEdge(VertexID from, VertexID to)
    {
        assert(isVertex(from) && isVertex(to));

        // Check if the edge already exists
        for (auto& edge : outgoing[from]) {
            if (edge.otherVertex == to) {
                return; // Edge already exists
            }
        }

        outgoing[from].emplace_back(to, INFTY);

        ++numEdges;

        // also add reverse edge
        incoming[to].emplace_back(from, INFTY);
    }

    // Loops over the outgoing edges of vertex and applies the function to
    // fromVertex == vertex, toVertex, edge weight
    template <typename FUNC>
    inline void doForAllOutgoingEdges(VertexID vertex, FUNC function)
    {
        assert(isVertex(vertex));

        for (auto& edge : outgoing[vertex]) {
            function(vertex, edge.otherVertex, static_cast<size_t>(edge.weight));
        }
    }

    // Loops over the outgoing edges of vertex and applies the function to
    // fromVertex == vertex, toVertex, edge weight, edge index
    template <typename FUNC>
    inline void doForAllOutgoingEdgeIDs(VertexID vertex, FUNC function)
    {
        assert(isVertex(vertex));

        for (auto& edge : outgoing[vertex]) {
            function(vertex, edge.otherVertex, static_cast<size_t>(edge.weight),
                edge.edgeIndex);
        }
    }

    // Loops over the incoming edges of vertex and applies the function to
    // fromVertex, toVertex == vertex, edge weight
    template <typename FUNC>
    inline void doForAllIncomingEdges(VertexID vertex, FUNC function)
    {
        assert(isVertex(vertex));

        for (auto& edge : incoming[vertex]) {
            function(edge.otherVertex, vertex, static_cast<size_t>(edge.weight));
        }
    }

    // Loops over the incoming edges of vertex and applies the function to
    // fromVertex, toVertex == vertex, edge weight, edge index
    template <typename FUNC>
    inline void doForAllIncomingEdgeIDs(VertexID vertex, FUNC function)
    {
        assert(isVertex(vertex));

        for (auto& edge : incoming[vertex]) {
            function(edge.otherVertex, vertex, static_cast<size_t>(edge.weight),
                edge.edgeIndex);
        }
    }

    // isolates this vertex
    inline void isolateVertex(VertexID vertex)
    {
        assert(isVertex(vertex));

        for (auto& edge : outgoing[vertex]) {
            --numEdges;
            for (size_t i(0); i < incoming[edge.otherVertex].size(); ++i) {
                Edge& edgeIncom = incoming[edge.otherVertex][i];
                if (edgeIncom.otherVertex == vertex) {
                    std::swap(incoming[edge.otherVertex].back(),
                        incoming[edge.otherVertex][i]);
                    incoming[edge.otherVertex].pop_back();
                    break;
                }
            }
        }

        outgoing[vertex].clear();

        for (auto& edge : incoming[vertex]) {
            removeEdge(edge.otherVertex, vertex);
        }

        incoming[vertex].clear();
    }

    // sort all edges (inside each adjacency vector)
    inline void sortEdges() noexcept
    {
        for (VertexID vertex(0); vertex < numberOfNodes(); ++vertex) {
            sortEdges(vertex);
        }
    }

    // sort all edges for this particular vertex
    inline void sortEdges(VertexID vertex) noexcept
    {
        assert(isVertex(vertex));

        std::sort(outgoing[vertex].begin(), outgoing[vertex].end());
        std::sort(incoming[vertex].begin(), incoming[vertex].end());
    }

    // sort all edges (inside each adjacency vector)
    inline void sortEdgesParallel() noexcept
    {
#pragma omp parallel for schedule(dynamic)
        for (VertexID vertex = 0; vertex < numberOfNodes(); ++vertex) {
            sortEdges(vertex);
        }
    }
    // the next methods compute the degree of the given vertex
    inline size_t getOutDegree(VertexID vertex) const
    {
        assert(isVertex(vertex));
        return outgoing[vertex].size();
    }

    inline size_t getInDegree(VertexID vertex) const
    {
        assert(isVertex(vertex));

        return incoming[vertex].size();
    }

    inline size_t getDegree(VertexID vertex) const
    {
        assert(isVertex(vertex));

        return getOutDegree(vertex) + getInDegree(vertex);
    }

    inline std::set<VertexID> neighbours(VertexID vertex) const
    {
        assert(isVertex(vertex));

        std::set<VertexID> result;

        for (auto& edge : outgoing[vertex])
            result.insert(edge.otherVertex);

        for (auto& edge : incoming[vertex])
            result.insert(edge.otherVertex);

        return result;
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

        std::vector<std::vector<Edge>> newOutgoing(numberOfNodes());
        std::vector<std::vector<Edge>> newIncoming(numberOfNodes());

        // to maintain the prefix for toAdj
        size_t runningEdgeIndex(0);

        for (size_t i = 0; i < order.size(); ++i) {
            auto& vertex = order[i];

            // now for all edges, create new ones
            doForAllOutgoingEdges(
                vertex, [&](auto /* from */, auto toVertexOfEdge, auto weightOfEdge) {
                    auto mappedToVertex = mapping[toVertexOfEdge];
                    newOutgoing[i].emplace_back(mappedToVertex, weightOfEdge,
                        runningEdgeIndex);
                    newIncoming[mappedToVertex].emplace_back(i, weightOfEdge,
                        runningEdgeIndex);

                    ++runningEdgeIndex;
                });
        }

        // assign the newly created fields to attributes
        outgoing = newOutgoing;
        incoming = newIncoming;
    }

    // asserts
    inline bool isVertex(VertexID vertex) const
    {
        return vertex < outgoing.size();
    }

    inline bool isEdge(Edge edge) const { return isVertex(edge.otherVertex); }

    inline size_t numberOfNodes() const { return outgoing.size(); }
    inline size_t numberOfEdges() const { return numEdges; }

    // show some stats about the graph
    inline void showGraphStatistics() const
    {
        int numVertices = outgoing.size();
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
                degree0Vertices++;
            } else if (degree == 1) {
                degree1Vertices++;
            } else if (degree == 2) {
                degree2Vertices++;
            }
        }

        double averageDegree = static_cast<double>(numEdges) / numVertices;

        std::cout << "**** DynamicGraph ****" << std::endl;
        std::cout << "   Number of vertices:           " << numVertices
                  << std::endl;
        std::cout << "   Number of edges:              " << numEdges << std::endl;
        std::cout << "   Average degree:               " << averageDegree
                  << std::endl;
        std::cout << "   Minimum degree:               " << minDegree << std::endl;
        std::cout << "   Maximum degree:               " << maxDegree << std::endl;
        std::cout << "   Vertices with degree 0:       " << degree0Vertices
                  << std::endl;
        std::cout << "   Vertices with degree 1:       " << degree1Vertices
                  << std::endl;
        std::cout << "   Vertices with degree 2:       " << degree2Vertices
                  << std::endl;

        size_t vertexVectorSize = outgoing.capacity() * sizeof(std::vector<Edge>);
        vertexVectorSize <<= 1; // for incoming edges
        size_t edgeStructSize = numEdges * sizeof(Edge);
        size_t totalMemoryConsumption = vertexVectorSize + edgeStructSize;

        std::cout << "   ~ total memory consumption:   " << totalMemoryConsumption
                  << std::endl;
    }

    // copies the graph from another dynamicgraph
    inline void copy(DynamicGraph& other)
    {
        outgoing.clear();
        incoming.clear();

        outgoing = other.outgoing;
        incoming = other.incoming;
        numEdges = other.numEdges;
    }

    inline void readDimacsFile(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        clear();

        std::string line;
        int numVertices = 0, numEdges = 0;

        while (std::getline(file, line)) {
            if (line.compare(0, 5, "p sp ") == 0) {
                std::istringstream iss(line.substr(5));
                iss >> numVertices >> numEdges;
                break;
            }
        }

        outgoing.resize(numVertices);
        incoming.resize(numVertices);

        for (int i(0); i < numVertices; ++i) {
            outgoing[i].reserve((int)(numEdges / (float)numVertices));
            incoming[i].reserve((int)(numEdges / (float)numVertices));
        }

        // Read edge lines and add edges to the graph
        while (std::getline(file, line)) {
            if (line[0] == 'a') {
                std::istringstream iss(line);
                char edgeType;
                VertexID source, destination;
                size_t weight;
                iss >> edgeType >> source >> destination >> weight;
                addEdge(source - 1, destination - 1, weight);
            }
        }

        file.close();
    }

    // read from binary (the RoutingKIT format)
    inline void readBinary(const std::string& fileName,
        const std::string& weightName) noexcept
    {
        std::cout << CLR_RED << "Reading Dynamic Graph" << CLR_RESET << " from "
                  << fileName << std::endl;
        clear();

        std::vector<unsigned> first_out = load_vector<unsigned>(fileName + "first_out");
        std::vector<unsigned> head = load_vector<unsigned>(fileName + "head");
        std::vector<unsigned> weight = load_vector<unsigned>(fileName + weightName);

        size_t numEdges = head.size();
        size_t numVertices = first_out.size() - 1;

        addVertices(numVertices);

        for (size_t i(0); i < numVertices; ++i) {
            reserveOutEdges(i, first_out[i + 1] - first_out[i]);

            for (Index j(first_out[i]); j < first_out[i + 1]; ++j) {
                assert(j < numEdges);
                addEdge(VertexID(i), VertexID(head[j]), weight[j], j);
            }
        }
        sortEdges();
    }

    inline void reset() noexcept
    {
        outgoing.clear();
        incoming.clear();
        numEdges = 0;
    }

private:
    std::vector<std::vector<Edge>> outgoing;
    std::vector<std::vector<Edge>> incoming;
    size_t numEdges;
};
} // namespace Graph
