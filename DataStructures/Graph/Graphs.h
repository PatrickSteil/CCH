/*
 * File: ./DataStructures/Graph/Graphs.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once
#include <tuple>

#include "../TypeDefs.h"

namespace Graph {
// This struct holds the information about one edge in a graph.
// Note that we do not define toVertex/fromVertex, but only otherVertex.
// This way, we can use this struct as both a forward and outgoing edge.

struct Edge {
    // Constructor for initializing Edge with optional parameters.
    Edge(VertexID otherVertex = 0, size_t weight = 0, Index edgeIndex = noIndex)
        : otherVertex(otherVertex)
        , weight(weight)
        , edgeIndex(edgeIndex)
    {
    }

    // The vertex connected by this edge.
    VertexID otherVertex;

    // The weight or cost associated with this edge.
    size_t weight;

    // Index of the edge in the graph structure. Default is noIndex.
    Index edgeIndex;

    // Comparison operator for sorting edges based on otherVertex and weight.
    bool operator<(const Edge& other) const
    {
        return std::tie(otherVertex, weight) < std::tie(other.otherVertex, other.weight);
    }
};
} // namespace Graph

#include "DiStaticGraph.h"
#include "DynamicGraph.h"
#include "StaticGraph.h"