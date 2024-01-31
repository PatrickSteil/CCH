/*
 * File: ./DataStructures/TypeDefs.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <stdint.h>

#include <limits>

// to not always write size_t and it is clearer what it means
typedef size_t Index;

constexpr size_t noIndex = size_t(-1);

// Vertex as for a node id in a graph
typedef uint32_t VertexID;

constexpr uint32_t noVertexID = uint32_t(-1);

// Edge Constants like infinity or never
constexpr uint32_t INFTY = std::numeric_limits<uint32_t>::max() / 2;
constexpr uint32_t never = INFTY;

// NOOP as a template parameter
struct NO_OPERATION {
    template <typename... ARGS>
    constexpr inline bool operator()(ARGS...) const noexcept
    {
        return false;
    }
};

NO_OPERATION NoOperation;
