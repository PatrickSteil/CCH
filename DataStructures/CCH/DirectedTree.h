/*
 * File: ./DataStructures/CCH/DirectedTree.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <functional>
#include <iostream>
#include <queue>
#include <stack>
#include <vector>

#include "../../DataStructures/TypeDefs.h"

class DirectedTree {
public:
    // Constructor that takes the parent information of each vertex
    DirectedTree(const std::vector<VertexID>& parents)
        : parents_(parents)
        , children_(parents.size())
    {
        // Build children information based on parent-child relationships
        for (VertexID v = 0; v < parents_.size(); ++v) {
            if (parents_[v] != v) {
                children_[parents_[v]].push_back(v);
            }
        }
    }

    // Function to compute the topological ordering using DFS
    std::vector<VertexID> computeTopologicalOrdering()
    {
        std::vector<bool> visited(parents_.size(), false);
        std::stack<VertexID> stack;
        std::vector<VertexID> topOrder;

        // DFS function to explore the tree and fill the stack
        std::function<void(VertexID)> dfs = [&](VertexID v) {
            visited[v] = true;

            // Explore all children of the current vertex
            for (VertexID child : children_[v]) {
                if (!visited[child]) {
                    dfs(child);
                }
            }

            // Push the current vertex to the stack after processing its children
            stack.push(v);
        };

        // Perform DFS starting from each unvisited vertex
        for (VertexID v = 0; v < parents_.size(); ++v) {
            if (!visited[v]) {
                dfs(v);
            }
        }

        // Pop elements from the stack to obtain the topological ordering
        while (!stack.empty()) {
            topOrder.push_back(stack.top());
            stack.pop();
        }

        return topOrder;
    }

    // Function to compute levels in the tree
    std::vector<int> computeLevels()
    {
        std::vector<int> levels(parents_.size(),
            -1); // Initialize levels to -1 (not visited)
        std::queue<VertexID> q;

        // Start BFS from the root (vertex with parent pointing to itself)
        for (VertexID v = 0; v < parents_.size(); ++v) {
            if (parents_[v] == v) {
                q.push(v);
                levels[v] = 0; // Root is at level 0
                break;
            }
        }

        // Perform BFS to assign levels
        while (!q.empty()) {
            VertexID current = q.front();
            q.pop();

            for (VertexID child : children_[current]) {
                if (levels[child] == -1) {
                    // If the child has not been visited, assign the next level and
                    // enqueue
                    levels[child] = levels[current] + 1;
                    q.push(child);
                }
            }
        }

        return levels;
    }

private:
    const std::vector<VertexID>& parents_; // Parent information for each vertex
    std::vector<std::vector<VertexID>>
        children_; // Children information for each vertex
};
