/*
 * File: ./Runnables/CCHTest.cpp
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#include <iostream>
#include <vector>

#include "../Algorithms/CCH/Customization.h"
#include "../Algorithms/CCH/DijkstraQuery.h"
#include "../Algorithms/CCH/EliminationTreeQuery.h"
#include "../Algorithms/CCH/Preprocessing.h"
#include "../Algorithms/CCH/Verifier.h"
#include "../DataStructures/CCH/Data.h"
#include "../DataStructures/Graph/Graphs.h"
#include "../Helpers/Colors.h"

// Example from here
// https://i11www.iti.kit.edu/_media/teaching/sommer2023/routenplanung/chap1-cch.pdf

int main()
{
    Graph::DynamicGraph graph(7);
    graph.addEdge(0, 6, 10, 0);
    graph.addEdge(6, 0, 10, 1);
    graph.addEdge(0, 5, 1, 2);
    graph.addEdge(5, 0, 1, 3);
    graph.addEdge(0, 3, 1, 4);
    graph.addEdge(3, 0, 1, 5);

    graph.addEdge(2, 6, 200, 6);
    graph.addEdge(6, 2, 2, 7);

    graph.addEdge(1, 6, 3, 8);
    graph.addEdge(6, 1, 3, 9);
    graph.addEdge(1, 4, 1, 10);
    graph.addEdge(4, 1, 1, 11);

    graph.addEdge(3, 4, 1, 12);
    graph.addEdge(4, 3, 100, 13);

    graph.showGraphStatistics();

    std::vector<VertexID> order = { 0, 1, 2, 3, 4, 5, 6 };

    graph.reorderNodes(order);

    CCH::Data data(graph);

    CCH::Preprocessing<Graph::DynamicGraph, true, true> builder(data);
    builder.runContraction();

    builder.showStats();

    std::cout << "Chordal Graph:" << std::endl;
    data.getChordalGraph().showGraphStatistics();

    data.computeEdgeMapping(graph);

    for (Index i(0); i < data.getChordalGraph().numberOfEdges(); ++i) {
        std::cout << i << "," << (data.getChordalGraph().getFromVertex(i) + 1)
                  << "," << (data.getChordalGraph().getToVertex(i) + 1) << ","
                  << data.getChordalGraph().getUpWeight(i) << ","
                  << data.getChordalGraph().getDownWeight(i) << "\t"
                  << (int)data.getRemovalFlags()[i] << std::endl;
    }

    std::cout << "After applying metric:\n";

    CCH::Customization<Graph::DynamicGraph, true, true> custo(data);
    custo.applyMetric(graph);

    custo.runCustomization();

    custo.runPerfectCustomization();

    custo.showStats();

    // data.reorderChrodalByOrder();

    for (Index i(0); i < data.getChordalGraph().numberOfEdges(); ++i) {
        std::cout << i << "," << (data.getChordalGraph().getFromVertex(i) + 1)
                  << "," << (data.getChordalGraph().getToVertex(i) + 1) << ","
                  << data.getChordalGraph().getUpWeight(i) << ","
                  << data.getChordalGraph().getDownWeight(i) << "\t"
                  << (int)data.getRemovalFlags()[i] << std::endl;
    }

    CCH::Verifier verifier(data, graph);
    verifier.verify();
}
