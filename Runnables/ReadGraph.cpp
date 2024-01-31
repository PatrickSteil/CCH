/*
 * File: ./Runnables/ReadGraph.cpp
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#include <iostream>
#include <vector>

#include "../Algorithms/CCH/Customization.h"
#include "../Algorithms/CCH/EliminationTreeQuery.h"
#include "../Algorithms/CCH/Preprocessing.h"
#include "../Algorithms/CCH/Verifier.h"
#include "../DataStructures/CCH/Data.h"
#include "../DataStructures/Graph/Graphs.h"
#include "../Helpers/Colors.h"

int main()
{
    Graph::DynamicGraph graph;

    graph.readBinary("../DataSets/stupferich/", "travel_time");
    /* graph.readBinary("../DataSets/algoDaten/praktikum/graph/karlsruhe/", */
    /*                  "travel_time"); */
    graph.showGraphStatistics();

    std::vector<unsigned> order = load_vector<unsigned>("../DataSets/stupferich/flowcutter_order");
    graph.reorderNodes(order);

    CCH::Data data(graph, 6);

    CCH::Preprocessing<Graph::DynamicGraph, true, true> builder(data);
    builder.runContraction();

    // std::cout << "Stats for the chordal graph building" << std::endl;
    builder.getStatsCollecter().writeStatsToCSV(
        "../DataSets/stupferich/stats.chordalGraphBuilding");

    std::cout << "Stats of the chordal graph" << std::endl;
    data.getChordalGraph().showGraphStatistics();

    data.buildLevelDataStructure();

    data.computeEdgeMapping(graph);

    CCH::Customization custo(data);
    custo.applyMetric(graph);

    custo.runCustomization();

    custo.runPerfectCustomization();

    // std::cout << "Stats for the customization" << std::endl;
    custo.getStatsCollecter().writeStatsToCSV(
        "../DataSets/stupferich/stats.customization");

    CCH::Verifier verifier(data, graph);
    verifier.verify();

    return 0;
}
