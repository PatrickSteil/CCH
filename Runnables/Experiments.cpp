/*
 * File: ./Runnables/Experiments.cpp
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#include <iostream>
#include <thread>
#include <vector>

#include "../Algorithms/CCH/Customization.h"
#include "../Algorithms/CCH/EliminationTreeQuery.h"
#include "../Algorithms/CCH/Preprocessing.h"
#include "../Algorithms/CCH/Verifier.h"
#include "../DataStructures/CCH/Data.h"
#include "../DataStructures/Graph/Graphs.h"
#include "../Helpers/Colors.h"

int main(int argc, char* argv[])
{
    std::string graph_name;
    std::string metric;
    std::string order_name;
    if (argc != 2) {
        std::cout << "Wrong argument count." << std::endl;
        return 1;
    } else {
        graph_name = argv[1];
        // metric = argv[2];
        metric = "travel_time";
        // order_name = argv[3];
        order_name = "flowcutter_order";
    }
    Graph::DynamicGraph graph;

    graph.readBinary("../../Datasets/" + graph_name + "/", metric);

    graph.showGraphStatistics();

    std::vector<unsigned> order = load_vector<unsigned>("../../Datasets/" + graph_name + "/" + order_name);
    graph.reorderNodes(order);

    int NUM_THREADS = std::thread::hardware_concurrency();
    CCH::Data data(graph, NUM_THREADS);

    CCH::Preprocessing builder(data);
    builder.runContraction();

    // std::cout << "Stats for the chordal graph building" << std::endl;
    builder.getStatsCollecter().writeStatsToCSV("../../Datasets/" + graph_name + "/stats.chordalGraphBuilding");

    std::cout << "Stats of the chordal graph" << std::endl;
    data.getChordalGraph().showGraphStatistics();

    data.buildLevelDataStructure();
    data.computeEdgeMapping(graph);

    CCH::Customization custo(data);
    custo.applyMetric(graph);

    custo.runCustomizationParallel();

    custo.runPerfectCustomizationParallel();

    // std::cout << "Stats for the customization" << std::endl;
    custo.getStatsCollecter().writeStatsToCSV("../../Datasets/" + graph_name + "/stats.customization");

    CCH::Verifier verifier(data, graph);
    verifier.verify();

    verifier.getEliminationTreeQuery().getStatsCollecter().writeStatsToCSV(
        "../../Datasets/" + graph_name + "/stats.elimination");
    verifier.getCCHDijkstraQuery().getStatsCollecter().writeStatsToCSV(
        "../../Datasets/" + graph_name + "/stats.bidijkstra");
    // verifier.getDijkstraQuery().getStatsCollecter().writeStatsToCSV(
    //   "../Datasets/" + graph_name + "/stats.dijkstra");

    return 0;
}
