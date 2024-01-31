/*
 * File: ./Algorithms/CCH/Verifier.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <random>
#include <vector>

#include "../../Extern/compare_vector.h"
#include "../../Helpers/ProgressBar.h"
#include "../Dijkstra/Dijkstra.h"
#include "DijkstraQuery.h"
#include "EliminationTreeQuery.h"

namespace CCH {

const size_t queryAmount = 500;

template <typename GRAPH>
class Verifier {
public:
    Verifier(Data<GRAPH>& data, GRAPH& originalGraph)
        : cchData(data)
        , originalGraph(originalGraph)
        , elimTreeQuery(data)
        , cchDijkstraQuery(data)
        , dijkstraQuery(originalGraph)
        , sources(queryAmount)
        , targets(queryAmount)
        , cchResults(queryAmount)
        , dijkstraResults(queryAmount)
    {
    }

    inline void generateSourcesAndTargets()
    {
        std::random_device dev;
        std::mt19937 randomGenerator(42);
        std::uniform_int_distribution<std::mt19937::result_type>
            sourcesDistribution(0, cchData.getChordalGraph().numberOfNodes() - 1);
        std::uniform_int_distribution<std::mt19937::result_type>
            targetsDistribution(0, cchData.getChordalGraph().numberOfNodes() - 1);
        for (size_t i = 0; i < queryAmount; i++) {
            sources[i] = sourcesDistribution(randomGenerator);
            targets[i] = targetsDistribution(randomGenerator);
        }
    }

    inline void verify()
    {
        std::cout << "Verifiying the CCH Query against Dijkstra" << std::endl;
        generateSourcesAndTargets();
        progressbar progressElim(queryAmount);
        progressbar progressCCHDijkstra(queryAmount);
        progressbar progressDijkstra(queryAmount);
        std::cout << CLR_BLUE << "[START]" << CLR_RESET << " Query Verifier"
                  << std::endl;

        std::cout << CLR_BLUE << "[START]" << CLR_RESET << " Elim Tree Query"
                  << std::endl;
        for (size_t i = 0; i < queryAmount; i++) {
            cchResults[i] = elimTreeQuery.runWithPathUnpacking(sources[i], targets[i]);
            progressElim.update();
        }
        std::cout << std::endl
                  << CLR_RED << "[DONE]" << CLR_RESET << " Elim Tree Query"
                  << std::endl;

        std::cout << CLR_BLUE << "[START]" << CLR_RESET << " Dijkstra Query"
                  << std::endl;

        for (size_t i = 0; i < queryAmount; i++) {
            cchDijkstraQuery.run(sources[i], targets[i]);
            dijkstraResults[i] = cchDijkstraQuery.getDistance();
            progressCCHDijkstra.update();
        }
        /* for (size_t i = 0; i < queryAmount; i++) {
          dijkstraQuery.run(sources[i], targets[i]);
          dijkstraResults[i] = dijkstraQuery.getDistance(targets[i]);
          progressDijkstra.update();
        } */
        std::cout << std::endl
                  << CLR_RED << "[DONE]" << CLR_RESET << " Dijkstra Query"
                  << std::endl;

        std::cout << CLR_RED << "[DONE]" << CLR_RESET << " Query Verifier"
                  << std::endl;
        // elimTreeQuery.showStats();
        // cchDijkstraQuery.showStats();
        // dijkstraQuery.showStats();

        for (size_t i = 0; i < queryAmount; i++) {
            if (cchResults[i] == dijkstraResults[i])
                continue;

            std::cout << CLR_RED;
            std::cout << sources[i] << " -> " << targets[i] << ": " << cchResults[i]
                      << ", " << dijkstraResults[i] << std::endl;
            std::cout << CLR_RESET;
        }

        compare_num_data(cchResults, dijkstraResults);
    }

    // to access the different query objects
    inline EliminationTreeQuery<GRAPH>& getEliminationTreeQuery() noexcept
    {
        return elimTreeQuery;
    }

    inline DijkstraQuery<GRAPH>& getCCHDijkstraQuery() noexcept
    {
        return cchDijkstraQuery;
    }

    inline Dijkstra::Dijkstra<GRAPH>& getDijkstraQuery() noexcept
    {
        return dijkstraQuery;
    }

private:
    Data<GRAPH>& cchData;
    GRAPH& originalGraph;

    EliminationTreeQuery<GRAPH> elimTreeQuery;
    DijkstraQuery<GRAPH> cchDijkstraQuery;
    Dijkstra::Dijkstra<GRAPH> dijkstraQuery;

    std::vector<unsigned> sources;
    std::vector<unsigned> targets;

    std::vector<unsigned> cchResults;
    std::vector<unsigned> dijkstraResults;
};
} // namespace CCH
