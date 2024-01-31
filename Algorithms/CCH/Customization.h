/*
 * File: ./Algorithms/CCH/Customization.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <cassert>

#include "../../DataStructures/CCH/Data.h"
#include "../../DataStructures/Graph/Graphs.h"
#include "../../Helpers/ProgressBar.h"
#include "../../Helpers/StatisticsCollecter.h"

namespace CCH {

// ******************************
// The next two definitions are used in the StatisticsCollecter
enum CUSTO_METRICS {
    NUM_OF_LOWER_TRIANGLES,
    NUM_OF_MIDDLE_TRIANGLES,
    NUM_OF_UPPER_TRIANGLES,
    NUM_OF_UPDATED_EDGES,
    NUM_OF_REMOVED_EDGES
};

enum CUSTO_PHASES {
    BASIC_CUSTO,
    PERFECT_CUSTO,
    BASIC_CUSTO_LOWER_TRIANGLE_ENUM,
    PERFECT_CUSTO_TRIANGLES,
    APPLY_METRIC
};

static std::vector<std::string> CUSTO_METRIC_NAMES = {
    "\"# of lower triangles enumerated:\"",
    "\"# of middle triangles enumerated:\"",
    "\"# of upper triangles enumerated:\"",
    "\"# of updated edge weights:\"",
    "\"# of removed edges:\"",
};

static std::vector<std::string> CUSTO_PHASE_NAMES = {
    "\"Basic Customization\"", "\"Perfect Customization\"",
    "\"Lower Triangle\"", "\"Upper & Middle Triangle\"", "\"Apply Metric\""
};
// ******************************

template <typename GRAPH, bool VERBOSE = true, bool DEBUG = false>
class Customization {
public:
    Customization(Data<GRAPH>& data)
        : data(data)
        , profiler(CUSTO_METRIC_NAMES, CUSTO_PHASE_NAMES) {};

    // performs the (basic) customization
    // first, apply the edge weights to "original" edges
    // second, loop over all lower triangles from bottom to top
    inline void runCustomization() noexcept
    {
        profiler.newRound();

        if (VERBOSE)
            std::cout << CLR_BLUE << "[START]" << CLR_RESET << " basic customization"
                      << std::endl;

        // progressbar progress(data.getChordalGraph().numberOfNodes());

        profiler.startPhase(BASIC_CUSTO);

        data.doForAllVerticesBottomToTop([&](VertexID vertex) {
            // loop over the outgoing edges of the current vertex (edge leads to y)
            // and loop over all lower triangles of {x, y}
            data.getChordalGraph().doForAllOutgoingEdgeIDs(
                vertex, [&](auto x, auto y, auto& xyUpWeight, auto& xyDownWeight, auto /* xyEdgeID */) {
                    profiler.startPhase(BASIC_CUSTO_LOWER_TRIANGLE_ENUM);
                    data.doForAllLowerTriangles(
                        x, y,
                        [&](auto /* x */, auto /* y */, auto /* z */, auto xzEdgeID,
                            auto yzEdgeID) {
                            profiler.count(NUM_OF_LOWER_TRIANGLES);

                            // get the sum of both up weights of the edge {x, z} and {y,
                            // z}
                            size_t newWeight = data.getChordalGraph().getDownWeight(xzEdgeID) + data.getChordalGraph().getUpWeight(yzEdgeID);

                            // can we update the {x, y} up weight? if yes, count it
                            if (newWeight < xyUpWeight) {
                                xyUpWeight = newWeight;
                                profiler.count(NUM_OF_UPDATED_EDGES);
                            }

                            // same for the down weight
                            newWeight = data.getChordalGraph().getUpWeight(xzEdgeID) + data.getChordalGraph().getDownWeight(yzEdgeID);

                            if (newWeight < xyDownWeight) {
                                xyDownWeight = newWeight;
                                profiler.count(NUM_OF_UPDATED_EDGES);
                            }
                        });
                    profiler.stopPhase(BASIC_CUSTO_LOWER_TRIANGLE_ENUM);
                });
            // progress.update();
        });

        profiler.stopPhase(BASIC_CUSTO);

        if (VERBOSE)
            std::cout << CLR_RED << "[DONE]" << CLR_RESET << " basic customization"
                      << std::endl;
        if (DEBUG)
            doAllTrianglesFulfillTriangleInequality();
    }

    // performs the (basic) customization in parallel
    // first, apply the edge weights to "original" edges
    // second, loop over all lower triangles from bottom to top
    inline void runCustomizationParallel() noexcept
    {
        profiler.newRound();

        if (VERBOSE)
            std::cout << CLR_BLUE << "[START]" << CLR_RESET
                      << " basic parallel customization" << std::endl;

        profiler.startPhase(BASIC_CUSTO);

        for (size_t level(0); level < data.numberOfLevels(); ++level) {
// inside one level, you can perform the basic customization in parallel
#pragma omp parallel for schedule(dynamic)
            for (Index i = data.getToAdjOfLevel(level);
                 i < data.getToAdjOfLevel(level + 1); ++i) {
                VertexID vertex = data.getVertexOfLevel(i);
                // loop over the outgoing edges of the current vertex (edge leads to y)
                // and loop over all lower triangles of {x, y}
                data.getChordalGraph().doForAllOutgoingEdgeIDs(
                    vertex, [&](auto x, auto y, auto& xyUpWeight, auto& xyDownWeight, auto /* xyEdgeID */) {
                        // profiler.startPhase(BASIC_CUSTO_LOWER_TRIANGLE_ENUM);
                        data.doForAllLowerTriangles(
                            x, y,
                            [&](auto /* x */, auto /* y */, auto /* z */, auto xzEdgeID,
                                auto yzEdgeID) {
                                // profiler.count(NUM_OF_LOWER_TRIANGLES);

                                // get the sum of both up weights of the edge {x, z} and {y,
                                // z}
                                size_t newWeight = data.getChordalGraph().getDownWeight(xzEdgeID) + data.getChordalGraph().getUpWeight(yzEdgeID);

                                // can we update the {x, y} up weight? if yes, count it
                                if (newWeight < xyUpWeight) {
                                    xyUpWeight = newWeight;
                                    // profiler.count(NUM_OF_UPDATED_EDGES);
                                }

                                // same for the down weight
                                newWeight = data.getChordalGraph().getUpWeight(xzEdgeID) + data.getChordalGraph().getDownWeight(yzEdgeID);

                                if (newWeight < xyDownWeight) {
                                    xyDownWeight = newWeight;
                                    // profiler.count(NUM_OF_UPDATED_EDGES);
                                }
                            });
                        // profiler.stopPhase(BASIC_CUSTO_LOWER_TRIANGLE_ENUM);
                    });
            }
        }

        profiler.stopPhase(BASIC_CUSTO);

        if (VERBOSE)
            std::cout << CLR_RED << "[DONE]" << CLR_RESET << " basic customization"
                      << std::endl;
        if (DEBUG)
            doAllTrianglesFulfillTriangleInequality();
    }

    // performs the perfect customization
    // loop over all middle & upper triangles from top to bottom
    // the edges (to removed) are marked inside the removalFlags vector in the
    // data class

    // loop over the outgoing edges of the current vertex (edge leads to y)
    // and loop over all middle & upper triangles of {x, y} and set: w(x,y) :=
    // min{w(x, y), w(x,z) + w(z,y)};

    // we only need to loop over upper triangles, since we 'construct' a
    // middle triangle but just 'flipping' the roles of the middle and
    // bottom edge ... recall that an upper triangle of {x, y} is x < y < z
    // and a middle triangle of {x, y} is x < z < y

    // TODO maybe code clean up, looks horrible

    inline void runPerfectCustomization() noexcept
    {
        if (VERBOSE)
            std::cout << CLR_BLUE << "[START]" << CLR_RESET
                      << " perfect customization" << std::endl;
        // progressbar progress(data.getChordalGraph().numberOfNodes());

        profiler.startPhase(PERFECT_CUSTO);

        data.doForAllVerticesTopToBottom([&](VertexID vertex) {
            data.getChordalGraph().doForAllOutgoingEdgeIDs(
                vertex, [&](auto x, auto y, auto& xyUpWeight, auto& xyDownWeight, auto xyEdgeID) {
                    profiler.startPhase(PERFECT_CUSTO_TRIANGLES);
                    data.doForAllUpperTriangles(
                        x, y,
                        [&](auto /* x */, auto /* y */, auto /* z */, auto xzEdgeID,
                            auto yzEdgeID) {
                            profiler.count(NUM_OF_UPPER_TRIANGLES);

                            // (x -> z) [UP] + (z -> y) [DOWN] < (x -> y) [UP] ?

                            size_t newWeight = data.getChordalGraph().getUpWeight(xzEdgeID) + data.getChordalGraph().getDownWeight(yzEdgeID);

                            if (newWeight < xyUpWeight) {
                                xyUpWeight = newWeight;

                                data.setForRemoval(xyEdgeID, UPWARD);
                                profiler.count(NUM_OF_REMOVED_EDGES);
                            }

                            // (z -> x) [DOWN] + (y -> z) [UP] < (y -> x) [DOWN] ?

                            newWeight = data.getChordalGraph().getDownWeight(xzEdgeID) + data.getChordalGraph().getUpWeight(yzEdgeID);

                            if (newWeight < xyDownWeight) {
                                xyDownWeight = newWeight;

                                data.setForRemoval(xyEdgeID, DOWNWARD);
                                profiler.count(NUM_OF_REMOVED_EDGES);
                            }

                            // now work the other triangle
                            profiler.count(NUM_OF_MIDDLE_TRIANGLES);

                            // (x -> y) [UP] + (y -> z) [UP] < (x -> z) [UP] ?
                            newWeight = xyUpWeight + data.getChordalGraph().getUpWeight(yzEdgeID);

                            if (newWeight < data.getChordalGraph().getUpWeight(xzEdgeID)) {
                                data.getChordalGraph().getUpWeight(xzEdgeID) = newWeight;

                                data.setForRemoval(xzEdgeID, UPWARD);
                                profiler.count(NUM_OF_REMOVED_EDGES);
                            }

                            // (y -> x) [DOWN] + (z -> y) [DOWN] < (z -> x) [DOWN] ?
                            newWeight = xyDownWeight + data.getChordalGraph().getDownWeight(yzEdgeID);

                            if (newWeight < data.getChordalGraph().getDownWeight(xzEdgeID)) {
                                data.getChordalGraph().getDownWeight(xzEdgeID) = newWeight;

                                data.setForRemoval(xzEdgeID, DOWNWARD);
                                profiler.count(NUM_OF_REMOVED_EDGES);
                            }
                        });
                    profiler.stopPhase(PERFECT_CUSTO_TRIANGLES);
                });
            // progress.update();
        });

        profiler.stopPhase(PERFECT_CUSTO);
        if (VERBOSE)
            std::cout << CLR_RED << "[DONE]" << CLR_RESET << " perfect customization"
                      << std::endl;
    }

    // perfect paralell customization
    inline void runPerfectCustomizationParallel() noexcept
    {
        if (VERBOSE)
            std::cout << CLR_BLUE << "[START]" << CLR_RESET
                      << " perfect parallel customization" << std::endl;

        profiler.startPhase(PERFECT_CUSTO);

        for (int level(data.numberOfLevels() - 1); level >= 0; --level) {
// inside one level, you can perform the perfect customization in parallel
#pragma omp parallel for schedule(dynamic)
            for (Index i = data.getToAdjOfLevel(level);
                 i < data.getToAdjOfLevel(level + 1); ++i) {
                VertexID vertex = data.getVertexOfLevel(i);
                data.getChordalGraph().doForAllOutgoingEdgeIDs(
                    vertex, [&](auto x, auto y, auto& xyUpWeight, auto& xyDownWeight, auto xyEdgeID) {
                        data.doForAllUpperTriangles(
                            x, y,
                            [&](auto /* x */, auto /* y */, auto /* z */, auto xzEdgeID,
                                auto yzEdgeID) {
                                // (x -> z) [UP] + (z -> y) [DOWN] < (x -> y) [UP] ?

                                size_t newWeight = data.getChordalGraph().getUpWeight(xzEdgeID) + data.getChordalGraph().getDownWeight(yzEdgeID);

                                if (newWeight < xyUpWeight) {
                                    xyUpWeight = newWeight;

                                    data.setForRemoval(xyEdgeID, UPWARD);
                                }

                                // (z -> x) [DOWN] + (y -> z) [UP] < (y -> x) [DOWN] ?

                                newWeight = data.getChordalGraph().getDownWeight(xzEdgeID) + data.getChordalGraph().getUpWeight(yzEdgeID);

                                if (newWeight < xyDownWeight) {
                                    xyDownWeight = newWeight;

                                    data.setForRemoval(xyEdgeID, DOWNWARD);
                                }

                                // now work the other triangle

                                // (x -> y) [UP] + (y -> z) [UP] < (x -> z) [UP] ?
                                newWeight = xyUpWeight + data.getChordalGraph().getUpWeight(yzEdgeID);

                                if (newWeight < data.getChordalGraph().getUpWeight(xzEdgeID)) {
                                    data.getChordalGraph().getUpWeight(xzEdgeID) = newWeight;

                                    data.setForRemoval(xzEdgeID, UPWARD);
                                }

                                // (y -> x) [DOWN] + (z -> y) [DOWN] < (z -> x) [DOWN] ?
                                newWeight = xyDownWeight + data.getChordalGraph().getDownWeight(yzEdgeID);

                                if (newWeight < data.getChordalGraph().getDownWeight(xzEdgeID)) {
                                    data.getChordalGraph().getDownWeight(xzEdgeID) = newWeight;

                                    data.setForRemoval(xzEdgeID, DOWNWARD);
                                }
                            });
                    });
            }
        }

        profiler.stopPhase(PERFECT_CUSTO);
        if (VERBOSE)
            std::cout << CLR_RED << "[DONE]" << CLR_RESET
                      << " perfect parallel customization" << std::endl;
    }

    // this method applies all edge weights to the edges in the chordal graph
    // all shortcuts get INFTY
    inline void applyMetric(Graph::DynamicGraph& graph) noexcept
    {
        profiler.startPhase(APPLY_METRIC);

        data.doForAllVerticesTopToBottom([&](auto node) {
            graph.doForAllOutgoingEdgeIDs(
                node, [&](auto /* from */, auto /* to */, auto weight, auto i) {
                    auto& mapping = data.getEdgeMapping(i);

                    // TODO check if that ever happens
                    if (mapping.first == noIndex) [[unlikely]]
                        return;

                    // is it an up- or down edge?
                    if (mapping.second) {
                        data.getChordalGraph().getUpWeight(mapping.first) = weight;
                    } else {
                        data.getChordalGraph().getDownWeight(mapping.first) = weight;
                    }
                });
        });
        profiler.stopPhase(APPLY_METRIC);
    }

    // Print all the collected stats
    inline void showStats() const noexcept { profiler.printStats(); }

    // to assert
    inline void doAllTrianglesFulfillTriangleInequality() noexcept
    {
        data.doForAllVerticesBottomToTop([&](VertexID vertex) {
            // loop over the outgoing edges of the current vertex (edge leads to y)
            // and loop over all lower triangles of {x, y}
            data.getChordalGraph().doForAllOutgoingEdgeIDs(
                vertex, [&](auto x, auto y, auto& xyUpWeight, auto& xyDownWeight, auto /* xyEdgeID */) {
                    data.doForAllLowerTriangles(
                        x, y,
                        [&](auto /* x */, auto /* y */, auto /* z */, auto xzEdgeID,
                            auto yzEdgeID) {
                            // assert that the weight of the xy edge is smaller than the
                            // detour
                            assert(xyUpWeight <= data.getChordalGraph().getDownWeight(xzEdgeID) + data.getChordalGraph().getUpWeight(yzEdgeID));
                            assert(xyDownWeight <= data.getChordalGraph().getUpWeight(xzEdgeID) + data.getChordalGraph().getDownWeight(yzEdgeID));
                        });
                });
        });
    }

    inline FullStatisticsCollecter& getStatsCollecter() noexcept
    {
        return profiler;
    }

private:
    Data<GRAPH>& data;

    // StatisticsCollecter, which keeps track of all sort of metrics.
    FullStatisticsCollecter profiler;
};

} // namespace CCH
