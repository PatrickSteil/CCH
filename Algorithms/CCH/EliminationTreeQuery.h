/*
 * File: ./Algorithms/CCH/EliminationTreeQuery.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
#pragma once

#include <cassert>

#include "../../DataStructures/CCH/Data.h"
#include "../../DataStructures/Graph/Graphs.h"
#include "../../DataStructures/TypeDefs.h"
#include "../../Helpers/Colors.h"
#include "../../Helpers/StatisticsCollecter.h"

namespace CCH {

// The next two definitions are used in the StatisticsCollecter
enum ELIM_TREE_METRICS {
    NUM_OF_RELAXED_FWD_EDGES,
    NUM_OF_RELAXED_BWD_EDGES,
    NUM_EXPLORED_FWD_VERTICES,
    NUM_EXPLORED_BWD_VERTICES,
    NUM_UPDATED_FWD_LABELS,
    NUM_UPDATED_BWD_LABELS,
    NUM_EDGES_IN_PATH,
};

enum ELIM_TREE_PHASES {
    INIT_DS,
    CLEAR,
    QUERY,
    EXTRACT_PATH,
    RELAX_UPWARD,
    RELAX_DOWNWARD
};

static std::vector<std::string> ELIM_TREE_METRIC_NAMES = {
    "\"# of relaxed forward edges:\"",
    "\"# of relaxed backward edges:\"",
    "\"# of explored forward vertices:\"",
    "\"# of explored backward vertices:\"",
    "\"# of updated forward labels:\"",
    "\"# of updated backward labels:\"",
    "\"# of edges in unpacked path:\"",
};

static std::vector<std::string> ELIM_TREE_PHASE_NAMES = {
    "\"Init datastructures\"",
    "\"Clear datastructures\"",
    "\"Run Query\"",
    "\"Extract Path\"",
    "\"Relax Upwards\"",
    "\"Relax Downwards\"",
};

template <typename GRAPH>
class EliminationTreeQuery {
public:
    struct VertexLabel {
        VertexLabel()
            : distance(INFTY)
            , parent(-1)
            , timeStamp(-1)
        {
        }
        inline void reset(int time)
        {
            distance = INFTY;
            parent = -1;
            timeStamp = time;
        }
        size_t distance;
        VertexID parent;
        int timeStamp;
    };

    // Simple struct for a Path Edge with weight and direction
    // because edgeId is not easily accessible for the given vertices
    struct PathEdge {
        PathEdge()
            : from(noVertexID)
            , to(noVertexID)
            , weight(INFTY)
            , upward(false)
        {
        }
        PathEdge(VertexID from, VertexID to, size_t weight, bool upward)
            : from(from)
            , to(to)
            , weight(weight)
            , upward(upward) {};

        VertexID from;
        VertexID to;
        size_t weight;
        // indicator if this is an upward edge or downward
        bool upward;
    };

    EliminationTreeQuery(Data<GRAPH>& data)
        : data(data)
        , parent(data.getChordalGraph().numberOfNodes(), noVertexID)
        , upLabels(data.getChordalGraph().numberOfNodes())
        , downLabels(data.getChordalGraph().numberOfNodes())
        , upDownPathRoot(noVertexID)
        , unpackedPath(data.getChordalGraph().numberOfNodes())
        , rightIndex(0)
        , timeStamp(0)
        , statsCounter(ELIM_TREE_METRIC_NAMES, ELIM_TREE_PHASE_NAMES)
    {
        buildStaticGraphs();
    };

    // Build static up and down graphs
    inline void buildStaticGraphs()
    {
        // Reserve space for static graphs
        upGraph.reserve(data.getChordalGraph().numberOfNodes(),
            data.getChordalGraph().numberOfEdges());
        downGraph.reserve(data.getChordalGraph().numberOfNodes(),
            data.getChordalGraph().numberOfEdges());

        int upPrefix = 0;
        int downPrefix = 0;
        for (VertexID vertex = 0; vertex < data.getChordalGraph().numberOfNodes();
             ++vertex) {
            upGraph.getToAdj().push_back(upPrefix);
            downGraph.getToAdj().push_back(downPrefix);

            data.getChordalGraph().doForAllOutgoingEdgeIDs(
                vertex, [&](auto /* vertex */, auto to, auto out_weight, auto in_weight, auto edgeId) {
                    // get edges for upgraph
                    // check if edge is removed or not
                    if (!(data.getRemovalFlags()[edgeId] & UPWARD)) {
                        upPrefix++;
                        upGraph.getToVertex().push_back(to);
                        upGraph.getWeight().push_back(out_weight);
                    }

                    // get edges for downgraph
                    // check if edge is removed or not
                    if (!(data.getRemovalFlags()[edgeId] & DOWNWARD)) {
                        downPrefix++;
                        downGraph.getToVertex().push_back(to);
                        downGraph.getWeight().push_back(in_weight);
                    }

                    // build elimination tree
                    if (parent[vertex] == noVertexID || to < parent[vertex]) {
                        parent[vertex] = to;
                    }
                });
        }
        upGraph.getToAdj().push_back(upPrefix);
        downGraph.getToAdj().push_back(downPrefix);
    }

    // Query Methods

    inline void clear()
    {
        statsCounter.startPhase(CLEAR);
        ++timeStamp;
        if (timeStamp == 0) {
            // reset distance arrays
            VertexLabel emptyLabel = VertexLabel();
            std::fill(upLabels.begin(), upLabels.end(), emptyLabel);
            std::fill(downLabels.begin(), downLabels.end(), emptyLabel);
        }
        statsCounter.stopPhase(CLEAR);
    }

    inline void initSourceAndTarget(const VertexID source,
        const VertexID target)
    {
        statsCounter.startPhase(INIT_DS);
        statsCounter.count(NUM_UPDATED_FWD_LABELS);
        VertexLabel& sourcelabel = getUpLabel(source);
        sourcelabel.distance = 0;
        sourcelabel.parent = source;

        statsCounter.count(NUM_UPDATED_BWD_LABELS);
        VertexLabel& targetlabel = getDownLabel(target);
        targetlabel.distance = 0;
        targetlabel.parent = target;
        statsCounter.stopPhase(INIT_DS);
    }

    inline unsigned run(const VertexID s, const VertexID t)
    {
        statsCounter.startPhase(QUERY);
        statsCounter.newRound();

        clear();
        initSourceAndTarget(s, t);
        // run
        VertexID source = s;
        VertexID target = t;

        // relax until lowest common ancestor (LCA)
        while (source != target) {
            if (source < target) {
                relaxUpward(source);
                source = parent[source];
            } else {
                relaxDownward(target);
                target = parent[target];
            }
        }
        // find shortest distance and path from LCA to root and continue to relax
        size_t minDistance = INFTY;
        upDownPathRoot = noVertexID;
        // relax from LCA to root (source == target)
        do {
            // distf(z) + distb(z) min?
            if (getUpDistance(source) + getDownDistance(source) < minDistance) {
                minDistance = getUpDistance(source) + getDownDistance(source);
                upDownPathRoot = source;
            }
            relaxUpward(source);
            relaxDownward(source);
            source = parent[source];
        } while (source != noVertexID);

        // std::cout << "source: " << s + 1 << "\t"
        //           << "target: " << t + 1 << "\t";
        // std::cout << "upDownPathRoot: " << upDownPathRoot + 1 << "\t";
        // std::cout << "minDistance: " << minDistance << std::endl;

        statsCounter.stopPhase(QUERY);

        return minDistance;
    }

    inline unsigned runWithPathUnpacking(const VertexID s, const VertexID t)
    {
        size_t minDistance = run(s, t);

        // does all the path unpacking things
        unpackPath();

        // in case the distance should be returned after the query
        return minDistance;
    }

    inline void unpackPath()
    {
        statsCounter.startPhase(EXTRACT_PATH);
        // first get shortest updown Path
        std::vector<PathEdge> path;

        // get up path from root to source
        VertexID source = upDownPathRoot;
        path.reserve(data.getChordalGraph().numberOfNodes());
        while (getUpParent(source) != source) {
            path.emplace_back(
                getUpParent(source), source,
                getUpDistance(source) - getUpDistance(getUpParent(source)), true);
            source = getUpParent(source);
        }
        // reverse it => source -> root
        std::reverse(path.begin(), path.end());

        // get down path from root to target
        VertexID target = upDownPathRoot;
        while (getDownParent(target) != target) {
            path.emplace_back(
                target, getDownParent(target),
                getDownDistance(target) - getDownDistance(getDownParent(target)),
                false);
            target = getDownParent(target);
        }

        for (size_t i = 0; i < path.size(); i++) {
            unpackShortcut(path[i]);
        }
        statsCounter.count(NUM_EDGES_IN_PATH, rightIndex);
        // printPath();
        rightIndex = 0;
        statsCounter.stopPhase(EXTRACT_PATH);
    }

    inline void printPath()
    {
        std::cout << "Unpacked Path: " << std::endl;
        for (size_t i = 0; i < rightIndex; i++) {
            std::cout << "From: " << unpackedPath[i].from + 1
                      << " to: " << unpackedPath[i].to + 1
                      << " with weight: " << unpackedPath[i].weight
                      << " and direction: "
                      << (unpackedPath[i].upward ? "upward" : "downward")
                      << std::endl;
        }
    }

    // unpacks the shortcut from a vertex x to y by finding the matching lower
    // triangle.
    inline void unpackShortcut(PathEdge shortcut)
    {
        // check if shortcut is upward or downward and handle it accordingly
        if (shortcut.upward) {
            // if matching triangle found ==> shortcut and unpack it recursively
            // otherwise ==> original edge, insert it into final path.
            if (!data.doForAllLowerTrianglesAndBreak(
                    shortcut.from, shortcut.to,
                    [&](auto /* x */, auto /* y */, auto z, auto xzEdgeID,
                        auto yzEdgeID) {
                        if (shortcut.weight == data.getChordalGraph().getDownWeight(xzEdgeID) + data.getChordalGraph().getUpWeight(yzEdgeID)) {
                            // matching triangle found

                            // unpack shortcut x to z
                            unpackShortcut(PathEdge(
                                shortcut.from, z,
                                data.getChordalGraph().getDownWeight(xzEdgeID), false));

                            // unpack shortcut z to y
                            unpackShortcut(PathEdge(
                                z, shortcut.to,
                                data.getChordalGraph().getUpWeight(yzEdgeID), true));

                            // stop after first lower triangle was found
                            return true;
                        }
                        return false;
                    })) {
                // no matching triangle was found
                unpackedPath[rightIndex] = shortcut;
                rightIndex++;
            }
        } else {
            // downward shortcut

            // doForAllLowerTrianglesAndBreak requires x < y (rank)
            if (!data.doForAllLowerTrianglesAndBreak(
                    shortcut.to, shortcut.from,
                    [&](auto /* x */, auto /* y */, auto z, auto xzEdgeID,
                        auto yzEdgeID) {
                        if (shortcut.weight == data.getChordalGraph().getDownWeight(yzEdgeID) + data.getChordalGraph().getUpWeight(xzEdgeID)) {
                            // matching triangle found

                            // unpack shortcut y to z
                            unpackShortcut(PathEdge(
                                shortcut.from, z,
                                data.getChordalGraph().getDownWeight(yzEdgeID), false));

                            // unpack shortcut z to x
                            unpackShortcut(PathEdge(
                                z, shortcut.to,
                                data.getChordalGraph().getUpWeight(xzEdgeID), true));

                            // stop after first lower triangle was found
                            return true;
                        }
                        return false;
                    })) {
                // no matching triangle was found
                unpackedPath[rightIndex] = shortcut;
                rightIndex++;
            }
        }
    }

    inline void relaxUpward(const VertexID vertex)
    {
        statsCounter.startPhase(RELAX_UPWARD);
        statsCounter.count(NUM_EXPLORED_FWD_VERTICES);

        VertexLabel& currentLabel = getUpLabel(vertex);
        upGraph.doForAllOutgoingEdges(vertex, [&](auto /* vertex */, auto to, auto weight, auto /* edgeId */) {
            statsCounter.count(NUM_OF_RELAXED_FWD_EDGES);
            VertexLabel& toLabel = getUpLabel(to);
            if (currentLabel.distance + weight < toLabel.distance) {
                statsCounter.count(NUM_UPDATED_FWD_LABELS);
                toLabel.distance = currentLabel.distance + weight;
                toLabel.parent = vertex;
            }
        });
        statsCounter.stopPhase(RELAX_UPWARD);
    }

    inline void relaxDownward(const VertexID vertex)
    {
        statsCounter.startPhase(RELAX_DOWNWARD);
        statsCounter.count(NUM_EXPLORED_BWD_VERTICES);
        VertexLabel& currentLabel = getDownLabel(vertex);
        downGraph.doForAllOutgoingEdges(
            vertex,
            [&](auto /* vertex */, auto to, auto weight, auto /* edgeId */) {
                statsCounter.count(NUM_OF_RELAXED_BWD_EDGES);
                VertexLabel& toLabel = getDownLabel(to);
                if (currentLabel.distance + weight < toLabel.distance) {
                    statsCounter.count(NUM_UPDATED_BWD_LABELS);

                    toLabel.distance = currentLabel.distance + weight;
                    toLabel.parent = vertex;
                }
            });
        statsCounter.stopPhase(RELAX_DOWNWARD);
    }

    // Access

    inline VertexLabel& getUpLabel(const VertexID vertex) noexcept
    {
        assert((size_t)vertex < upLabels.size());
        VertexLabel& result = upLabels[vertex];
        if (result.timeStamp != timeStamp)
            result.reset(timeStamp);
        return result;
    }

    inline size_t getUpDistance(const VertexID vertex)
    {
        VertexLabel& label = getUpLabel(vertex);
        return label.distance;
    }

    inline VertexID getUpParent(const VertexID vertex)
    {
        VertexLabel& label = getUpLabel(vertex);
        return label.parent;
    }

    inline VertexLabel& getDownLabel(const VertexID vertex) noexcept
    {
        assert((size_t)vertex < downLabels.size());
        VertexLabel& result = downLabels[vertex];
        if (result.timeStamp != timeStamp)
            result.reset(timeStamp);
        return result;
    }

    inline size_t getDownDistance(const VertexID vertex)
    {
        VertexLabel& label = getDownLabel(vertex);
        return label.distance;
    }

    inline VertexID getDownParent(const VertexID vertex)
    {
        VertexLabel& label = getDownLabel(vertex);
        return label.parent;
    }

    // Print all the collected stats
    inline void showStats() const noexcept { statsCounter.printStats(); }

    inline FullStatisticsCollecter& getStatsCollecter() noexcept
    {
        return statsCounter;
    }

private:
    Data<GRAPH>& data;

    std::vector<VertexID> parent;

    Graph::StaticGraph upGraph;
    Graph::StaticGraph downGraph;

    std::vector<VertexLabel> upLabels;
    std::vector<VertexLabel> downLabels;

    VertexID upDownPathRoot;

    std::vector<PathEdge> unpackedPath;
    size_t rightIndex;

    int timeStamp;

    // StatisticsCollecter, which keeps track of all sort of metrics.
    FullStatisticsCollecter statsCounter;
};

} // namespace CCH