/*
 * File: ./Helpers/StatisticsCollecter.h
 * Author: Patrick Steil
 * Author: Daniel-Delong Zhang 
 * Year: 2023/2024
 */
// Include guards to prevent multiple inclusion of the header file
#pragma once

// Standard C++ library headers
#include <iostream>
#include <string>
#include <vector>

// Custom headers for additional functionality
#include "Colors.h"
#include "Timer.h"

// Class for collecting statistics on various metrics
class StatisticsCollecter {
public:
    // Constructor taking a vector of metric names
    StatisticsCollecter(std::vector<std::string>& metricNames)
        : metricNames(metricNames)
        , counters(metricNames.size(), 0)
        , // Renamed to 'counters' for clarity
        numRounds(0) {};

    // Start a new round of data collection
    void newRound() noexcept { ++numRounds; }

    // Increase the counter associated with the specified metric
    void count(size_t i) noexcept
    {
        assert(i < metricNames.size());
        ++counters[i];
    }

    // Reset the counter for a specific metric
    void reset(size_t i) noexcept
    {
        assert(i < metricNames.size());
        counters[i] = 0;
    }

    // Clear all counters, keeping metric names
    void clear() noexcept { counters.assign(metricNames.size(), 0); }

    // Hard reset, assign new metrics
    void reset(std::vector<std::string>& newMetricNames)
    {
        metricNames = newMetricNames;
        clear();
    }

    // Print collected statistics along with metric names
    void printStats() const noexcept
    {
        if (numRounds == 0) {
            std::cout << "0 Rounds have been tracked, maybe something is wrong?"
                      << std::endl;
            return;
        }
        std::cout << "********** STATISTICS COLLECTED **********" << std::endl;
        std::cout << "Metric Name, Avg. counter, Total Counter" << std::endl;

        for (size_t i(0); i < metricNames.size(); ++i) {
            std::cout << metricNames[i] << ","
                      << (long double)(counters[i] / (long double)numRounds) << ","
                      << counters[i] << std::endl;
        }

        std::cout << "******************************************" << std::endl;
    }

    // Write collected statistics to a CSV file
    void writeStatsToCSV(const std::string& fileName) const
    {
        std::ofstream outputFile(fileName + ".csv");

        if (!outputFile.is_open()) {
            std::cerr << "Error: Unable to open the file " << fileName
                      << ".csv for writing." << std::endl;
            return;
        }

        outputFile << "Metric Name, Avg. counter, Total Counter" << std::endl;

        for (size_t i(0); i < metricNames.size(); ++i) {
            outputFile << metricNames[i] << ","
                       << (long double)(counters[i] / (long double)numRounds) << ","
                       << counters[i] << std::endl;
        }

        outputFile.close();
        std::cout << "Statistics written to " << fileName << ".csv successfully."
                  << std::endl;
    }

private:
    std::vector<std::string> metricNames;
    std::vector<size_t> counters; // Renamed to 'counters' for clarity
    size_t numRounds;
};

// Extension of StatisticsCollecter that also measures execution time of phases
class FullStatisticsCollecter {
public:
    // Constructor taking vectors of metric names and phase names
    FullStatisticsCollecter(std::vector<std::string>& newMetricNames,
        std::vector<std::string>& newPhasesNames)
        : metricNames(newMetricNames)
        , phasesNames(newPhasesNames)
        , counters(metricNames.size(), 0)
        , // Renamed to 'counters' for clarity
        timersPerPhase(phasesNames.size())
        , totalTimePerPhase(phasesNames.size(), 0)
        , totalCallsPerPhase(phasesNames.size(), 0)
        , numRounds(0) {};

    // Start a new round of data collection
    void newRound() noexcept { ++numRounds; }

    // Increase the counter associated with the specified metric
    void count(size_t i) noexcept
    {
        assert(i < metricNames.size());
        ++counters[i];
    }

    // Increase the counter associated with the specified metric by a specified
    // amount
    void count(size_t i, size_t amount) noexcept
    {
        assert(i < metricNames.size());
        counters[i] = counters[i] + amount;
    }

    // Start the timer for a specific phase
    void startPhase(size_t i) noexcept
    {
        assert(i < phasesNames.size());
        timersPerPhase[i].start();
        ++totalCallsPerPhase[i];
    }

    // Stop the timer for a specific phase
    void stopPhase(size_t i) noexcept
    {
        assert(i < phasesNames.size());
        timersPerPhase[i].stop();
        totalTimePerPhase[i] += timersPerPhase[i].elapsedMilliseconds();
    }

    // Reset the counter for a specific metric
    void resetMetric(size_t i) noexcept
    {
        assert(i < metricNames.size());
        counters[i] = 0;
    }

    // Reset the timer for a specific phase
    void resetPhase(size_t i) noexcept
    {
        assert(i < phasesNames.size());
        timersPerPhase[i].reset();
    }

    // Clear all metric counters, keeping metric names
    void clearMetric() noexcept { counters.assign(metricNames.size(), 0); }

    // Clear all phase timers, keeping phase names
    void clearPhases() noexcept
    {
        timersPerPhase.clear();
        timersPerPhase.reserve(phasesNames.size());
        totalTimePerPhase.assign(phasesNames.size(), 0);
    }

    // Hard reset, assign new metrics
    void resetMetric(std::vector<std::string>& newMetricNames)
    {
        metricNames = newMetricNames;
        clearMetric();
    }

    // Hard reset, assign new phases
    void resetPhase(std::vector<std::string>& newPhasesNames)
    {
        phasesNames = newPhasesNames;
        clearPhases();
    }

    // Print collected statistics along with metric and phase names
    void printStats() const noexcept
    {
        if (numRounds == 0) {
            std::cout << "0 Rounds have been tracked, maybe something is wrong?"
                      << std::endl;
            return;
        }
        std::cout << "********** STATISTICS COLLECTED **********" << std::endl;
        std::cout << "Metric Name, Avg. counter, Total Counter" << std::endl;

        for (size_t i(0); i < metricNames.size(); ++i) {
            std::cout << metricNames[i] << ","
                      << (long double)(counters[i] / (long double)numRounds) << ","
                      << counters[i] << std::endl;
        }

        std::cout << "******************************************" << std::endl;
        std::cout << "Phase Name, Avg. Time [ms]" << std::endl;

        for (size_t i(0); i < phasesNames.size(); ++i) {
            std::cout << phasesNames[i] << ","
                      << (long double)(totalTimePerPhase[i] / (long double)totalCallsPerPhase[i])
                      << std::endl;
        }

        std::cout << "******************************************" << std::endl;
    }

    // Write collected statistics to a CSV file
    void writeStatsToCSV(const std::string& fileName) const
    {
        std::ofstream outputFileMetric(fileName + ".metrics.csv");

        if (!outputFileMetric.is_open()) {
            std::cerr << "Error: Unable to open the file " << fileName
                      << ".metrics.csv for writing." << std::endl;
            return;
        }

        outputFileMetric << "Metric Name, Avg. counter, Total Counter" << std::endl;

        for (size_t i(0); i < metricNames.size(); ++i) {
            outputFileMetric << metricNames[i] << ","
                             << (long double)(counters[i] / (long double)numRounds)
                             << "," << counters[i] << std::endl;
        }

        outputFileMetric.close();

        std::ofstream outputFilePhases(fileName + ".phases.csv");

        if (!outputFilePhases.is_open()) {
            std::cerr << "Error: Unable to open the file " << fileName
                      << ".phases.csv for writing." << std::endl;
            return;
        }

        outputFilePhases << "Phase Name, Avg. Time [ms]" << std::endl;

        for (size_t i(0); i < phasesNames.size(); ++i) {
            outputFilePhases << phasesNames[i] << ","
                             << (long double)(totalTimePerPhase[i] / (long double)totalCallsPerPhase[i])
                             << std::endl;
        }

        outputFilePhases.close();
        std::cout << "Statistics written to " << fileName
                  << ".phases.csv successfully." << std::endl;
    }

private:
    std::vector<std::string> metricNames;
    std::vector<std::string> phasesNames;
    std::vector<size_t> counters; // Renamed to 'counters' for clarity
    std::vector<Timer> timersPerPhase;
    std::vector<long double> totalTimePerPhase;
    std::vector<long double> totalCallsPerPhase;
    size_t numRounds;
};
