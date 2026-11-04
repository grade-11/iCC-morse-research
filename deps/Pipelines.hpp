#pragma once

#include <cstddef>
#include <optional>

#include "DependencySupportPostProcessor.hpp"
#include "PipelineCommon.hpp"
#include "WindowState.hpp"

class PersistentHomologyPipeline
{
public:
    PersistentHomologyPipeline(
        DistanceMatrix distance_matrix,
        PipelineConfig config);

    void run();

private:
    PipelineRuntime runtime_;
};

class PwphPipeline
{
public:
    PwphPipeline(
        DistanceMatrix distance_matrix,
        PipelineConfig config);

    void run();

private:
    void runGuardExtendLadder();
    std::optional<std::size_t> runLiveGuardExtendPrefix();
    std::optional<std::size_t> runExtendSegment(std::size_t first_cut, double eps_lo);
    void runExtendAll();
    // --extend-rule: commit at an admitted cut only if a commit followed by one window to the endpoint
    // is predicted to cost less than extending the current window to the endpoint.
    bool ruleCommits(double epsilon, std::size_t cut);
    struct RuleCalibration
    {
        double build_s = 0, build_w = 0, extend_s = 0, extend_w = 0;
        void add(bool extend, double seconds, double work)
        {
            if (work <= 0) return;
            (extend ? extend_s : build_s) += seconds;
            (extend ? extend_w : build_w) += work;
        }
    };
    RuleCalibration rule_calibration_;

    PipelineRuntime runtime_;
    WindowState window_state_;
    DependencySupportPostProcessor support_processor_;
};

class ApparentPairPipeline
{
public:
    ApparentPairPipeline(
        DistanceMatrix distance_matrix,
        PipelineConfig config);

    void run();

private:
    PipelineRuntime runtime_;
};
