#pragma once

#include <cstddef>
#include <optional>
#include <vector>

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
    // is predicted to cost less than extending the current window to the endpoint; otherwise stop
    // committing for the rest of the run.
    bool ruleCommits(double epsilon, std::size_t cut);
    bool ruleDecide(const DependencySupportPostProcessor::Worth& w, const char* stage);
    struct RuleCalibration
    {
        double build_s = 0, build_w = 0, extend_s = 0, extend_w = 0, prep_s = 0;
        // The latest window's time over its prediction from the rate before it (1 before any prediction).
        double last_ratio = 1.0;
        // Every window so far: predicted work at its lower and upper scale (lower 0 for a build), time, build?
        struct Window { double w_lo, w_hi, seconds; bool build; };
        std::vector<Window> windows;
        double rate() const
        {
            const double work = build_w + extend_w;
            return work > 0 ? (build_s + extend_s) / work : 0.0;
        }
        void add(bool extend, double seconds, double w_lo, double w_hi)
        {
            windows.push_back({extend ? w_lo : 0.0, w_hi, seconds, !extend});
            const double work = w_hi - (extend ? w_lo : 0.0);
            if (work <= 0) return;
            if (rate() > 0) last_ratio = seconds / (rate() * work);
            (extend ? extend_s : build_s) += seconds;
            (extend ? extend_w : build_w) += work;
        }
        // How fast the time per unit of work grows with the window: the first window's build and its own
        // extensions give a rate (time over work) at the middle of each step; the slope of log rate against
        // log work, plus one, is beta, kept in [1, 3]; 1 with fewer than three steps.
        double beta() const;
        // Time to cover work w_from..w_to when the rate grows as W^(beta-1), scaled to the latest window
        // (an extension step has its mean rate at its middle; a build from zero, the rate at its top over beta).
        // Falls back to the pooled rate when no window has work.
        double cost(double b, double w_from, double w_to) const;
    };
    RuleCalibration rule_calibration_;
    std::size_t commits_done_ = 0;
    bool no_more_commits_ = false;

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
