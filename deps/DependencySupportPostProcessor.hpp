#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include "DependencySupport.hpp"
#include "PipelineCommon.hpp"
#include "WindowState.hpp"

class DependencySupportPostProcessor
{
public:
    // Current frozen-PV policy. A future absorbing policy belongs beside this
    // method and would need to flatten PV-containing supports before commit.
    std::size_t processPwph(
        const PipelineRuntime& runtime,
        WindowState& window,
        double eps_hi,
        DependencySupportBatch&& batch) const;

    // EXTEND gain rule: predicted work removed by the latest commit, and its cut (-1: none yet).
    [[nodiscard]] double lastGain() const noexcept { return last_gain_; }
    [[nodiscard]] double lastGainEpsilon() const noexcept { return last_gain_epsilon_; }
    // Whether the latest commit predicts that one more cut pays (see extend_cost_rule).
    [[nodiscard]] bool lastKeepCutting() const noexcept { return last_keep_cutting_; }

    // Predicted top-dimension work at the latest processed cut (--extend-worth-log, --extend-rule):
    // W0 without and Wsel with the committed supports, at the cut, the next cut and the endpoint.
    struct Worth
    {
        double epsilon = -1, next = 0, end = 0, w0_r = 0, w0_next = 0, wsel_next = 0, w0_end = 0, wsel_end = 0;
        double F_s = 0, C_s = 0;
        std::size_t committed = 0;
    };
    [[nodiscard]] const Worth& lastWorth() const noexcept { return last_worth_; }
    // --extend-rule: a decision taken on the eligible supports, before the risk filters run. The filters
    // only drop or shrink supports, so these supports overstate what a commit removes; if even they do not
    // pay, the processor commits nothing and skips the filters. Returns true to go on (commit if it pays).
    void setPreFilterRule(std::function<bool(const Worth&)> rule) { pre_filter_rule_ = std::move(rule); }
    // Whether the latest processPwph stopped at the pre-filter decision.
    [[nodiscard]] bool lastStoppedBeforeFilters() const noexcept { return last_stopped_before_filters_; }

private:
    mutable Worth last_worth_;
    std::function<bool(const Worth&)> pre_filter_rule_;
    mutable bool last_stopped_before_filters_ = false;
    mutable double last_gain_ = 0.0;
    mutable double last_gain_epsilon_ = -1.0;
    mutable bool last_keep_cutting_ = true;

    static constexpr std::size_t MAX_PV_CARDINALITY = 64;

    static double maxPairwiseDistance(
        const DistanceMatrix& distance_matrix,
        const std::vector<std::size_t>& vertices);

    static double minCrossDistance(
        const DistanceMatrix& distance_matrix,
        const std::vector<std::size_t>& lhs,
        const std::vector<std::size_t>& rhs);
};
