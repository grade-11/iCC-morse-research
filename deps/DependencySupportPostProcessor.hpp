#pragma once

#include <cstddef>
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

private:
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
