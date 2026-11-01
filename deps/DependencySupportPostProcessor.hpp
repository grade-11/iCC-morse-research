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

private:
    static constexpr std::size_t MAX_PV_CARDINALITY = 64;

    static double maxPairwiseDistance(
        const DistanceMatrix& distance_matrix,
        const std::vector<std::size_t>& vertices);

    static double minCrossDistance(
        const DistanceMatrix& distance_matrix,
        const std::vector<std::size_t>& lhs,
        const std::vector<std::size_t>& rhs);
};
