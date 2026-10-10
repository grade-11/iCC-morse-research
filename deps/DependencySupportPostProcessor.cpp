#include "DependencySupportPostProcessor.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <utility>

void DependencySupportPostProcessor::processPwph(
    const PipelineRuntime& runtime,
    WindowState& window,
    const double eps_hi,
    DependencySupportBatch&& batch) const
{
    const std::size_t raw_support_count = batch.supports.size();
    const std::size_t original_vertex_count = window.originalVertexCount();

    const double diameter_limit =
        runtime.config().pv_cap_scale * eps_hi;
    const double min_separation = runtime.config().absoluteMinSeparation();
    const auto& distance_matrix = runtime.distanceMatrix();

    std::vector<SelectedPV> selected_pvs;
    std::unordered_set<std::size_t> new_absorbed_labels;
    std::size_t eligible_support_count = 0;
    std::size_t rejected_overlap_count = 0;

    const auto is_separated = [&](const std::vector<std::size_t>& vertices)
    {
        if (min_separation <= 0.0)
            return true;

        // PwPH never deactivates an existing PV. Iterating active PV labels is
        // equivalent today and is the correct lifetime boundary for a future
        // absorbing policy, where historical inactive records must be ignored.
        for (const std::size_t label : window.activeLabels())
        {
            if (label < original_vertex_count)
                continue;

            const auto& existing =
                window.pseudoVertices()[label - original_vertex_count].representatives;
            if (minCrossDistance(distance_matrix, vertices, existing) <
                min_separation)
            {
                return false;
            }
        }

        for (const auto& selected : selected_pvs)
        {
            if (minCrossDistance(
                    distance_matrix,
                    vertices,
                    selected.representatives) < min_separation)
            {
                return false;
            }
        }

        return true;
    };

    // Matching discovers supports while traversing facets from high to low
    // rank. PwPH applies its static filters while traversing that discovery
    // vector in reverse, preserving the existing overlap-selection order
    // without allocating a second vector of support containers.
    for (auto iter = batch.supports.rbegin();
         iter != batch.supports.rend();
         ++iter)
    {
        auto& labels = iter->labels;

        const bool contains_pv = std::any_of(
            labels.begin(),
            labels.end(),
            [original_vertex_count](const std::size_t label)
            {
                return label >= original_vertex_count;
            });

        // Preserve the frozen-PV policy and its ordering. This test precedes
        // protection, so the additional protected PV labels retained by the
        // generic frame cannot affect current PwPH selection.
        if (contains_pv)
            continue;

        const bool contains_protected_label = std::any_of(
            labels.begin(),
            labels.end(),
            [&](const std::size_t label)
            {
                return batch.protected_labels.contains(label);
            });
        if (contains_protected_label)
            continue;

        ++eligible_support_count;

        bool overlaps_new_absorbed_label = false;
        for (const std::size_t label : labels)
        {
            if (new_absorbed_labels.contains(label))
            {
                overlaps_new_absorbed_label = true;
                break;
            }
        }
        if (overlaps_new_absorbed_label)
        {
            ++rejected_overlap_count;
            continue;
        }

        if (labels.empty() || labels.size() >= MAX_PV_CARDINALITY)
            continue;

        const double diameter = maxPairwiseDistance(
            distance_matrix, labels);
        if (diameter > diameter_limit)
            continue;
        if (!is_separated(labels))
            continue;

        new_absorbed_labels.insert(labels.begin(), labels.end());
        selected_pvs.push_back(
            SelectedPV{std::move(labels), diameter});
    }

    if (runtime.config().verbose && runtime.config().maxdim >= 2)
    {
        std::cout << "raw pv support cofacets size = " << raw_support_count << '\n';
        std::cout << "pv support label sets with no pv/protected contents size = "
                  << eligible_support_count << '\n';
        std::cout << "supports rejected by overlap with newly selected PVs = "
                  << rejected_overlap_count << '\n';
    }

    window.commitSelectedPVs(
        std::move(selected_pvs),
        new_absorbed_labels);
}

double DependencySupportPostProcessor::maxPairwiseDistance(
    const DistanceMatrix& distance_matrix,
    const std::vector<std::size_t>& vertices)
{
    double maximum = 0.0;
    for (auto lhs = vertices.begin(); lhs != vertices.end(); ++lhs)
    {
        auto rhs = lhs;
        ++rhs;
        for (; rhs != vertices.end(); ++rhs)
            maximum = std::max(maximum, distance_matrix.distance(*lhs, *rhs));
    }
    return maximum;
}

double DependencySupportPostProcessor::minCrossDistance(
    const DistanceMatrix& distance_matrix,
    const std::vector<std::size_t>& lhs,
    const std::vector<std::size_t>& rhs)
{
    double minimum = std::numeric_limits<double>::infinity();
    for (const std::size_t lhs_vertex : lhs)
    {
        for (const std::size_t rhs_vertex : rhs)
        {
            minimum = std::min(
                minimum,
                distance_matrix.distance(lhs_vertex, rhs_vertex));
        }
    }
    return minimum;
}
