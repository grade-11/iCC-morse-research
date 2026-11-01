#include "WindowState.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

#include <omp.h>

#include "SimplexEnumerator.hpp"
#include "SimplexUtility.hpp"

namespace
{
LabelPairKey makeLabelPairKey(size_t lhs, size_t rhs)
{
    if (lhs > rhs)
        std::swap(lhs, rhs);

    constexpr size_t MAX_PACKED_LABEL = std::numeric_limits<uint32_t>::max();
    if (lhs > MAX_PACKED_LABEL || rhs > MAX_PACKED_LABEL)
        throw std::overflow_error("PV label-distance keys support labels up to uint32_t max.");

    return (static_cast<LabelPairKey>(lhs) << 32) |
           static_cast<LabelPairKey>(rhs);
}

double computeLabelDistance(
    const PipelineRuntime& runtime,
    const std::vector<SelectedPV>& pseudo_vertices,
    const size_t lhs,
    const size_t rhs)
{
    const auto& distance_matrix = runtime.distanceMatrix();
    const size_t original_vertex_count = distance_matrix.vertexCount();

    if (lhs < original_vertex_count && rhs < original_vertex_count)
        return distance_matrix.distance(lhs, rhs);

    double minimum_distance = std::numeric_limits<double>::infinity();

    if (lhs < original_vertex_count)
    {
        for (const size_t rhs_representative :
             pseudo_vertices[rhs - original_vertex_count].representatives)
        {
            minimum_distance = std::min(
                minimum_distance,
                distance_matrix.distance(lhs, rhs_representative));
        }
        return minimum_distance;
    }

    if (rhs < original_vertex_count)
    {
        for (const size_t lhs_representative :
             pseudo_vertices[lhs - original_vertex_count].representatives)
        {
            minimum_distance = std::min(
                minimum_distance,
                distance_matrix.distance(lhs_representative, rhs));
        }
        return minimum_distance;
    }

    for (const size_t lhs_representative :
         pseudo_vertices[lhs - original_vertex_count].representatives)
    {
        for (const size_t rhs_representative :
             pseudo_vertices[rhs - original_vertex_count].representatives)
        {
            minimum_distance = std::min(
                minimum_distance,
                distance_matrix.distance(lhs_representative, rhs_representative));
        }
    }

    return minimum_distance;
}

struct QuotientEdgeData
{
    robin_hood::unordered_map<LabelPairKey, double> pv_label_distances;
    SimplexList sorted_edges;
};

QuotientEdgeData buildQuotientEdges(
    const PipelineRuntime& runtime,
    const std::vector<size_t>& active_labels,
    const std::vector<SelectedPV>& pseudo_vertices,
    const double eps_hi)
{
    const int worker_count = runtime.config().threads;
    const size_t original_vertex_count = runtime.distanceMatrix().vertexCount();

    std::vector<robin_hood::unordered_map<LabelPairKey, double>>
        thread_pv_label_distance_hashes(static_cast<size_t>(worker_count));
    std::vector<SimplexList> thread_edge_workspaces(static_cast<size_t>(worker_count));

    omp_set_num_threads(worker_count);

#pragma omp parallel for schedule(dynamic) num_threads(worker_count)
    for (size_t i = 0; i < active_labels.size(); ++i)
    {
        const int thread_id = omp_get_thread_num();
        auto& thread_distances = thread_pv_label_distance_hashes[thread_id];
        auto& thread_edges = thread_edge_workspaces[thread_id];

        for (size_t j = i + 1; j < active_labels.size(); ++j)
        {
            size_t lhs = active_labels[i];
            size_t rhs = active_labels[j];
            if (lhs > rhs)
                std::swap(lhs, rhs);

            const double weight = computeLabelDistance(
                runtime,
                pseudo_vertices,
                lhs,
                rhs);

            if (lhs >= original_vertex_count || rhs >= original_vertex_count)
                thread_distances.emplace(makeLabelPairKey(lhs, rhs), weight);

            // Preserve quotient-edge semantics: zero-distance pairs are not edges.
            if (weight > 0.0 && weight < eps_hi)
            {
                const SimplexBindex bindex = SimplexUtility::getEdgeBinomialIndex(
                    runtime.binomialTable(),
                    rhs,
                    lhs);
                thread_edges.emplace_back(bindex, weight);
            }
        }
    }

    size_t pv_label_pair_count = 0;
    for (const auto& thread_distances : thread_pv_label_distance_hashes)
        pv_label_pair_count += thread_distances.size();

    robin_hood::unordered_map<LabelPairKey, double> pv_label_distances;
    pv_label_distances.reserve(pv_label_pair_count);
    for (const auto& thread_distances : thread_pv_label_distance_hashes)
    {
        pv_label_distances.insert(
            thread_distances.begin(),
            thread_distances.end());
    }

    auto sorted_edges = SimplexUtility::sortAndMergeSimplexChunks(
        thread_edge_workspaces,
        worker_count);

    return QuotientEdgeData{
        std::move(pv_label_distances),
        std::move(sorted_edges)};
}
}

WindowState::WindowState(const size_t original_vertex_count)
    : original_vertex_count_(original_vertex_count),
      active_label_list_(original_vertex_count)
{
    std::iota(active_label_list_.begin(), active_label_list_.end(), size_t{0});
}

SimplexList WindowState::prepareWindow(
    const PipelineRuntime& runtime,
    const WindowBounds& bounds)
{
    std::vector<uint8_t> active_label_mask(totalLabelCount(), 0);
    for (const size_t label : active_label_list_)
        active_label_mask[label] = 1;

    robin_hood::unordered_map<LabelPairKey, double> pv_label_distances;
    SimplexList sorted_edges;

    switch (runtime.mode())
    {
    case PipelineMode::RegVRPH:
    {
        SimplexEnumerator simplex_enumerator(
            runtime.distanceMatrix(),
            runtime.binomialTable());
        // The ordinary enumerator intentionally retains zero-distance edges.
        sorted_edges = simplex_enumerator.getSortedVREdges(bounds.eps_hi);
        break;
    }

    case PipelineMode::PwPH:
    {
        auto quotient_edges = buildQuotientEdges(
            runtime,
            active_label_list_,
            pv_list_,
            bounds.eps_hi);
        pv_label_distances = std::move(quotient_edges.pv_label_distances);
        sorted_edges = std::move(quotient_edges.sorted_edges);
        break;
    }

    }

    active_label_mask_ = std::move(active_label_mask);
    pv_label_distance_hash_ = std::move(pv_label_distances);
    bounds_ = bounds;

    return sorted_edges;
}

void WindowState::extendFirstWindowTo(const double eps_hi)
{
    if (bounds_.eps_lo != 0.0 || !pv_list_.empty() ||
        active_label_list_.size() != original_vertex_count_ ||
        active_label_mask_.size() != original_vertex_count_)
    {
        throw std::logic_error("EXTEND requires prepared first-window geometry without PVs.");
    }
    if (!std::isfinite(eps_hi) || !(eps_hi > bounds_.eps_hi))
        throw std::invalid_argument("EXTEND requires a strictly larger finite cut.");
    bounds_.eps_hi = eps_hi;
}

void WindowState::commitSelectedPVs(
    std::vector<SelectedPV>&& selected_pvs,
    const std::unordered_set<size_t>& new_absorbed_labels)
{
    // The pipeline ends the prepared-window lifetime before support
    // postprocessing mutates persistent label state.
    assert(active_label_mask_.empty());
    assert(pv_label_distance_hash_.empty());

    std::vector<size_t> next_active_labels;
    next_active_labels.reserve(active_label_list_.size() + selected_pvs.size());

    for (const size_t label : active_label_list_)
    {
        if (!new_absorbed_labels.contains(label))
            next_active_labels.push_back(label);
    }

    size_t next_pv_label = totalLabelCount();
    pv_list_.reserve(pv_list_.size() + selected_pvs.size());
    for (auto& selected : selected_pvs)
    {
        pv_list_.push_back(std::move(selected));
        next_active_labels.push_back(next_pv_label++);
    }

    active_label_list_ = std::move(next_active_labels);
}

void WindowState::invalidateCurrentWindow() noexcept
{
    // These caches can be proportional to the complete active-label graph.
    // Empty swaps release their backing storage at the window boundary.
    decltype(active_label_mask_){}.swap(active_label_mask_);
    decltype(pv_label_distance_hash_){}.swap(pv_label_distance_hash_);
    bounds_ = WindowBounds{};
}
