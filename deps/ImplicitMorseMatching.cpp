#include "ImplicitMorseMatching.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include "DimensionFrame.hpp"
#include "SimplexUtility.hpp"

namespace
{
template <typename Comparator>
class ReusableIndexQueue
    : public std::priority_queue<
          std::size_t,
          std::vector<std::size_t>,
          Comparator>
{
public:
    using Base = std::priority_queue<
        std::size_t,
        std::vector<std::size_t>,
        Comparator>;
    using Base::Base;

    std::vector<std::size_t>& container() noexcept
    {
        return this->c;
    }
};

using MinIndexQueue = ReusableIndexQueue<std::greater<std::size_t>>;
using MaxIndexQueue = ReusableIndexQueue<std::less<std::size_t>>;

template <typename Queue>
void recoverQueueStorage(
    Queue& queue,
    std::vector<std::size_t>& workspace)
{
    auto& container = queue.container();
    container.clear();
    workspace = std::move(container);
}
}

std::int64_t ImplicitMorseMatching::findOptimizedApparentPairCofacet(
    DimensionFrame& frame,
    const std::int64_t facet_bindex,
    const double facet_weight)
{
    auto& workspace = *frame.interface_workspace_;
    const auto& distance_matrix = frame.runtime_.distanceMatrix();
    const auto& active_label_mask = frame.window_.activeLabelMask();
    const auto& pv_label_distance_hash = frame.window_.pvLabelDistanceHash();
    const std::size_t original_vertex_count = distance_matrix.vertexCount();

    std::int64_t tied_cofacet_rank = -1;
    const auto probe_exact_tie =
        [&](const std::int64_t candidate_bindex) -> bool
    {
        std::size_t cofacet_rank = 0;
        if (!frame.tryFindCofacetRank(candidate_bindex, cofacet_rank))
            return false;

        if (frame.cofacet_list_[cofacet_rank].second != facet_weight)
            return false;

        // Candidates arrive in increasing bindex order. Since no cofacet can
        // weigh less than its facet, the first stored exact tie is the minimum
        // cofacet under the complete (weight, bindex) list order.
        tied_cofacet_rank = static_cast<std::int64_t>(cofacet_rank);
        return true;
    };

    const auto check_apparent_pair_capability =
        [&](const std::size_t covertex,
            const std::int64_t candidate_bindex) -> bool
    {
        if (active_label_mask[covertex] == 0)
            return false;

        for (const std::size_t facet_label : workspace.vertex_workspace)
        {
            double label_distance = 0.0;
            if (covertex < original_vertex_count &&
                facet_label < original_vertex_count)
            {
                label_distance =
                    distance_matrix.distance(covertex, facet_label);
            }
            else
            {
                label_distance = SimplexUtility::getPVLabelDistance(
                    pv_label_distance_hash,
                    covertex,
                    facet_label);

                // Both labels are active, so WindowState must have memoized
                // every PV-incident pair during window preparation.
                assert(label_distance >= 0.0);
            }

            if (label_distance > facet_weight)
                return false;
        }

        return probe_exact_tie(candidate_bindex);
    };

    SimplexUtility::forEachImmediateCofacetInBindexOrder(
        frame.binomial_table_,
        workspace.vertex_workspace,
        facet_bindex,
        frame.window_.totalLabelCount(),
        frame.dimension_ - 1,
        check_apparent_pair_capability);

    return tied_cofacet_rank;
}

bool ImplicitMorseMatching::tryInstallApparentPair(
    DimensionFrame& frame,
    const std::size_t facet_list_index)
{
    auto& workspace = *frame.interface_workspace_;
    const auto& facet = frame.facet_list_[facet_list_index];

    const auto try_install_candidate =
        [&](const std::size_t cofacet_rank) -> bool
    {
        if (workspace.cofacet_mates[cofacet_rank] != workspace.UNMATCHED)
            return false;

        SimplexUtility::getFacetListIndicesInPlace(
            frame.binomial_table_,
            frame.active_facet_index_,
            workspace.facet_indices,
            workspace.vertex_workspace,
            frame.cofacet_list_[cofacet_rank].first,
            frame.window_.totalLabelCount(),
            frame.dimension_);

        if (workspace.facet_indices.empty() ||
            facet_list_index != *std::max_element(
                workspace.facet_indices.begin(),
                workspace.facet_indices.end()))
        {
            return false;
        }

        workspace.facet_mates[facet_list_index] =
            static_cast<DimensionFrame::RankIndex>(cofacet_rank);
        workspace.cofacet_mates[cofacet_rank] =
            static_cast<DimensionFrame::RankIndex>(facet_list_index);
        return true;
    };

    const std::int64_t optimized_candidate =
        findOptimizedApparentPairCofacet(
            frame,
            facet.first,
            facet.second);

    if (optimized_candidate >= 0 &&
        try_install_candidate(
            static_cast<std::size_t>(optimized_candidate)))
    {
        return true;
    }

    // The optimized search is only a filter. Rebuild the complete immediate
    // cofacet list before applying the ordinary apparent-pair test.
    collectImmediateCofacetRanks(frame, facet.first);
    if (workspace.cofacet_indices.empty())
        return false;

    const std::size_t minimum_cofacet_rank =
        *std::min_element(
            workspace.cofacet_indices.begin(),
            workspace.cofacet_indices.end());

    if (frame.cofacet_list_[minimum_cofacet_rank].second != facet.second)
        return false;

    return try_install_candidate(minimum_cofacet_rank);
}

void ImplicitMorseMatching::collectImmediateCofacetRanks(
    DimensionFrame& frame,
    const std::int64_t facet_bindex)
{
    auto& workspace = *frame.interface_workspace_;
    workspace.cofacet_indices.clear();

    const auto collect_existing_cofacet =
        [&](const std::size_t,
            const std::int64_t candidate_bindex) -> bool
    {
        std::size_t cofacet_rank = 0;
        if (frame.tryFindCofacetRank(candidate_bindex, cofacet_rank))
            workspace.cofacet_indices.push_back(cofacet_rank);
        return false;
    };

    SimplexUtility::forEachImmediateCofacetInBindexOrder(
        frame.binomial_table_,
        workspace.vertex_workspace,
        facet_bindex,
        frame.window_.totalLabelCount(),
        frame.dimension_ - 1,
        collect_existing_cofacet);
}

void ImplicitMorseMatching::matchPersistence(
    DimensionFrame& frame,
    const bool collect_dependency_support,
    const bool resume_matching) const
{
    auto& workspace = *frame.interface_workspace_;
    const std::size_t cofacet_count = frame.cofacet_list_.size();
    const std::size_t facet_count = frame.facet_list_.size();
    const bool collect_raw_support =
        collect_dependency_support && frame.isTopDimension();
    const bool trace_live_reductions = frame.retain_for_extension_;

    workspace.cofacet_indices.reserve(frame.dimension_ < 5 ? 32 : 64);
    workspace.facet_indices.reserve(frame.dimension_ + 1);
    workspace.vertex_workspace.reserve(frame.dimension_ + 1);
    workspace.priority_queue_workspace.reserve(cofacet_count / 2);

    workspace.raw_support_info.support_cofacet_indices.clear();
    workspace.raw_support_info.protected_facet_indices.clear();
    if (collect_raw_support)
    {
        workspace.raw_support_info.support_cofacet_indices.reserve(
            frame.dimension_ < 5 ? 64 : 128);
    }
    if (collect_dependency_support)
    {
        workspace.raw_support_info.protected_facet_indices.reserve(
            frame.dimension_ < 5 ? 32 : 64);
    }

    // Processing order is part of persistence and support-selection behavior.
    for (std::size_t reverse_index = facet_count;
         reverse_index > 0;
         --reverse_index)
    {
        const std::size_t facet_list_index = reverse_index - 1;
        const auto& facet = frame.facet_list_[facet_list_index];

        if (frame.active_facet_index_.find(facet.first) ==
            frame.active_facet_index_.end())
        {
            continue;
        }

        // An extension preserves every old finite pair. Only new facets and
        // old unmatched facets need another pass through the same matcher.
        if (resume_matching &&
            workspace.facet_mates[facet_list_index] != workspace.UNMATCHED)
        {
            continue;
        }

        // A live EXTEND reduction resumes on the appended cofacets only. Its
        // facet already failed the apparent test on unchanged older cofacets.
        const bool live_resume =
            trace_live_reductions && resume_matching &&
            frame.live_reductions_.contains(facet_list_index);
        if (live_resume)
        {
            collectLiveReductionTail(frame, facet_list_index);
            workspace.reduction_trace.clear();
        }
        else
        {
            if (tryInstallApparentPair(frame, facet_list_index))
                continue;

            // tryInstallApparentPair leaves the complete immediate-cofacet ranks
            // in the reusable workspace whenever no pair was installed.
            if (workspace.cofacet_indices.empty())
            {
                frame.persistent_pairs_.push_back(
                    {facet.second, -1.0, facet.first, -1});
                if (collect_dependency_support)
                {
                    workspace.raw_support_info.protected_facet_indices.push_back(
                        facet_list_index);
                }
                continue;
            }

            if (trace_live_reductions)
            {
                workspace.reduction_trace.assign(
                    1, static_cast<DimensionFrame::RankIndex>(facet_list_index));
            }
        }

        workspace.record_reduction_trace = trace_live_reductions;
        const std::int64_t terminal_cofacet =
            findCompressedTerminalCofacet(
                frame,
                workspace.cofacet_indices);
        workspace.record_reduction_trace = false;

        if (trace_live_reductions)
        {
            if (terminal_cofacet < 0)
                saveLiveReduction(frame, facet_list_index);
            else
                frame.live_reductions_.erase(facet_list_index);
        }

        if (terminal_cofacet < 0)
        {
            frame.persistent_pairs_.push_back(
                {facet.second, -1.0, facet.first, -1});
            if (collect_dependency_support)
            {
                workspace.raw_support_info.protected_facet_indices.push_back(
                    facet_list_index);
            }
            continue;
        }

        const std::size_t terminal_cofacet_rank =
            static_cast<std::size_t>(terminal_cofacet);

        if (collect_raw_support)
        {
            // Replay must see the pre-install matching: the current facet is
            // the unmatched terminal row of the reduced cofacet boundary.
            workspace.raw_support_info.support_cofacet_indices.push_back(
                collectReducedColumnSupport(
                    frame,
                    terminal_cofacet_rank,
                    facet_list_index));
        }

        const auto& cofacet = frame.cofacet_list_[terminal_cofacet_rank];

        workspace.facet_mates[facet_list_index] =
            static_cast<DimensionFrame::RankIndex>(terminal_cofacet_rank);
        workspace.cofacet_mates[terminal_cofacet_rank] =
            static_cast<DimensionFrame::RankIndex>(facet_list_index);

        if (facet.second != cofacet.second)
        {
            frame.persistent_pairs_.push_back(
                {facet.second,
                 cofacet.second,
                 facet.first,
                 cofacet.first});
        }
    }
}

void ImplicitMorseMatching::matchApparentPairsOnly(
    DimensionFrame& frame) const
{
    auto& workspace = *frame.interface_workspace_;
    const std::size_t cofacet_count = frame.cofacet_list_.size();
    const std::size_t facet_count = frame.facet_list_.size();

    workspace.cofacet_indices.reserve(frame.dimension_ < 5 ? 32 : 64);
    workspace.facet_indices.reserve(frame.dimension_ + 1);
    workspace.vertex_workspace.reserve(frame.dimension_ + 1);

    for (std::size_t reverse_index = facet_count;
         reverse_index > 0;
         --reverse_index)
    {
        const std::size_t facet_list_index = reverse_index - 1;
        const auto& facet = frame.facet_list_[facet_list_index];

        if (frame.active_facet_index_.find(facet.first) ==
            frame.active_facet_index_.end())
        {
            continue;
        }

        // This performs both the optimized probe and the complete ordinary
        // apparent fallback. It never enters compressed reduction.
        tryInstallApparentPair(frame, facet_list_index);
    }

    // Report active unmatched facets in ascending filtration-list order.
    for (std::size_t facet_list_index = 0;
         facet_list_index < facet_count;
         ++facet_list_index)
    {
        const auto& facet = frame.facet_list_[facet_list_index];
        if (frame.active_facet_index_.find(facet.first) ==
            frame.active_facet_index_.end())
        {
            continue;
        }

        if (workspace.facet_mates[facet_list_index] == workspace.UNMATCHED)
            frame.unmatched_facet_weights_.push_back(facet.second);
    }

    if (frame.isTopDimension())
    {
        // The top cofacet list is never promoted to a facet list, so expose
        // its unmatched cells explicitly in ascending filtration-list order.
        for (std::size_t cofacet_rank = 0;
             cofacet_rank < cofacet_count;
             ++cofacet_rank)
        {
            if (workspace.cofacet_mates[cofacet_rank] == workspace.UNMATCHED)
            {
                frame.unmatched_top_cofacet_weights_.push_back(
                    frame.cofacet_list_[cofacet_rank].second);
            }
        }
    }
}

void ImplicitMorseMatching::collectLiveReductionTail(
    DimensionFrame& frame,
    const std::size_t facet_list_index)
{
    auto& workspace = *frame.interface_workspace_;
    const auto& live = frame.live_reductions_.at(facet_list_index);
    const auto& distance_matrix = frame.runtime_.distanceMatrix();
    const double eps_hi = frame.window_.bounds().eps_hi;
    workspace.cofacet_indices.clear();

    for (const auto reducer : live.reducer_facets)
    {
        const auto& vertices = workspace.vertex_workspace;
        // In a first window without PVs a facet weighs its longest edge.
        const double facet_weight = frame.facet_list_[reducer].second;
        // Only appended cofacets, weight in [eps_lo, eps_hi), can be new.
        // The weight filter is exact for this frame without PVs; the rank
        // lookup still decides membership, as in collectImmediateCofacetRanks.
        const auto collect_appended_cofacet =
            [&](const std::size_t covertex,
                const std::int64_t candidate_bindex) -> bool
        {
            double weight = facet_weight;
            for (const std::size_t vertex : vertices)
                weight = std::max(weight, distance_matrix.distance(vertex, covertex));
            if (!(weight >= live.eps_lo && weight < eps_hi))
                return false;

            std::size_t cofacet_rank = 0;
            if (frame.tryFindCofacetRank(candidate_bindex, cofacet_rank))
            {
                if (cofacet_rank < live.cofacet_frontier)
                {
                    throw std::logic_error(
                        "collectLiveReductionTail: appended-weight cofacet has an old rank");
                }
                workspace.cofacet_indices.push_back(cofacet_rank);
            }
            return false;
        };

        SimplexUtility::forEachImmediateCofacetInBindexOrder(
            frame.binomial_table_,
            workspace.vertex_workspace,
            frame.facet_list_[reducer].first,
            frame.window_.totalLabelCount(),
            frame.dimension_ - 1,
            collect_appended_cofacet);
    }

    // The compressed queue cancels equal ranks in pairs, so mod 2 suffices.
    auto& ranks = workspace.cofacet_indices;
    std::sort(ranks.begin(), ranks.end());
    std::size_t kept = 0;
    for (std::size_t i = 0; i < ranks.size();)
    {
        std::size_t j = i;
        while (j < ranks.size() && ranks[j] == ranks[i])
            ++j;
        if (((j - i) & 1U) != 0U)
            ranks[kept++] = ranks[i];
        i = j;
    }
    ranks.resize(kept);
}

void ImplicitMorseMatching::saveLiveReduction(
    DimensionFrame& frame,
    const std::size_t facet_list_index)
{
    auto& workspace = *frame.interface_workspace_;
    auto& live = frame.live_reductions_[facet_list_index];
    auto& reducers = live.reducer_facets;
    reducers.insert(
        reducers.end(),
        workspace.reduction_trace.begin(),
        workspace.reduction_trace.end());
    std::sort(reducers.begin(), reducers.end());
    std::size_t kept = 0;
    for (std::size_t i = 0; i < reducers.size();)
    {
        std::size_t j = i;
        while (j < reducers.size() && reducers[j] == reducers[i])
            ++j;
        if (((j - i) & 1U) != 0U)
            reducers[kept++] = reducers[i];
        i = j;
    }
    reducers.resize(kept);
    live.cofacet_frontier = frame.cofacet_list_.size();
    live.eps_lo = frame.window_.bounds().eps_hi;
    workspace.reduction_trace.clear();
}

std::int64_t ImplicitMorseMatching::findCompressedTerminalCofacet(
    DimensionFrame& frame,
    const std::vector<std::size_t>& start_cofacet_indices)
{
    auto& workspace = *frame.interface_workspace_;
    workspace.priority_queue_workspace.clear();

    MinIndexQueue cofacet_queue(
        std::greater<std::size_t>{},
        std::move(workspace.priority_queue_workspace));

    for (const std::size_t cofacet_rank : start_cofacet_indices)
        cofacet_queue.push(cofacet_rank);

    while (!cofacet_queue.empty())
    {
        const std::size_t top_cofacet = cofacet_queue.top();
        cofacet_queue.pop();

        std::size_t multiplicity = 1;
        while (!cofacet_queue.empty() &&
               cofacet_queue.top() == top_cofacet)
        {
            cofacet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL)
            continue;

        if (workspace.cofacet_mates[top_cofacet] == workspace.UNMATCHED)
        {
            recoverQueueStorage(
                cofacet_queue,
                workspace.priority_queue_workspace);
            return static_cast<std::int64_t>(top_cofacet);
        }

        const std::size_t next_facet_list_index =
            static_cast<std::size_t>(workspace.cofacet_mates[top_cofacet]);

        enqueueReducedCompressedColumnTail(
            frame,
            next_facet_list_index,
            top_cofacet,
            cofacet_queue);
    }

    recoverQueueStorage(
        cofacet_queue,
        workspace.priority_queue_workspace);
    return -1;
}

template <typename MinQueue>
void ImplicitMorseMatching::enqueueReducedCompressedColumnTail(
    DimensionFrame& frame,
    const std::size_t facet_list_index,
    const std::size_t expected_pivot_cofacet,
    MinQueue& target_queue)
{
    auto& workspace = *frame.interface_workspace_;
    if (workspace.record_reduction_trace)
    {
        workspace.reduction_trace.push_back(
            static_cast<DimensionFrame::RankIndex>(facet_list_index));
    }

    collectImmediateCofacetRanks(
        frame,
        frame.facet_list_[facet_list_index].first);

    const auto minimum_cofacet = std::min_element(
        workspace.cofacet_indices.begin(),
        workspace.cofacet_indices.end());
    if (minimum_cofacet != workspace.cofacet_indices.end() &&
        *minimum_cofacet == expected_pivot_cofacet)
    {
        for (const std::size_t cofacet_rank : workspace.cofacet_indices)
        {
            if (cofacet_rank != expected_pivot_cofacet)
                target_queue.push(cofacet_rank);
        }
        return;
    }

    MinQueue column_queue(std::greater<std::size_t>{});
    for (const std::size_t cofacet_rank : workspace.cofacet_indices)
        column_queue.push(cofacet_rank);

    while (!column_queue.empty())
    {
        const std::size_t top_cofacet = column_queue.top();
        column_queue.pop();

        std::size_t multiplicity = 1;
        while (!column_queue.empty() &&
               column_queue.top() == top_cofacet)
        {
            column_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL)
            continue;

        if (top_cofacet == expected_pivot_cofacet)
        {
            // This is the term canceled by the matched column.
            break;
        }

        if (top_cofacet > expected_pivot_cofacet)
        {
            throw std::logic_error(
                "enqueueReducedCompressedColumnTail: passed expected pivot");
        }

        const auto matched_facet =
            workspace.cofacet_mates[top_cofacet];
        if (matched_facet == workspace.UNMATCHED)
        {
            throw std::logic_error(
                "enqueueReducedCompressedColumnTail: unmatched term before expected pivot");
        }

        const std::size_t matched_facet_list_index =
            static_cast<std::size_t>(matched_facet);
        if (matched_facet_list_index <= facet_list_index)
        {
            throw std::logic_error(
                "enqueueReducedCompressedColumnTail: reducer is not available in processing order");
        }

        enqueueReducedCompressedColumnTail(
            frame,
            matched_facet_list_index,
            top_cofacet,
            column_queue);
    }

    while (!column_queue.empty())
    {
        const std::size_t queued_cofacet = column_queue.top();
        column_queue.pop();

        std::size_t multiplicity = 1;
        while (!column_queue.empty() &&
               column_queue.top() == queued_cofacet)
        {
            column_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) != 0ULL &&
            queued_cofacet != expected_pivot_cofacet)
        {
            target_queue.push(queued_cofacet);
        }
    }
}

std::vector<std::size_t>
ImplicitMorseMatching::collectReducedColumnSupport(
    DimensionFrame& frame,
    const std::size_t terminal_cofacet,
    const std::size_t expected_facet_list_index)
{
    static_cast<void>(expected_facet_list_index);

    auto& workspace = *frame.interface_workspace_;

    workspace.facet_indices.clear();
    workspace.priority_queue_workspace.clear();

    MaxIndexQueue facet_queue(
        std::less<std::size_t>{},
        std::move(workspace.priority_queue_workspace));

    std::vector<std::size_t> cofacet_trace;
    cofacet_trace.reserve(frame.dimension_ < 5 ? 16 : 32);
    cofacet_trace.push_back(terminal_cofacet);

    SimplexUtility::getFacetListIndicesInPlace(
        frame.binomial_table_,
        frame.active_facet_index_,
        workspace.facet_indices,
        workspace.vertex_workspace,
        frame.cofacet_list_[terminal_cofacet].first,
        frame.window_.totalLabelCount(),
        frame.dimension_);

    for (const std::size_t facet_rank : workspace.facet_indices)
        facet_queue.push(facet_rank);

    while (!facet_queue.empty())
    {
        const std::size_t top_facet_list_index = facet_queue.top();
        facet_queue.pop();

        std::size_t multiplicity = 1;
        while (!facet_queue.empty() &&
               facet_queue.top() == top_facet_list_index)
        {
            facet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL)
            continue;

        const auto next_cofacet =
            workspace.facet_mates[top_facet_list_index];
        if (next_cofacet == workspace.UNMATCHED)
        {
            assert(expected_facet_list_index == top_facet_list_index);
            break;
        }

        const std::size_t cofacet_rank =
            static_cast<std::size_t>(next_cofacet);
        cofacet_trace.push_back(cofacet_rank);

        enqueueReducedCofacetBoundaryTail(
            frame,
            cofacet_rank,
            top_facet_list_index,
            facet_queue,
            cofacet_trace);
    }

    recoverQueueStorage(
        facet_queue,
        workspace.priority_queue_workspace);
    return cofacet_trace;
}

template <typename MaxQueue>
void ImplicitMorseMatching::enqueueReducedCofacetBoundaryTail(
    DimensionFrame& frame,
    const std::size_t cofacet_list_index,
    const std::size_t pivot_facet_list_index,
    MaxQueue& target_queue,
    std::vector<std::size_t>& cofacet_trace)
{
    auto& workspace = *frame.interface_workspace_;

    SimplexUtility::getFacetListIndicesInPlace(
        frame.binomial_table_,
        frame.active_facet_index_,
        workspace.facet_indices,
        workspace.vertex_workspace,
        frame.cofacet_list_[cofacet_list_index].first,
        frame.window_.totalLabelCount(),
        frame.dimension_);

    const auto maximum_facet = std::max_element(
        workspace.facet_indices.begin(),
        workspace.facet_indices.end());
    if (maximum_facet != workspace.facet_indices.end() &&
        *maximum_facet == pivot_facet_list_index)
    {
        for (const std::size_t facet_rank : workspace.facet_indices)
        {
            if (facet_rank != pivot_facet_list_index)
                target_queue.push(facet_rank);
        }
        return;
    }

    MaxQueue facet_queue(std::less<std::size_t>{});
    for (const std::size_t facet_rank : workspace.facet_indices)
        facet_queue.push(facet_rank);

    while (!facet_queue.empty())
    {
        const std::size_t top_facet_list_index = facet_queue.top();
        facet_queue.pop();

        std::size_t multiplicity = 1;
        while (!facet_queue.empty() &&
               facet_queue.top() == top_facet_list_index)
        {
            facet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) == 0ULL)
            continue;

        const auto next_cofacet =
            workspace.facet_mates[top_facet_list_index];
        if (next_cofacet != workspace.UNMATCHED &&
            top_facet_list_index > pivot_facet_list_index)
        {
            const std::size_t cofacet_rank =
                static_cast<std::size_t>(next_cofacet);
            cofacet_trace.push_back(cofacet_rank);

            enqueueReducedCofacetBoundaryTail(
                frame,
                cofacet_rank,
                top_facet_list_index,
                facet_queue,
                cofacet_trace);
            continue;
        }

        if (top_facet_list_index != pivot_facet_list_index)
        {
            assert(
                false &&
                "cofacet boundary replay did not stop at its pivot facet");
            target_queue.push(top_facet_list_index);
        }

        break;
    }

    while (!facet_queue.empty())
    {
        const std::size_t queued_facet = facet_queue.top();
        facet_queue.pop();

        std::size_t multiplicity = 1;
        while (!facet_queue.empty() &&
               facet_queue.top() == queued_facet)
        {
            facet_queue.pop();
            ++multiplicity;
        }

        if ((multiplicity & 1ULL) != 0ULL &&
            queued_facet != pivot_facet_list_index)
        {
            target_queue.push(queued_facet);
        }
    }
}
