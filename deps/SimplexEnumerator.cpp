#include <algorithm>

#include "omp.h"
#include <chrono>
#include <iostream>
#include "SimplexEnumerator.hpp"
#include "SimplexUtility.hpp"
#include "WindowState.hpp"

std::vector<std::pair<int64_t, double>> SimplexEnumerator::getSortedVREdges(const double maxeps)
{
    //********************use openmp later************************//

    std::vector<std::pair<int64_t, double>> sorted_edge;

    size_t npt = distance_matrix_.vertexCount();

    for (size_t i = 0; i < npt - 1; i++)
    {
        for (size_t j = i + 1; j < npt; j++)
        {
            double weight = distance_matrix_.distance(i, j);
            if (weight < maxeps)
            {
                int64_t bindex = SimplexUtility::getEdgeBinomialIndex(binomial_table_, j, i);
                sorted_edge.emplace_back(bindex, weight);
            }
        }
    }

    SimplexUtility::sortSimplexByWeightThenIndex(sorted_edge);

    return sorted_edge;
}

std::vector<std::pair<int64_t, double>> SimplexEnumerator::getSortedVRCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex_list, const size_t dim, const double maxeps, const int threadnum)
{
    // dim == simplex dimension == cofacet dimension - 1
    const size_t npts = distance_matrix_.vertexCount();

    const int worker_count = threadnum > 0 ? threadnum : 1;
    std::vector<std::vector<std::pair<int64_t, double>>> thread_workspace(static_cast<size_t>(worker_count));

    omp_set_num_threads(worker_count);

#pragma omp parallel for schedule(dynamic) num_threads(worker_count)
    for (size_t i = 0; i < sorted_simplex_list.size(); ++i)
    {
        int threadid = omp_get_thread_num();
        auto& thread_cofacets = thread_workspace[threadid];

        const int64_t bindex = sorted_simplex_list[i].first;
        const double weight = sorted_simplex_list[i].second;

        std::vector<size_t> simplex_vertices = SimplexUtility::getSimplexVertices(binomial_table_, bindex, npts, dim);

        const size_t minvt = simplex_vertices.back(); // sorted in descending order
        for (size_t covt = 0; covt < minvt; ++covt)
        {
            double newweight = 0.0;
            for (const auto& vt : simplex_vertices)
            {
                const double distance = distance_matrix_.distance(covt, vt);
                if (distance > newweight)
                    newweight = distance;
            }

            double cofacetweight = newweight > weight ? newweight : weight;

            if (newweight < maxeps)
            {
                int64_t shiftedbindex = SimplexUtility::getBinomialIndex(binomial_table_, simplex_vertices, 1);
                int64_t cofacetbindex = shiftedbindex + covt;
                thread_cofacets.emplace_back(cofacetbindex, cofacetweight);
            }
        }
    }

    return SimplexUtility::sortAndMergeSimplexChunks(thread_workspace, threadnum);
}

// ======================= Geometric PV enumeration =======================

void SimplexEnumerator::prepareFacetWitnessContext(WitnessWorkspace& ws,
                                                                const std::vector<size_t>& facet_labels,
                                                                const std::vector<SelectedPV>& pseudo_vertices,
                                                                const size_t originalvtnum, const double maxeps) const
{
    const size_t K = facet_labels.size();

    if (ws.singleton_slots.size() < K + 1)
        ws.singleton_slots.resize(K + 1);
    if (ws.rep_ptrs.size() < K + 1)
        ws.rep_ptrs.resize(K + 1);

    for (size_t i = 0; i < K; ++i)
    {
        const size_t label = facet_labels[i];

        if (label < originalvtnum)
        {
            ws.singleton_slots[i].assign(1, label);
            ws.rep_ptrs[i] = &ws.singleton_slots[i];
        }
        else
        {
            ws.rep_ptrs[i] = &pseudo_vertices[label - originalvtnum].representatives;
        }
    }

    ws.facet_edges.clear();

    for (size_t i = 0; i + 1 < K; ++i)
    {
        const auto& rep_i = *ws.rep_ptrs[i];

        for (size_t j = i + 1; j < K; ++j)
        {
            const auto& rep_j = *ws.rep_ptrs[j];

            for (size_t a = 0; a < rep_i.size(); ++a)
            {
                for (size_t b = 0; b < rep_j.size(); ++b)
                {
                    const double w = distance_matrix_.distance(rep_i[a], rep_j[b]);
                    if (w < maxeps)
                    {
                        const size_t u = std::max(rep_i[a], rep_j[b]);
                        const size_t v = std::min(rep_i[a], rep_j[b]);
                        const int64_t edge_bindex = binomial_table_[u][2] + static_cast<int64_t>(v);
                        ws.facet_edges.push_back({w, edge_bindex, static_cast<uint8_t>(i), static_cast<uint8_t>(j),
                                                  static_cast<uint8_t>(a), static_cast<uint8_t>(b)});
                    }
                }
            }
        }
    }

    std::sort(ws.facet_edges.begin(), ws.facet_edges.end());
}

void SimplexEnumerator::prepareCovtWitnessGroup(WitnessWorkspace& ws, const size_t covt, const size_t facet_label_count,
                                                             const std::vector<SelectedPV>& pseudo_vertices,
                                                             const size_t originalvtnum, const double maxeps) const
{
    const size_t K = facet_label_count;

    if (covt < originalvtnum)
    {
        ws.singleton_slots[K].assign(1, covt);
        ws.rep_ptrs[K] = &ws.singleton_slots[K];
    }
    else
    {
        ws.rep_ptrs[K] = &pseudo_vertices[covt - originalvtnum].representatives;
    }

    const auto& rep_c = *ws.rep_ptrs[K];

    ws.covt_edges.clear();

    for (size_t i = 0; i < K; ++i)
    {
        const auto& rep_i = *ws.rep_ptrs[i];

        for (size_t a = 0; a < rep_i.size(); ++a)
        {
            for (size_t b = 0; b < rep_c.size(); ++b)
            {
                const double w = distance_matrix_.distance(rep_i[a], rep_c[b]);
                if (w < maxeps)
                {
                    const size_t u = std::max(rep_i[a], rep_c[b]);
                    const size_t v = std::min(rep_i[a], rep_c[b]);
                    const int64_t edge_bindex = binomial_table_[u][2] + static_cast<int64_t>(v);
                    ws.covt_edges.push_back({w, edge_bindex, static_cast<uint8_t>(i), static_cast<uint8_t>(K),
                                             static_cast<uint8_t>(a), static_cast<uint8_t>(b)});
                }
            }
        }
    }

    std::sort(ws.covt_edges.begin(), ws.covt_edges.end());
}

template <bool RecordRealization>
auto SimplexEnumerator::getGeometricPVSimplexWeight(SelectedWitnessWorkspace<RecordRealization>& ws,
                                                                 const size_t target_simplex_label_count,
                                                                 const double lower_bound) const
    -> WitnessWeightResult<RecordRealization>
{
    const size_t flattened_adjacency_mask_words = target_simplex_label_count * target_simplex_label_count * MAX_PV_CARDINALITY_;
    if (ws.flattened_adjacency_mask.size() < flattened_adjacency_mask_words)
        ws.flattened_adjacency_mask.resize(flattened_adjacency_mask_words);
    std::fill_n(ws.flattened_adjacency_mask.begin(), flattened_adjacency_mask_words, 0ULL);

    if (ws.current_local_indices.size() < target_simplex_label_count)
        ws.current_local_indices.resize(target_simplex_label_count);
    if (ws.recursion_candidate_local_index_mask_stack.size() < target_simplex_label_count)
        ws.recursion_candidate_local_index_mask_stack.resize(target_simplex_label_count);
    for (auto& level : ws.recursion_candidate_local_index_mask_stack)
    {
        if (level.size() < target_simplex_label_count)
            level.resize(target_simplex_label_count);
    }

    auto& initial_candidate_local_index_mask = ws.recursion_candidate_local_index_mask_stack[0];

    uint64_t* flattened_adjacency_mask = ws.flattened_adjacency_mask.data();

    const auto adjacencyMask = [flattened_adjacency_mask, target_simplex_label_count](const size_t g0, const size_t g1, const size_t local0) -> uint64_t &
    {
        return flattened_adjacency_mask[(g0 * target_simplex_label_count + g1) * MAX_PV_CARDINALITY_ + local0];
    };

    const auto fullMask = [&ws](const size_t g) -> uint64_t
    {
        const size_t repsize = ws.rep_ptrs[g]->size();
        return (repsize >= 64) ? ~0ULL : ((1ULL << repsize) - 1ULL);
    };

    uint64_t coverage = 0ULL;
    bool covered = false;

    const auto& fedges = ws.facet_edges;
    const auto& cedges = ws.covt_edges;
    size_t fi = 0;
    size_t ci = 0;

    while (fi < fedges.size() || ci < cedges.size())
    {
        bool take_facet_edge;
        if (fi >= fedges.size())
            take_facet_edge = false;
        else if (ci >= cedges.size())
            take_facet_edge = true;
        else
            take_facet_edge = !(cedges[ci] < fedges[fi]);

        const auto& edge = take_facet_edge ? fedges[fi++] : cedges[ci++];

        adjacencyMask(edge.groupidx0, edge.groupidx1, edge.localidx0) |= (1ULL << edge.localidx1);
        adjacencyMask(edge.groupidx1, edge.groupidx0, edge.localidx1) |= (1ULL << edge.localidx0);

        if (!covered)
        {
            coverage |= (1ULL << edge.groupidx0) | (1ULL << edge.groupidx1);
            if (static_cast<size_t>(__builtin_popcountll(coverage)) != target_simplex_label_count)
                continue;
            covered = true;
        }

        if (edge.weight < lower_bound)
            continue;

        const size_t g0 = edge.groupidx0;
        const size_t g1 = edge.groupidx1;
        const size_t l0 = edge.localidx0;
        const size_t l1 = edge.localidx1;

        // Anchored search: before inserting this edge, no witness clique had been found.
        // Therefore any clique that first appears now must include this newest edge.
        // Pin its two local representatives and pre-intersect every other label's
        // candidate mask with both endpoint adjacency rows before recursing.
        bool feasible = true;
        for (size_t g = 0; g < target_simplex_label_count; ++g)
        {
            if (g == g0)
            {
                initial_candidate_local_index_mask[g] = (1ULL << l0);
                ws.current_local_indices[g] = l0;
                continue;
            }
            if (g == g1)
            {
                initial_candidate_local_index_mask[g] = (1ULL << l1);
                ws.current_local_indices[g] = l1;
                continue;
            }

            ws.current_local_indices[g] = UNCHOSEN_;
            initial_candidate_local_index_mask[g] = fullMask(g) & adjacencyMask(g0, g, l0) & adjacencyMask(g1, g, l1);

            if (initial_candidate_local_index_mask[g] == 0ULL)
            {
                feasible = false;
                break;
            }
        }

        if (feasible && findCliqueRecursive(flattened_adjacency_mask, target_simplex_label_count, ws, 2))
        {
            if constexpr (RecordRealization)
            {
                uint64_t packed_realization = 0ULL;
                for (size_t g = 0; g < target_simplex_label_count; ++g)
                    packed_realization |= (static_cast<uint64_t>(ws.current_local_indices[g] & 0xFFULL) << (8 * g));

                return {edge.weight, packed_realization};
            }
            else
            {
                return edge.weight;
            }
        }
    }

    if constexpr (RecordRealization)
        return {std::numeric_limits<double>::infinity(), 0ULL};
    else
        return std::numeric_limits<double>::infinity();
}

bool SimplexEnumerator::findCliqueRecursive(const uint64_t* flattened_adjacency_mask, const size_t target_simplex_label_count, WitnessWorkspace& ws,
                                                         const size_t current_local_index_count) const
{
    if (current_local_index_count == target_simplex_label_count)
        return true;

    const size_t current_mask_index = current_local_index_count - 2;
    const size_t next_mask_index = current_local_index_count - 1;
    const auto& candidate_local_index_mask = ws.recursion_candidate_local_index_mask_stack[current_mask_index];

    size_t pivot = UNCHOSEN_;
    int min_candidates = std::numeric_limits<int>::max();

    for (size_t g = 0; g < target_simplex_label_count; ++g)
    {
        if (ws.current_local_indices[g] == UNCHOSEN_)
        {
            const int candnum = __builtin_popcountll(candidate_local_index_mask[g]);
            if (candnum < min_candidates)
            {
                min_candidates = candnum;
                pivot = g;
            }
        }
    }

    if (pivot == UNCHOSEN_)
        return false;

    auto& next_candidate_local_index_mask = ws.recursion_candidate_local_index_mask_stack[next_mask_index];

    uint64_t opts = candidate_local_index_mask[pivot];
    while (opts != 0ULL)
    {
        const size_t local = static_cast<size_t>(__builtin_ctzll(opts));
        opts &= (opts - 1);

        ws.current_local_indices[pivot] = local;

        bool feasible = true;
        for (size_t g = 0; g < target_simplex_label_count; ++g)
        {
            next_candidate_local_index_mask[g] = candidate_local_index_mask[g];

            if (g == pivot || ws.current_local_indices[g] != UNCHOSEN_)
                continue;

            next_candidate_local_index_mask[g] &= flattened_adjacency_mask[(pivot * target_simplex_label_count + g) * MAX_PV_CARDINALITY_ + local];

            if (next_candidate_local_index_mask[g] == 0ULL)
            {
                feasible = false;
                break;
            }
        }

        if (feasible && findCliqueRecursive(flattened_adjacency_mask, target_simplex_label_count, ws, current_local_index_count + 1))
            return true;
    }

    ws.current_local_indices[pivot] = UNCHOSEN_;
    return false;
}

template <bool RecordRealization>
std::vector<std::pair<int64_t, double>> SimplexEnumerator::enumerateGeometricCofacets(
    const std::vector<std::pair<int64_t, double>>& sorted_quotient_simplex_list,
    const std::vector<size_t>& active_labels,
    const std::vector<SelectedPV>& pseudo_vertices,
    const robin_hood::unordered_map<uint64_t, double>& pv_label_distance_hash,
    const size_t dim,
    const double maxeps,
    const int threadnum,
    robin_hood::unordered_map<int64_t, uint64_t>* pv_realization_out)
{
    const int worker_count = threadnum;
    std::vector<std::vector<std::pair<int64_t, double>>> thread_workspace(static_cast<size_t>(worker_count));

    const size_t originalvtnum = distance_matrix_.vertexCount();
    const size_t npts = originalvtnum + pseudo_vertices.size();

    // Avoid duplicating original-original distances in the PV-specific hash.
    const auto getQuotientLabelDistance =
        [this, originalvtnum, &pv_label_distance_hash](const size_t label_i, const size_t label_j)
    {
        if (label_i < originalvtnum && label_j < originalvtnum)
            return distance_matrix_.distance(label_i, label_j);
        return SimplexUtility::getPVLabelDistance(pv_label_distance_hash, label_i, label_j);
    };

    std::vector<SelectedWitnessWorkspace<RecordRealization>> witness_workspaces(static_cast<size_t>(worker_count));

    using Clock = std::chrono::steady_clock;

    struct WitnessTiming
    {
        size_t calls = 0;
        double search_ms = 0.0;
        double facet_prepare_ms = 0.0;
        double covt_prepare_ms = 0.0;
    };

    std::vector<WitnessTiming> witness_timings(
        static_cast<size_t>(worker_count));

    omp_set_num_threads(worker_count);

#pragma omp parallel for schedule(dynamic) num_threads(worker_count)
    for (size_t i = 0; i < sorted_quotient_simplex_list.size(); ++i)
    {
        const int threadid = omp_get_thread_num();
        auto& thread_cofacets = thread_workspace[threadid];
        auto& ws = witness_workspaces[threadid];

        const int64_t bindex = sorted_quotient_simplex_list[i].first;
        const double weight = sorted_quotient_simplex_list[i].second;

        SimplexUtility::getSimplexVerticesInPlace(binomial_table_, ws.facet_vertices, bindex, npts, dim);
        const auto& simplex_vertices = ws.facet_vertices;

        const size_t minfacetvt = simplex_vertices.back();
        const auto iter = std::lower_bound(active_labels.begin(), active_labels.end(), minfacetvt);
        const size_t vtpos = static_cast<size_t>(std::distance(active_labels.begin(), iter));
        if (vtpos == 0)
            continue;

        const bool has_pv_label = (simplex_vertices.front() >= originalvtnum);
        const size_t facet_label_count = dim + 1;
        bool witness_context_ready = false;

        for (size_t j = 0; j < vtpos; ++j)
        {
            const size_t covt = active_labels[j];

            double cofacetweight = 0.0;

            if (has_pv_label)
            {
                double lower_bound = weight;
                for (const auto& vt : simplex_vertices)
                {
                    const double d = getQuotientLabelDistance(covt, vt);
                    if (d > lower_bound)
                        lower_bound = d;
                }

                if (lower_bound >= maxeps)
                    continue;

                if (!witness_context_ready)
                {
                    const auto facet_prepare_start =
                        RecordRealization ? Clock::now() : Clock::time_point{};

                    prepareFacetWitnessContext(
                        ws, simplex_vertices, pseudo_vertices,
                        originalvtnum, maxeps);

                    if constexpr (RecordRealization)
                    {
                        const auto facet_prepare_end = Clock::now();
                        witness_timings[threadid].facet_prepare_ms +=
                            std::chrono::duration<double, std::milli>(
                                facet_prepare_end - facet_prepare_start).count();
                    }

                    witness_context_ready = true;
                }

                const auto covt_prepare_start =
                    RecordRealization ? Clock::now() : Clock::time_point{};

                prepareCovtWitnessGroup(
                    ws, covt, facet_label_count, pseudo_vertices,
                    originalvtnum, maxeps);

                if constexpr (RecordRealization)
                {
                    const auto covt_prepare_end = Clock::now();
                    witness_timings[threadid].covt_prepare_ms +=
                        std::chrono::duration<double, std::milli>(
                            covt_prepare_end - covt_prepare_start).count();
                }
                if constexpr (RecordRealization)
                {
                    const auto search_start = Clock::now();
                    const auto pv_witness =
                        getGeometricPVSimplexWeight<RecordRealization>(
                            ws, facet_label_count + 1, lower_bound);
                    const auto search_end = Clock::now();

                    auto& timing = witness_timings[threadid];
                    ++timing.calls;
                    timing.search_ms +=
                        std::chrono::duration<double, std::milli>(
                            search_end - search_start).count();

                    cofacetweight = pv_witness.first;
                    ws.packed_realization = pv_witness.second;
                }
                else
                {
                    cofacetweight =
                        getGeometricPVSimplexWeight<RecordRealization>(ws, facet_label_count + 1, lower_bound);
                }

                if (!(cofacetweight < maxeps))
                    continue;
            }
            else
            {
                cofacetweight = weight;
                for (const auto& vt : simplex_vertices)
                {
                    const double distance = distance_matrix_.distance(covt, vt);
                    if (distance > cofacetweight)
                        cofacetweight = distance;
                }
            }

            if (cofacetweight > 0.0 && cofacetweight < maxeps)
            {
                ws.cofacet_vertices.assign(simplex_vertices.begin(), simplex_vertices.end());
                ws.cofacet_vertices.push_back(covt);

                const int64_t cofacetbindex =
                    SimplexUtility::getBinomialIndex(binomial_table_, ws.cofacet_vertices, 0);
                thread_cofacets.emplace_back(cofacetbindex, cofacetweight);

                if constexpr (RecordRealization)
                {
                    if (has_pv_label)
                        ws.realizations.emplace_back(cofacetbindex, ws.packed_realization);
                }
            }
        }
    }

    if constexpr (RecordRealization)
    {
        size_t search_calls = 0;
        double search_worker_ms = 0.0;
        double facet_prepare_worker_ms = 0.0;
        double covt_prepare_worker_ms = 0.0;

        for (const auto& timing : witness_timings)
        {
            search_calls += timing.calls;
            search_worker_ms += timing.search_ms;
            facet_prepare_worker_ms += timing.facet_prepare_ms;
            covt_prepare_worker_ms += timing.covt_prepare_ms;
        }

        std::cout << "[witness-timing]"
                  << " eps_hi=" << maxeps
                  << " dim=" << (dim + 1)
                  << " search_calls=" << search_calls
                  << " search_worker_ms=" << search_worker_ms
                  << " facet_prepare_worker_ms=" << facet_prepare_worker_ms
                  << " covt_prepare_worker_ms=" << covt_prepare_worker_ms << '\n';
        size_t realization_count = 0;
        for (const auto& ws : witness_workspaces)
            realization_count += ws.realizations.size();

        pv_realization_out->clear();
        pv_realization_out->reserve(realization_count);
        for (const auto& ws : witness_workspaces)
        {
            for (const auto& [cofacet_bindex, packed] : ws.realizations)
                pv_realization_out->emplace(cofacet_bindex, packed);
        }
    }

    return SimplexUtility::sortAndMergeSimplexChunks(thread_workspace, threadnum);
}

std::vector<std::pair<int64_t, double>> SimplexEnumerator::getGeometricCofacetList(
    const std::vector<std::pair<int64_t, double>>& sorted_quotient_simplex_list,
    const std::vector<size_t>& active_labels,
    const std::vector<SelectedPV>& pseudo_vertices,
    const robin_hood::unordered_map<uint64_t, double>& pv_label_distance_hash,
    const size_t dim,
    const double maxeps,
    const int threadnum)
{
    return enumerateGeometricCofacets<false>(sorted_quotient_simplex_list, active_labels,
                                             pseudo_vertices, pv_label_distance_hash,
                                             dim, maxeps, threadnum, nullptr);
}

std::vector<std::pair<int64_t, double>> SimplexEnumerator::getGeometricCofacetListWithRealizations(
    const std::vector<std::pair<int64_t, double>>& sorted_quotient_simplex_list,
    const std::vector<size_t>& active_labels,
    const std::vector<SelectedPV>& pseudo_vertices,
    const robin_hood::unordered_map<uint64_t, double>& pv_label_distance_hash,
    const size_t dim,
    const double maxeps,
    const int threadnum,
    robin_hood::unordered_map<int64_t, uint64_t>& pv_realization_out)
{
    return enumerateGeometricCofacets<true>(sorted_quotient_simplex_list, active_labels,
                                            pseudo_vertices, pv_label_distance_hash,
                                            dim, maxeps, threadnum, &pv_realization_out);
}
