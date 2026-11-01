#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "robin_hood.h"

#include "PipelineCommon.hpp"
#include "SimplexUtility.hpp"

struct SelectedPV;

class SimplexEnumerator
{
public:
    SimplexEnumerator(
        const DistanceMatrix& distance_matrix,
        const std::vector<std::vector<int64_t>>& binomial_table)
        : distance_matrix_(distance_matrix), binomial_table_(binomial_table) {}

    std::vector<std::pair<int64_t, double>> getSortedVREdges(const double maxeps);

    // Native first-window PwPH tails (no pseudo-vertices), with the same strict
    // upper bound and zero-edge exclusion as WindowState::prepareWindow.
    SimplexList getFirstWindowEdgeTail(double eps_lo, double eps_hi, int threads);
    SimplexList getFirstWindowTriangleTail(
        const SimplexList& new_edges, double eps_lo, double eps_hi, int threads);

    std::vector<std::pair<int64_t, double>> getSortedVRCofacets(const std::vector<std::pair<int64_t, double>>& sorted_simplex, const size_t dim, const double maxeps, const int threadnum);

    // Geometric enumeration for quotient complexes with pseudo-vertices.
    // Each PV contains sorted original-vertex representatives, and simplex
    // weight is the smallest clique realization weight among representative choices.
    std::vector<std::pair<int64_t, double>> getGeometricCofacetList(const std::vector<std::pair<int64_t, double>>& sorted_quotient_simplex_list,
                                        const std::vector<size_t>& active_labels,
                                        const std::vector<SelectedPV>& pseudo_vertices,
                                        const robin_hood::unordered_map<uint64_t, double>& pv_label_distance_hash,
                                        const size_t dim, const double maxeps, const int threadnum);

    // Verbose diagnostics use a separate enumerator so the regular path does not allocate
    // realization storage, pack witnesses, or impose canonical edge tie-breaking.
    std::vector<std::pair<int64_t, double>> getGeometricCofacetListWithRealizations(
                                        const std::vector<std::pair<int64_t, double>>& sorted_quotient_simplex_list,
                                        const std::vector<size_t>& active_labels,
                                        const std::vector<SelectedPV>& pseudo_vertices,
                                        const robin_hood::unordered_map<uint64_t, double>& pv_label_distance_hash,
                                        const size_t dim, const double maxeps, const int threadnum,
                                        robin_hood::unordered_map<int64_t, uint64_t>& pv_realization_out);

private:
    const DistanceMatrix& distance_matrix_;
    const std::vector<std::vector<int64_t>>& binomial_table_;

    // Local representative choices are encoded in uint64_t masks; PV cardinality must stay below this cap.
    static constexpr size_t MAX_PV_CARDINALITY_ = 64;
    static constexpr size_t UNCHOSEN_ = std::numeric_limits<size_t>::max();

    // Keep edge ordering consistent with simplex ordering: weight first, then
    // binomial index for deterministic handling of equal-weight edges.
    struct EdgeRecord
    {
        double weight;
        int64_t edge_bindex;
        uint8_t groupidx0, groupidx1;
        uint8_t localidx0, localidx1;

        bool operator<(const EdgeRecord& edge) const
        {
            if (weight != edge.weight)
                return weight < edge.weight;
            return edge_bindex < edge.edge_bindex;
        }
    };

    struct WitnessWorkspace
    {
        std::vector<size_t> facet_vertices;
        std::vector<size_t> cofacet_vertices;
        std::vector<std::vector<size_t>> singleton_slots;
        std::vector<const std::vector<size_t> *> rep_ptrs;
        std::vector<EdgeRecord> facet_edges;
        std::vector<EdgeRecord> covt_edges;
        std::vector<uint64_t> flattened_adjacency_mask;                  // (group0, group1, local0) -> local1 bitmask
        std::vector<std::vector<uint64_t>> recursion_candidate_local_index_mask_stack; // recursion depth -> group candidate masks
        std::vector<size_t> current_local_indices;                       // group -> selected local representative
    };

    struct RecordingWitnessWorkspace : WitnessWorkspace
    {
        uint64_t packed_realization = 0ULL;
        std::vector<std::pair<int64_t, uint64_t>> realizations;
    };

    template <bool RecordRealization>
    using SelectedWitnessWorkspace = std::conditional_t<RecordRealization, RecordingWitnessWorkspace, WitnessWorkspace>;

    template <bool RecordRealization>
    using WitnessWeightResult = std::conditional_t<RecordRealization, std::pair<double, uint64_t>, double>;

    void prepareFacetWitnessContext(WitnessWorkspace& ws,
                                    const std::vector<size_t>& facet_labels,
                                    const std::vector<SelectedPV>& pseudo_vertices,
                                    const size_t originalvtnum, const double maxeps) const;

    void prepareCovtWitnessGroup(WitnessWorkspace& ws, const size_t covt, const size_t facet_label_count,
                                 const std::vector<SelectedPV>& pseudo_vertices,
                                 const size_t originalvtnum, const double maxeps) const;

    template <bool RecordRealization>
    WitnessWeightResult<RecordRealization> getGeometricPVSimplexWeight(SelectedWitnessWorkspace<RecordRealization>& ws,
                                                                       const size_t target_simplex_label_count,
                                                                       const double lower_bound) const;

    bool findCliqueRecursive(const uint64_t* flattened_adjacency_mask, const size_t target_simplex_label_count, WitnessWorkspace& ws,
                             const size_t current_local_index_count) const;

    template <bool RecordRealization>
    std::vector<std::pair<int64_t, double>> enumerateGeometricCofacets(
                                        const std::vector<std::pair<int64_t, double>>& sorted_quotient_simplex_list,
                                        const std::vector<size_t>& active_labels,
                                        const std::vector<SelectedPV>& pseudo_vertices,
                                        const robin_hood::unordered_map<uint64_t, double>& pv_label_distance_hash,
                                        const size_t dim, const double maxeps, const int threadnum,
                                        robin_hood::unordered_map<int64_t, uint64_t>* pv_realization_out);
};
