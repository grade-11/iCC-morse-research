#pragma once

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DependencySupport.hpp"
#include "PipelineCommon.hpp"
#include "SimplexEnumerator.hpp"
#include "WindowState.hpp"
#include "robin_hood.h"

class ImplicitMorseMatching;

class DimensionFrame
{
public:
    DimensionFrame(
        const PipelineRuntime& runtime,
        const WindowState& window,
        SimplexList prepared_edges,
        bool retain_for_extension = false);

    DimensionFrame(const DimensionFrame&) = delete;
    DimensionFrame& operator=(const DimensionFrame&) = delete;
    DimensionFrame(DimensionFrame&&) = delete;
    DimensionFrame& operator=(DimensionFrame&&) = delete;

    // EXTEND at max dimension 3: the retained triangle->tetrahedron frame over a
    // matched retained edge->triangle frame. Extend `lower` first, then this frame.
    [[nodiscard]]
    static std::unique_ptr<DimensionFrame> retainedUpperFrame(const DimensionFrame& lower);

    // Match once per dimension, consume the results, then call advance().
    // Continue matching only if advance() returns true.
    void matchPersistence(bool collect_dependency_support);
    // First PwPH window only. Retain the matching workspace across larger cuts.
    void matchPersistenceForExtension(bool collect_dependency_support = true);
    void extendTo(double eps_hi);
    void extendWith(const DimensionFrame& lower);
    [[nodiscard]] DependencySupportBatch snapshotDependencySupportBatch() const;
    void traceComplex() const;
    void matchApparentPairsOnly();

    [[nodiscard]]
    bool advance();

    [[nodiscard]]
    std::size_t dimension() const noexcept
    {
        return dimension_;
    }

    [[nodiscard]]
    std::size_t facetDimension() const noexcept
    {
        return dimension_ == 1 ? 1 : dimension_ - 1;
    }

    [[nodiscard]]
    bool isTopDimension() const noexcept
    {
        return dimension_ == runtime_.config().maxdim;
    }

    [[nodiscard]]
    std::size_t facetCount() const noexcept
    {
        return facet_list_.size();
    }

    [[nodiscard]]
    std::size_t cofacetCount() const noexcept
    {
        return cofacet_list_.size();
    }

    // Gate-only views of the filtration-ordered lists.
    [[nodiscard]]
    const SimplexList& facetList() const noexcept
    {
        return facet_list_;
    }

    [[nodiscard]]
    const SimplexList& cofacetList() const noexcept
    {
        return cofacet_list_;
    }

    [[nodiscard]]
    const std::vector<PersistentPairInfo>& persistentPairs() const noexcept
    {
        return persistent_pairs_;
    }

    [[nodiscard]]
    const std::vector<double>& unmatchedFacetWeights() const noexcept
    {
        return unmatched_facet_weights_;
    }

    [[nodiscard]]
    const std::vector<double>& unmatchedTopCofacetWeights() const noexcept
    {
        return unmatched_top_cofacet_weights_;
    }

    DependencySupportBatch takeDependencySupportBatch();

private:
    friend class ImplicitMorseMatching;

    using RankIndex = std::uint32_t;

    class FullCofacetIndex
    {
    public:
        explicit FullCofacetIndex(const SimplexList& cofacets);
        void appendRanks(std::size_t first_rank);
        // EXTEND tetrahedra: index an appended tail, every cofacet of weight
        // >= weight_lo, in a table of its own. One growing table would be
        // rehashed at almost every cut.
        void appendSegment(std::size_t first_rank, double weight_lo);

        [[nodiscard]]
        bool segmented() const noexcept
        {
            return !segment_lo_.empty();
        }

        [[nodiscard]]
        bool tryFindRank(
            SimplexBindex bindex,
            std::size_t& rank_out) const noexcept;

        // `weight` is the candidate's filtration weight; it selects the table.
        [[nodiscard]]
        bool tryFindRank(
            SimplexBindex bindex,
            double weight,
            std::size_t& rank_out) const noexcept;

    private:
        struct IndirectHash
        {
            using is_transparent = void;

            const SimplexList* cofacets;

            std::size_t operator()(const RankIndex rank) const noexcept
            {
                return robin_hood::hash<SimplexBindex>{}((*cofacets)[rank].first);
            }

            std::size_t operator()(const SimplexBindex bindex) const noexcept
            {
                return robin_hood::hash<SimplexBindex>{}(bindex);
            }
        };

        struct IndirectEqual
        {
            using is_transparent = void;

            const SimplexList* cofacets;

            bool operator()(const RankIndex lhs, const RankIndex rhs) const noexcept
            {
                return (*cofacets)[lhs].first == (*cofacets)[rhs].first;
            }

            bool operator()(const SimplexBindex bindex, const RankIndex rank) const noexcept
            {
                return bindex == (*cofacets)[rank].first;
            }

            bool operator()(const RankIndex rank, const SimplexBindex bindex) const noexcept
            {
                return (*cofacets)[rank].first == bindex;
            }
        };

        using RankSet = robin_hood::unordered_flat_set<RankIndex, IndirectHash, IndirectEqual>;

        // Borrow the vector object, not data(): append may reallocate its storage.
        const SimplexList* cofacets_;
        RankSet ranks_;
        // Segment s holds the ranks appended with weights in
        // [segment_lo_[s], segment_lo_[s + 1]); ranks_ holds the weights below.
        std::vector<double> segment_lo_;
        std::vector<RankSet> segments_;
    };

    struct RawMatchSupportInfo
    {
        std::vector<std::vector<std::size_t>> support_cofacet_indices;
        std::vector<std::size_t> protected_facet_indices;
    };

    struct InterfaceWorkspace
    {
        static constexpr RankIndex UNMATCHED =
            std::numeric_limits<RankIndex>::max();

        InterfaceWorkspace(
            const SimplexList& cofacets,
            std::size_t facet_count,
            std::size_t facet_capacity_hint = 0,
            std::size_t cofacet_capacity_hint = 0);

        FullCofacetIndex cofacet_lookup;
        // facet_mates[f] = c and cofacet_mates[c] = f; otherwise UNMATCHED.
        std::vector<RankIndex> facet_mates;
        std::vector<RankIndex> cofacet_mates;

        std::vector<std::size_t> cofacet_indices;
        std::vector<std::size_t> facet_indices;
        std::vector<std::size_t> vertex_workspace;
        std::vector<std::size_t> priority_queue_workspace;
        // EXTEND only: facets whose coboundaries a compressed reduction added.
        std::vector<RankIndex> reduction_trace;
        bool record_reduction_trace = false;

        RawMatchSupportInfo raw_support_info;
    };

    // EXTEND only: an active facet whose compressed reduction emptied in the
    // current frame. After appending cofacets of weight >= eps_lo, a fresh
    // reduction replays every older pop exactly; what reaches the new cofacets
    // is, mod 2, the coboundary of the reducer set on those cofacets alone.
    struct LiveReduction
    {
        std::vector<RankIndex> reducer_facets; // sorted; each present once (mod 2)
        std::size_t cofacet_frontier = 0;
        double eps_lo = 0.0;
    };

    struct FfiRealizationState
    {
        robin_hood::unordered_map<SimplexBindex, std::uint64_t> facet_realizations;
        robin_hood::unordered_map<SimplexBindex, std::uint64_t> cofacet_realizations;
    };

    const PipelineRuntime& runtime_;
    const WindowState& window_;
    const BinomialTable& binomial_table_;

    std::size_t dimension_ = 2;
    bool retain_for_extension_ = false;
    // EXTEND of a window with pseudo vertices: tails come from the geometric enumerator, and the
    // tetrahedron lookup is not segmented by weight.
    bool extension_geometric_ = false;
    bool extension_matched_ = false;
    double extension_hi_ = 0.0;
    double last_enum_ms_ = 0.0;  // enumeration part of the latest extension tail (timing log only)
    // Capacity reserved for later tails, as a multiple of the first window.
    std::size_t extension_reserve_factor_ = 1;
    // First facet appended by the latest extendTo(); the edge tail for tetrahedra.
    std::size_t extension_facet_tail_ = 0;

    SimplexEnumerator simplex_enumerator_;

    SimplexList facet_list_;
    SimplexList cofacet_list_;

    robin_hood::unordered_map<SimplexBindex, std::size_t> active_facet_index_;
    std::vector<std::size_t> extension_kruskal_parent_;

    std::vector<std::uint8_t> next_active_facet_mask_;
    std::size_t next_active_facet_count_ = 0;

    std::vector<PersistentPairInfo> persistent_pairs_;
    std::vector<std::optional<PersistentPairInfo>> extension_pairs_by_facet_;
    std::vector<std::optional<DependencySupport>> extension_supports_by_facet_;
    std::unordered_map<std::size_t, LiveReduction> live_reductions_;
    std::vector<double> unmatched_facet_weights_;
    std::vector<double> unmatched_top_cofacet_weights_;

    DependencySupportBatch dependency_support_batch_;
    std::unordered_set<std::size_t> protected_labels_;

    std::optional<FfiRealizationState> ffi_realizations_;
    std::optional<InterfaceWorkspace> interface_workspace_;

    struct UpperFrameTag
    {
    };
    DimensionFrame(UpperFrameTag, const DimensionFrame& lower);

    void resumeExtension(
        std::size_t old_cofacet_count, double eps_hi,
        std::chrono::steady_clock::time_point tail_start);
    void buildInitialActiveFacetIndex();
    void enumerateInitialCofacets();
    SimplexList enumerateNextCofacets(
        robin_hood::unordered_map<SimplexBindex, std::uint64_t>* realization_output);

    void materializeMatchSupport(RawMatchSupportInfo raw_support_info);
    void mergeExtensionResults();
    void finalizeNextActiveFacets();
    void releaseActiveFacetIndex();
    void initializeInterfaceWorkspace();
    RawMatchSupportInfo releaseInterfaceWorkspace(
        bool preserve_raw_support_info);
    void clearCurrentInterfaceResults();

    [[nodiscard]]
    bool tryFindCofacetRank(
        SimplexBindex bindex,
        std::size_t& rank_out) const noexcept;

    [[nodiscard]]
    bool tryFindCofacetRank(
        SimplexBindex bindex,
        double weight,
        std::size_t& rank_out) const noexcept;

    // True when cofacet lookups must carry the candidate's weight.
    [[nodiscard]]
    bool segmentedCofacetIndex() const noexcept
    {
        return interface_workspace_.has_value() && interface_workspace_->cofacet_lookup.segmented();
    }

    void reportFalseFacetIdentificationStats() const;
};
