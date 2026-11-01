#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

class DimensionFrame;

// Stateless algorithm driver. DimensionFrame owns every list, lookup, mate
// array, and reusable interface workspace touched by these operations.
class ImplicitMorseMatching
{
private:
    friend class DimensionFrame;

    void matchPersistence(
        DimensionFrame& frame,
        bool collect_dependency_support,
        bool resume_matching = false) const;

    void matchApparentPairsOnly(DimensionFrame& frame) const;

private:
    static std::int64_t findOptimizedApparentPairCofacet(
        DimensionFrame& frame,
        std::int64_t facet_bindex,
        double facet_weight);

    static bool tryInstallApparentPair(
        DimensionFrame& frame,
        std::size_t facet_list_index);

    static void collectImmediateCofacetRanks(
        DimensionFrame& frame,
        std::int64_t facet_bindex);

    static std::int64_t findCompressedTerminalCofacet(
        DimensionFrame& frame,
        const std::vector<std::size_t>& start_cofacet_indices);

    template <typename MinQueue>
    static void enqueueReducedCompressedColumnTail(
        DimensionFrame& frame,
        std::size_t facet_list_index,
        std::size_t expected_pivot_cofacet,
        MinQueue& target_queue);

    static std::vector<std::size_t> collectReducedColumnSupport(
        DimensionFrame& frame,
        std::size_t terminal_cofacet,
        std::size_t expected_facet_list_index);

    template <typename MaxQueue>
    static void enqueueReducedCofacetBoundaryTail(
        DimensionFrame& frame,
        std::size_t cofacet_list_index,
        std::size_t pivot_facet_list_index,
        MaxQueue& target_queue,
        std::vector<std::size_t>& cofacet_trace);
};
