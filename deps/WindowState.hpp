#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "PipelineCommon.hpp"
#include "robin_hood.h"

struct SelectedPV
{
    // Canonical sorted original-vertex indices represented by this PV.
    std::vector<size_t> representatives;
    double diameter = 0.0;
};

class WindowState
{
public:
    explicit WindowState(size_t original_vertex_count);

    [[nodiscard]]
    size_t originalVertexCount() const noexcept
    {
        return original_vertex_count_;
    }

    [[nodiscard]]
    size_t totalLabelCount() const noexcept
    {
        return original_vertex_count_ + pv_list_.size();
    }

    [[nodiscard]]
    const std::vector<SelectedPV>& pseudoVertices() const noexcept
    {
        return pv_list_;
    }

    [[nodiscard]]
    const std::vector<size_t>& activeLabels() const noexcept
    {
        return active_label_list_;
    }

    [[nodiscard]]
    const std::vector<uint8_t>& activeLabelMask() const noexcept
    {
        return active_label_mask_;
    }

    [[nodiscard]]
    const robin_hood::unordered_map<LabelPairKey, double>&
    pvLabelDistanceHash() const noexcept
    {
        return pv_label_distance_hash_;
    }

    [[nodiscard]]
    const WindowBounds& bounds() const noexcept
    {
        return bounds_;
    }

    SimplexList prepareWindow(
        const PipelineRuntime& runtime,
        const WindowBounds& bounds);

    // No quotient has occurred: widening changes no label geometry or caches.
    void extendFirstWindowTo(double eps_hi);

    // EXTEND of any prepared PwPH window, with or without pseudo vertices: prepareWindow caches the
    // distance of every active pair involving a PV, so widening changes no label geometry.
    void extendPreparedWindowTo(double eps_hi, bool is_final);

    // Quotient edges of the prepared window weighing in [eps_lo, eps_hi), in prepareWindow's order.
    [[nodiscard]]
    SimplexList quotientEdgeTail(const PipelineRuntime& runtime, double eps_lo, double eps_hi) const;

    // End the prepared-window geometry lifetime after its DimensionFrame dies.
    // Persistent active-label and pseudo-vertex state remains available for
    // support postprocessing and the next window.
    void invalidateCurrentWindow() noexcept;

    void commitSelectedPVs(
        std::vector<SelectedPV>&& selected_pvs,
        const std::unordered_set<size_t>& new_absorbed_labels);

private:
    size_t original_vertex_count_ = 0;

    // Append-only for the complete pipeline run.
    std::vector<SelectedPV> pv_list_;

    // Sorted labels active in the current quotient complex.
    std::vector<size_t> active_label_list_;

    // Current-window caches, all defined over the complete historical label range.
    std::vector<uint8_t> active_label_mask_;
    robin_hood::unordered_map<LabelPairKey, double> pv_label_distance_hash_;

    WindowBounds bounds_;
};
