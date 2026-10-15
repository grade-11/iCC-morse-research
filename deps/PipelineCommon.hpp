#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "DistanceMatrix.hpp"

using SimplexBindex = int64_t;
using LabelPairKey = uint64_t;
using BinomialTable = std::vector<std::vector<int64_t>>;
using SimplexRecord = std::pair<SimplexBindex, double>;
using SimplexList = std::vector<SimplexRecord>;

struct PersistentPairInfo
{
    double facet_weight = -1.0;
    double cofacet_weight = -1.0;
    SimplexBindex facet_bindex = -1;
    SimplexBindex cofacet_bindex = -1;
};

struct PipelineConfig
{
    std::vector<double> eps_breaks;
    size_t maxdim = 0;
    int threads = 1;
    double pv_cap_scale = 1.0;
    double pv_min_separation = 0.0;
    bool verbose = false;
    bool adaptive = false;
    double adaptive_target = 0.35;
    bool trace_pwph = false;

    void validate() const;

    [[nodiscard]]
    double finalEpsilon() const;

    [[nodiscard]]
    double absoluteMinSeparation() const;
};

enum class PipelineMode : uint8_t
{
    RegVRPH,
    PwPH
};

class PipelineRuntime
{
public:
    PipelineRuntime(
        DistanceMatrix distance_matrix,
        PipelineConfig config,
        PipelineMode mode);

    [[nodiscard]]
    const DistanceMatrix& distanceMatrix() const noexcept
    {
        return distance_matrix_;
    }

    [[nodiscard]]
    const PipelineConfig& config() const noexcept
    {
        return config_;
    }

    [[nodiscard]]
    PipelineMode mode() const noexcept
    {
        return mode_;
    }

    [[nodiscard]]
    const BinomialTable& binomialTable() const noexcept
    {
        return binomial_table_;
    }

    void ensureBinomialCapacity(size_t total_label_count);

private:
    DistanceMatrix distance_matrix_;
    PipelineConfig config_;
    PipelineMode mode_;
    BinomialTable binomial_table_;
};

struct WindowBounds
{
    double eps_lo = 0.0;
    double eps_hi = 0.0;
    bool is_final = false;
};
