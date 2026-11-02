#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include <string>

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
    std::string g_model;
    double g_margin = .05;
    double g_memory_mib = 4096;
    bool cost_telemetry = false;
    std::vector<int> quotient_mask;
    size_t stop_after_window = 0;
    std::string ac_policy;
    std::string ac_model;
    std::string h1_risk;
    bool h1_realizations = false;
    std::string pack_mode;
    std::string pack_model;
    std::string payback_action = "auto";
    bool full_endpoint = false;
    bool calibration_targets = false;
    double payback_aL = 0, payback_aM = 0, payback_bC = 0;
    // Gate-only snapshots are written outside the compact persistence stream.
    std::string extend_dump;
    bool extend_force_chain = false;
    bool extend_direct_probe = false;
    // Keep the pre-quotient birth of a class alive across a quotient; veto a quotient
    // that two classes of one dimension would survive.
    bool birth_repair = false;
    // Veto every quotient while any class of any dimension is alive (implies birth_repair's checks).
    bool live_guard = false;
    // With live_guard: retain the first window across vetoed cuts and extend it.
    bool live_guard_extend = false;
    // Gate only: extend through every cut and compare with a direct window at each.
    bool extend_gate = false;
    // With live_guard_extend: after the first committed quotient, skip the remaining
    // cuts and close with one window up to the endpoint, as the H1 EXTEND ladder does.
    bool extend_commit_once = false;

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
