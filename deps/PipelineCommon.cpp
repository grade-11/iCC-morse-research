#include "PipelineCommon.hpp"

#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

#include "SimplexUtility.hpp"

namespace
{
constexpr size_t MAX_SIMPLEX_DIMENSION =
    std::numeric_limits<uint64_t>::digits - 1;
}

void PipelineConfig::validate() const
{
    if(h1_realizations && (maxdim!=2 || !verbose))
        throw std::invalid_argument("H1 realization diagnostics require --max-dim 2 --verbose");
    // Under the live guard the H1 filter also runs in H2 computations (max-dim 3).
    if (!h1_risk.empty() && ((maxdim!=2 && !(maxdim==3 && live_guard)) || (h1_risk!="edge-fiber" && h1_risk!="selective" && h1_risk!="selective-fast" && h1_risk!="selective-graph" && h1_risk!="observe")))
        throw std::invalid_argument("H1 risk mode must be edge-fiber, selective or observe, max-dim=2 (or 3 with --live-guard)");
    const bool growth_policy=ac_policy.starts_with("growth") || ac_policy.starts_with("payback");
    const bool defer_retry_policy=ac_policy=="guard-defer-retry" || ac_policy=="guard-defer-retry-control";
    const bool extend_policy=ac_policy=="guard-extend-ladder";
    if(extend_policy || extend_direct_probe || extend_force_chain || !extend_dump.empty())
    {
        if ((!extend_policy && !extend_direct_probe) ||
            (extend_direct_probe && !ac_policy.empty()) ||
            (extend_force_chain && !extend_policy) ||
            ((extend_force_chain || extend_direct_probe) && extend_dump.empty()) ||
            maxdim!=2 || h1_risk!="selective-fast" || pv_cap_scale!=.85 || pv_min_separation!=0 ||
            !ac_model.empty() || adaptive || !g_model.empty() || !quotient_mask.empty() ||
            !pack_mode.empty() || payback_action!="auto" || stop_after_window || cost_telemetry ||
            calibration_targets || h1_realizations ||
            (extend_policy && eps_breaks.size()<2))
            throw std::invalid_argument("EXTEND requires an explicit ladder, H1 selective-fast, alpha=.85, delta=0, no conflicting scheduler; probes require --extend-dump");
    }
    if(defer_retry_policy && (h1_risk!="selective-fast" || !pack_mode.empty() || pv_cap_scale!=.85 || pv_min_separation!=0))
        throw std::invalid_argument("DEFER_RETRY requires unchanged selective-fast, alpha=.85, delta=0 and no packing mode");
    if (!extend_policy && !ac_policy.empty() && ((ac_policy!="cost" && ac_policy!="rank" && ac_policy!="rank-base" && ac_policy!="rank-extra" && ac_policy!="guard" && ac_policy!="guard-base" && ac_policy!="guard50" && ac_policy!="guard70" && !growth_policy && !defer_retry_policy) || ac_model.empty() || adaptive || !g_model.empty() || !quotient_mask.empty() ||
        ((growth_policy || defer_retry_policy) ? eps_breaks.size()!=1 : ac_policy=="cost" ? eps_breaks!=std::vector<double>{.6,.9,1.2,1.5,2.} : eps_breaks.size()!=5) || maxdim!=2))
        throw std::invalid_argument("AC cost policy requires H1 candidate grid .6,.9,1.2,1.5,2 and independent model");
    if (!std::isfinite(g_margin) || g_margin<0 || g_margin>1 || !std::isfinite(g_memory_mib) || g_memory_mib<=0)
        throw std::invalid_argument("Invalid G margin or memory budget");
    if ((!g_model.empty() && (adaptive || !quotient_mask.empty())) || (adaptive && !quotient_mask.empty()))
        throw std::invalid_argument("G, v1 and forced masks are mutually exclusive");
    if (!quotient_mask.empty() && (quotient_mask.size()+1!=eps_breaks.size() ||
        std::any_of(quotient_mask.begin(),quotient_mask.end(),[](int v){return v!=0 && v!=1;})))
        throw std::invalid_argument("Quotient mask requires one 0/1 per nonfinal checkpoint");
    if (stop_after_window>eps_breaks.size() || (stop_after_window && !cost_telemetry))
        throw std::invalid_argument("Diagnostic stop requires cost telemetry and a valid window count");
    if (!std::isfinite(adaptive_target) || adaptive_target <= 0.0 || adaptive_target > 1.0)
        throw std::invalid_argument("Adaptive reduction target must be in (0, 1].");
    if (adaptive && eps_breaks != std::vector<double>{1.2,1.4,1.5,1.6,1.7,1.8,2.0})
        throw std::invalid_argument("Adaptive pilot requires grid 1.2,1.4,1.5,1.6,1.7,1.8,2.0.");
    if (eps_breaks.empty())
        throw std::invalid_argument("Pipeline epsilon breaks must not be empty.");

    for (size_t i = 0; i < eps_breaks.size(); ++i)
    {
        const double epsilon = eps_breaks[i];
        if (!std::isfinite(epsilon) || epsilon <= 0.0)
            throw std::invalid_argument("Pipeline epsilon breaks must be finite and positive.");

        if (i > 0 && epsilon <= eps_breaks[i - 1])
            throw std::invalid_argument("Pipeline epsilon breaks must be strictly increasing.");
    }

    if (maxdim == 0 || maxdim > MAX_SIMPLEX_DIMENSION)
        throw std::invalid_argument("Pipeline maximum dimension must be in the range [1, 63].");

    if (threads <= 0)
        throw std::invalid_argument("Pipeline thread count must be positive.");

    if (!std::isfinite(pv_cap_scale) || pv_cap_scale <= 0.0)
        throw std::invalid_argument("Pipeline PV cap scale must be finite and positive.");

    if (!std::isfinite(pv_min_separation) || pv_min_separation < 0.0)
        throw std::invalid_argument("Pipeline PV minimum separation must be finite and nonnegative.");

    if (!std::isfinite(pv_cap_scale * eps_breaks.back()))
        throw std::invalid_argument("Pipeline PV cap at the final epsilon must be finite.");

    if (!std::isfinite(pv_min_separation * eps_breaks.back()))
        throw std::invalid_argument("Pipeline absolute PV minimum separation must be finite.");
}

double PipelineConfig::finalEpsilon() const
{
    return eps_breaks.back();
}

double PipelineConfig::absoluteMinSeparation() const
{
    return pv_min_separation * finalEpsilon();
}

PipelineRuntime::PipelineRuntime(
    DistanceMatrix distance_matrix,
    PipelineConfig config,
    const PipelineMode mode)
    : distance_matrix_(std::move(distance_matrix)),
      config_(std::move(config)),
      mode_(mode)
{
    config_.validate();

    switch (mode_)
    {
    case PipelineMode::RegVRPH:
        config_.eps_breaks.assign(1, config_.finalEpsilon());
        break;
    case PipelineMode::PwPH:
        break;
    default:
        throw std::invalid_argument("Unknown pipeline mode.");
    }

    binomial_table_ = SimplexUtility::getBinomialTable(
        distance_matrix_.vertexCount(),
        config_.maxdim);
}

void PipelineRuntime::ensureBinomialCapacity(const size_t total_label_count)
{
    const size_t original_vertex_count = distance_matrix_.vertexCount();
    if (binomial_table_.size() > total_label_count)
        return;

    // Grow a staged copy so an overflow or allocation failure cannot leave the
    // runtime advertising partially initialized rows.
    BinomialTable expanded_table = binomial_table_;
    SimplexUtility::updateBinomialTable(
        expanded_table,
        original_vertex_count,
        total_label_count - original_vertex_count,
        config_.maxdim);
    binomial_table_.swap(expanded_table);
}
