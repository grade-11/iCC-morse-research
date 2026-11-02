#include "Pipelines.hpp"
#include <iomanip>
#include <iostream>
#include <utility>
#include <chrono>
#include <cstdio>
#include <exception>
#include "DimensionFrame.hpp"
#include "SimplexUtility.hpp"
#include "AdaptiveSchedule.hpp"
#include "GCostController.hpp"
#include "ACWindowController.hpp"
#include "ACRankController.hpp"
#include "DeferRetryController.hpp"
#include "GeometryShadow.hpp"
#include "PolicyCostLog.hpp"
#include "PackDmtPayback.hpp"
#include "GrowthController.hpp"
#include "PaybackAnalyticController.hpp"
#include <fstream>
#include <bit>
#include <algorithm>
#include <optional>
#include <stdexcept>


namespace
{
class PwphWindowResources
{
    using Clock = std::chrono::steady_clock;

    struct RssSnapshot
    {
        long rss_kib = -1;
        long peak_kib = -1;

        bool valid() const
        {
            return rss_kib >= 0 && peak_kib >= 0;
        }
    };

    static RssSnapshot readRss()
    {
        RssSnapshot result;
        std::FILE* file = std::fopen("/proc/self/status", "r");
        if (file == nullptr)
            return result;

        char line[512];
        while (std::fgets(line, sizeof(line), file) != nullptr)
        {
            long value = -1;
            if (std::sscanf(line, "VmRSS: %ld", &value) == 1)
                result.rss_kib = value;
            else if (std::sscanf(line, "VmHWM: %ld", &value) == 1)
                result.peak_kib = value;
        }

        const bool read_ok = std::ferror(file) == 0;
        std::fclose(file);
        return read_ok ? result : RssSnapshot{};
    }

    static bool resetPeak()
    {
        std::FILE* file = std::fopen("/proc/self/clear_refs", "w");
        if (file == nullptr)
            return false;

        const bool write_ok = std::fputs("5\n", file) >= 0;
        const bool close_ok = std::fclose(file) == 0;
        return write_ok && close_ok;
    }

public:
    PwphWindowResources(bool enabled, double eps_lo, double eps_hi)
        : enabled_(enabled), eps_lo_(eps_lo), eps_hi_(eps_hi)
    {
        if (!enabled_)
            return;

        initial_rss_ = readRss();
        if (initial_rss_.valid())
            reset_ok_ = resetPeak();

        start_ = Clock::now();
    }

    PwphWindowResources(const PwphWindowResources&) = delete;
    PwphWindowResources& operator=(const PwphWindowResources&) = delete;

    ~PwphWindowResources()
    {
        if (!enabled_ || std::uncaught_exceptions() > 0)
            return;

        const auto end = Clock::now();
        const auto final_rss = readRss();

        const double elapsed_ms =
            std::chrono::duration<double, std::milli>(
                end - start_).count();

        const bool memory_ok =
            initial_rss_.valid() && reset_ok_ && final_rss.valid();

        const char* memory_status =
            !initial_rss_.valid() ? "read-start-failed" :
            !reset_ok_ ? "reset-failed" :
            !final_rss.valid() ? "read-end-failed" : "ok";

        std::cout << "[window-resources]"
                  << " eps_lo=" << eps_lo_
                  << " eps_hi=" << eps_hi_
                  << " elapsed_ms=" << elapsed_ms
                  << " rss_start_kib=" << initial_rss_.rss_kib
                  << " rss_end_kib=" << final_rss.rss_kib
                  << " rss_peak_kernel_kib="
                  << (memory_ok ? final_rss.peak_kib : -1)
                  << " memory_status=" << memory_status << '\n';
    }

private:
    bool enabled_;
    double eps_lo_;
    double eps_hi_;
    RssSnapshot initial_rss_;
    bool reset_ok_ = false;
    Clock::time_point start_{};
};
void printSimplexLabels(
    const BinomialTable& binomial_table,
    const char* const prefix,
    const SimplexBindex simplex_bindex,
    const std::size_t total_label_count,
    const std::size_t simplex_dimension)
{
    std::cout << "    " << prefix << " labels: ";
    if (simplex_bindex < 0)
    {
        std::cout << "(none)\n";
        return;
    }

    const auto labels = SimplexUtility::getSimplexVertices(
        binomial_table,
        simplex_bindex,
        total_label_count,
        simplex_dimension);

    if (labels.empty())
    {
        std::cout << "(empty)";
    }
    else
    {
        for (std::size_t i = 0; i < labels.size(); ++i)
        {
            if (i > 0)
                std::cout << ' ';
            std::cout << labels[i];
        }
    }
    std::cout << '\n';
}

void printBirthDeath(const PersistentPairInfo& pair)
{
    const auto old_precision = std::cout.precision();
    std::cout << std::setprecision(17)
              << pair.facet_weight << ", " << pair.cofacet_weight;
    std::cout.precision(old_precision);
}

void printPersistencePair(
    const PersistentPairInfo& pair,
    const BinomialTable& binomial_table,
    const std::size_t dimension,
    const std::size_t total_label_count)
{
    std::cout << "  (";
    printBirthDeath(pair);
    std::cout << ")\n";

    printSimplexLabels(
        binomial_table,
        "birth facet",
        pair.facet_bindex,
        total_label_count,
        dimension - 1);
    printSimplexLabels(
        binomial_table,
        "death cofacet",
        pair.cofacet_bindex,
        total_label_count,
        dimension);
}
}

PersistentHomologyPipeline::PersistentHomologyPipeline(
    DistanceMatrix distance_matrix,
    PipelineConfig config)
    : runtime_(
          std::move(distance_matrix),
          std::move(config),
          PipelineMode::RegVRPH)
{
}

void PersistentHomologyPipeline::run()
{
    WindowState window(runtime_.distanceMatrix().vertexCount());

    const WindowBounds bounds{
        .eps_lo = 0.0,
        .eps_hi = runtime_.config().finalEpsilon(),
        .is_final = true};
    auto edges = window.prepareWindow(runtime_, bounds);

    std::cout << "total point number: "
              << runtime_.distanceMatrix().vertexCount() << std::endl;

    DimensionFrame frame(runtime_, window, std::move(edges));
    for (;;)
    {
        frame.traceComplex();
        frame.matchPersistence(false);

        if (frame.dimension() >= 2)
        {
            std::cout << "Processing dim " << frame.dimension()
                      << ", cofacet number: " << frame.cofacetCount()
                      << ", facet number: " << frame.facetCount()
                      << '\n'
                      << "dimensional persistent pairs:" << std::endl;

            if (frame.persistentPairs().empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (const auto& pair : frame.persistentPairs())
                {
                    if (runtime_.config().verbose)
                    {
                        printPersistencePair(
                            pair,
                            runtime_.binomialTable(),
                            frame.dimension(),
                            window.totalLabelCount());
                    }
                    else
                    {
                        std::cout << "  (";
                        printBirthDeath(pair);
                        std::cout << ")\n";
                    }
                }
            }
        }

        if (!frame.advance())
            break;
    }
}

PwphPipeline::PwphPipeline(
    DistanceMatrix distance_matrix,
    PipelineConfig config)
    : runtime_(
          std::move(distance_matrix),
          std::move(config),
          PipelineMode::PwPH),
      window_state_(runtime_.distanceMatrix().vertexCount())
{
}

void PwphPipeline::run()
{
    if(runtime_.config().ac_policy=="guard-extend-ladder" || runtime_.config().extend_direct_probe)
    {
        runGuardExtendLadder();
        return;
    }
    const auto& config = runtime_.config();
    auto policy_cost = PolicyCost::Logger::create(runtime_.distanceMatrix(),config);
    if(policy_cost)policy_cost->mark("GEO_PREPASS_BEGIN");
    auto geometry_shadow = GeoShadow::Logger::create(runtime_.distanceMatrix(),config);
    if(policy_cost){policy_cost->mark("GEO_PREPASS_END");policy_cost->mark("OPENING_READY");}
    const bool g_enabled=!config.g_model.empty();
    const bool ac_enabled=!config.ac_policy.empty();
    const bool retry_control_enabled=config.ac_policy=="guard-defer-retry-control";
    const bool defer_retry_enabled=config.ac_policy=="guard-defer-retry" || retry_control_enabled;
    bool retry_control_forced_defer=false;
    const bool payback_enabled=config.ac_policy.starts_with("payback");
    const bool growth_enabled=config.ac_policy.starts_with("growth") || config.ac_policy.starts_with("payback");
    const bool rank_enabled=ac_enabled && config.ac_policy!="cost" && !growth_enabled && !defer_retry_enabled;
    const bool guard_enabled=config.ac_policy.starts_with("guard") || growth_enabled;
    const bool measure=config.adaptive || config.cost_telemetry || g_enabled || payback_enabled || (ac_enabled && !rank_enabled && !growth_enabled);
    ACWindowController ac_controller;
    if(ac_enabled && !rank_enabled && !growth_enabled && !defer_retry_enabled)ac_controller.load(config.ac_model);
    DeferRetryController defer_retry_controller;
    if(defer_retry_enabled)defer_retry_controller.load(config.ac_model,runtime_.distanceMatrix(),config.finalEpsilon());
    GrowthController growth_controller;
    if(growth_enabled)growth_controller.load(config.ac_model,runtime_.distanceMatrix(),config.finalEpsilon());
    PaybackAnalyticController payback_controller;
    PaybackAnalyticController::Decision payback_decision;
    if(payback_enabled)payback_controller.load(config.ac_model);
    ACRankController rank_controller;
    if(rank_enabled)rank_controller.load(config.ac_model,runtime_.distanceMatrix(),config.finalEpsilon());
    ACWindowController::Decision ac_decision;
    GCostController controller;
    if (g_enabled) controller.load(config.g_model,config.maxdim);
    double previous_prediction=0;
    const auto peak_kib=[]()
    {
        std::ifstream in("/proc/self/status");std::string line;
        while(std::getline(in,line)) if(line.starts_with("VmHWM:"))
            return std::stod(line.substr(6));
        return 0.;
    };
    using Clock = std::chrono::steady_clock;
    struct Observation
    {
        double lo, hi, next = -1, reduction = 0;
        size_t active_before = 0, active_after = 0, top_cofacets = 0;
        double matching_ms = 0, rebuild_ms = 0, elapsed_ms = 0, controller_us = 0;
        double work_s=0,q_s=0,peak=0,native_s=0,feature_H=0,feature_s=0,commit_bookkeeping_s=0;
        size_t pv_before=0,pv_after=0,edges_built=0,facets_stored=0,raw=0,entries=0,protected_count=0;
        uint64_t state_hash=0;
        std::string action="FINAL";
        GCostController::Decision decision;
    };
    std::vector<Observation> observations;
    std::vector<bool> guard_quotients;
    // Birth repair. After a quotient at q the next window is rebuilt from scale zero, and the
    // witness weights preserve the filtration only from q on, so a class alive across q is
    // re-born in the quotient complex at a wrong scale below q. Its death is unaffected. With at
    // most one such class per dimension the class is identified by its new birth lying below q,
    // and the birth seen before the quotient is printed instead.
    std::vector<std::optional<double>> carried_birth(config.maxdim+1);
    std::vector<std::size_t> carried_seen(config.maxdim+1);
    std::vector<bool> carried_died(config.maxdim+1);
    std::vector<std::vector<double>> live_births(config.maxdim+1);
    double repair_cut=-1;
    double eps_lo = 0.0;
    double eps_hi = defer_retry_enabled?defer_retry_controller.first:growth_enabled?growth_controller.first:rank_enabled?rank_controller.first:config.eps_breaks.front();
    std::size_t first_window_index = 0;
    if(config.live_guard_extend)
    {
        const auto committed_cut = runLiveGuardExtendPrefix();
        if(!committed_cut)
            return;
        // Nothing was alive at the committed cut; the repair checks stay armed.
        first_window_index = config.extend_commit_once ? config.eps_breaks.size() - 1 : *committed_cut + 1;
        repair_cut = eps_lo = config.eps_breaks[*committed_cut];
        eps_hi = config.eps_breaks[first_window_index];
        // Gain rule: a commit that removes too little predicted work ends the cutting.
        if((config.extend_min_gain>0 || config.extend_cost_rule) && !support_processor_.lastKeepCutting())
            eps_hi = config.finalEpsilon();
    }
    for (std::size_t window_index = first_window_index;; ++window_index)
    {
        if(geometry_shadow)geometry_shadow->noteWindow();
        const WindowBounds bounds{
            .eps_lo = eps_lo,
            .eps_hi = eps_hi,
            .is_final = eps_hi == config.finalEpsilon()};
        if(policy_cost)policy_cost->setWindow(window_index,bounds.eps_lo,bounds.eps_hi);
        Observation observation{.lo=eps_lo, .hi=eps_hi};
        const auto feature_start=Clock::now();
        if(config.calibration_targets)observation.feature_H=PackDmtPayback::actualH(runtime_.distanceMatrix(),window_state_,{},bounds.eps_hi,
                                                     bounds.is_final && config.full_endpoint);
        observation.feature_s=std::chrono::duration<double>(Clock::now()-feature_start).count();
        const auto window_start = Clock::now();
        double interval_io_s=0;
        observation.active_before = window_state_.activeLabels().size();
        observation.pv_before = window_state_.pseudoVertices().size();

        PwphWindowResources window_resources(
            config.verbose, bounds.eps_lo, bounds.eps_hi);
            
        const auto header_io_start=Clock::now();
        if (!config.verbose)
        {
            std::cout << "window lower bound, window upper bound\n";
            std::cout << std::setprecision(17) << bounds.eps_lo << ", " << bounds.eps_hi << '\n' << std::setprecision(6);
            std::cout << "dimension, birth weight, death weight\n";
        }

        interval_io_s+=std::chrono::duration<double>(Clock::now()-header_io_start).count();
        using Clock = std::chrono::steady_clock;
        if(policy_cost)policy_cost->mark("EDGE_REBUILD_BEGIN");
        const auto edges_start = Clock::now();
        auto edges = window_state_.prepareWindow(runtime_, bounds);
        const auto edges_end = Clock::now();
        if(policy_cost)policy_cost->mark("EDGE_REBUILD_END");
        if (measure) observation.edges_built=edges.size();
        if (measure)
            observation.rebuild_ms += std::chrono::duration<double, std::milli>(edges_end-edges_start).count();

        if (config.verbose)
        {
            const double edges_ms =
                std::chrono::duration<double, std::milli>(
                    edges_end - edges_start).count();

            std::cout << "[rebuild-timing] stage=edges"
                      << " eps_lo=" << bounds.eps_lo
                      << " eps_hi=" << bounds.eps_hi
                      << " edges_ms=" << edges_ms << '\n';
        }

        if(config.birth_repair)
            for(std::size_t d=0;d<=config.maxdim;++d)
            {carried_seen[d]=0;carried_died[d]=false;live_births[d].clear();}

        DependencySupportBatch support_batch;
        {
            const auto frame_start = Clock::now();
            if(policy_cost)policy_cost->mark("FRAME_BUILD_BEGIN",2);
            DimensionFrame frame(runtime_, window_state_, std::move(edges));
            if(policy_cost)policy_cost->mark("FRAME_BUILD_END",2);
            if (measure)
                observation.rebuild_ms += std::chrono::duration<double, std::milli>(Clock::now()-frame_start).count();
            for (;;)
            {
                frame.traceComplex();
                if(policy_cost)policy_cost->mark("MATCH_BEGIN",2);
                const auto match_start = Clock::now();
                frame.matchPersistence(!bounds.is_final);
                if(policy_cost)policy_cost->mark("MATCH_END",2);
                if (measure)
                {
                    observation.matching_ms += std::chrono::duration<double, std::milli>(Clock::now()-match_start).count();
                    if (frame.isTopDimension()) {
                        observation.top_cofacets = frame.cofacetCount();
                        observation.facets_stored = frame.facetCount();
                    }
                }

                const auto interval_io_start=Clock::now();
                if (frame.dimension() >= 2)
                {
                    if (config.verbose)
                    {
                        std::cout << std::setprecision(17) << "in eps range " << bounds.eps_lo << "  "
                                  << bounds.eps_hi << "    dimension = "
                                  << frame.dimension()
                                  << "  cofacet num = " << frame.cofacetCount()
                                  << "  facet num = " << frame.facetCount() << '\n'
                                  << "   persistent pairs:" << std::endl;
                    }

                    bool printed_any = false;
                    for (const auto& pair : frame.persistentPairs())
                    {
                        if (!(pair.cofacet_weight >= bounds.eps_lo ||
                              pair.cofacet_weight < 0.0))
                        {
                            continue;
                        }

                        PersistentPairInfo shown = pair;
                        if (config.birth_repair)
                        {
                            const std::size_t d = frame.dimension();
                            if (repair_cut >= 0.0 && pair.facet_weight < repair_cut)
                            {
                                if (!carried_birth[d] || ++carried_seen[d] > 1)
                                    throw std::runtime_error("birth repair: unexpected class born before the last quotient");
                                shown.facet_weight = *carried_birth[d];
                                carried_died[d] = pair.cofacet_weight >= 0.0;
                            }
                            if (pair.cofacet_weight < 0.0)
                                live_births[d].push_back(shown.facet_weight);
                        }

                        if (config.verbose)
                        {
                            printPersistencePair(
                                shown,
                                runtime_.binomialTable(),
                                frame.dimension(),
                                window_state_.totalLabelCount());
                        }
                        else
                        {
                            std::cout << frame.dimension() << ", ";
                            printBirthDeath(shown);
                            std::cout << '\n';
                        }

                        printed_any = true;
                    }

                    if (config.verbose && !printed_any)
                    {
                        std::cout
                            << "  (no new interval or surviving interval from previous eps range)"
                            << std::endl;
                    }
                }

                interval_io_s+=std::chrono::duration<double>(Clock::now()-interval_io_start).count();
                if (frame.isTopDimension() && !bounds.is_final)
                    support_batch = frame.takeDependencySupportBatch();

                if(policy_cost)policy_cost->mark("ADVANCE_BEGIN",2);
                const auto rebuild_start = Clock::now();
                const bool advanced = frame.advance();
                if(policy_cost)policy_cost->mark("ADVANCE_END",2);
                if (measure)
                    observation.rebuild_ms += std::chrono::duration<double, std::milli>(Clock::now()-rebuild_start).count();
                if (!advanced)
                    break;
            }
        }

        observation.protected_count=support_batch.protected_labels.size();
        if(config.birth_repair)
            for(std::size_t d=0;d<=config.maxdim;++d)
            {
                if(carried_birth[d] && carried_seen[d]!=1)
                    throw std::runtime_error("birth repair: class alive across the last quotient is missing");
                if(carried_died[d])carried_birth[d].reset();
            }
        // DimensionFrame was the final consumer of this window's geometry.
        // Support postprocessing needs persistent labels/PVs and eps_hi only.
        window_state_.invalidateCurrentWindow();
        if(policy_cost)policy_cost->mark("WINDOW_CORE_END");

        const auto blank_io_start=Clock::now();
        if (!config.verbose)
            std::cout << '\n';
        interval_io_s+=std::chrono::duration<double>(Clock::now()-blank_io_start).count();

        if (measure)
        {
            observation.work_s=std::chrono::duration<double>(Clock::now()-window_start).count();
            observation.native_s=observation.work_s-interval_io_s;
            observation.raw=support_batch.supports.size();
            if(config.cost_telemetry || g_enabled || (ac_enabled && !rank_enabled && !growth_enabled))
                observation.peak=peak_kib();
            // Full support hashing is diagnostic bookkeeping, not a model feature.
            // Compact deployment still charges every statistic actually used.
            if(config.cost_telemetry || g_enabled)
            {
            uint64_t hash=1469598103934665603ULL;
            const auto add=[&](uint64_t value){hash^=value;hash*=1099511628211ULL;};
            add(std::bit_cast<uint64_t>(bounds.eps_hi));
            for(auto label:window_state_.activeLabels())add(label);
            for(const auto& pv:window_state_.pseudoVertices())
            {add(pv.representatives.size());for(auto label:pv.representatives)add(label);}
            for(const auto& support:support_batch.supports)
            {add(support.labels.size());for(auto label:support.labels)add(label);observation.entries+=support.labels.size();}
            uint64_t protected_hash=0;for(auto label:support_batch.protected_labels)protected_hash^=(label+1)*1099511628211ULL;
            add(protected_hash);add(support_batch.protected_labels.size());observation.state_hash=hash;
            }
            if (g_enabled && previous_prediction>0)controller.update(observation.work_s,previous_prediction);
        }

        if (bounds.is_final || (config.stop_after_window && window_index+1==config.stop_after_window))
        {
            if (measure)
            {
                observation.active_after = window_state_.activeLabels().size();
                observation.pv_after = window_state_.pseudoVertices().size();
                observation.elapsed_ms = std::chrono::duration<double, std::milli>(Clock::now()-window_start).count();
                observation.action=bounds.is_final?"FINAL":"STOP";
                observations.push_back(observation);
            }
            break;
        }
        const std::size_t active_labels_before =
        window_state_.activeLabels().size();
        const std::size_t previous_pv_count =
            window_state_.pseudoVertices().size();
        if(config.cost_telemetry)
        {
            // Diagnostic digest of selected semantic fields, not a heap snapshot.
            std::cout<<std::setprecision(17)<<"[retry-preaction] q="<<bounds.eps_hi
                <<" bits="<<std::bit_cast<uint64_t>(bounds.eps_hi)
                <<" active="<<window_state_.activeLabels().size()
                <<" pv="<<window_state_.pseudoVertices().size()
                <<" raw="<<support_batch.supports.size()
                <<" protected_count="<<support_batch.protected_labels.size()
                <<" state_hash="<<observation.state_hash<<'\n'<<std::setprecision(6);
            if(config.trace_pwph)
            {
                for(size_t id=0;id<support_batch.supports.size();++id)
                {
                    std::cout<<"[retry-preaction-support] order="<<id<<" members=";
                    for(auto label:support_batch.supports[id].labels)std::cout<<label<<',';
                    std::cout<<'\n';
                }
                std::vector<size_t> protected_ids(support_batch.protected_labels.begin(),support_batch.protected_labels.end());
                std::sort(protected_ids.begin(),protected_ids.end());
                std::cout<<"[retry-preaction-protected] members=";
                for(auto label:protected_ids)std::cout<<label<<',';
                std::cout<<'\n';
            }
        }
        bool do_quotient=config.quotient_mask.empty() || config.quotient_mask[window_index]!=0;
        if(guard_enabled)
            do_quotient=support_batch.protected_labels.empty() && !support_batch.supports.empty();
        if (g_enabled)
        {
            observation.decision=controller.decide(observation.active_before,observation.top_cofacets,
                observation.pv_before,observation.raw,observation.entries,bounds.eps_hi,
                config.eps_breaks[window_index+1],observation.peak,config.g_memory_mib*1024.,config.g_margin);
            do_quotient=observation.decision.quotient;
            previous_prediction=do_quotient?observation.decision.r+observation.decision.f:observation.decision.cont;
        }
        // Diagnostic causal arm: identical early path, planner and retry prefix;
        // suppress only the second nonfinal action after common guard checks.
        if(retry_control_enabled && window_index==1)
        {do_quotient=false;retry_control_forced_defer=true;}
        if(config.birth_repair)
        {
            bool crowded=false;
            const std::size_t allowed=config.live_guard?0:1;
            for(const auto& births:live_births)crowded=crowded || births.size()>allowed;
            if(crowded)do_quotient=false;
            std::cout<<"[birth-repair] epsilon="<<std::setprecision(17)<<bounds.eps_hi<<std::setprecision(6)<<" live=";
            for(std::size_t d=2;d<=config.maxdim;++d)std::cout<<(d>2?",":"")<<live_births[d].size();
            std::cout<<" action="<<(crowded?"VETO":"ALLOW")<<'\n';
        }
        observation.action=do_quotient?"QUOTIENT":"DEFER";
        if(guard_enabled)guard_quotients.push_back(do_quotient);
        if(geometry_shadow)geometry_shadow->firstDecision(bounds.eps_hi,do_quotient,
            window_state_.activeLabels().size(),window_state_.pseudoVertices().size(),window_state_.totalLabelCount());
        const auto quotient_start=Clock::now();
        if (do_quotient)
        {
            const auto committed=support_processor_.processPwph(runtime_,window_state_,bounds.eps_hi,std::move(support_batch));
            if(config.birth_repair && committed)
            {
                repair_cut=bounds.eps_hi;
                for(std::size_t d=0;d<=config.maxdim;++d)
                    carried_birth[d]=live_births[d].empty()?std::nullopt:std::optional<double>(live_births[d].front());
            }
            if(committed){
                if(policy_cost)policy_cost->mark("BINOMIAL_BEGIN");
                const auto bookkeeping_start=Clock::now();
                runtime_.ensureBinomialCapacity(window_state_.totalLabelCount());
                observation.commit_bookkeeping_s=std::chrono::duration<double>(Clock::now()-bookkeeping_start).count();
                if(policy_cost)policy_cost->mark("BINOMIAL_END");
            }
            if((!config.h1_risk.empty() || defer_retry_enabled) && committed==0)
            {
                do_quotient=false;observation.action="DEFER";
                if(guard_enabled)guard_quotients.back()=false;
            }

        }
        else support_batch=DependencySupportBatch{};
        observation.q_s=std::chrono::duration<double>(Clock::now()-quotient_start).count();
        if (config.trace_pwph)
        {
            for (size_t i=0;i<window_state_.pseudoVertices().size();++i)
            {
                std::cout << "[pv-members] eps_hi=" << bounds.eps_hi << " label="
                          << window_state_.originalVertexCount()+i << " representatives=";
                for (auto v : window_state_.pseudoVertices()[i].representatives) std::cout << v << ',';
                std::cout << '\n';
            }
        }
        const std::size_t new_pv_count =
            window_state_.pseudoVertices().size() - previous_pv_count;
        if(geometry_shadow)geometry_shadow->baselineResult(new_pv_count,
            window_state_.activeLabels().size(),window_state_.pseudoVertices().size(),window_state_.totalLabelCount());
        if(policy_cost)policy_cost->actual(new_pv_count,window_state_.activeLabels().size(),
            window_state_.pseudoVertices().size(),window_state_.totalLabelCount());
        if (measure) observation.pv_after=window_state_.pseudoVertices().size();

        eps_lo = eps_hi;
        if(defer_retry_enabled)
        {
            const auto controller_start=Clock::now();
            eps_hi=defer_retry_controller.decide(do_quotient,eps_lo,config.finalEpsilon());
            if(retry_control_forced_defer)defer_retry_controller.reason="retry-control-forced-DEFER/direct-final";
            defer_retry_controller.controller_s+=std::chrono::duration<double>(Clock::now()-controller_start).count();
        }
        else if (ac_enabled)
        {
            const auto controller_start=Clock::now();
            observation.active_after=window_state_.activeLabels().size();
            if(window_index==0)
            {
                if(payback_enabled)
                {
                    payback_decision=payback_controller.decide(active_labels_before,observation.active_after,
                        observation.pv_after,observation.top_cofacets,observation.work_s,observation.q_s,
                        eps_lo,growth_controller.middle,growth_controller.late,config.finalEpsilon());
                    eps_hi=payback_decision.next;
                    growth_controller.next=eps_hi;growth_controller.reduction=payback_decision.reduction;
                    growth_controller.deferred=!do_quotient;growth_controller.reason=payback_decision.reason;
                }
                else if(growth_enabled)
                    eps_hi=growth_controller.decide(config.ac_policy,active_labels_before,observation.active_after,!do_quotient,config.finalEpsilon());
                else if(rank_enabled)
                    eps_hi=rank_controller.decide(config.ac_policy,active_labels_before,observation.active_after,
                        window_state_.originalVertexCount(),peak_kib(),config.g_memory_mib*1024.,config.finalEpsilon(),!do_quotient,observation.protected_count);
                else
                {
                    ac_decision=ac_controller.decide(window_state_.originalVertexCount(),observation.top_cofacets,
                    observation.active_after,observation.pv_after,observation.raw,observation.edges_built,
                    observation.peak,config.g_memory_mib*1024.);
                    eps_hi=ac_decision.next;
                }
            }
            else eps_hi=config.finalEpsilon();
            ac_decision.controller_s+=std::chrono::duration<double>(Clock::now()-controller_start).count();
            rank_controller.controller_s=ac_decision.controller_s;
        }
        else if (config.adaptive)
        {
            const auto controller_start = Clock::now();
            observation.active_after = window_state_.activeLabels().size();
            observation.reduction = 1.0 - static_cast<double>(observation.active_after) / active_labels_before;
            eps_hi = chooseAdaptiveEpsilon(config.eps_breaks, eps_lo, observation.reduction, config.adaptive_target);
            observation.next = eps_hi;
            observation.controller_us = std::chrono::duration<double, std::micro>(Clock::now()-controller_start).count();
            observation.elapsed_ms = std::chrono::duration<double, std::milli>(Clock::now()-window_start).count();
            observations.push_back(observation);
        }
        else eps_hi = config.eps_breaks[window_index + 1];
        if((config.extend_min_gain>0 || config.extend_cost_rule) && support_processor_.lastGainEpsilon()==eps_lo &&
           !support_processor_.lastKeepCutting())
            eps_hi=config.finalEpsilon();
        if (measure && !config.adaptive)
        {
            observation.active_after=window_state_.activeLabels().size();
            observation.next=eps_hi;
            observation.elapsed_ms=std::chrono::duration<double,std::milli>(Clock::now()-window_start).count();
            observations.push_back(observation);
        }

        if (config.verbose)
        {
            std::cout << "after eps range " << bounds.eps_lo << "  "
                      << bounds.eps_hi << "  new pv number = " << new_pv_count
                      << "  total pv number = "
                      << window_state_.pseudoVertices().size()
                      << "  active labels before quotient = " << active_labels_before
                      << "  active labels after quotient = "
                      << window_state_.activeLabels().size() << std::endl;

            std::cout << "pv representative statistics:" << std::endl;
            if (window_state_.pseudoVertices().empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (std::size_t pv_index = 0;
                     pv_index < window_state_.pseudoVertices().size();
                     ++pv_index)
                {
                    const auto& pv = window_state_.pseudoVertices()[pv_index];
                    std::cout << "  [" << pv_index << "] size = "
                              << pv.representatives.size()
                              << "  max diameter = " << pv.diameter << '\n';
                }
                std::cout << std::endl;
            }
        }
    }
    if(policy_cost){policy_cost->global();policy_cost->mark("PIPELINE_SUMMARY_BEGIN");}
    // Deferred output: no logging in the timed controller, no verbose/FFI needed.
    if(defer_retry_enabled)
        std::cout<<std::setprecision(17)<<"[defer-retry] policy=DEFER_RETRY_V1 research_only=1"
            <<" diagnostic_control="<<retry_control_enabled<<" control_forced_defer="<<retry_control_forced_defer
            <<" degree="<<defer_retry_controller.degree<<" retry_fraction="<<defer_retry_controller.retry_fraction
            <<" first="<<defer_retry_controller.first<<" retry="<<defer_retry_controller.retry
            <<" first_bits="<<std::bit_cast<uint64_t>(defer_retry_controller.first)
            <<" retry_bits="<<std::bit_cast<uint64_t>(defer_retry_controller.retry)
            <<" early_available="<<defer_retry_controller.early_available<<" retry_available="<<defer_retry_controller.retry_available
            <<" retry_scheduled="<<defer_retry_controller.retry_scheduled
            <<" actual_last_action="<<(defer_retry_controller.decisions?(defer_retry_controller.last_committed?"COMMIT":"DEFER"):"NONE")
            <<" next="<<defer_retry_controller.next<<" reason="<<defer_retry_controller.reason
            <<" positive_pairs="<<defer_retry_controller.positive_pairs<<" distinct_pairs="<<defer_retry_controller.distinct_pairs
            <<" early_rank="<<defer_retry_controller.early_rank<<" planner_s="<<defer_retry_controller.initialization_s
            <<" controller_s="<<defer_retry_controller.controller_s<<'\n'<<std::setprecision(6);
    if(guard_enabled)
    {
        std::cout << "[ac-actions] count=" << guard_quotients.size();
        for(size_t i=0;i<guard_quotients.size();++i)std::cout << " q" << i << '=' << int(guard_quotients[i]);
        std::cout << '\n';
    }
    if(rank_enabled && config.cost_telemetry)
        std::cout << std::setprecision(17) << "[ac-state] policy=" << config.ac_policy
            << " first=" << rank_controller.first << " candidate=" << rank_controller.second
            << " late=" << rank_controller.late << " first_protected=" << rank_controller.first_protected
            << " first_action=" << (rank_controller.first_deferred?"DEFER":"QUOTIENT")
            << " next=" << rank_controller.next << " reduction=" << rank_controller.reduction
            << " threshold=" << rank_controller.threshold << " initialization_s=" << rank_controller.initialization_s
            << " controller_s=" << rank_controller.controller_s << " reason=" << rank_controller.reason << '\n' << std::setprecision(6);
    if(growth_enabled)
        std::cout << std::setprecision(17) << "[growth-state] policy=" << config.ac_policy
            << " first=" << growth_controller.first << " late=" << growth_controller.late
            << " first_quantile=" << growth_controller.first_fraction << " late_quantile=" << growth_controller.late_fraction
            << " samples=" << growth_controller.samples << " initialization_s=" << growth_controller.initialization_s
            << " controller_s=" << growth_controller.controller_s << " next=" << growth_controller.next
            << " reduction=" << growth_controller.reduction << " deferred=" << int(growth_controller.deferred)
            << " closure_threshold=" << growth_controller.closure_threshold
            << " coverage_threshold=" << growth_controller.coverage_threshold
            << " placement_reason=" << growth_controller.placement_reason
            << " reason=" << growth_controller.reason << '\n' << std::setprecision(6);
    for(const auto& candidate:growth_controller.candidates)
        std::cout << std::setprecision(17) << "[growth-candidate] fraction=" << candidate.fraction
            << " epsilon=" << candidate.epsilon << " wedges=" << candidate.wedges
            << " estimated_triangles=" << candidate.triangles << " closure=" << candidate.closure
            << " coverage=" << candidate.coverage << " accepted=" << int(candidate.accepted) << '\n' << std::setprecision(6);
    if(payback_enabled)
    {
        double actual_remaining_s=0;
        for(size_t i=1;i<observations.size();++i)actual_remaining_s+=observations[i].work_s+observations[i].q_s;
        // Future observations below are output only after all decisions/work end.
        std::cout << std::setprecision(17) << "[payback-state] policy=" << config.ac_policy
            << " next=" << payback_decision.next << " selected=" << payback_decision.selected
            << " predicted0_s=" << payback_decision.seconds[0]
            << " predicted1_s=" << payback_decision.seconds[1]
            << " predicted2_s=" << payback_decision.seconds[2]
            << " predicted_active1=" << payback_decision.predicted_active[0]
            << " predicted_active2=" << payback_decision.predicted_active[1]
            << " first_reduction=" << payback_decision.reduction << " margin=" << payback_controller.margin
            << " actual_remaining_s=" << actual_remaining_s
            << " first_work_s=" << (observations.empty()?0:observations[0].work_s)
            << " first_q_s=" << (observations.empty()?0:observations[0].q_s)
            << " first_A=" << (observations.empty()?0:observations[0].active_before)
            << " first_C=" << (observations.empty()?0:observations[0].top_cofacets)
            << " first_active_after=" << (observations.empty()?0:observations[0].active_after)
            << " first_PV_after=" << (observations.empty()?0:observations[0].pv_after)
            << " first_protected=" << (observations.empty()?0:observations[0].protected_count)
            << " first_supports=" << (observations.empty()?0:observations[0].raw)
            << " model_loading_s=" << payback_controller.loading_s
            << " controller_s=" << ac_decision.controller_s
            << " reason=" << payback_decision.reason << '\n' << std::setprecision(6);
    }
    if(ac_enabled && !rank_enabled && !growth_enabled && !defer_retry_enabled && config.cost_telemetry)
    {
        std::cout << std::setprecision(17) << "[ac-state] policy=" << config.ac_policy
            << " next=" << ac_decision.next << " selected=" << ac_decision.selected
            << " controller_s=" << ac_decision.controller_s << " reason=" << ac_decision.reason;
        for(size_t k=0;k<8;++k)std::cout << " x" << k << '=' << ac_decision.x[k];
        for(size_t j=0;j<4;++j)std::cout << " t" << j << '=' << ac_decision.seconds[j] << " m" << j << '=' << ac_decision.kib[j];
        std::cout << '\n' << std::setprecision(6);
    }
    for (const auto& o : observations)
    {
        if (config.adaptive)
        std::cout << "[adaptive] eps_lo=" << o.lo << " eps_hi=" << o.hi
                  << " next=" << o.next << " reason=" << (o.next<0 ? "final" : "reduction-scaled-width")
                  << " active_before=" << o.active_before << " active_after=" << o.active_after
                  << " reduction=" << o.reduction << " target=" << config.adaptive_target
                  << " top_cofacets=" << o.top_cofacets << " matching_ms=" << o.matching_ms
                  << " rebuild_ms=" << o.rebuild_ms << " elapsed_ms=" << o.elapsed_ms
                  << " controller_us=" << o.controller_us << '\n';
        if (config.cost_telemetry || g_enabled)
            std::cout << std::setprecision(17) << "[state] eps=" << o.hi << " lo=" << o.lo
                << " h=" << o.next << " A=" << o.active_before << " C=" << o.top_cofacets
                << " pv=" << o.pv_before << " raw=" << o.raw << " entries=" << o.entries
                << " state_hash=" << o.state_hash << " native_s=" << o.native_s
                << " native_valid=" << (config.calibration_targets && !config.verbose && !config.trace_pwph && !config.h1_realizations)
                << " feature_H=" << o.feature_H << " feature_s=" << o.feature_s << " work_s=" << o.work_s
                << " rebuild_s=" << o.rebuild_ms/1000. << " residual_s=" << o.work_s-o.rebuild_ms/1000.
                << " commit_bookkeeping_s=" << o.commit_bookkeeping_s << " q_s=" << o.q_s << " active_after=" << o.active_after << " peak_kib=" << o.peak
                << " edges_built=" << o.edges_built << " facets_stored=" << o.facets_stored
                << " pv_after=" << o.pv_after << " matching_s=" << o.matching_ms/1000.
                << " protected_count=" << o.protected_count
                << " elapsed_s=" << o.elapsed_ms/1000.
                << " action=" << o.action << '\n' << std::setprecision(6);
        if (g_enabled && o.action!="FINAL" && o.action!="STOP")
            std::cout << std::setprecision(17) << "[g] eps=" << o.hi << " h=" << o.decision.horizon
                << " continue_s=" << o.decision.cont << " select_quotient_s=" << o.decision.q
                << " rebuild_s=" << o.decision.r << " future_s=" << o.decision.f
                << " G_s=" << o.decision.g << " threshold_s=" << o.decision.threshold
                << " mem_continue_kib=" << o.decision.mem_continue << " mem_quotient_kib=" << o.decision.mem_quotient
                << " action=" << o.action << " reason=" << o.decision.reason << '\n' << std::setprecision(6);
    }
    if(policy_cost){policy_cost->mark("PIPELINE_SUMMARY_END");policy_cost->mark("GEO_FINISH_BEGIN");}
    if(geometry_shadow)geometry_shadow->finish();
    if(policy_cost){policy_cost->mark("GEO_FINISH_END");policy_cost->finish();}
}

ApparentPairPipeline::ApparentPairPipeline(
    DistanceMatrix distance_matrix,
    PipelineConfig config)
    : runtime_(
          std::move(distance_matrix),
          std::move(config),
          PipelineMode::RegVRPH)
{
}

void ApparentPairPipeline::run()
{
    WindowState window(runtime_.distanceMatrix().vertexCount());

    const WindowBounds bounds{
        .eps_lo = 0.0,
        .eps_hi = runtime_.config().finalEpsilon(),
        .is_final = true};
    auto edges = window.prepareWindow(runtime_, bounds);

    std::cout << "dimension, unmatched weight\n";

    DimensionFrame frame(runtime_, window, std::move(edges));
    for (;;)
    {
        frame.matchApparentPairsOnly();

        for (const double weight : frame.unmatchedFacetWeights())
            std::cout << frame.facetDimension() << ", " << weight << '\n';

        if (frame.isTopDimension() && frame.dimension() >= 2)
        {
            for (const double weight : frame.unmatchedTopCofacetWeights())
                std::cout << frame.dimension() << ", " << weight << '\n';
        }

        if (!frame.advance())
            break;
    }
}
