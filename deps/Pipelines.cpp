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
#include <fstream>
#include <bit>


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
    const auto& config = runtime_.config();
    const bool g_enabled=!config.g_model.empty();
    const bool measure=config.adaptive || config.cost_telemetry || g_enabled;
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
        double work_s=0,q_s=0,peak=0;
        size_t pv_before=0,raw=0,entries=0;
        uint64_t state_hash=0;
        std::string action="FINAL";
        GCostController::Decision decision;
    };
    std::vector<Observation> observations;
    double eps_lo = 0.0;
    double eps_hi = config.eps_breaks.front();
    for (std::size_t window_index = 0;; ++window_index)
    {
        const WindowBounds bounds{
            .eps_lo = eps_lo,
            .eps_hi = eps_hi,
            .is_final = eps_hi == config.finalEpsilon()};
        const auto window_start = Clock::now();
        Observation observation{.lo=eps_lo, .hi=eps_hi};
        observation.active_before = window_state_.activeLabels().size();
        observation.pv_before = window_state_.pseudoVertices().size();

        PwphWindowResources window_resources(
            config.verbose, bounds.eps_lo, bounds.eps_hi);
            
        if (!config.verbose)
        {
            std::cout << "window lower bound, window upper bound\n";
            std::cout << bounds.eps_lo << ", " << bounds.eps_hi << '\n';
            std::cout << "dimension, birth weight, death weight\n";
        }

        using Clock = std::chrono::steady_clock;
        const auto edges_start = Clock::now();
        auto edges = window_state_.prepareWindow(runtime_, bounds);
        const auto edges_end = Clock::now();
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

        DependencySupportBatch support_batch;
        {
            const auto frame_start = Clock::now();
            DimensionFrame frame(runtime_, window_state_, std::move(edges));
            if (measure)
                observation.rebuild_ms += std::chrono::duration<double, std::milli>(Clock::now()-frame_start).count();
            for (;;)
            {
                frame.traceComplex();
                const auto match_start = Clock::now();
                frame.matchPersistence(!bounds.is_final);
                if (measure)
                {
                    observation.matching_ms += std::chrono::duration<double, std::milli>(Clock::now()-match_start).count();
                    if (frame.isTopDimension()) observation.top_cofacets = frame.cofacetCount();
                }

                if (frame.dimension() >= 2)
                {
                    if (config.verbose)
                    {
                        std::cout << "in eps range " << bounds.eps_lo << "  "
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

                        if (config.verbose)
                        {
                            printPersistencePair(
                                pair,
                                runtime_.binomialTable(),
                                frame.dimension(),
                                window_state_.totalLabelCount());
                        }
                        else
                        {
                            std::cout << frame.dimension() << ", ";
                            printBirthDeath(pair);
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

                if (frame.isTopDimension() && !bounds.is_final)
                    support_batch = frame.takeDependencySupportBatch();

                const auto rebuild_start = Clock::now();
                const bool advanced = frame.advance();
                if (measure)
                    observation.rebuild_ms += std::chrono::duration<double, std::milli>(Clock::now()-rebuild_start).count();
                if (!advanced)
                    break;
            }
        }

        // DimensionFrame was the final consumer of this window's geometry.
        // Support postprocessing needs persistent labels/PVs and eps_hi only.
        window_state_.invalidateCurrentWindow();

        if (!config.verbose)
            std::cout << '\n';

        if (measure)
        {
            observation.work_s=std::chrono::duration<double>(Clock::now()-window_start).count();
            observation.peak=peak_kib();
            observation.raw=support_batch.supports.size();
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
            if (g_enabled && previous_prediction>0)controller.update(observation.work_s,previous_prediction);
        }

        if (bounds.is_final || (config.stop_after_window && window_index+1==config.stop_after_window))
        {
            if (measure)
            {
                observation.active_after = window_state_.activeLabels().size();
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
        bool do_quotient=config.quotient_mask.empty() || config.quotient_mask[window_index]!=0;
        if (g_enabled)
        {
            observation.decision=controller.decide(observation.active_before,observation.top_cofacets,
                observation.pv_before,observation.raw,observation.entries,bounds.eps_hi,
                config.eps_breaks[window_index+1],observation.peak,config.g_memory_mib*1024.,config.g_margin);
            do_quotient=observation.decision.quotient;
            previous_prediction=do_quotient?observation.decision.r+observation.decision.f:observation.decision.cont;
        }
        observation.action=do_quotient?"QUOTIENT":"DEFER";
        const auto quotient_start=Clock::now();
        if (do_quotient)
        {
            support_processor_.processPwph(runtime_,window_state_,bounds.eps_hi,std::move(support_batch));
            runtime_.ensureBinomialCapacity(window_state_.totalLabelCount());
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

        eps_lo = eps_hi;
        if (config.adaptive)
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
    // Deferred output: no logging in the timed controller, no verbose/FFI needed.
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
                << " state_hash=" << o.state_hash << " work_s=" << o.work_s
                << " rebuild_s=" << o.rebuild_ms/1000. << " residual_s=" << o.work_s-o.rebuild_ms/1000.
                << " q_s=" << o.q_s << " active_after=" << o.active_after << " peak_kib=" << o.peak
                << " action=" << o.action << '\n' << std::setprecision(6);
        if (g_enabled && o.action!="FINAL" && o.action!="STOP")
            std::cout << std::setprecision(17) << "[g] eps=" << o.hi << " h=" << o.decision.horizon
                << " continue_s=" << o.decision.cont << " select_quotient_s=" << o.decision.q
                << " rebuild_s=" << o.decision.r << " future_s=" << o.decision.f
                << " G_s=" << o.decision.g << " threshold_s=" << o.decision.threshold
                << " mem_continue_kib=" << o.decision.mem_continue << " mem_quotient_kib=" << o.decision.mem_quotient
                << " action=" << o.action << " reason=" << o.decision.reason << '\n' << std::setprecision(6);
    }
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
