#include "Pipelines.hpp"
#include <iomanip>
#include <iostream>
#include <utility>
#include <chrono>
#include <cstdio>
#include <exception>
#include "DimensionFrame.hpp"
#include "SimplexUtility.hpp"


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
    for (std::size_t window_index = 0;
         window_index < config.eps_breaks.size();
         ++window_index)
    {
        const WindowBounds bounds{
            .eps_lo = window_index == 0
                ? 0.0
                : config.eps_breaks[window_index - 1],
            .eps_hi = config.eps_breaks[window_index],
            .is_final = window_index + 1 == config.eps_breaks.size()};

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
            DimensionFrame frame(runtime_, window_state_, std::move(edges));
            for (;;)
            {
                frame.matchPersistence(!bounds.is_final);

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

                if (!frame.advance())
                    break;
            }
        }

        // DimensionFrame was the final consumer of this window's geometry.
        // Support postprocessing needs persistent labels/PVs and eps_hi only.
        window_state_.invalidateCurrentWindow();

        if (!config.verbose)
            std::cout << '\n';

        if (bounds.is_final)
            break;
        const std::size_t active_labels_before =
        window_state_.activeLabels().size();
        const std::size_t previous_pv_count =
            window_state_.pseudoVertices().size();
        support_processor_.processPwph(
            runtime_, window_state_, bounds.eps_hi, std::move(support_batch));
        runtime_.ensureBinomialCapacity(window_state_.totalLabelCount());
        const std::size_t new_pv_count =
            window_state_.pseudoVertices().size() - previous_pv_count;

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
