#include "Pipelines.hpp"
#include <iomanip>
#include <iostream>
#include <utility>

#include "DimensionFrame.hpp"
#include "SimplexUtility.hpp"

namespace
{
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

        if (!config.verbose)
        {
            std::cout << "window lower bound, window upper bound\n";
            std::cout << bounds.eps_lo << ", " << bounds.eps_hi << '\n';
            std::cout << "dimension, birth weight, death weight\n";
        }

        auto edges = window_state_.prepareWindow(runtime_, bounds);

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
