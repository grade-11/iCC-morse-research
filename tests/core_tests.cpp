#include "DependencySupportPostProcessor.hpp"
#include "DimensionFrame.hpp"
#include "PipelineCommon.hpp"
#include "Pipelines.hpp"
#include "SimplexUtility.hpp"
#include "WindowState.hpp"
#include "AdaptiveSchedule.hpp"
#include "GCostController.hpp"
#include "FiberRisk.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
int failure_count = 0;

void recordFailure(
    const char* expression,
    const char* file,
    const int line,
    const std::string& detail = {})
{
    ++failure_count;
    std::cerr << file << ':' << line << ": CHECK(" << expression << ") failed";
    if (!detail.empty())
        std::cerr << ": " << detail;
    std::cerr << '\n';
}

#define CHECK(expression)                                                        \
    do                                                                           \
    {                                                                            \
        if (!(expression))                                                       \
            recordFailure(#expression, __FILE__, __LINE__);                      \
    } while (false)

template <typename Exception, typename Function>
void checkThrows(
    Function&& function,
    const char* expression,
    const char* file,
    const int line)
{
    try
    {
        std::invoke(std::forward<Function>(function));
    }
    catch (const Exception&)
    {
        return;
    }
    catch (const std::exception& error)
    {
        recordFailure(
            expression,
            file,
            line,
            std::string("unexpected exception: ") + error.what());
        return;
    }
    catch (...)
    {
        recordFailure(expression, file, line, "unexpected non-standard exception");
        return;
    }

    recordFailure(expression, file, line, "expected exception was not thrown");
}

#define CHECK_THROWS_AS(expression, exception_type)                              \
    checkThrows<exception_type>(                                                 \
        [&]() { static_cast<void>(expression); },                                \
        #expression,                                                             \
        __FILE__,                                                                \
        __LINE__)

PipelineConfig makeConfig(
    std::vector<double> eps_breaks = {1.0},
    const std::size_t maxdim = 2,
    const int threads = 1,
    const double cap_scale = 1.0,
    const double separation_scale = 0.0)
{
    PipelineConfig config;
    config.eps_breaks = std::move(eps_breaks);
    config.maxdim = maxdim;
    config.threads = threads;
    config.pv_cap_scale = cap_scale;
    config.pv_min_separation = separation_scale;
    return config;
}

WindowBounds boundsFor(
    const PipelineConfig& config,
    const std::size_t index)
{
    return WindowBounds{
        .eps_lo = index == 0 ? 0.0 : config.eps_breaks[index - 1],
        .eps_hi = config.eps_breaks[index],
        .is_final = index + 1 == config.eps_breaks.size()};
}

DependencySupport support(std::initializer_list<std::size_t> labels)
{
    std::vector<std::size_t> sorted_labels(labels);
    std::sort(sorted_labels.begin(), sorted_labels.end());
    sorted_labels.erase(
        std::unique(sorted_labels.begin(), sorted_labels.end()),
        sorted_labels.end());
    return DependencySupport{.labels = std::move(sorted_labels)};
}

void processCompletedPwphWindow(
    const DependencySupportPostProcessor& processor,
    const PipelineRuntime& runtime,
    WindowState& window,
    DependencySupportBatch&& batch)
{
    const double eps_hi = window.bounds().eps_hi;
    window.invalidateCurrentWindow();
    processor.processPwph(
        runtime, window, eps_hi, std::move(batch));
}

void testPipelineConfig()
{
    const auto valid = makeConfig({1.0, 2.0}, 3, 2, 0.75, 0.25);
    valid.validate();
    CHECK(valid.finalEpsilon() == 2.0);
    CHECK(valid.absoluteMinSeparation() == 0.5);

    auto invalid = valid;
    invalid.eps_breaks.clear();
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.eps_breaks = {1.0, 1.0};
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.eps_breaks = {2.0, 1.0};
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.eps_breaks = {1.0, std::numeric_limits<double>::infinity()};
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.eps_breaks = {std::numeric_limits<double>::quiet_NaN()};
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.maxdim = 0;
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.maxdim = 64;
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.threads = 0;
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.pv_cap_scale = 0.0;
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.pv_cap_scale = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.pv_min_separation = -0.01;
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    invalid = valid;
    invalid.pv_min_separation = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(invalid.validate(), std::invalid_argument);

    PipelineRuntime ordinary_runtime(
        DistanceMatrix({{0.0}}),
        valid,
        PipelineMode::RegVRPH);
    CHECK(ordinary_runtime.mode() == PipelineMode::RegVRPH);
    CHECK((ordinary_runtime.config().eps_breaks == std::vector<double>{2.0}));

    PipelineRuntime pwph_runtime(
        DistanceMatrix({{0.0}}),
        valid,
        PipelineMode::PwPH);
    CHECK(pwph_runtime.mode() == PipelineMode::PwPH);
    CHECK((pwph_runtime.config().eps_breaks == std::vector<double>{1.0, 2.0}));

    const auto overflow_config = makeConfig({1.0}, 32);
    PipelineRuntime overflow_runtime(
        DistanceMatrix({{0.0}, {1.0}}),
        overflow_config,
        PipelineMode::RegVRPH);
    const std::size_t original_row_count =
        overflow_runtime.binomialTable().size();
    CHECK_THROWS_AS(
        overflow_runtime.ensureBinomialCapacity(70),
        std::overflow_error);
    CHECK(overflow_runtime.binomialTable().size() == original_row_count);
    CHECK_THROWS_AS(
        overflow_runtime.ensureBinomialCapacity(70),
        std::overflow_error);
    CHECK(overflow_runtime.binomialTable().size() == original_row_count);
}

void testDistanceMatrixValidation()
{
    CHECK_THROWS_AS(
        DistanceMatrix({{0.0}, {std::numeric_limits<double>::quiet_NaN()}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        DistanceMatrix({{0.0}, {std::numeric_limits<double>::infinity()}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        DistanceMatrix({{-std::numeric_limits<double>::max()},
                        {std::numeric_limits<double>::max()}}),
        std::invalid_argument);
}

void testWindowStatePreparation()
{
    const std::vector<std::vector<double>> points{
        {0.0, 0.0},
        {0.0, 0.0},
        {1.0, 0.0}};
    const auto config = makeConfig({2.0});
    PipelineRuntime ordinary_runtime(
        DistanceMatrix(points), config, PipelineMode::RegVRPH);

    WindowState ordinary(points.size());
    CHECK(ordinary.originalVertexCount() == 3);
    CHECK(ordinary.totalLabelCount() == 3);
    CHECK((ordinary.activeLabels() == std::vector<std::size_t>{0, 1, 2}));
    CHECK(ordinary.activeLabelMask().empty());
    CHECK(ordinary.pseudoVertices().empty());

    const auto ordinary_edges = ordinary.prepareWindow(
        ordinary_runtime,
        boundsFor(ordinary_runtime.config(), 0));
    CHECK(ordinary_runtime.mode() == PipelineMode::RegVRPH);
    CHECK((ordinary.activeLabelMask() == std::vector<std::uint8_t>{1, 1, 1}));
    CHECK(ordinary_edges.size() == 3);
    CHECK(std::count_if(
              ordinary_edges.begin(),
              ordinary_edges.end(),
              [](const SimplexRecord& edge) { return edge.second == 0.0; }) == 1);
    PipelineRuntime pwph_runtime(
        DistanceMatrix(points), config, PipelineMode::PwPH);
    WindowState quotient(points.size());
    const auto quotient_edges = quotient.prepareWindow(
        pwph_runtime,
        boundsFor(pwph_runtime.config(), 0));
    CHECK(pwph_runtime.mode() == PipelineMode::PwPH);
    CHECK(quotient_edges.size() == 2);
    CHECK(std::none_of(
        quotient_edges.begin(),
        quotient_edges.end(),
        [](const SimplexRecord& edge) { return edge.second == 0.0; }));

    DependencySupportPostProcessor processor;
    processCompletedPwphWindow(
        processor, pwph_runtime, quotient, DependencySupportBatch{});
    CHECK(quotient.activeLabelMask().empty());
    CHECK(quotient.pvLabelDistanceHash().empty());
}

struct FrameRun
{
    std::size_t initial_facet_count = 0;
    std::size_t initial_cofacet_count = 0;
    std::size_t promoted_facet_count = 0;
    std::size_t top_cofacet_count = 0;
    std::vector<double> first_unmatched_facets;
    std::vector<double> top_unmatched_facets;
    std::vector<double> top_unmatched_cofacets;
};

FrameRun runApparentFrame(const int threads)
{
    const std::vector<std::vector<double>> points{
        {0.0, 0.0},
        {1.0, 0.0},
        {0.2, 0.9},
        {1.3, 1.1},
        {0.6, 1.8}};
    const auto config = makeConfig({3.0}, 3, threads);
    PipelineRuntime runtime(
        DistanceMatrix(points), config, PipelineMode::RegVRPH);
    WindowState window(points.size());
    auto edges = window.prepareWindow(
        runtime,
        boundsFor(runtime.config(), 0));
    DimensionFrame frame(runtime, window, std::move(edges));

    FrameRun result;
    result.initial_facet_count = frame.facetCount();
    result.initial_cofacet_count = frame.cofacetCount();

    CHECK(frame.dimension() == 2);
    CHECK(frame.facetDimension() == 1);
    CHECK(!frame.isTopDimension());

    frame.matchApparentPairsOnly();
    result.first_unmatched_facets = frame.unmatchedFacetWeights();

    CHECK(frame.advance());
    result.promoted_facet_count = frame.facetCount();
    result.top_cofacet_count = frame.cofacetCount();
    CHECK(frame.dimension() == 3);
    CHECK(frame.facetDimension() == 2);
    CHECK(frame.isTopDimension());
    CHECK(result.promoted_facet_count == result.initial_cofacet_count);

    frame.matchApparentPairsOnly();
    result.top_unmatched_facets = frame.unmatchedFacetWeights();
    result.top_unmatched_cofacets = frame.unmatchedTopCofacetWeights();
    CHECK(!frame.advance());

    return result;
}

void testDimensionFrameAdvancementAndDeterminism()
{
    const auto single_thread = runApparentFrame(1);
    const auto multi_thread = runApparentFrame(2);

    CHECK(single_thread.initial_facet_count == multi_thread.initial_facet_count);
    CHECK(single_thread.initial_cofacet_count == multi_thread.initial_cofacet_count);
    CHECK(single_thread.promoted_facet_count == multi_thread.promoted_facet_count);
    CHECK(single_thread.top_cofacet_count == multi_thread.top_cofacet_count);
    CHECK(single_thread.first_unmatched_facets == multi_thread.first_unmatched_facets);
    CHECK(single_thread.top_unmatched_facets == multi_thread.top_unmatched_facets);
    CHECK(single_thread.top_unmatched_cofacets == multi_thread.top_unmatched_cofacets);

    CHECK(single_thread.first_unmatched_facets.size() == 1);
    if (single_thread.first_unmatched_facets.size() == 1)
    {
        CHECK(std::abs(single_thread.first_unmatched_facets.front() -
                       std::sqrt(1.3)) < 1e-12);
    }
    CHECK(single_thread.top_unmatched_facets.size() == 1);
    if (single_thread.top_unmatched_facets.size() == 1)
    {
        CHECK(std::abs(single_thread.top_unmatched_facets.front() -
                       std::sqrt(1.45)) < 1e-12);
    }
    CHECK(single_thread.top_unmatched_cofacets.size() == 1);
    if (single_thread.top_unmatched_cofacets.size() == 1)
    {
        CHECK(std::abs(single_thread.top_unmatched_cofacets.front() -
                       std::sqrt(3.6)) < 1e-12);
    }
}

void testDimensionFrameWideCofacetIndices()
{
    struct InterfaceResult
    {
        std::vector<std::pair<double, double>> persistent_pair_weights;
        std::vector<double> unmatched_facet_weights;
        std::vector<double> unmatched_top_cofacet_weights;

        bool operator==(const InterfaceResult&) const = default;
    };

    const auto run_frame = [](
        const std::size_t isolated_vertex_count,
        const bool persistence)
    {
        // Isolated vertices shift the clique's binomial indices without changing
        // its cofacet counts, filtration order, or positive-dimensional matching.
        std::vector<std::vector<double>> points;
        for (std::size_t i = 0; i < isolated_vertex_count; ++i)
            points.push_back({100.0 * static_cast<double>(i + 1), 0.0});
        points.insert(points.end(), {
            {0.0, 0.0},
            {1.0, 0.0},
            {0.2, 0.9},
            {1.3, 1.1},
            {0.6, 1.8}});

        const auto config = makeConfig({3.0}, 5);
        PipelineRuntime runtime(
            DistanceMatrix(points), config, PipelineMode::RegVRPH);
        if (isolated_vertex_count > 0)
        {
            const auto smallest_tetrahedron_bindex =
                SimplexUtility::getBinomialIndex(
                    runtime.binomialTable(),
                    {isolated_vertex_count + 3, isolated_vertex_count + 2,
                     isolated_vertex_count + 1, isolated_vertex_count},
                    0);
            CHECK(smallest_tetrahedron_bindex >
                  std::numeric_limits<std::uint32_t>::max());
        }

        WindowState window(points.size());
        auto edges = window.prepareWindow(
            runtime, boundsFor(runtime.config(), 0));
        DimensionFrame frame(runtime, window, std::move(edges));

        const std::vector<std::size_t> expected_facet_counts{10, 10, 5, 1};
        const std::vector<std::size_t> expected_cofacet_counts{10, 5, 1, 0};
        std::vector<InterfaceResult> results;
        do
        {
            const std::size_t index = frame.dimension() - 2;
            CHECK(frame.facetCount() == expected_facet_counts.at(index));
            CHECK(frame.cofacetCount() == expected_cofacet_counts.at(index));
            if (persistence)
                frame.matchPersistence(false);
            else
                frame.matchApparentPairsOnly();

            InterfaceResult result;
            for (const auto& pair : frame.persistentPairs())
            {
                result.persistent_pair_weights.emplace_back(
                    pair.facet_weight, pair.cofacet_weight);
            }
            result.unmatched_facet_weights = frame.unmatchedFacetWeights();
            result.unmatched_top_cofacet_weights =
                frame.unmatchedTopCofacetWeights();
            results.push_back(std::move(result));
        } while (frame.advance());
        CHECK(results.size() == expected_cofacet_counts.size());
        return results;
    };

    for (const bool persistence : {false, true})
        CHECK(run_frame(0, persistence) == run_frame(600, persistence));
}

void testPersistenceMatchingAndSupportExtraction()
{
    const std::vector<std::vector<double>> points{
        {0.0, 0.0},
        {1.0, 0.0},
        {0.2, 0.9},
        {1.3, 1.1},
        {0.6, 1.8}};
    const auto config = makeConfig({3.0}, 3);
    PipelineRuntime runtime(
        DistanceMatrix(points), config, PipelineMode::RegVRPH);
    WindowState window(points.size());
    auto edges = window.prepareWindow(
        runtime,
        boundsFor(runtime.config(), 0));
    DimensionFrame frame(runtime, window, std::move(edges));

    frame.matchPersistence(true);
    CHECK(frame.persistentPairs().size() == 1);
    if (frame.persistentPairs().size() == 1)
    {
        const auto& pair = frame.persistentPairs().front();
        CHECK(std::abs(pair.facet_weight - std::sqrt(1.3)) < 1e-12);
        CHECK(std::abs(pair.cofacet_weight - std::sqrt(1.45)) < 1e-12);
        CHECK(pair.facet_bindex >= 0);
        CHECK(pair.cofacet_bindex >= 0);
    }
    CHECK(frame.advance());
    frame.matchPersistence(true);
    CHECK(frame.persistentPairs().empty());
    auto batch = frame.takeDependencySupportBatch();
    CHECK(batch.supports.empty());
    CHECK(batch.protected_labels.empty());
    CHECK(!frame.advance());
}

void testMatchingDerivedPwphSupportAndPvGeometry()
{
    const std::vector<std::vector<double>> points{
        {0.0, 0.0},
        {1.0, 0.0},
        {0.5, 0.8660254037844386},
        {-0.5, 0.8660254037844386},
        {-1.0, 0.0},
        {-0.5, -0.8660254037844386},
        {0.5, -0.8660254037844386},
        {2.2, 0.2},
        {2.7, 0.9},
        {2.8, -0.7}};
    const auto config = makeConfig({1.1, 1.8, 3.5}, 2, 1, 2.0);
    PipelineRuntime runtime(
        DistanceMatrix(points), config, PipelineMode::PwPH);
    WindowState window(points.size());

    DependencySupportBatch first_batch;
    {
        auto edges = window.prepareWindow(
            runtime,
            boundsFor(runtime.config(), 0));
        DimensionFrame frame(runtime, window, std::move(edges));
        CHECK(frame.facetCount() == 14);
        CHECK(frame.cofacetCount() == 6);
        frame.matchPersistence(true);
        CHECK(frame.persistentPairs().size() == 2);
        if (frame.persistentPairs().size() == 2)
        {
            CHECK(frame.persistentPairs()[0].facet_bindex == 16);
            CHECK(frame.persistentPairs()[0].cofacet_bindex == 20);
            CHECK(frame.persistentPairs()[1].facet_bindex == 14);
            CHECK(frame.persistentPairs()[1].cofacet_bindex == 16);
        }
        first_batch = frame.takeDependencySupportBatch();
        // Each four-label support requires replay through multiple triangles.
        // Their order also controls which overlapping support PwPH consumes first.
        CHECK(first_batch.supports.size() == 2);
        if (first_batch.supports.size() == 2)
        {
            CHECK((first_batch.supports[0].labels ==
                   std::vector<std::size_t>{0, 1, 2, 6}));
            CHECK((first_batch.supports[1].labels ==
                   std::vector<std::size_t>{0, 3, 4, 5}));
        }
        CHECK(first_batch.protected_labels.empty());
        CHECK(!frame.advance());
    }

    DependencySupportPostProcessor processor;
    processCompletedPwphWindow(
        processor, runtime, window, std::move(first_batch));
    CHECK(window.pseudoVertices().size() == 1);
    if (window.pseudoVertices().size() == 1)
        CHECK(window.pseudoVertices().front().representatives.size() == 4);
    CHECK(window.activeLabels().size() == 7);

    runtime.ensureBinomialCapacity(window.totalLabelCount());
    auto pv_edges = window.prepareWindow(
        runtime,
        boundsFor(runtime.config(), 1));
    DimensionFrame pv_frame(runtime, window, std::move(pv_edges));
    CHECK(pv_frame.facetCount() == 10);
    CHECK(pv_frame.cofacetCount() == 5);
    pv_frame.matchPersistence(true);
    CHECK(pv_frame.persistentPairs().size() == 1);
    if (pv_frame.persistentPairs().size() == 1)
    {
        CHECK(pv_frame.persistentPairs()[0].facet_bindex == 51);
        CHECK(pv_frame.persistentPairs()[0].cofacet_bindex == 136);
    }
    const auto pv_batch = pv_frame.takeDependencySupportBatch();
    CHECK(pv_batch.supports.size() == 1);
    if (pv_batch.supports.size() == 1)
    {
        CHECK((pv_batch.supports[0].labels ==
               std::vector<std::size_t>{1, 2, 6, 10}));
    }
    CHECK(pv_batch.protected_labels.empty());
    CHECK(!pv_frame.advance());
}

void testEdgeOnlyDimensionFrame()
{
    const std::vector<std::vector<double>> points{
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0}};
    const auto config = makeConfig({2.0}, 1);
    PipelineRuntime runtime(
        DistanceMatrix(points), config, PipelineMode::RegVRPH);
    WindowState window(points.size());
    auto edges = window.prepareWindow(
        runtime,
        boundsFor(runtime.config(), 0));
    DimensionFrame frame(runtime, window, std::move(edges));

    CHECK(frame.dimension() == 1);
    CHECK(frame.facetDimension() == 1);
    CHECK(frame.isTopDimension());
    CHECK(frame.cofacetCount() == 0);

    frame.matchApparentPairsOnly();
    CHECK(frame.unmatchedFacetWeights().size() == 1);
    CHECK(frame.unmatchedTopCofacetWeights().empty());
    CHECK(!frame.advance());
}

void testSingleWindowPipelinesUseFinalEpsilon()
{
    const std::vector<std::vector<double>> point{{0.0}};
    const auto config = makeConfig({0.5, 1.0}, 1);

    PersistentHomologyPipeline direct(DistanceMatrix(point), config);
    direct.run();

    ApparentPairPipeline apparent(DistanceMatrix(point), config);
    apparent.run();
}

void testPwphProtectionAndDiameter()
{
    DependencySupportPostProcessor processor;

    {
        const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.2, 0.0}};
        const auto config = makeConfig({1.0});
        PipelineRuntime runtime(
            DistanceMatrix(points), config, PipelineMode::PwPH);
        WindowState window(points.size());
        static_cast<void>(window.prepareWindow(
            runtime,
            boundsFor(runtime.config(), 0)));

        DependencySupportBatch batch;
        batch.supports.push_back(support({0, 1}));
        batch.protected_labels.insert(1);
        processCompletedPwphWindow(
            processor, runtime, window, std::move(batch));

        CHECK(window.pseudoVertices().empty());
        CHECK((window.activeLabels() == std::vector<std::size_t>{0, 1}));
    }

    {
        const std::vector<std::vector<double>> points{{0.0, 0.0}, {2.0, 0.0}};
        const auto config = makeConfig({3.0}, 2, 1, 0.5);
        PipelineRuntime runtime(
            DistanceMatrix(points), config, PipelineMode::PwPH);
        WindowState window(points.size());
        static_cast<void>(window.prepareWindow(
            runtime,
            boundsFor(runtime.config(), 0)));

        DependencySupportBatch batch;
        batch.supports.push_back(support({0, 1}));
        processCompletedPwphWindow(
            processor, runtime, window, std::move(batch));

        CHECK(window.pseudoVertices().empty());
        CHECK((window.activeLabels() == std::vector<std::size_t>{0, 1}));
    }
}

void testPwphOverlapOrderAndSeparation()
{
    DependencySupportPostProcessor processor;

    {
        const std::vector<std::vector<double>> points{
            {0.0, 0.0},
            {0.1, 0.0},
            {0.2, 0.0}};
        const auto config = makeConfig({1.0});
        PipelineRuntime runtime(
            DistanceMatrix(points), config, PipelineMode::PwPH);
        WindowState window(points.size());
        static_cast<void>(window.prepareWindow(
            runtime,
            boundsFor(runtime.config(), 0)));

        DependencySupportBatch batch;
        batch.supports.push_back(support({0, 1}));
        batch.supports.push_back(support({1, 2}));
        processCompletedPwphWindow(
            processor, runtime, window, std::move(batch));

        CHECK(window.pseudoVertices().size() == 1);
        CHECK((window.pseudoVertices()[0].representatives ==
               std::vector<std::size_t>{1, 2}));
        CHECK((window.activeLabels() == std::vector<std::size_t>{0, 3}));
    }

    {
        const std::vector<std::vector<double>> points{
            {0.0, 0.0},
            {0.0, 0.2},
            {0.5, 0.0},
            {0.5, 0.2}};
        const auto config = makeConfig({2.0}, 2, 1, 1.0, 0.5);
        PipelineRuntime runtime(
            DistanceMatrix(points), config, PipelineMode::PwPH);
        WindowState window(points.size());
        static_cast<void>(window.prepareWindow(
            runtime,
            boundsFor(runtime.config(), 0)));

        DependencySupportBatch batch;
        batch.supports.push_back(support({0, 1}));
        batch.supports.push_back(support({2, 3}));
        processCompletedPwphWindow(
            processor, runtime, window, std::move(batch));

        CHECK(window.pseudoVertices().size() == 1);
        CHECK((window.pseudoVertices()[0].representatives ==
               std::vector<std::size_t>{2, 3}));
        CHECK((window.activeLabels() == std::vector<std::size_t>{0, 1, 4}));
    }
}

void testSequentialFilterDiagnostics()
{
    std::vector<std::vector<double>> points;
    for (std::size_t i = 0; i < 64; ++i)
        points.push_back({0.1 * i, 0.0});
    auto config = makeConfig({1.0}, 2, 1, 1.0, 0.5);
    config.verbose = true;
    PipelineRuntime runtime(DistanceMatrix(points), config, PipelineMode::PwPH);
    WindowState window(points.size());
    static_cast<void>(window.prepareWindow(runtime, boundsFor(config, 0)));
    DependencySupportBatch batch;
    // Reverse discovery traversal: select, reject separation, diameter,
    // cardinality, overlap, protection, and existing PV in that order.
    batch.supports = {support({64}), support({3, 4}), support({1, 2}),
                      support({}), support({4, 63}), support({2, 4}), support({0, 1})};
    batch.protected_labels.insert(3);
    std::ostringstream output;
    auto* previous = std::cout.rdbuf(output.rdbuf());
    try
    {
        processCompletedPwphWindow(DependencySupportPostProcessor{}, runtime,
                                   window, std::move(batch));
    }
    catch (...)
    {
        std::cout.rdbuf(previous);
        throw;
    }
    std::cout.rdbuf(previous);
    CHECK(output.str().find("raw=7 rejected_pv=1 rejected_protected=1 label_clean=5 "
        "rejected_overlap=1 rejected_cardinality=1 rejected_diameter=1 "
        "rejected_separation=1 eligible_at_visit=1 selected=1") != std::string::npos);
    CHECK(window.pseudoVertices().size() == 1);
    CHECK((window.pseudoVertices()[0].representatives == std::vector<std::size_t>{0, 1}));
    CHECK(window.activeLabels().size() == 63);
}

void testPwphFrozenPvsAndAppendOnlyLabels()
{
    const std::vector<std::vector<double>> points{
        {0.0, 0.0},
        {0.1, 0.0},
        {2.0, 0.0},
        {2.1, 0.0},
        {5.0, 0.0}};
    const auto config = makeConfig({1.0, 3.0});
    PipelineRuntime runtime(
        DistanceMatrix(points), config, PipelineMode::PwPH);
    WindowState window(points.size());
    DependencySupportPostProcessor processor;

    static_cast<void>(window.prepareWindow(
        runtime,
        boundsFor(runtime.config(), 0)));
    DependencySupportBatch first_batch;
    first_batch.supports.push_back(support({0, 1}));
    processCompletedPwphWindow(
        processor, runtime, window, std::move(first_batch));

    CHECK(window.pseudoVertices().size() == 1);
    CHECK(window.totalLabelCount() == 6);
    CHECK((window.pseudoVertices()[0].representatives ==
           std::vector<std::size_t>{0, 1}));
    CHECK((window.activeLabels() == std::vector<std::size_t>{2, 3, 4, 5}));

    runtime.ensureBinomialCapacity(window.totalLabelCount());
    static_cast<void>(window.prepareWindow(
        runtime,
        boundsFor(runtime.config(), 1)));
    CHECK((window.pseudoVertices()[0].representatives ==
           std::vector<std::size_t>{0, 1}));
    CHECK(window.activeLabelMask().size() == 6);
    CHECK(window.activeLabelMask()[0] == 0);
    CHECK(window.activeLabelMask()[1] == 0);
    CHECK(window.activeLabelMask()[5] == 1);
    CHECK(!window.pvLabelDistanceHash().empty());

    DependencySupportBatch second_batch;
    second_batch.supports.push_back(support({4, 5}));
    second_batch.supports.push_back(support({2, 3}));
    processCompletedPwphWindow(
        processor, runtime, window, std::move(second_batch));

    CHECK(window.pseudoVertices().size() == 2);
    CHECK(window.totalLabelCount() == 7);
    CHECK((window.pseudoVertices()[0].representatives ==
           std::vector<std::size_t>{0, 1}));
    CHECK((window.pseudoVertices()[1].representatives ==
           std::vector<std::size_t>{2, 3}));
    CHECK((window.activeLabels() == std::vector<std::size_t>{4, 5, 6}));
    CHECK(window.activeLabelMask().empty());
    CHECK(window.pvLabelDistanceHash().empty());
    CHECK(window.bounds().eps_lo == 0.0);
    CHECK(window.bounds().eps_hi == 0.0);
}

void testAdaptiveSchedule()
{
    const std::vector<double> grid{1.2,1.4,1.5,1.6,1.7,1.8,2.0};
    CHECK(chooseAdaptiveEpsilon(grid,1.2,0.4,0.35)==1.7);
    CHECK(chooseAdaptiveEpsilon(grid,1.2,0.0,0.35)==1.4);
    CHECK(chooseAdaptiveEpsilon(grid,1.5,0.0,0.35)==1.7);
    CHECK(chooseAdaptiveEpsilon(grid,1.8,0.0,0.35)==2.0);
    for (double reduction : {0.0,0.1,0.4,1.0})
    {
        double current=1.2;
        size_t count=0;
        while (current<2.0 && count++<grid.size())
        {
            const double next=chooseAdaptiveEpsilon(grid,current,reduction,0.35);
            CHECK(next>current);
            CHECK(std::find(grid.begin(),grid.end(),next)!=grid.end());
            current=next;
        }
        CHECK(current==2.0);
    }
    PipelineConfig invalid{.eps_breaks=grid,.maxdim=3,.adaptive=true,.adaptive_target=0.0};
    bool rejected=false;
    try { invalid.validate(); } catch (const std::invalid_argument&) { rejected=true; }
    CHECK(rejected);
}

void testGCostDecision()
{
    GCostController g;g.dimension=3;g.max_labels=200;g.max_cofacets=1e6;
    g.rebuild[0]=.01;g.residual[0]=.02;g.residual[7]=-.015;g.select[0]=.001;
    g.memory[0]=20000;g.memory[7]=-10000;
    const auto positive=g.decide(100,1e5,0,10,40,1.2,1.7,5000,1e6,.05);
    CHECK(positive.quotient);
    CHECK(std::abs(positive.g-.014)<1e-12);
    CHECK(positive.horizon==1.7);
    g.select[0]=.02;
    CHECK(!g.decide(100,1e5,0,10,40,1.2,1.7,5000,1e6,.05).quotient);
    const auto memory=g.decide(100,1e5,0,10,40,1.2,1.7,5000,20000,.05);
    CHECK(memory.quotient && memory.reason=="memory-constraint");
    g.scale=2.;g.update(.04,.04);CHECK(g.scale==2.);
    CHECK(!g.decide(100,1e5,0,0,0,1.2,1.7,5000,1e6,.05).quotient);
}

void testPwphAbsorbsPvsWhenEnabled()
{
    const std::vector<std::vector<double>> points{
        {0.0, 0.0},
        {0.1, 0.0},
        {2.0, 0.0},
        {2.1, 0.0},
        {5.0, 0.0}};
    auto config = makeConfig({1.0, 3.0});
    config.absorb_pv = true;
    PipelineRuntime runtime(
        DistanceMatrix(points), config, PipelineMode::PwPH);
    WindowState window(points.size());
    DependencySupportPostProcessor processor;

    static_cast<void>(window.prepareWindow(
        runtime,
        boundsFor(runtime.config(), 0)));
    DependencySupportBatch first_batch;
    first_batch.supports.push_back(support({0, 1}));
    processCompletedPwphWindow(
        processor, runtime, window, std::move(first_batch));
    CHECK((window.activeLabels() == std::vector<std::size_t>{2, 3, 4, 5}));

    runtime.ensureBinomialCapacity(window.totalLabelCount());
    static_cast<void>(window.prepareWindow(
        runtime,
        boundsFor(runtime.config(), 1)));
    // The support {2, PV 5} flattens to {0, 1, 2} (diameter 2.0 <= 3.0): PV 5 goes inactive,
    // its record stays, and the new PV gets the next label.
    DependencySupportBatch second_batch;
    second_batch.supports.push_back(support({2, 5}));
    processCompletedPwphWindow(
        processor, runtime, window, std::move(second_batch));
    CHECK(window.pseudoVertices().size() == 2);
    CHECK((window.pseudoVertices()[0].representatives ==
           std::vector<std::size_t>{0, 1}));
    CHECK((window.pseudoVertices()[1].representatives ==
           std::vector<std::size_t>{0, 1, 2}));
    CHECK((window.activeLabels() == std::vector<std::size_t>{3, 4, 6}));
}

void testFiberDominationKeepsLateVoids()
{
    // Random clouds near a sphere (voids) and in three clusters: removing dominated vertices must
    // keep, at every scale >= from, whether a void is alive and the first scale where one is.
    uint64_t state=20261103;
    const auto uniform=[&state]{state=state*6364136223846793005ULL+1442695040888963407ULL;
                                 return static_cast<double>(state>>11)*0x1.0p-53;};
    size_t voids=0,pruned=0;
    for(int trial=0;trial<300;++trial)
    {
        const size_t m=6+static_cast<size_t>(uniform()*18);
        std::vector<std::vector<double>> points;
        for(size_t i=0;i<m;++i)
        {
            if(trial%2==0)
            {
                const double z=2*uniform()-1,t=6.283185307179586*uniform(),r=std::sqrt(1-z*z)*(0.9+0.2*uniform());
                points.push_back({r*std::cos(t),r*std::sin(t),z});
            }
            else
            {
                const double c=static_cast<double>(i%3);
                points.push_back({c+0.6*uniform(),0.5*c*c+0.6*uniform(),0.6*uniform()});
            }
        }
        const DistanceMatrix d(points);
        std::vector<size_t> vs(m);for(size_t i=0;i<m;++i)vs[i]=i;
        const double from=0.2+1.2*uniform(),end=from+0.05+2*uniform();
        const auto kept=FiberRisk::undominated(d,vs,from,end);
        const double full=FiberRisk::firstVoid(d,vs,from,end),fast=FiberRisk::firstVoid(d,kept,from,end);
        CHECK((full<0)==(fast<0));
        if(full>=0)CHECK(std::max(full,from)==std::max(fast,from));
        voids+=full>=0;pruned+=kept.size()<m;
    }
    CHECK(voids>=20);
    CHECK(pruned>=100);
}

void testPwphShrinkRestoresSafePart()
{
    // Support {0,1,2} meets {3,4,5} and {6,7,8} each in two separate places, so the pair filter
    // drops it whole (the cheapest drop); --shrink-risky gives back {0,1}, which passes.
    const std::vector<std::vector<double>> points{
        {1.0, 1.39}, {1.38, 1.08}, {1.27, 0.8},
        {0.5, 1.18}, {0.05, 0.8}, {0.62, 0.68},
        {0.02, 0.65}, {0.32, 0.84}, {0.51, 0.39}};
    for (const bool shrink : {false, true})
    {
        auto config = makeConfig({0.8, 3.0});
        config.h1_risk = "selective-fast";
        config.shrink_risky = shrink;
        PipelineRuntime runtime(
            DistanceMatrix(points), config, PipelineMode::PwPH);
        WindowState window(points.size());
        DependencySupportPostProcessor processor;
        static_cast<void>(window.prepareWindow(
            runtime,
            boundsFor(runtime.config(), 0)));
        DependencySupportBatch batch;
        batch.supports.push_back(support({0, 1, 2}));
        batch.supports.push_back(support({3, 4, 5}));
        batch.supports.push_back(support({6, 7, 8}));
        processCompletedPwphWindow(
            processor, runtime, window, std::move(batch));
        std::vector<std::vector<std::size_t>> pvs;
        for (const auto& pv : window.pseudoVertices())
            pvs.push_back(pv.representatives);
        std::sort(pvs.begin(), pvs.end());
        const std::vector<std::vector<std::size_t>> expected = shrink
            ? std::vector<std::vector<std::size_t>>{{0, 1}, {3, 4, 5}, {6, 7, 8}}
            : std::vector<std::vector<std::size_t>>{{3, 4, 5}, {6, 7, 8}};
        if (pvs != expected)
            for (const auto& pv : pvs)
            {
                std::cout << "pv:";
                for (auto v : pv) std::cout << ' ' << v;
                std::cout << '\n';
            }
        CHECK(pvs == expected);
    }
}

void runTest(const char* name, void (*test)())
{
    const int failures_before = failure_count;
    try
    {
        test();
    }
    catch (const std::exception& error)
    {
        recordFailure(name, __FILE__, __LINE__, std::string("uncaught exception: ") + error.what());
    }
    catch (...)
    {
        recordFailure(name, __FILE__, __LINE__, "uncaught non-standard exception");
    }

    std::cout << (failure_count == failures_before ? "[PASS] " : "[FAIL] ")
              << name << '\n';
}
}

int main()
{
    runTest("PipelineConfig", testPipelineConfig);
    runTest("Adaptive finite-grid progress and validation", testAdaptiveSchedule);
    runTest("G common-horizon cost, margin, memory and update",testGCostDecision);
    runTest("DistanceMatrix validation", testDistanceMatrixValidation);
    runTest("WindowState preparation", testWindowStatePreparation);
    runTest("DimensionFrame advancement and determinism", testDimensionFrameAdvancementAndDeterminism);
    runTest("DimensionFrame wide cofacet indices", testDimensionFrameWideCofacetIndices);
    runTest("Persistence matching and support extraction", testPersistenceMatchingAndSupportExtraction);
    runTest("Matching-derived PwPH support and PV geometry", testMatchingDerivedPwphSupportAndPvGeometry);
    runTest("DimensionFrame edge-only mode", testEdgeOnlyDimensionFrame);
    runTest("Single-window pipelines use final epsilon", testSingleWindowPipelinesUseFinalEpsilon);
    runTest("PwPH protection and diameter", testPwphProtectionAndDiameter);
    runTest("PwPH overlap order and separation", testPwphOverlapOrderAndSeparation);
    runTest("Sequential filter diagnostics", testSequentialFilterDiagnostics);
    runTest("PwPH frozen PVs and append-only labels", testPwphFrozenPvsAndAppendOnlyLabels);
    runTest("PwPH absorbs PVs with --absorb-pv", testPwphAbsorbsPvsWhenEnabled);
    runTest("Fiber domination keeps late voids", testFiberDominationKeepsLateVoids);
    runTest("PwPH --shrink-risky restores the safe part", testPwphShrinkRestoresSafePart);

    if (failure_count != 0)
    {
        std::cerr << failure_count << " check(s) failed\n";
        return 1;
    }

    std::cout << "All core tests passed\n";
    return 0;
}
