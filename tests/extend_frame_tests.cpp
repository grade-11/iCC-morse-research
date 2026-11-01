#include "DimensionFrame.hpp"
#include "PipelineCommon.hpp"
#include "WindowState.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
std::size_t comparisons = 0;

void require(const bool value, const std::string& detail)
{
    if (!value)
        throw std::runtime_error(detail);
}

std::uint64_t bits(const double value)
{
    return std::bit_cast<std::uint64_t>(value);
}

void compare(const DimensionFrame& extended, const DimensionFrame& direct,
             const std::string& label)
{
    require(extended.facetCount() == direct.facetCount(), label + ": facets");
    require(extended.cofacetCount() == direct.cofacetCount(), label + ": cofacets");
    const auto& lhs = extended.persistentPairs();
    const auto& rhs = direct.persistentPairs();
    require(lhs.size() == rhs.size(), label + ": pair count");
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        require(bits(lhs[i].facet_weight) == bits(rhs[i].facet_weight) &&
                bits(lhs[i].cofacet_weight) == bits(rhs[i].cofacet_weight) &&
                lhs[i].facet_bindex == rhs[i].facet_bindex &&
                lhs[i].cofacet_bindex == rhs[i].cofacet_bindex,
                label + ": ordered pair " + std::to_string(i));
    }
    const auto a = extended.snapshotDependencySupportBatch();
    const auto b = direct.snapshotDependencySupportBatch();
    require(a.protected_labels == b.protected_labels, label + ": protected labels");
    require(a.supports.size() == b.supports.size(), label + ": support count");
    for (std::size_t i = 0; i < a.supports.size(); ++i)
        require(a.supports[i].labels == b.supports[i].labels,
                label + ": ordered support " + std::to_string(i));
    ++comparisons;
}

void exercise(const std::vector<std::vector<double>>& points,
              const std::vector<double>& cuts, const std::string& name,
              const int threads = 1)
{
    PipelineConfig config;
    config.eps_breaks = cuts;
    config.maxdim = 2;
    config.threads = threads;
    PipelineRuntime runtime(DistanceMatrix(points), config, PipelineMode::PwPH);
    WindowState extended_window(points.size());
    auto edges = extended_window.prepareWindow(runtime, {0.0, cuts.front(), false});
    DimensionFrame extended(runtime, extended_window, std::move(edges), true);
    extended.matchPersistenceForExtension(true);
    for (std::size_t i = 0; i < cuts.size(); ++i)
    {
        if (i > 0)
        {
            extended_window.extendFirstWindowTo(cuts[i]);
            extended.extendTo(cuts[i]);
        }
        WindowState direct_window(points.size());
        auto direct_edges = direct_window.prepareWindow(runtime, {0.0, cuts[i], false});
        DimensionFrame direct(runtime, direct_window, std::move(direct_edges));
        direct.matchPersistence(true);
        compare(extended, direct, name + ":cut=" + std::to_string(i));
    }
}
}

int main()
{
    try
    {
        // Old H1 infinity/protection must be revoked after the diagonals arrive.
        exercise({{0,0}, {1,0}, {1,1}, {0,1}}, {.5,1.01,1.42,2.0}, "square-ties");
        // The direct canonical owner (2,1) is old; a different edge is maximal.
        exercise({{0,0}, {2,0}, {2,.1}}, {.2,1.9,2.01}, "old-owner-new-maximum");
        // Preserve direct zero-edge semantics in both canonical-owner cases.
        exercise({{0,0}, {1,0}, {1,0}}, {.1,.8,1.1}, "zero-owner");
        exercise({{0,0}, {0,0}, {1,0}}, {.1,.8,1.1}, "zero-other-edge");

        for (const std::size_t n : {6U, 9U, 14U, 24U})
        {
            for (std::uint64_t seed = 1; seed <= 12; ++seed)
            {
                std::mt19937_64 generator(seed * 7919 + n);
                std::uniform_real_distribution<double> coordinate(-1.0, 1.0);
                std::vector<std::vector<double>> points(n, std::vector<double>(3));
                for (auto& point : points)
                    for (auto& value : point)
                        value = coordinate(generator);
                DistanceMatrix distances(points);
                std::vector<double> values;
                for (std::size_t i = 0; i < n; ++i)
                    for (std::size_t j = i + 1; j < n; ++j)
                        values.push_back(distances.distance(i, j));
                std::sort(values.begin(), values.end());
                std::vector<double> cuts;
                for (const std::size_t numerator : {1U, 2U, 4U, 6U, 8U})
                    cuts.push_back(values[(values.size() - 1) * numerator / 10]);
                cuts.push_back(std::nextafter(values.back(), std::numeric_limits<double>::infinity()));
                cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
                exercise(points, cuts, "random-n=" + std::to_string(n) +
                         ":seed=" + std::to_string(seed), seed % 2 ? 1 : 4);
            }
        }
        std::cout << "PASS synthetic frame comparisons=" << comparisons << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
