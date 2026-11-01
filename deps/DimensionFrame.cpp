#include "DimensionFrame.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <chrono>
#include <omp.h>

#include "ImplicitMorseMatching.hpp"
#include "SimplexUtility.hpp"
#include "H1RealizationDiagnostics.hpp"

namespace
{
constexpr std::size_t MAX_FFI_PACKED_LABELS = sizeof(std::uint64_t);
}

DimensionFrame::FullCofacetIndex::FullCofacetIndex(const SimplexList& cofacets)
    : cofacets_(&cofacets),
      ranks_(0, IndirectHash{&cofacets}, IndirectEqual{&cofacets})
{
    appendRanks(0);
}

void DimensionFrame::FullCofacetIndex::appendRanks(const std::size_t first_rank)
{
    const auto& cofacets = *cofacets_;
    if (cofacets.empty())
        return;

    if (cofacets.size() - 1 > std::numeric_limits<RankIndex>::max())
    {
        throw std::overflow_error("Cofacet lookup ranks exceed the uint32_t range.");
    }

    ranks_.reserve(cofacets.size());
    for (std::size_t rank = first_rank; rank < cofacets.size(); ++rank)
        ranks_.emplace(static_cast<RankIndex>(rank));
}

bool DimensionFrame::FullCofacetIndex::tryFindRank(
    const SimplexBindex bindex,
    std::size_t& rank_out) const noexcept
{
    const auto it = ranks_.find(bindex);
    if (it == ranks_.end())
        return false;

    rank_out = *it;
    return true;
}

DimensionFrame::InterfaceWorkspace::InterfaceWorkspace(
    const SimplexList& cofacets,
    const std::size_t facet_count)
    : cofacet_lookup(cofacets),
      facet_mates(facet_count, UNMATCHED),
      cofacet_mates(cofacets.size(), UNMATCHED)
{
}

DimensionFrame::DimensionFrame(
    const PipelineRuntime& runtime,
    const WindowState& window,
    SimplexList prepared_edges,
    const bool retain_for_extension)
    : runtime_(runtime),
      window_(window),
      binomial_table_(runtime.binomialTable()),
      dimension_(runtime.config().maxdim == 1 ? 1 : 2),
      retain_for_extension_(retain_for_extension),
      extension_hi_(window.bounds().eps_hi),
      simplex_enumerator_(runtime.distanceMatrix(), binomial_table_),
      facet_list_(std::move(prepared_edges))
{
    if (retain_for_extension_ &&
        (runtime.mode() != PipelineMode::PwPH || runtime.config().maxdim != 2 ||
         window.bounds().eps_lo != 0.0 || !window.pseudoVertices().empty() ||
         window.activeLabels().size() != window.originalVertexCount()))
    {
        throw std::invalid_argument("EXTEND is restricted to the first PwPH H1 frame without PVs.");
    }
    buildInitialActiveFacetIndex();

    if (dimension_ >= 2)
    {
        using Clock = std::chrono::steady_clock;
        const auto cofacets_start = Clock::now();
        enumerateInitialCofacets();
        const auto cofacets_end = Clock::now();

        if (runtime_.mode() == PipelineMode::PwPH &&
            runtime_.config().verbose)
        {
            const double cofacets_ms =
                std::chrono::duration<double, std::milli>(
                    cofacets_end - cofacets_start).count();

            std::cout << "[rebuild-timing] stage=initial-cofacets"
                      << " eps_lo=" << window_.bounds().eps_lo
                      << " eps_hi=" << window_.bounds().eps_hi
                      << " dim=" << dimension_
                      << " cofacets_ms=" << cofacets_ms << '\n';
        }
    }
}

void DimensionFrame::buildInitialActiveFacetIndex()
{
    using Clock = std::chrono::steady_clock;
    const auto active_facets_start = Clock::now();

    if (!retain_for_extension_)
    {
        active_facet_index_ = SimplexUtility::getActiveEdgeIndexHashTable(
            binomial_table_, facet_list_, window_.totalLabelCount());
    }
    else
    {
        extension_kruskal_parent_.resize(window_.totalLabelCount());
        std::iota(extension_kruskal_parent_.begin(), extension_kruskal_parent_.end(), std::size_t{0});
        active_facet_index_.reserve(facet_list_.size());
        for (std::size_t i = 0; i < facet_list_.size(); ++i)
        {
            const auto labels = SimplexUtility::getSimplexVertices(
                binomial_table_, facet_list_[i].first, window_.totalLabelCount(), 1);
            if (SimplexUtility::mstFindRoot(extension_kruskal_parent_, labels[0]) !=
                SimplexUtility::mstFindRoot(extension_kruskal_parent_, labels[1]))
            {
                SimplexUtility::mstSetUnion(extension_kruskal_parent_, labels[0], labels[1]);
            }
            else
            {
                active_facet_index_.emplace(facet_list_[i].first, i);
            }
        }
    }

    const auto active_facets_end = Clock::now();

    if (runtime_.mode() == PipelineMode::PwPH &&
        runtime_.config().verbose)
    {
        const double active_facets_ms =
            std::chrono::duration<double, std::milli>(
                active_facets_end - active_facets_start).count();

        std::cout << "[rebuild-timing] stage=initial-active-facets"
                  << " eps_lo=" << window_.bounds().eps_lo
                  << " eps_hi=" << window_.bounds().eps_hi
                  << " dim=" << dimension_
                  << " active_facets_ms=" << active_facets_ms << '\n';
    }

    if (runtime_.config().verbose)
    {
        // Read the complement of the actual MSF clearing mask. Do not run a
        // second MST, and do not count historical inactive labels as components.
        // These are window-local graph diagnostics, not stitched PwPH H0 bars.
        const auto merges = facet_list_.size() - active_facet_index_.size();
        const auto components = window_.activeLabels().size() - merges;
        const auto old_precision = std::cout.precision();
        std::cout << std::setprecision(17)
                  << "[h0-msf] eps_lo=" << window_.bounds().eps_lo
                  << " eps_hi=" << window_.bounds().eps_hi
                  << " vertices=" << window_.activeLabels().size()
                  << " merges=" << merges << " components=" << components
                  << " cutoff_rule=strict_lt zero_edges="
                  << (runtime_.mode() == PipelineMode::PwPH ? "excluded" : "included")
                  << '\n';
        for (const auto& edge : facet_list_)
        {
            if (!active_facet_index_.contains(edge.first))
                std::cout << "[h0-msf-pair] birth=0 death=" << edge.second << '\n';
        }
        for (std::size_t i = 0; i < components; ++i)
            std::cout << "[h0-msf-pair] birth=0 death=-1\n";
        std::cout.precision(old_precision);
    }
}

void DimensionFrame::enumerateInitialCofacets()
{
    const auto& config = runtime_.config();
    const bool ffi_requested =
        config.verbose && !window_.pseudoVertices().empty() && config.maxdim >= 3;
    const bool collect_h1 = config.h1_realizations && config.maxdim==2 && !window_.pseudoVertices().empty();
    const bool collect_ffi = (ffi_requested && config.maxdim < MAX_FFI_PACKED_LABELS) || collect_h1;

    if (ffi_requested && !collect_ffi)
    {
        std::cout << "  [ffi] diagnostics skipped: packed witness realizations support max dimension "
                  << (MAX_FFI_PACKED_LABELS - 1) << '\n';
    }

    if (runtime_.mode() == PipelineMode::RegVRPH)
    {
        cofacet_list_ = simplex_enumerator_.getSortedVRCofacets(
            facet_list_, 1, window_.bounds().eps_hi, config.threads);
        return;
    }

    if (collect_ffi)
    {
        ffi_realizations_.emplace();
        cofacet_list_ = simplex_enumerator_.getGeometricCofacetListWithRealizations(
            facet_list_,
            window_.activeLabels(),
            window_.pseudoVertices(),
            window_.pvLabelDistanceHash(),
            1,
            window_.bounds().eps_hi,
            config.threads,
            ffi_realizations_->cofacet_realizations);
        if(collect_h1)
            H1RealizationDiagnostics::report(runtime_,window_,facet_list_,cofacet_list_,ffi_realizations_->cofacet_realizations);
        return;
    }

    cofacet_list_ = simplex_enumerator_.getGeometricCofacetList(
        facet_list_,
        window_.activeLabels(),
        window_.pseudoVertices(),
        window_.pvLabelDistanceHash(),
        1,
        window_.bounds().eps_hi,
        config.threads);
}

void DimensionFrame::traceComplex() const
{
    if (!runtime_.config().trace_pwph) return;
    const auto emit = [&](const SimplexList& simplices, size_t dim)
    {
        if (simplices.size()>10000) return;
        for (const auto& [index,weight] : simplices)
        {
            std::cout << "[simplex] eps_hi=" << window_.bounds().eps_hi
                      << " dimension=" << dim << " labels=";
            const auto labels=SimplexUtility::getSimplexVertices(binomial_table_,index,window_.totalLabelCount(),dim);
            for (auto label : labels) std::cout << label << ',';
            std::cout << " weight=" << std::setprecision(17) << weight << std::setprecision(6) << '\n';
        }
    };
    if (dimension_==2) emit(facet_list_,1);
    emit(cofacet_list_,dimension_);
}

void DimensionFrame::matchPersistence(const bool collect_dependency_support)
{
    using Clock = std::chrono::steady_clock;
    const bool measure_interface =
        runtime_.mode() == PipelineMode::PwPH &&
        runtime_.config().verbose;

    const auto log_interface_time =
        [&](const char* stage, Clock::time_point start, Clock::time_point end)
    {
        if (!measure_interface)
            return;

        const double elapsed_ms =
            std::chrono::duration<double, std::milli>(end - start).count();

        std::cout << "[interface-timing] stage=" << stage
                  << " eps_lo=" << window_.bounds().eps_lo
                  << " eps_hi=" << window_.bounds().eps_hi
                  << " dim=" << dimension_
                  << " elapsed_ms=" << elapsed_ms << '\n';
    };

    const auto timed_step = [&](const char* stage, auto&& operation)
    {
        if (!measure_interface)
        {
            operation();
            return;
        }

        const auto start = Clock::now();
        operation();
        const auto end = Clock::now();
        log_interface_time(stage, start, end);
    };

    timed_step("clear-results", [&]()
    {
        clearCurrentInterfaceResults();
    });

    // maxdim=1 is a valid edge-only run. H0 is resolved by the initial MSF and
    // the direct/PwPH output contracts intentionally print no H0 intervals.
    if (dimension_ == 1)
    {
        timed_step("active-index-release", [&]()
        {
            releaseActiveFacetIndex();
        });
        return;
    }

    if (isTopDimension())
    {
        timed_step("ffi-diagnostics-call", [&]()
        {
            reportFalseFacetIdentificationStats();
        });
    }

    timed_step("workspace-init", [&]()
    {
        initializeInterfaceWorkspace();
    });

    ImplicitMorseMatching matcher;

    const auto matching_start = Clock::now();
    matcher.matchPersistence(*this, collect_dependency_support);
    const auto matching_end = Clock::now();

    if (measure_interface)
    {
        const double matching_ms =
            std::chrono::duration<double, std::milli>(
                matching_end - matching_start).count();

        std::cout << "[matching-timing]"
                  << " eps_lo=" << window_.bounds().eps_lo
                  << " eps_hi=" << window_.bounds().eps_hi
                  << " dim=" << dimension_
                  << " matching_ms=" << matching_ms << '\n';
    }

    timed_step("active-index-release", [&]()
    {
        releaseActiveFacetIndex();
    });

    // Only a non-top interface needs to communicate unmatched cofacets to
    // advance(). A top-dimensional frame finishes without consuming this mask.
    if (!isTopDimension())
    {
        timed_step("next-active-mask", [&]()
        {
            finalizeNextActiveFacets();
        });
    }

    const auto workspace_release_start =
        measure_interface ? Clock::now() : Clock::time_point{};

    auto raw_support_info =
        releaseInterfaceWorkspace(collect_dependency_support);

    const auto workspace_release_end =
        measure_interface ? Clock::now() : Clock::time_point{};

    log_interface_time(
        "workspace-release", workspace_release_start, workspace_release_end);

    if (collect_dependency_support)
    {
        timed_step("support-materialize", [&]()
        {
            materializeMatchSupport(std::move(raw_support_info));
        });
    }
}

void DimensionFrame::matchPersistenceForExtension(const bool collect_dependency_support)
{
    if (!retain_for_extension_ || extension_matched_)
        throw std::logic_error("EXTEND initial matching requires a fresh retained frame.");
    if (!collect_dependency_support)
        throw std::invalid_argument("EXTEND must retain first-pass dependency supports.");

    clearCurrentInterfaceResults();
    initializeInterfaceWorkspace();
    extension_pairs_by_facet_.resize(facet_list_.size());
    extension_supports_by_facet_.resize(facet_list_.size());
    ImplicitMorseMatching matcher;
    matcher.matchPersistence(*this, true);
    mergeExtensionResults();
    extension_matched_ = true;
}

void DimensionFrame::extendTo(const double eps_hi)
{
    if (!retain_for_extension_ || !extension_matched_ || !interface_workspace_)
        throw std::logic_error("EXTEND requires a retained matched first-window frame.");
    if (!std::isfinite(eps_hi) || !(eps_hi > extension_hi_) ||
        window_.bounds().eps_hi != eps_hi || !window_.pseudoVertices().empty())
    {
        throw std::invalid_argument("Widen WindowState to the strictly larger EXTEND cut first.");
    }

    auto new_edges = simplex_enumerator_.getFirstWindowEdgeTail(
        extension_hi_, eps_hi, runtime_.config().threads);
    auto new_triangles = simplex_enumerator_.getFirstWindowTriangleTail(
        new_edges, extension_hi_, eps_hi, runtime_.config().threads);
    const std::size_t old_facet_count = facet_list_.size();
    const std::size_t old_cofacet_count = cofacet_list_.size();
    constexpr std::size_t max_mate_count = InterfaceWorkspace::UNMATCHED;
    if (new_edges.size() > max_mate_count - old_facet_count ||
        new_triangles.size() > max_mate_count - old_cofacet_count)
    {
        throw std::overflow_error("EXTEND exceeds the uint32_t mate capacity.");
    }

    // Old lists are strict filtration prefixes, so every old rank stays fixed.
    facet_list_.insert(facet_list_.end(), new_edges.begin(), new_edges.end());
    cofacet_list_.insert(cofacet_list_.end(), new_triangles.begin(), new_triangles.end());
    active_facet_index_.reserve(facet_list_.size());
    for (std::size_t i = old_facet_count; i < facet_list_.size(); ++i)
    {
        const auto labels = SimplexUtility::getSimplexVertices(
            binomial_table_, facet_list_[i].first, window_.totalLabelCount(), 1);
        if (SimplexUtility::mstFindRoot(extension_kruskal_parent_, labels[0]) !=
            SimplexUtility::mstFindRoot(extension_kruskal_parent_, labels[1]))
        {
            SimplexUtility::mstSetUnion(extension_kruskal_parent_, labels[0], labels[1]);
        }
        else
        {
            active_facet_index_.emplace(facet_list_[i].first, i);
        }
    }

    auto& workspace = *interface_workspace_;
    workspace.facet_mates.resize(facet_list_.size(), InterfaceWorkspace::UNMATCHED);
    workspace.cofacet_mates.resize(cofacet_list_.size(), InterfaceWorkspace::UNMATCHED);
    workspace.cofacet_lookup.appendRanks(old_cofacet_count);
    extension_pairs_by_facet_.resize(facet_list_.size());
    extension_supports_by_facet_.resize(facet_list_.size());
    for (auto& pair : extension_pairs_by_facet_)
    {
        // These classes may die in the tail; their old -1 output is provisional.
        if (pair && pair->cofacet_weight < 0.0)
            pair.reset();
    }
    persistent_pairs_.clear();
    ImplicitMorseMatching matcher;
    matcher.matchPersistence(*this, true, true);
    mergeExtensionResults();
    extension_hi_ = eps_hi;
}

void DimensionFrame::mergeExtensionResults()
{
    auto& workspace = *interface_workspace_;
    for (const auto& pair : persistent_pairs_)
    {
        const auto owner = active_facet_index_.find(pair.facet_bindex);
        if (owner == active_facet_index_.end())
            throw std::logic_error("EXTEND pair has no active owning facet.");
        extension_pairs_by_facet_[owner->second] = pair;
    }

    std::vector<std::size_t> support_owners;
    support_owners.reserve(workspace.raw_support_info.support_cofacet_indices.size());
    for (const auto& ranks : workspace.raw_support_info.support_cofacet_indices)
    {
        if (ranks.empty() || workspace.cofacet_mates[ranks.front()] == InterfaceWorkspace::UNMATCHED)
            throw std::logic_error("EXTEND support has no matched terminal cofacet.");
        support_owners.push_back(workspace.cofacet_mates[ranks.front()]);
    }

    // The current unmatched set can shrink: do not union protection across cuts.
    protected_labels_.clear();
    materializeMatchSupport(std::move(workspace.raw_support_info));
    if (dependency_support_batch_.supports.size() != support_owners.size())
        throw std::logic_error("EXTEND support ownership count differs from materialization.");
    for (std::size_t i = 0; i < support_owners.size(); ++i)
    {
        extension_supports_by_facet_[support_owners[i]] =
            std::move(dependency_support_batch_.supports[i]);
    }

    persistent_pairs_.clear();
    dependency_support_batch_.supports.clear();
    for (std::size_t i = facet_list_.size(); i > 0; --i)
    {
        const std::size_t owner = i - 1;
        if (extension_pairs_by_facet_[owner])
            persistent_pairs_.push_back(*extension_pairs_by_facet_[owner]);
        if (extension_supports_by_facet_[owner])
            dependency_support_batch_.supports.push_back(*extension_supports_by_facet_[owner]);
    }
}

DependencySupportBatch DimensionFrame::snapshotDependencySupportBatch() const
{
    DependencySupportBatch batch;
    batch.supports = dependency_support_batch_.supports;
    batch.protected_labels = protected_labels_;
    return batch;
}

void DimensionFrame::matchApparentPairsOnly()
{
    clearCurrentInterfaceResults();

    if (dimension_ == 1)
    {
        for (const auto& edge : facet_list_)
        {
            if (active_facet_index_.contains(edge.first))
                unmatched_facet_weights_.push_back(edge.second);
        }
        releaseActiveFacetIndex();
        return;
    }

    initializeInterfaceWorkspace();

    ImplicitMorseMatching matcher;
    matcher.matchApparentPairsOnly(*this);
    releaseActiveFacetIndex();

    if (!isTopDimension())
        finalizeNextActiveFacets();

    releaseInterfaceWorkspace(false);
}

void DimensionFrame::materializeMatchSupport(
    RawMatchSupportInfo raw_support_info)
{
    const std::size_t total_label_count = window_.totalLabelCount();
    const std::size_t facet_dimension = dimension_ - 1;

    for (const std::size_t facet_index :
         raw_support_info.protected_facet_indices)
    {
        const auto labels = SimplexUtility::getSimplexVertices(
            binomial_table_, facet_list_[facet_index].first,
            total_label_count, facet_dimension);
        protected_labels_.insert(labels.begin(), labels.end());
    }

    if (!isTopDimension())
        return;

    std::vector<DependencySupport> materialized_supports;
    materialized_supports.reserve(
        raw_support_info.support_cofacet_indices.size());
    std::vector<std::uint8_t> label_seen(total_label_count, 0);

    for (const auto& support_ranks :
         raw_support_info.support_cofacet_indices)
    {
        DependencySupport support;
        for (const std::size_t cofacet_index : support_ranks)
        {
            const auto labels = SimplexUtility::getSimplexVertices(
                binomial_table_, cofacet_list_[cofacet_index].first,
                total_label_count, dimension_);
            for (const std::size_t label : labels)
            {
                if (label_seen[label] != 0)
                    continue;

                label_seen[label] = 1;
                support.labels.push_back(label);
            }
        }

        for (const std::size_t label : support.labels)
            label_seen[label] = 0;

        std::sort(support.labels.begin(), support.labels.end());
        materialized_supports.push_back(std::move(support));
    }

    // Replace rather than append so retrying after a later allocation failure
    // cannot duplicate supports already materialized by the first attempt.
    dependency_support_batch_.supports = std::move(materialized_supports);
}

void DimensionFrame::finalizeNextActiveFacets()
{
    assert(interface_workspace_.has_value());
    const auto& cofacet_mates = interface_workspace_->cofacet_mates;

    next_active_facet_mask_.assign(cofacet_list_.size(), 0);
    next_active_facet_count_ = 0;
    for (std::size_t i = 0; i < cofacet_list_.size(); ++i)
    {
        if (cofacet_mates[i] != InterfaceWorkspace::UNMATCHED)
            continue;

        next_active_facet_mask_[i] = 1;
        ++next_active_facet_count_;
    }
}

void DimensionFrame::releaseActiveFacetIndex()
{
    decltype(active_facet_index_){}.swap(active_facet_index_);
}

void DimensionFrame::initializeInterfaceWorkspace()
{
    // Reserve the maximum value for UNMATCHED before allocating either side.
    constexpr std::size_t max_mate_count = InterfaceWorkspace::UNMATCHED;
    if (facet_list_.size() > max_mate_count ||
        cofacet_list_.size() > max_mate_count)
    {
        throw std::overflow_error("Facet or cofacet count exceeds the uint32_t mate capacity.");
    }

    interface_workspace_.emplace(cofacet_list_, facet_list_.size());
}

DimensionFrame::RawMatchSupportInfo
DimensionFrame::releaseInterfaceWorkspace(
    const bool preserve_raw_support_info)
{
    assert(interface_workspace_.has_value());

    RawMatchSupportInfo raw_support_info;
    if (preserve_raw_support_info)
    {
        raw_support_info =
            std::move(interface_workspace_->raw_support_info);
    }

    interface_workspace_.reset();
    return raw_support_info;
}

bool DimensionFrame::tryFindCofacetRank(
    const SimplexBindex bindex,
    std::size_t& rank_out) const noexcept
{
    // This is a complete membership contract: false means the candidate
    // cofacet does not exist in the current complex, not merely that it lacks
    // a stored pivot. A pivot-only index cannot implement this interface.
    return interface_workspace_.has_value() &&
           interface_workspace_->cofacet_lookup.tryFindRank(bindex, rank_out);
}

SimplexList DimensionFrame::enumerateNextCofacets(
    robin_hood::unordered_map<SimplexBindex, std::uint64_t>* const realization_output)
{
    const auto& config = runtime_.config();
    if (runtime_.mode() == PipelineMode::RegVRPH)
    {
        return simplex_enumerator_.getSortedVRCofacets(
            cofacet_list_, dimension_, window_.bounds().eps_hi, config.threads);
    }

    if (realization_output != nullptr)
    {
        return simplex_enumerator_.getGeometricCofacetListWithRealizations(
            cofacet_list_,
            window_.activeLabels(),
            window_.pseudoVertices(),
            window_.pvLabelDistanceHash(),
            dimension_,
            window_.bounds().eps_hi,
            config.threads,
            *realization_output);
    }

    return simplex_enumerator_.getGeometricCofacetList(
        cofacet_list_,
        window_.activeLabels(),
        window_.pseudoVertices(),
        window_.pvLabelDistanceHash(),
        dimension_,
        window_.bounds().eps_hi,
        config.threads);
}

bool DimensionFrame::advance()
{
    if (isTopDimension())
        return false;

    // The old active-facet index was released immediately after matching.
    // Release the previous facet list before enumeration as well. Building the
    // next active-facet index afterward is an intentional peak-RSS improvement
    // over the old Q&E sequence.
    SimplexList{}.swap(facet_list_);
    std::optional<robin_hood::unordered_map<SimplexBindex, std::uint64_t>>
        next_cofacet_realizations;
    if (ffi_realizations_.has_value())
    {
        decltype(ffi_realizations_->facet_realizations){}.swap(
            ffi_realizations_->facet_realizations);
        next_cofacet_realizations.emplace();
    }

    using Clock = std::chrono::steady_clock;
    const auto cofacets_start = Clock::now();
    auto next_cofacets = enumerateNextCofacets(
        next_cofacet_realizations.has_value() ? &*next_cofacet_realizations : nullptr);
    const auto cofacets_end = Clock::now();

    if (runtime_.mode() == PipelineMode::PwPH &&
        runtime_.config().verbose)
    {
        const double cofacets_ms =
            std::chrono::duration<double, std::milli>(
                cofacets_end - cofacets_start).count();

        std::cout << "[rebuild-timing] stage=next-cofacets"
                  << " eps_lo=" << window_.bounds().eps_lo
                  << " eps_hi=" << window_.bounds().eps_hi
                  << " dim=" << (dimension_ + 1)
                  << " cofacets_ms=" << cofacets_ms << '\n';
    }

    const auto active_facets_start = Clock::now();

    robin_hood::unordered_map<SimplexBindex, std::size_t> next_active_facets;
    next_active_facets.reserve(next_active_facet_count_);
    for (std::size_t i = 0; i < cofacet_list_.size(); ++i)
    {
        if (next_active_facet_mask_[i] != 0)
            next_active_facets.emplace(cofacet_list_[i].first, i);
    }

    const auto active_facets_end = Clock::now();

    if (runtime_.mode() == PipelineMode::PwPH &&
        runtime_.config().verbose)
    {
        const double active_facets_ms =
            std::chrono::duration<double, std::milli>(
                active_facets_end - active_facets_start).count();

        std::cout << "[rebuild-timing] stage=next-active-facets"
                  << " eps_lo=" << window_.bounds().eps_lo
                  << " eps_hi=" << window_.bounds().eps_hi
                  << " dim=" << (dimension_ + 1)
                  << " active_facets_ms=" << active_facets_ms << '\n';
    }

    facet_list_ = std::move(cofacet_list_);
    cofacet_list_ = std::move(next_cofacets);
    active_facet_index_ = std::move(next_active_facets);

    if (ffi_realizations_.has_value())
    {
        ffi_realizations_->facet_realizations =
            std::move(ffi_realizations_->cofacet_realizations);
        ffi_realizations_->cofacet_realizations =
            std::move(*next_cofacet_realizations);
    }

    ++dimension_;
    clearCurrentInterfaceResults();
    return true;
}

void DimensionFrame::clearCurrentInterfaceResults()
{
    decltype(persistent_pairs_){}.swap(persistent_pairs_);
    decltype(unmatched_facet_weights_){}.swap(unmatched_facet_weights_);
    decltype(unmatched_top_cofacet_weights_){}.swap(
        unmatched_top_cofacet_weights_);
    decltype(next_active_facet_mask_){}.swap(next_active_facet_mask_);
    next_active_facet_count_ = 0;
}

DependencySupportBatch DimensionFrame::takeDependencySupportBatch()
{
    dependency_support_batch_.protected_labels = std::move(protected_labels_);
    return std::move(dependency_support_batch_);
}

void DimensionFrame::reportFalseFacetIdentificationStats() const
{
    if (!ffi_realizations_.has_value())
        return;

    const std::size_t original_vertex_count = window_.originalVertexCount();
    const std::size_t pv_count = window_.pseudoVertices().size();
    if (pv_count == 0 || dimension_ < 3 || dimension_ >= MAX_FFI_PACKED_LABELS)
        return;

    const auto& pseudo_vertices = window_.pseudoVertices();
    const std::size_t total_label_count = window_.totalLabelCount();
    const std::size_t cofacet_label_count = dimension_ + 1;
    const std::size_t facet_label_count = dimension_;
    const auto facet_hash = SimplexUtility::getSimplexIndexHashTable(facet_list_);

    struct FfiStats
    {
        std::size_t cofacets_with_pv = 0;
        std::size_t incidence_total = 0;
        std::size_t facet_all_original = 0;
        std::size_t missing_facet_realization = 0;
        std::size_t diff_multi_safe = 0;
        std::size_t diff_multi_flagged = 0;
        std::size_t gap_ratio_zero = 0;

        std::vector<float> gap_ratio_samples;

        double worst_overshoot = 0.0;
        double worst_overshoot_cofacet_weight = 0.0;
        double worst_overshoot_facet_weight = 0.0;

        void accumulate(FfiStats& other)
        {
            cofacets_with_pv += other.cofacets_with_pv;
            incidence_total += other.incidence_total;
            facet_all_original += other.facet_all_original;
            missing_facet_realization += other.missing_facet_realization;
            diff_multi_safe += other.diff_multi_safe;
            diff_multi_flagged += other.diff_multi_flagged;
            gap_ratio_zero += other.gap_ratio_zero;
            gap_ratio_samples.insert(
                gap_ratio_samples.end(),
                other.gap_ratio_samples.begin(),
                other.gap_ratio_samples.end());
            std::vector<float>().swap(other.gap_ratio_samples);

            if (other.worst_overshoot > worst_overshoot)
            {
                worst_overshoot = other.worst_overshoot;
                worst_overshoot_cofacet_weight =
                    other.worst_overshoot_cofacet_weight;
                worst_overshoot_facet_weight =
                    other.worst_overshoot_facet_weight;
            }
        }
    };

    const int worker_count = runtime_.config().threads;
    std::vector<FfiStats> thread_stats(static_cast<std::size_t>(worker_count));

    struct FfiWorkspace
    {
        std::vector<std::size_t> cofacet_labels;
        std::vector<std::size_t> cofacet_witnesses;
        std::vector<std::size_t> facet_witnesses;
        std::vector<std::size_t> restricted_witnesses;
        std::vector<std::size_t> differing_positions;
    };
    std::vector<FfiWorkspace> thread_workspaces(
        static_cast<std::size_t>(worker_count));

    omp_set_num_threads(worker_count);

#pragma omp parallel for schedule(dynamic, 64) num_threads(worker_count)
    for (std::size_t ci = 0; ci < cofacet_list_.size(); ++ci)
    {
        const int thread_id = omp_get_thread_num();
        auto& stats = thread_stats[static_cast<std::size_t>(thread_id)];
        auto& workspace = thread_workspaces[static_cast<std::size_t>(thread_id)];

        const SimplexBindex cofacet_bindex = cofacet_list_[ci].first;
        const double cofacet_weight = cofacet_list_[ci].second;

        const auto realization_iter =
            ffi_realizations_->cofacet_realizations.find(cofacet_bindex);
        if (realization_iter == ffi_realizations_->cofacet_realizations.end())
            continue;

        ++stats.cofacets_with_pv;
        const std::uint64_t packed_cofacet_realization = realization_iter->second;

        SimplexUtility::getSimplexVerticesInPlace(
            binomial_table_, workspace.cofacet_labels,
            cofacet_bindex, total_label_count, dimension_);

        workspace.cofacet_witnesses.resize(cofacet_label_count);
        for (std::size_t group = 0; group < cofacet_label_count; ++group)
        {
            const std::size_t label = workspace.cofacet_labels[group];
            const std::size_t local_index =
                (packed_cofacet_realization >> (8 * group)) & 0xFFULL;
            workspace.cofacet_witnesses[group] =
                label < original_vertex_count
                    ? label
                    : pseudo_vertices[label - original_vertex_count]
                          .representatives[local_index];
        }

        SimplexBindex above = 0;
        SimplexBindex below = cofacet_bindex;
        std::size_t k = dimension_;

        for (std::size_t dropped_position = 0;
             dropped_position < cofacet_label_count;
             ++dropped_position)
        {
            const std::size_t dropped_label =
                workspace.cofacet_labels[dropped_position];
            below -= binomial_table_[dropped_label][k + 1];
            const SimplexBindex facet_bindex = above + below;
            above += binomial_table_[dropped_label][k];
            if (dropped_position + 1 < cofacet_label_count)
                --k;

            ++stats.incidence_total;

            std::size_t facet_pv_label_count = 0;
            for (std::size_t q = 0; q < facet_label_count; ++q)
            {
                const std::size_t label = workspace.cofacet_labels[
                    q < dropped_position ? q : q + 1];
                if (label < original_vertex_count)
                    break;
                ++facet_pv_label_count;
            }

            if (facet_pv_label_count == 0)
            {
                ++stats.facet_all_original;
                continue;
            }

            const auto facet_iter = facet_hash.find(facet_bindex);
            const auto facet_realization_iter =
                ffi_realizations_->facet_realizations.find(facet_bindex);
            if (facet_iter == facet_hash.end() ||
                facet_realization_iter ==
                    ffi_realizations_->facet_realizations.end())
            {
                ++stats.missing_facet_realization;
                continue;
            }

            const double facet_weight = facet_list_[facet_iter->second].second;
            const std::uint64_t packed_facet_realization =
                facet_realization_iter->second;

            workspace.facet_witnesses.resize(facet_label_count);
            workspace.restricted_witnesses.resize(facet_label_count);
            workspace.differing_positions.clear();

            for (std::size_t q = 0; q < facet_label_count; ++q)
            {
                const std::size_t cofacet_position =
                    q < dropped_position ? q : q + 1;
                const std::size_t label =
                    workspace.cofacet_labels[cofacet_position];
                const std::size_t facet_local_index =
                    (packed_facet_realization >> (8 * q)) & 0xFFULL;

                workspace.facet_witnesses[q] =
                    label < original_vertex_count
                        ? label
                        : pseudo_vertices[label - original_vertex_count]
                              .representatives[facet_local_index];
                workspace.restricted_witnesses[q] =
                    workspace.cofacet_witnesses[cofacet_position];

                if (workspace.facet_witnesses[q] !=
                    workspace.restricted_witnesses[q])
                {
                    workspace.differing_positions.push_back(q);
                }
            }

            double restricted_weight = 0.0;
            for (std::size_t a = 0; a + 1 < facet_label_count; ++a)
            {
                for (std::size_t b = a + 1; b < facet_label_count; ++b)
                {
                    restricted_weight = std::max(
                        restricted_weight,
                        runtime_.distanceMatrix().distance(
                            workspace.restricted_witnesses[a],
                            workspace.restricted_witnesses[b]));
                }
            }

            const double realization_gap = restricted_weight - facet_weight;
            const double realization_gap_ratio = realization_gap / facet_weight;
            stats.gap_ratio_samples.push_back(
                static_cast<float>(realization_gap_ratio));

            if (realization_gap_ratio == 0.0)
                ++stats.gap_ratio_zero;

            if (workspace.differing_positions.size() > 1)
            {
                double max_cross_distance = 0.0;
                for (const std::size_t q : workspace.differing_positions)
                {
                    for (std::size_t r = 0; r < facet_label_count; ++r)
                    {
                        max_cross_distance = std::max(
                            max_cross_distance,
                            runtime_.distanceMatrix().distance(
                                workspace.restricted_witnesses[q],
                                workspace.facet_witnesses[r]));
                    }
                }

                if (max_cross_distance <= cofacet_weight)
                {
                    ++stats.diff_multi_safe;
                }
                else
                {
                    ++stats.diff_multi_flagged;
                    const double overshoot =
                        max_cross_distance - cofacet_weight;
                    if (overshoot > stats.worst_overshoot)
                    {
                        stats.worst_overshoot = overshoot;
                        stats.worst_overshoot_cofacet_weight = cofacet_weight;
                        stats.worst_overshoot_facet_weight = facet_weight;
                    }
                }
            }
        }
    }

    FfiStats total;
    for (auto& stats : thread_stats)
        total.accumulate(stats);

    std::sort(total.gap_ratio_samples.begin(), total.gap_ratio_samples.end());
    const auto quantile = [&](const double q) -> double
    {
        if (total.gap_ratio_samples.empty())
            return 0.0;
        const std::size_t position = std::min(
            total.gap_ratio_samples.size() - 1,
            static_cast<std::size_t>(
                q * (total.gap_ratio_samples.size() - 1) + 0.5));
        return total.gap_ratio_samples[position];
    };

    const std::size_t pv_incident =
        total.incidence_total - total.facet_all_original;
    const std::size_t measured_pv_incident = total.gap_ratio_samples.size();
    const std::size_t multi_incidence =
        total.diff_multi_safe + total.diff_multi_flagged;
    const auto percentage = [&](const std::size_t count) -> double
    {
        if (measured_pv_incident == 0)
            return 0.0;
        return 100.0 * static_cast<double>(count) /
               static_cast<double>(measured_pv_incident);
    };

    if (pv_incident == 0)
        return;

    if (measured_pv_incident > 0)
    {
        std::cout << std::endl << "FFI stats \n";
        std::cout << "  [ffi2] realization gap ratio (w(Y|F) - w(F)) / w(F): zero = "
                  << percentage(total.gap_ratio_zero) << "%"
                  << "  p50 = " << quantile(0.50)
                  << "  p90 = " << quantile(0.90)
                  << "  p99 = " << quantile(0.99) << '\n';

        std::cout << "  [ffi2] differing PV coordinates:"
                  << "  multi indcidence = " << percentage(multi_incidence) << "%"
                  << "  clique-certified = " << percentage(total.diff_multi_safe) << "%"
                  << "  void capable = "
                  << percentage(total.diff_multi_flagged) << "%\n";
    }

    if (total.diff_multi_flagged > 0)
    {
        std::cout << "  [ffi2] worst no-clique-certificate incidence: carrier overshoot = "
                  << total.worst_overshoot
                  << "  w(facet) = " << total.worst_overshoot_facet_weight
                  << "  w(cofacet) = " << total.worst_overshoot_cofacet_weight
                  << '\n';
    }

    if (total.missing_facet_realization > 0)
    {
        std::cout << "  [ffi2] WARNING: missing facet realizations = "
                  << total.missing_facet_realization << '\n';
    }

    std::cout << "  [ffi2-line]"
              << " eps_lo=" << window_.bounds().eps_lo
              << " eps_hi=" << window_.bounds().eps_hi
              << " dim=" << dimension_
              << " pv_cofacets =" << total.cofacets_with_pv
              << " diffm_safe_pct=" << percentage(total.diff_multi_safe) << "%"
              << " diffm_flagged_pct=" << percentage(total.diff_multi_flagged) << "%"
              << "\nEnd of FFI stats"
              << std::endl;
}
