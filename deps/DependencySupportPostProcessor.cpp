#include "DependencySupportPostProcessor.hpp"
#include "H1EdgeRisk.hpp"
#include "FiberRisk.hpp"
#include "PackDmtPayback.hpp"
#include "GeometryShadow.hpp"
#include "PolicyCostLog.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <utility>
#include <chrono>

namespace
{
// Labels that committing `selected` absorbs: their original vertices and every
// active PV whose vertices they contain (only with --absorb-pv).
std::unordered_set<std::size_t> absorbedLabels(const WindowState& window,const std::vector<SelectedPV>& selected)
{
    const std::size_t n=window.originalVertexCount();
    std::unordered_set<std::size_t> out;std::vector<char> taken(n,0);
    for(const auto& pv:selected)for(auto v:pv.representatives){out.insert(v);taken[v]=1;}
    for(std::size_t label:window.activeLabels())
        if(label>=n && taken[window.pseudoVertices()[label-n].representatives.front()])out.insert(label);
    return out;
}

// --shrink-risky. First drops that leave a clean set: the usual greedy drops (pair filter, then
// fiber filter), or with `joint` (--joint-drops) one greedy pass over the warnings of both
// filters, cheapest (size-1)/warnings first, then putting back any drop no warning needs. A
// dropped support's vertices become singletons (which only --fiber-full can flag) and the old PVs
// it absorbed become blocks again, so the set is inspected again until clean. Then each dropped
// support is offered back, in drop order, as large a part as stays clean. Every new warning
// involves the offered part, so only its fibers are inspected. Pair warning: the
// part keeps its crossing-graph component with the most of its vertices. Fiber warning: the
// part loses its core vertices. With `search` (--shrink-search), each step first tries
// removing a single unit, smallest first, and keeps the first part that passes. Absorbed PVs
// go whole; a part with fewer than two units is given up. A last full inspection confirms the
// result.
struct ShrinkStats
{
    size_t pair_warnings=0,fiber_warnings=0,late_drops=0,offered=0,offered_vertices=0,restored=0,restored_vertices=0,checks=0;
    double seconds=0;H1EdgeRisk::Result first_pair;std::vector<std::vector<size_t>> dropped;
};

ShrinkStats shrinkRisky(const DistanceMatrix& d,const WindowState& window,std::vector<SelectedPV>& selected,
                        double q,double endpoint,bool fiber,bool search,bool joint,bool full,bool h2)
{
    const auto start=std::chrono::steady_clock::now();ShrinkStats s;
    const size_t n=window.originalVertexCount();
    std::vector<size_t> unit(n);
    for(size_t v=0;v<n;++v)unit[v]=v;
    for(size_t label:window.activeLabels())if(label>=n){
        const auto& reps=window.pseudoVertices()[label-n].representatives;
        for(auto v:reps)unit[v]=reps.front();
    }
    std::vector<SelectedPV> dropped;
    const auto drop=[&](const std::vector<size_t>& ids){
        std::vector<char> gone(selected.size(),0);for(size_t id:ids)gone[id]=1;
        std::vector<SelectedPV> kept;
        for(size_t id=0;id<selected.size();++id)(gone[id]?dropped:kept).push_back(std::move(selected[id]));
        selected=std::move(kept);
    };
    // Proposed supports (indices into `selected`) of every warning of both filters.
    const auto warnings=[&](const H1EdgeRisk::Result& pr,const FiberRisk::Result& fr){
        std::vector<std::vector<size_t>> sets;
        for(const auto& w:pr.warnings){
            sets.emplace_back();
            for(size_t b:{w.lhs,w.rhs})if(pr.blocks[b].proposed)sets.back().push_back(pr.blocks[b].label-window.totalLabelCount());
        }
        for(const auto& w:fr.warnings){
            sets.emplace_back();
            for(size_t b:{w.a,w.b,w.c})if(b!=FiberRisk::NONE && fr.blocks[b].proposed)sets.back().push_back(fr.blocks[b].proposed_id);
        }
        return sets;
    };
    const auto jointDrops=[&](const std::vector<std::vector<size_t>>& sets){
        std::vector<size_t> w(selected.size());
        for(size_t p=0;p<selected.size();++p)w[p]=selected[p].representatives.size()-1;
        std::vector<char> gone(selected.size(),0),live(sets.size(),1);size_t left=sets.size();
        while(left){
            std::vector<size_t> inc(selected.size(),0);
            for(size_t k=0;k<sets.size();++k)if(live[k])for(size_t p:sets[k])++inc[p];
            size_t best=selected.size();
            for(size_t p=0;p<selected.size();++p)
                if(inc[p] && (best==selected.size() || w[p]*inc[best]<w[best]*inc[p]))best=p;
            if(best==selected.size())throw std::logic_error("Warning without a proposed support");
            gone[best]=1;
            for(size_t k=0;k<sets.size();++k)if(live[k] && std::count(sets[k].begin(),sets[k].end(),best)){live[k]=0;--left;}
        }
        std::vector<size_t> ids;for(size_t p=0;p<selected.size();++p)if(gone[p])ids.push_back(p);
        std::stable_sort(ids.begin(),ids.end(),[&](size_t x,size_t y){return w[x]>w[y];});
        for(size_t p:ids){
            gone[p]=0;
            for(const auto& set:sets)
                if(std::count(set.begin(),set.end(),p) && std::none_of(set.begin(),set.end(),[&](size_t x){return gone[x];})){gone[p]=1;break;}
        }
        ids.clear();for(size_t p=0;p<selected.size();++p)if(gone[p])ids.push_back(p);
        return ids;
    };
    const auto settle=[&]{
        for(;;){
            const auto pr=H1EdgeRisk::inspect(d,window,selected,q,endpoint);
            const auto fr=fiber?FiberRisk::inspect(d,window,selected,q,endpoint,FiberRisk::NONE,full,h2):FiberRisk::Result{};
            const auto sets=warnings(pr,fr);
            if(sets.empty())return;
            const auto ids=jointDrops(sets);s.late_drops+=ids.size();drop(ids);
        }
    };
    s.first_pair=H1EdgeRisk::inspect(d,window,selected,q,endpoint);
    s.pair_warnings=s.first_pair.warnings.size();
    if(joint){
        const auto fr=fiber?FiberRisk::inspect(d,window,selected,q,endpoint,FiberRisk::NONE,full,h2):FiberRisk::Result{};
        s.fiber_warnings=fr.warnings.size();
        drop(jointDrops(warnings(s.first_pair,fr)));
    }else{
        if(!s.first_pair.warnings.empty())drop(H1EdgeRisk::greedyDrops(s.first_pair,selected,window.totalLabelCount()));
        if(fiber && !selected.empty()){
            const auto f=FiberRisk::inspect(d,window,selected,q,endpoint,FiberRisk::NONE,full,h2);
            s.fiber_warnings=f.warnings.size();
            if(!f.warnings.empty())drop(FiberRisk::greedyDrops(f,selected.size()));
        }
    }
    settle();
    const auto units=[&](const std::vector<size_t>& vs){
        std::vector<size_t> u;for(auto v:vs)u.push_back(unit[v]);
        std::sort(u.begin(),u.end());return static_cast<size_t>(std::unique(u.begin(),u.end())-u.begin());
    };
    const auto without=[&](const std::vector<size_t>& part,const std::vector<size_t>& out){
        std::vector<size_t> gone;for(auto v:out)gone.push_back(unit[v]);
        std::sort(gone.begin(),gone.end());
        std::vector<size_t> next;
        for(auto v:part)if(!std::binary_search(gone.begin(),gone.end(),unit[v]))next.push_back(v);
        return next;
    };
    // The vertices the rule removes from `part`; empty when the part passes both filters.
    const auto risky=[&](const std::vector<size_t>& part){
        selected.push_back({part,0.0});
        const size_t focus=selected.size()-1;
        std::vector<size_t> out;
        ++s.checks;
        const auto pair=H1EdgeRisk::inspect(d,window,selected,q,endpoint,focus);
        if(!pair.warnings.empty()){
            const auto& w=pair.warnings.front();
            const bool lhs=pair.blocks[w.lhs].proposed && pair.blocks[w.lhs].label==window.totalLabelCount()+focus;
            const auto& comp=lhs?w.lhs_component:w.rhs_component;
            std::vector<size_t> count;
            for(size_t c:comp)if(c!=H1EdgeRisk::NONE){if(c>=count.size())count.resize(c+1,0);++count[c];}
            const size_t keep=std::max_element(count.begin(),count.end())-count.begin();
            for(size_t v=0;v<comp.size();++v)if(comp[v]!=H1EdgeRisk::NONE && comp[v]!=keep)out.push_back(part[v]);
        }else if(fiber){
            const auto f=FiberRisk::inspect(d,window,selected,q,endpoint,focus,full,h2);
            if(!f.warnings.empty()){
                for(auto v:f.warnings.front().core)if(std::binary_search(part.begin(),part.end(),v))out.push_back(v);
                if(out.empty())out=part;
            }
        }
        selected.pop_back();
        return out;
    };
    const auto offered=std::move(dropped);dropped.clear();
    for(const auto& pv:offered){
        ++s.offered;s.offered_vertices+=pv.representatives.size();
        std::vector<size_t> part=pv.representatives;bool restored=false;
        while(!restored && units(part)>=2){
            auto out=risky(part);
            if(!out.empty() && search){
                // Units of the part, smallest first (singletons before absorbed PVs), then by vertex.
                std::vector<size_t> ids;for(auto v:part)ids.push_back(unit[v]);
                std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
                std::vector<std::pair<size_t,size_t>> order;
                for(auto u:ids)order.push_back({static_cast<size_t>(std::count_if(part.begin(),part.end(),[&](size_t v){return unit[v]==u;})),u});
                std::sort(order.begin(),order.end());
                for(const auto& [size,u]:order){
                    auto trial=without(part,{u});
                    if(units(trial)>=2 && risky(trial).empty()){part=std::move(trial);out.clear();break;}
                }
            }
            if(out.empty()){
                double diameter=0;
                for(size_t i=0;i<part.size();++i)for(size_t j=i+1;j<part.size();++j)diameter=std::max(diameter,d.distance(part[i],part[j]));
                ++s.restored;s.restored_vertices+=part.size();selected.push_back({part,diameter});
                restored=true;continue;
            }
            part=without(part,out);
        }
        if(!restored)s.dropped.push_back(pv.representatives);
    }
    const size_t before=dropped.size();settle();
    for(size_t k=before;k<dropped.size();++k)s.dropped.push_back(dropped[k].representatives);
    s.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return s;
}
}

std::size_t DependencySupportPostProcessor::processPwph(
    const PipelineRuntime& runtime,
    WindowState& window,
    const double eps_hi,
    DependencySupportBatch&& batch) const
{
    using Clock = std::chrono::steady_clock;
    const auto selection_start = Clock::now();
    if(auto* cost=PolicyCost::active)cost->mark("SELECTOR_BEGIN");

    const std::size_t raw_support_count = batch.supports.size();
    const std::size_t original_vertex_count = window.originalVertexCount();

    const double diameter_limit =
        runtime.config().pv_cap_scale * eps_hi;
    const double min_separation = runtime.config().absoluteMinSeparation();
    const auto& distance_matrix = runtime.distanceMatrix();
    if (runtime.config().trace_pwph)
    {
        std::vector<size_t> protected_labels(batch.protected_labels.begin(),batch.protected_labels.end());
        std::sort(protected_labels.begin(),protected_labels.end());
        std::cout << "[support-trace] eps_hi=" << eps_hi << " protected=";
        for (auto label:protected_labels) std::cout << label << ',';
        std::cout << '\n';
        for (const auto& support:batch.supports)
        {
            std::cout << "[support-trace] eps_hi=" << eps_hi << " support=";
            for (auto label:support.labels) std::cout << label << ',';
            std::cout << '\n';
        }
    }

    std::vector<SelectedPV> selected_pvs;
    std::unordered_set<std::size_t> new_absorbed_labels;
    std::size_t label_clean_support_count = 0;
    std::size_t rejected_pv_count = 0;
    std::size_t rejected_protected_count = 0;
    std::size_t rejected_overlap_count = 0;
    std::size_t rejected_cardinality_count = 0;
    std::size_t rejected_diameter_count = 0;
    std::size_t rejected_separation_count = 0;

    const auto is_separated = [&](const std::vector<std::size_t>& vertices)
    {
        if (min_separation <= 0.0)
            return true;

        // PwPH never deactivates an existing PV. Iterating active PV labels is
        // equivalent today and is the correct lifetime boundary for a future
        // absorbing policy, where historical inactive records must be ignored.
        for (const std::size_t label : window.activeLabels())
        {
            if (label < original_vertex_count)
                continue;

            const auto& existing =
                window.pseudoVertices()[label - original_vertex_count].representatives;
            if (minCrossDistance(distance_matrix, vertices, existing) <
                min_separation)
            {
                return false;
            }
        }

        for (const auto& selected : selected_pvs)
        {
            if (minCrossDistance(
                    distance_matrix,
                    vertices,
                    selected.representatives) < min_separation)
            {
                return false;
            }
        }

        return true;
    };

    // Matching discovers supports while traversing facets from high to low
    // rank. PwPH applies its static filters while traversing that discovery
    // vector in reverse, preserving the existing overlap-selection order
    // without allocating a second vector of support containers.
    for (auto iter = batch.supports.rbegin();
         iter != batch.supports.rend();
         ++iter)
    {
        auto& labels = iter->labels;

        const bool contains_pv = std::any_of(
            labels.begin(),
            labels.end(),
            [original_vertex_count](const std::size_t label)
            {
                return label >= original_vertex_count;
            });

        // Preserve the frozen-PV policy and its ordering. This test precedes
        // protection, so the additional protected PV labels retained by the
        // generic frame cannot affect current PwPH selection. With --absorb-pv
        // a support may contain PV labels and is flattened below.
        if (contains_pv && !runtime.config().absorb_pv)
        {
            ++rejected_pv_count;
            continue;
        }

        const bool contains_protected_label = std::any_of(
            labels.begin(),
            labels.end(),
            [&](const std::size_t label)
            {
                return batch.protected_labels.contains(label);
            });
        if (contains_protected_label)
        {
            ++rejected_protected_count;
            continue;
        }

        ++label_clean_support_count;

        bool overlaps_new_absorbed_label = false;
        for (const std::size_t label : labels)
        {
            if (new_absorbed_labels.contains(label))
            {
                overlaps_new_absorbed_label = true;
                break;
            }
        }
        if (overlaps_new_absorbed_label)
        {
            ++rejected_overlap_count;
            continue;
        }

        // The flattened set of original vertices: the labels themselves unless
        // the support absorbs existing PVs. The caps apply to this set.
        std::vector<std::size_t> flattened;
        if (contains_pv)
        {
            for (const std::size_t label : labels)
            {
                if (label < original_vertex_count)
                    flattened.push_back(label);
                else
                {
                    const auto& representatives =
                        window.pseudoVertices()[label - original_vertex_count].representatives;
                    flattened.insert(flattened.end(), representatives.begin(), representatives.end());
                }
            }
            std::sort(flattened.begin(), flattened.end());
        }
        const auto& vertices = contains_pv ? flattened : labels;

        if (vertices.empty() || vertices.size() >= MAX_PV_CARDINALITY)
        {
            ++rejected_cardinality_count;
            continue;
        }

        const double diameter = maxPairwiseDistance(
            distance_matrix, vertices);
        if (diameter > diameter_limit)
        {
            ++rejected_diameter_count;
            continue;
        }
        if (!is_separated(vertices))
        {
            ++rejected_separation_count;
            continue;
        }

        new_absorbed_labels.insert(labels.begin(), labels.end());
        selected_pvs.push_back(
            SelectedPV{contains_pv ? std::move(flattened) : std::move(labels), diameter});
    }
    
    const auto selection_end = Clock::now();
    if(auto* cost=PolicyCost::active)cost->mark("FILTER_END");
    
    if (runtime.config().verbose && runtime.config().maxdim >= 2)
    {
        std::cout << "raw pv support cofacets size = " << raw_support_count << '\n';
        std::cout << "pv support label sets with no pv/protected contents size = "
                  << label_clean_support_count << '\n';
        std::cout << "supports rejected by overlap with newly selected PVs = "
                  << rejected_overlap_count << '\n';
        // Sequential first-rejection counts, in reverse discovery order.
        // Eligibility is evaluated against the PVs selected so far; no claim
        // is made about supports skipped before later predicates are reached.
        std::cout << "[pv-filters] eps_hi=" << eps_hi
                  << " raw=" << raw_support_count
                  << " rejected_pv=" << rejected_pv_count
                  << " rejected_protected=" << rejected_protected_count
                  << " label_clean=" << label_clean_support_count
                  << " rejected_overlap=" << rejected_overlap_count
                  << " rejected_cardinality=" << rejected_cardinality_count
                  << " rejected_diameter=" << rejected_diameter_count
                  << " rejected_separation=" << rejected_separation_count
                  << " eligible_at_visit=" << selected_pvs.size()
                  << " selected=" << selected_pvs.size() << '\n';
    }

    const auto& cfg=runtime.config();
    const auto& risk_mode=cfg.h1_risk;
    const bool selective=risk_mode=="selective" || risk_mode=="selective-fast" || risk_mode=="selective-graph";
    const bool fast=risk_mode=="selective-fast" || risk_mode=="selective-graph";
    const size_t eligible_count=selected_pvs.size();
    H1EdgeRisk::Result risk;
    size_t risk_scans=0,remaining_warnings=0;
    bool skipped_live=risk_mode=="selective-fast" && !batch.protected_labels.empty();
    std::vector<std::vector<size_t>> removed_supports;
    if(cfg.shrink_risky && batch.protected_labels.empty() && !selected_pvs.empty())
    {
        auto shrink=shrinkRisky(distance_matrix,window,selected_pvs,eps_hi,cfg.finalEpsilon(),cfg.fiber_risk,cfg.shrink_search,cfg.joint_drops,cfg.fiber_full,cfg.maxdim>=3);
        risk=std::move(shrink.first_pair);risk.seconds=shrink.seconds;
        risk_scans=shrink.checks+1;remaining_warnings=0;removed_supports=std::move(shrink.dropped);
        new_absorbed_labels=absorbedLabels(window,selected_pvs);
        std::cout<<std::setprecision(17)<<"[shrink] epsilon="<<eps_hi<<" eligible="<<eligible_count
            <<" pair_warnings="<<shrink.pair_warnings<<" fiber_warnings="<<shrink.fiber_warnings<<" late_drops="<<shrink.late_drops
            <<" offered="<<shrink.offered<<" offered_vertices="<<shrink.offered_vertices<<" restored="<<shrink.restored
            <<" restored_vertices="<<shrink.restored_vertices
            <<" checks="<<shrink.checks<<" kept="<<selected_pvs.size()<<" seconds="<<shrink.seconds<<'\n'<<std::setprecision(6);
    }
    else if(!risk_mode.empty())
    {
        if(!skipped_live){
            risk=H1EdgeRisk::inspect(distance_matrix,window,selected_pvs,eps_hi,cfg.finalEpsilon());
            risk_scans=1;remaining_warnings=risk.warnings.size();
        }
        if((risk_mode=="edge-fiber" && (!batch.protected_labels.empty() || !risk.warnings.empty())) ||
           (selective && !batch.protected_labels.empty()))
        {
            selected_pvs.clear();new_absorbed_labels.clear();
        }
        else if(selective && !risk.warnings.empty())
        {
            if(fast){
                const auto drops=H1EdgeRisk::greedyDrops(risk,selected_pvs,window.totalLabelCount());
                std::vector<bool> kept(selected_pvs.size(),true);
                for(size_t id:drops){kept[id]=false;removed_supports.push_back(selected_pvs[id].representatives);}
                std::vector<SelectedPV> remaining;
                for(size_t id=0;id<selected_pvs.size();++id)if(kept[id])remaining.push_back(std::move(selected_pvs[id]));
                selected_pvs=std::move(remaining);remaining_warnings=0;
            }else{
                auto current=risk;
                while(!current.warnings.empty() && !selected_pvs.empty())
                {
                    std::vector<size_t> incident(selected_pvs.size(),0);
                    for(const auto& warning:current.warnings)
                        for(size_t block:{warning.lhs,warning.rhs})
                            if(current.blocks[block].proposed)
                                ++incident.at(current.blocks[block].label-window.totalLabelCount());
                    size_t remove=selected_pvs.size();double best=std::numeric_limits<double>::infinity();
                    for(size_t p=0;p<selected_pvs.size();++p)if(incident[p])
                    {
                        const double cost=double(selected_pvs[p].representatives.size()-1)/incident[p];
                        if(cost<best){best=cost;remove=p;}
                    }
                    if(remove==selected_pvs.size())throw std::logic_error("Risk pair has no proposed support");
                    removed_supports.push_back(selected_pvs[remove].representatives);
                    selected_pvs.erase(selected_pvs.begin()+remove);
                    current=H1EdgeRisk::inspect(distance_matrix,window,selected_pvs,eps_hi,cfg.finalEpsilon());
                    ++risk_scans;risk.seconds+=current.seconds;
                }
                remaining_warnings=current.warnings.size();
            }
            new_absorbed_labels=absorbedLabels(window,selected_pvs);
        }
    }
    // Fiber check (three blocks; with --fiber-full also pairs and pairs with a vertex), after the pair filter above.
    if(cfg.fiber_risk && !cfg.shrink_risky && !selected_pvs.empty() && batch.protected_labels.empty())
    {
        const auto fiber=FiberRisk::inspect(distance_matrix,window,selected_pvs,eps_hi,cfg.finalEpsilon(),FiberRisk::NONE,cfg.fiber_full,cfg.maxdim>=3);
        const auto drops=FiberRisk::greedyDrops(fiber,selected_pvs.size());
        std::vector<bool> kept(selected_pvs.size(),true);for(size_t id:drops)kept[id]=false;
        std::vector<SelectedPV> remaining;
        for(size_t id=0;id<selected_pvs.size();++id)if(kept[id])remaining.push_back(std::move(selected_pvs[id]));
        selected_pvs=std::move(remaining);
        new_absorbed_labels=absorbedLabels(window,selected_pvs);
        std::cout<<std::setprecision(17)<<"[fiber-risk] epsilon="<<eps_hi<<" blocks="<<fiber.blocks.size()
            <<" triples="<<fiber.triples;
        if(cfg.fiber_full)std::cout<<" pairs="<<fiber.pairs<<" pair_vertex="<<fiber.pair_vertex;
        std::cout<<" warnings="<<fiber.warnings.size()<<" removed="<<drops.size()
            <<" kept="<<selected_pvs.size()<<" vertices="<<fiber.vertices<<" undominated="<<fiber.undominated
            <<" seconds="<<fiber.seconds;
        for(const auto& w:fiber.warnings)std::cout<<" void@"<<w.epsilon;
        std::cout<<'\n'<<std::setprecision(6);
    }
    const size_t retained_count=selected_pvs.size();
    uint64_t retained_hash=1469598103934665603ULL;
    std::vector<std::vector<size_t>> retained_members;
    for(const auto& pv:selected_pvs){
        retained_hash^=pv.representatives.size();retained_hash*=1099511628211ULL;
        for(auto v:pv.representatives){retained_hash^=v;retained_hash*=1099511628211ULL;}
        if(cfg.verbose)retained_members.push_back(pv.representatives);
    }
    if(auto* cost=PolicyCost::active)
        cost->retained(raw_support_count,batch.protected_labels.size(),eligible_count,remaining_warnings,
                       skipped_live,retained_hash,selected_pvs,window.activeLabels().size(),
                       window.pseudoVertices().size(),window.totalLabelCount());
    if(auto* shadow=GeoShadow::active;shadow && shadow->awaitingCommit())
        shadow->retained(eps_hi,raw_support_count,batch.protected_labels.size(),eligible_count,
                         remaining_warnings,skipped_live,retained_hash,selected_pvs);
    if(auto* cost=PolicyCost::active)cost->mark("VALUATION_BEGIN");
    const auto valuation_start=Clock::now();
    // Calibration valuations use the actual partition; they never modify it.
    const double Hbefore=PackDmtPayback::actualH(distance_matrix,window,{},cfg.finalEpsilon(),cfg.full_endpoint);
    const double Hafter=PackDmtPayback::actualH(distance_matrix,window,selected_pvs,cfg.finalEpsilon(),cfg.full_endpoint);
    const double deltaH=Hbefore-Hafter;
    const double C_hat=cfg.payback_bC*double(original_vertex_count)*original_vertex_count;
    const bool gated=cfg.pack_mode=="plan" || cfg.pack_mode=="replay";
    const bool payback=cfg.payback_aL*deltaH>1.25*C_hat;
    std::string reason=!batch.protected_labels.empty() && selective?"live-H1-history":eligible_count==0?"no-useful-support":retained_count==0?"selective-empty":"retained-supports";
    if(auto* cost=PolicyCost::active){cost->mark("VALUATION_END");cost->mark("DECISION_READY");}
    const bool force_defer=cfg.payback_action=="defer";
    if(retained_count && (force_defer || (gated && cfg.payback_action=="auto" && !payback))){
        selected_pvs.clear();new_absorbed_labels.clear();
        reason=force_defer?"development-forced-DEFER":"payback-not-confirmed";
    }
    const auto valuation_end=Clock::now();
    const size_t committed_count=selected_pvs.size();
    // Membership output is diagnostic I/O, outside the precommit and commit targets.
    const auto f_end=Clock::now();
    if(cfg.verbose)for(size_t id=0;id<retained_members.size();++id){
        std::cout<<std::setprecision(17)<<"[pv-retained] epsilon="<<eps_hi<<" order="<<id<<" members=";
        for(auto v:retained_members[id])std::cout<<v<<',';std::cout<<'\n'<<std::setprecision(6);
    }
    if(cfg.verbose)for(size_t id=0;id<selected_pvs.size();++id){
        std::cout<<std::setprecision(17)<<"[pv-selected] epsilon="<<eps_hi<<" order="<<id<<" members=";
        for(auto v:selected_pvs[id].representatives)std::cout<<v<<',';std::cout<<'\n'<<std::setprecision(6);
    }
    // EXTEND gain rule: the share of the next window's predicted top-dimension work this
    // commit removes. The pipeline keeps cutting only while it reaches the threshold.
    // The cost rule assumes the next commit removes the same share g of the final window's work, and
    // weighs that saving against the work of the extra window up to the next cut.
    if((cfg.extend_min_gain>0 || cfg.extend_cost_rule) && committed_count)
    {
        const auto next=std::upper_bound(cfg.eps_breaks.begin(),cfg.eps_breaks.end(),eps_hi);
        const double x=next==cfg.eps_breaks.end()?cfg.finalEpsilon():*next;
        const double before=PackDmtPayback::topWork(distance_matrix,window,{},x,cfg.maxdim);
        const double after=PackDmtPayback::topWork(distance_matrix,window,selected_pvs,x,cfg.maxdim);
        const double after_end=PackDmtPayback::topWork(distance_matrix,window,selected_pvs,cfg.finalEpsilon(),cfg.maxdim);
        last_gain_=before>0?1.0-after/before:0.0;
        last_gain_epsilon_=eps_hi;
        last_keep_cutting_=cfg.extend_cost_rule?last_gain_*after_end>after:last_gain_>=cfg.extend_min_gain;
        std::cout<<std::setprecision(17)<<"[extend-gain] epsilon="<<eps_hi<<" next="<<x<<" work_before="<<before
            <<" work_after="<<after<<" work_after_end="<<after_end<<" gain="<<last_gain_
            <<" predicted_saving="<<last_gain_*after_end<<" rule="<<(cfg.extend_cost_rule?"cost":"threshold")
            <<" threshold="<<cfg.extend_min_gain<<" action="<<(last_keep_cutting_?"KEEP-CUTTING":"CLOSE")<<'\n'<<std::setprecision(6);
    }
    if(auto* cost=PolicyCost::active)cost->mark("MUTATION_BEGIN");
    const auto commit_start=Clock::now();
    if(cfg.absorb_pv && committed_count){
        size_t absorbed_pvs=0;for(auto label:new_absorbed_labels)absorbed_pvs+=label>=original_vertex_count;
        size_t largest=0;for(const auto& pv:selected_pvs)largest=std::max(largest,pv.representatives.size());
        std::cout<<std::setprecision(17)<<"[absorb] epsilon="<<eps_hi<<" absorbed_pvs="<<absorbed_pvs
            <<" committed="<<committed_count<<" largest="<<largest<<'\n'<<std::setprecision(6);
    }
    if(committed_count)window.commitSelectedPVs(std::move(selected_pvs),new_absorbed_labels);
    const auto commit_end=Clock::now();
    if(auto* cost=PolicyCost::active)cost->mark("MUTATION_END");
    const double F=std::chrono::duration<double>(f_end-selection_start).count();
    const double C=std::chrono::duration<double>(commit_end-commit_start).count();
    if(!risk_mode.empty()){
        if(skipped_live)std::cout<<std::setprecision(17)<<"[h1-risk] epsilon="<<eps_hi<<" horizon="<<cfg.finalEpsilon()
            <<" mode="<<risk_mode<<" eligible="<<eligible_count<<" committed="<<committed_count
            <<" live_protected="<<batch.protected_labels.size()<<" risky_pairs=NA pairs_checked=NA edges_checked=NA event_batches=NA seconds=0 risk_scan=SKIPPED/live-H1 reason=live-H1-history\n"<<std::setprecision(6);
        else H1EdgeRisk::report(risk,eps_hi,cfg.finalEpsilon(),batch.protected_labels.size(),eligible_count,committed_count,risk_mode,cfg.verbose);
    }
    if(selective){
        std::cout<<std::setprecision(17)<<"[h1-risk-selection] epsilon="<<eps_hi
            <<" scans="<<risk_scans<<" removed="<<removed_supports.size()
            <<" initial_warnings="<<(skipped_live?"NA":std::to_string(risk.warnings.size()))
            <<" remaining_warnings="<<(skipped_live?"NA":std::to_string(remaining_warnings))
            <<" total_scan_s="<<risk.seconds<<" committed="<<committed_count<<'\n'<<std::setprecision(6);
        if(cfg.verbose)for(const auto& members:removed_supports){
            std::cout<<std::setprecision(17)<<"[h1-risk-removed] epsilon="<<eps_hi<<" members=";
            for(auto v:members)std::cout<<v<<',';std::cout<<'\n'<<std::setprecision(6);
        }
    }
    uint64_t post=1469598103934665603ULL;
    const auto add=[&](uint64_t value){post^=value;post*=1099511628211ULL;};
    for(auto label:window.activeLabels())add(label);
    for(const auto& pv:window.pseudoVertices()){add(pv.representatives.size());for(auto v:pv.representatives)add(v);}
    std::cout<<std::setprecision(17)<<"[postcommit] epsilon="<<eps_hi<<" eps_bits="<<std::bit_cast<uint64_t>(eps_hi)
        <<" eligible="<<eligible_count<<" retained="<<retained_count<<" committed="<<committed_count
        <<" retained_hash="<<retained_hash<<" Hbefore="<<Hbefore<<" Hafter="<<Hafter<<" deltaH="<<deltaH
        <<" marginal_M="<<cfg.payback_aM*deltaH-C_hat<<" C_hat="<<C_hat<<" payback="<<payback
        <<" F_s="<<F<<" C_s="<<C<<" C_observed="<<(committed_count>0)
        <<" valuation_s="<<std::chrono::duration<double>(valuation_end-valuation_start).count()
        <<" diagnostic="<<cfg.verbose<<" post_state_hash="<<post<<" active="<<window.activeLabels().size()
        <<" reason="<<reason<<'\n'<<std::setprecision(6);
    if (runtime.config().verbose && runtime.config().maxdim >= 2)
    {
        const double selection_ms =
            std::chrono::duration<double, std::milli>(
                selection_end - selection_start).count();
        const double commit_ms =
            std::chrono::duration<double, std::milli>(
                commit_end - commit_start).count();

        std::cout << "[pv-timing] eps_hi=" << eps_hi
                  << " selection_ms=" << selection_ms
                  << " commit_ms=" << commit_ms << '\n';
    }
    if(auto* cost=PolicyCost::active)cost->mark("DSP_RETURN");
    return committed_count;
}

double DependencySupportPostProcessor::maxPairwiseDistance(
    const DistanceMatrix& distance_matrix,
    const std::vector<std::size_t>& vertices)
{
    double maximum = 0.0;
    for (auto lhs = vertices.begin(); lhs != vertices.end(); ++lhs)
    {
        auto rhs = lhs;
        ++rhs;
        for (; rhs != vertices.end(); ++rhs)
            maximum = std::max(maximum, distance_matrix.distance(*lhs, *rhs));
    }
    return maximum;
}

double DependencySupportPostProcessor::minCrossDistance(
    const DistanceMatrix& distance_matrix,
    const std::vector<std::size_t>& lhs,
    const std::vector<std::size_t>& rhs)
{
    double minimum = std::numeric_limits<double>::infinity();
    for (const std::size_t lhs_vertex : lhs)
    {
        for (const std::size_t rhs_vertex : rhs)
        {
            minimum = std::min(
                minimum,
                distance_matrix.distance(lhs_vertex, rhs_vertex));
        }
    }
    return minimum;
}
