#include "Pipelines.hpp"
#include "DimensionFrame.hpp"
#include "PackDmtPayback.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace
{
// The frozen postprocessor is intentionally reused without algorithm changes.
// A copied window is invalidated before selection, so even a positive preview
// cannot invalidate the geometry borrowed by the live incremental frame.
class CaptureStdout
{
public:
    explicit CaptureStdout(std::ostream& destination)
        : flags_(std::cout.flags()), precision_(std::cout.precision()),
          fill_(std::cout.fill()), original_(std::cout.rdbuf(destination.rdbuf()))
    {
    }

    ~CaptureStdout()
    {
        std::cout.rdbuf(original_);
        std::cout.flags(flags_);
        std::cout.precision(precision_);
        std::cout.fill(fill_);
    }

    CaptureStdout(const CaptureStdout&) = delete;
    CaptureStdout& operator=(const CaptureStdout&) = delete;

private:
    std::ios_base::fmtflags flags_;
    std::streamsize precision_;
    char fill_;
    std::streambuf* original_;
};

struct PreviewSelection
{
    WindowState selected_window;
    std::string output;
    std::size_t eligible_count = 0;
    std::size_t retained_count = 0;
    std::uint64_t retained_hash = 0;
    std::size_t committed_count = 0;
    std::string reason;
};

PreviewSelection previewSelection(
    const PipelineRuntime& runtime,
    const WindowState& original_window,
    const DependencySupportPostProcessor& processor,
    const double epsilon,
    const DependencySupportBatch& batch)
{
    PreviewSelection result{.selected_window=original_window};
    result.selected_window.invalidateCurrentWindow();
    auto copied_batch=batch;
    std::ostringstream captured;
    {
        CaptureStdout redirect(captured);
        result.committed_count=processor.processPwph(
            runtime,result.selected_window,epsilon,std::move(copied_batch));
    }
    result.output=captured.str();
    std::istringstream lines(result.output);
    std::string line;
    std::unordered_map<std::string,std::string> fields;
    bool found=false;
    while(std::getline(lines,line))
    {
        if(!line.starts_with("[postcommit] "))
            continue;
        if(found)
            throw std::logic_error("EXTEND preview emitted repeated postcommit records");
        found=true;
        std::istringstream tokens(line);
        std::string token;
        tokens>>token;
        while(tokens>>token)
        {
            const auto equals=token.find('=');
            if(equals!=std::string::npos)
                fields.emplace(token.substr(0,equals),token.substr(equals+1));
        }
    }
    if(!found)
        throw std::logic_error("EXTEND preview is missing its postcommit receipt");
    result.eligible_count=std::stoull(fields.at("eligible"));
    result.retained_count=std::stoull(fields.at("retained"));
    result.retained_hash=std::stoull(fields.at("retained_hash"));
    if(result.committed_count!=std::stoull(fields.at("committed")))
        throw std::logic_error("EXTEND preview commit count disagrees with receipt");
    result.reason=fields.at("reason");
    return result;
}

void emitWindow(const DimensionFrame& frame,const double lo,const double hi,
                const DimensionFrame* upper=nullptr)
{
    const auto precision=std::cout.precision();
    std::cout<<"window lower bound, window upper bound\n"
             <<std::setprecision(17)<<lo<<", "<<hi<<'\n'
             <<"dimension, birth weight, death weight\n";
    for(const auto* shown:{&frame,upper})
    {
        if(!shown)
            continue;
        for(const auto& pair:shown->persistentPairs())
        {
            if(pair.cofacet_weight>=lo || pair.cofacet_weight<0.0)
                std::cout<<shown->dimension()<<", "<<pair.facet_weight<<", "
                         <<pair.cofacet_weight<<'\n';
        }
    }
    std::cout<<'\n';
    std::cout.precision(precision);
}

std::size_t liveClassCount(const DimensionFrame& frame)
{
    return static_cast<std::size_t>(std::count_if(
        frame.persistentPairs().begin(),frame.persistentPairs().end(),
        [](const PersistentPairInfo& pair){return pair.cofacet_weight<0.0;}));
}

bool samePairs(const std::vector<PersistentPairInfo>& lhs,const std::vector<PersistentPairInfo>& rhs)
{
    return std::equal(lhs.begin(),lhs.end(),rhs.begin(),rhs.end(),
        [](const PersistentPairInfo& a,const PersistentPairInfo& b)
        {
            return a.facet_bindex==b.facet_bindex && a.cofacet_bindex==b.cofacet_bindex &&
                   std::bit_cast<std::uint64_t>(a.facet_weight)==std::bit_cast<std::uint64_t>(b.facet_weight) &&
                   std::bit_cast<std::uint64_t>(a.cofacet_weight)==std::bit_cast<std::uint64_t>(b.cofacet_weight);
        });
}

// Gate only: the retained frames at this cut must equal a fresh direct window
// [0, epsilon) on its own WindowState, built by the ordinary non-retained path
// with the geometric enumerators for both triangles and tetrahedra.
void checkAgainstDirect(const PipelineRuntime& runtime,const WindowState& state,const double eps_lo,
                        const double epsilon,const bool final,
                        const DimensionFrame& lower,const DimensionFrame* upper)
{
    // A copy of the current quotient (no PVs before the first commit), prepared from scratch.
    WindowState direct_state=state;
    direct_state.invalidateCurrentWindow();
    auto edges=direct_state.prepareWindow(runtime,{.eps_lo=eps_lo,.eps_hi=epsilon,.is_final=final});
    std::string failure;
    const auto check=[&](const bool same,const char* what){if(!same && failure.empty())failure=what;};
    std::size_t direct_lower_pairs=0;
    DependencySupportBatch direct_batch;
    {
        DimensionFrame direct(runtime,direct_state,std::move(edges));
        check(direct.facetList()==lower.facetList(),"edges");
        check(direct.cofacetList()==lower.cofacetList(),"triangles");
        direct.matchPersistence(true);
        check(samePairs(direct.persistentPairs(),lower.persistentPairs()),"H1-pairs");
        direct_lower_pairs=direct.persistentPairs().size();
        if(upper)
        {
            check(direct.advance(),"advance");
            check(direct.facetList()==upper->facetList(),"triangle-facets");
            check(direct.cofacetList()==upper->cofacetList(),"tetrahedra");
            direct.matchPersistence(true);
            check(samePairs(direct.persistentPairs(),upper->persistentPairs()),"H2-pairs");
        }
        direct_batch=direct.takeDependencySupportBatch();
    }
    auto chain_batch=(upper?*upper:lower).snapshotDependencySupportBatch();
    if(upper)
    {
        const auto lower_batch=lower.snapshotDependencySupportBatch();
        chain_batch.protected_labels.insert(lower_batch.protected_labels.begin(),lower_batch.protected_labels.end());
    }
    check(direct_batch.protected_labels==chain_batch.protected_labels,"protected-labels");
    check(std::equal(direct_batch.supports.begin(),direct_batch.supports.end(),
                     chain_batch.supports.begin(),chain_batch.supports.end(),
                     [](const DependencySupport& a,const DependencySupport& b){return a.labels==b.labels;}),
          "supports");
    const auto precision=std::cout.precision();
    std::cout<<std::setprecision(17)<<"[extend-gate] epsilon="<<epsilon
             <<" edges="<<lower.facetCount()<<" triangles="<<lower.cofacetCount()
             <<" tetrahedra="<<(upper?upper->cofacetCount():0)
             <<" pairs="<<direct_lower_pairs<<','<<(upper?upper->persistentPairs().size():0)
             <<" supports="<<chain_batch.supports.size()<<" protected="<<chain_batch.protected_labels.size()
             <<" status="<<(failure.empty()?"MATCH":"MISMATCH:"+failure)<<'\n';
    std::cout.precision(precision);
    if(!failure.empty())
        throw std::runtime_error("EXTEND gate: retained frames differ from the direct window ("+failure+")");
}

void dumpSnapshot(
    std::ostream& out,
    const char* kind,
    const double epsilon,
    const double previous_epsilon,
    const bool final,
    const DimensionFrame& frame,
    const DependencySupportBatch& batch,
    const PreviewSelection& preview)
{
    const bool quotient=batch.protected_labels.empty() && preview.committed_count>0;
    out<<"{\"schema\":\"EXTEND_GATE_V1\",\"kind\":\""<<kind
       <<"\",\"eps_bits\":"<<std::bit_cast<std::uint64_t>(epsilon)
       <<",\"previous_eps_bits\":"<<std::bit_cast<std::uint64_t>(previous_epsilon)
       <<",\"facet_count\":"<<frame.facetCount()
       <<",\"cofacet_count\":"<<frame.cofacetCount()
       <<",\"pair_count\":"<<frame.persistentPairs().size()
       <<",\"pairs\":[";
    bool first=true;
    for(const auto& pair:frame.persistentPairs())
    {
        if(!first)out<<',';
        first=false;
        out<<'['<<pair.facet_bindex<<','<<pair.cofacet_bindex<<','
           <<std::bit_cast<std::uint64_t>(pair.facet_weight)<<','
           <<std::bit_cast<std::uint64_t>(pair.cofacet_weight)<<']';
    }
    std::vector<std::size_t> protected_labels(
        batch.protected_labels.begin(),batch.protected_labels.end());
    std::sort(protected_labels.begin(),protected_labels.end());
    out<<"],\"protected_count\":"<<protected_labels.size()
       <<",\"protected_labels\":[";
    for(std::size_t i=0;i<protected_labels.size();++i)
    {
        if(i)out<<',';
        out<<protected_labels[i];
    }
    out<<"],\"support_count\":"<<batch.supports.size()<<",\"supports\":[";
    for(std::size_t i=0;i<batch.supports.size();++i)
    {
        if(i)out<<',';
        out<<'[';
        for(std::size_t j=0;j<batch.supports[i].labels.size();++j)
        {
            if(j)out<<',';
            out<<batch.supports[i].labels[j];
        }
        out<<']';
    }
    out<<"],\"eligible_count\":"<<preview.eligible_count
       <<",\"retained_count\":"<<preview.retained_count
       <<",\"retained_hash\":"<<preview.retained_hash
       <<",\"committed_count\":"<<preview.committed_count
       <<",\"action\":\""<<(final?"FINAL":quotient?"QUOTIENT":"DEFER")
       <<"\",\"reason\":\""<<(final?"final-window":preview.reason)
       <<"\",\"preview_action\":\""<<(preview.committed_count?"QUOTIENT":"DEFER")
       <<"\",\"preview_reason\":\""<<preview.reason<<"\"}\n";
    out.flush();
    if(!out)
        throw std::runtime_error("Could not write EXTEND diagnostic snapshot");
}
}

void PwphPipeline::runGuardExtendLadder()
{
    const auto& config=runtime_.config();
    std::ofstream dump;
    if(!config.extend_dump.empty())
    {
        dump.open(config.extend_dump,std::ios::out|std::ios::trunc);
        if(!dump)
            throw std::runtime_error("Could not open EXTEND diagnostic file");
    }

    const double endpoint=config.finalEpsilon();
    const double first=config.eps_breaks.front();
    const WindowBounds first_bounds{.eps_lo=0.0,.eps_hi=first,.is_final=first==endpoint};
    auto first_edges=window_state_.prepareWindow(runtime_,first_bounds);

    if(config.extend_direct_probe)
    {
        // Independent oracle: use the original non-retained constructor and
        // ordinary matching, never the incremental initialization path.
        DimensionFrame direct(runtime_,window_state_,std::move(first_edges));
        direct.matchPersistence(true);
        const auto batch=direct.snapshotDependencySupportBatch();
        const auto preview=previewSelection(runtime_,window_state_,support_processor_,first,batch);
        dumpSnapshot(dump,"DIRECT",first,0.0,first==endpoint,direct,batch,preview);
        return;
    }

    auto frame=std::make_unique<DimensionFrame>(
        runtime_,window_state_,std::move(first_edges),true);
    frame->matchPersistenceForExtension(true);
    double previous=0.0;
    for(std::size_t cut=0;cut<config.eps_breaks.size();++cut)
    {
        const double epsilon=config.eps_breaks[cut];
        const bool final=epsilon==endpoint;
        if(cut)
        {
            window_state_.extendFirstWindowTo(epsilon);
            frame->extendTo(epsilon);
        }
        if(final && !dump.is_open())
        {
            emitWindow(*frame,0.0,endpoint);
            frame.reset();
            window_state_.invalidateCurrentWindow();
            return;
        }
        const auto batch=frame->snapshotDependencySupportBatch();
        auto preview=previewSelection(runtime_,window_state_,support_processor_,epsilon,batch);
        if(dump.is_open())
            dumpSnapshot(dump,config.extend_force_chain?"CHAIN":"POLICY",epsilon,
                         previous,final,*frame,batch,preview);
        if(final)
        {
            emitWindow(*frame,0.0,endpoint);
            frame.reset();
            window_state_.invalidateCurrentWindow();
            return;
        }

        const bool quotient=batch.protected_labels.empty() && preview.committed_count>0;
        if(quotient && !config.extend_force_chain)
        {
            // Close the first window once. The clone already contains the
            // exact frozen selection, so no second selector/risk pass occurs.
            emitWindow(*frame,0.0,epsilon);
            frame.reset();
            window_state_.invalidateCurrentWindow();
            window_state_=std::move(preview.selected_window);
            std::cout<<preview.output;
            runtime_.ensureBinomialCapacity(window_state_.totalLabelCount());

            const WindowBounds final_bounds{
                .eps_lo=epsilon,.eps_hi=endpoint,.is_final=true};
            auto final_edges=window_state_.prepareWindow(runtime_,final_bounds);
            {
                DimensionFrame final_frame(runtime_,window_state_,std::move(final_edges));
                final_frame.matchPersistence(false);
                emitWindow(final_frame,epsilon,endpoint);
            }
            window_state_.invalidateCurrentWindow();
            return;
        }
        // A zero-commit preview only consumed copies. Original labels, window
        // caches and the live matching frame remain available for extension.
        previous=epsilon;
    }
    throw std::logic_error("EXTEND ladder ended without its final endpoint");
}

std::optional<std::size_t> PwphPipeline::runLiveGuardExtendPrefix()
{
    // Under the live guard a cut is admissible only when no class of any
    // dimension is alive. Instead of rebuilding the first window from scale zero
    // at every vetoed cut, retain its frames and extend them to the next cut.
    // Returns the cut of the first committed quotient, after which the plain
    // ladder continues, or nothing when the first window reached the endpoint.
    const auto& config=runtime_.config();
    const double endpoint=config.finalEpsilon();
    const double first=config.eps_breaks.front();
    // --extend-worth-log (dry run): time and predicted top-dimension work of the build and of each extension.
    const auto worth_window=[&](const char* mode,double lo,double hi,std::chrono::steady_clock::time_point start)
    {
        if(!config.extend_worth_log)return;
        const double core_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const auto& d=runtime_.distanceMatrix();
        std::cout<<std::setprecision(17)<<"[worth-window] mode="<<mode<<" lo="<<lo<<" hi="<<hi<<" core_s="<<core_s
            <<" W_lo="<<(lo>0?PackDmtPayback::topWork(d,window_state_,{},lo,config.maxdim):0.0)
            <<" W_hi="<<PackDmtPayback::topWork(d,window_state_,{},hi,config.maxdim)
            <<" active="<<window_state_.activeLabels().size()<<'\n'<<std::setprecision(6);
    };
    const auto build_start=std::chrono::steady_clock::now();
    auto first_edges=window_state_.prepareWindow(
        runtime_,{.eps_lo=0.0,.eps_hi=first,.is_final=first==endpoint});
    auto lower=std::make_unique<DimensionFrame>(runtime_,window_state_,std::move(first_edges),true);
    lower->matchPersistenceForExtension(true);
    std::unique_ptr<DimensionFrame> upper;
    if(config.maxdim==3)
    {
        upper=DimensionFrame::retainedUpperFrame(*lower);
        upper->matchPersistenceForExtension(true);
    }
    worth_window("build",0.0,first,build_start);
    for(std::size_t cut=0;cut<config.eps_breaks.size();++cut)
    {
        const double epsilon=config.eps_breaks[cut];
        const bool final=epsilon==endpoint;
        if(cut)
        {
            const auto extend_start=std::chrono::steady_clock::now();
            window_state_.extendFirstWindowTo(epsilon);
            lower->extendTo(epsilon);
            if(upper)
                upper->extendWith(*lower);
            worth_window("extend",config.eps_breaks[cut-1],epsilon,extend_start);
        }
        if(config.extend_gate)
            checkAgainstDirect(runtime_,WindowState(runtime_.distanceMatrix().vertexCount()),0.0,
                               epsilon,final,*lower,upper.get());
        if(final)
        {
            emitWindow(*lower,0.0,endpoint,upper.get());
            upper.reset();
            lower.reset();
            window_state_.invalidateCurrentWindow();
            return std::nullopt;
        }

        const std::size_t live_h1=liveClassCount(*lower);
        const std::size_t live_h2=upper?liveClassCount(*upper):0;
        const auto report=[&](const char* action)
        {
            std::cout<<"[birth-repair] epsilon="<<std::setprecision(17)<<epsilon<<std::setprecision(6)
                     <<" live="<<live_h1;
            if(upper)
                std::cout<<','<<live_h2;
            std::cout<<" action="<<action<<'\n';
        };
        if(config.extend_gate || live_h1+live_h2>0)
        {
            report(config.extend_gate?"CHAIN":"VETO");
            continue;
        }
        // Every class is dead, so no label is protected; the supports are the
        // top frame's, as in the plain ladder.
        const auto batch=(upper?*upper:*lower).snapshotDependencySupportBatch();
        auto preview=previewSelection(runtime_,window_state_,support_processor_,epsilon,batch);
        if(preview.committed_count==0)
        {
            // Nothing to collapse: the frames stay valid and keep extending.
            report("ALLOW");
            std::cout<<preview.output;
            continue;
        }
        emitWindow(*lower,0.0,epsilon,upper.get());
        report("ALLOW");
        upper.reset();
        lower.reset();
        window_state_.invalidateCurrentWindow();
        window_state_=std::move(preview.selected_window);
        std::cout<<preview.output;
        runtime_.ensureBinomialCapacity(window_state_.totalLabelCount());
        return cut;
    }
    throw std::logic_error("Live-guard EXTEND ended without its final endpoint");
}

// --extend-all: every window, not only the first, is retained and extended through the cuts where it
// does not commit. A segment builds the window [eps_lo, cut) on the current quotient from scale zero and
// extends it cut by cut, exactly as the first window is extended before the first commit. At an admitted
// cut with supports left after the filters, --extend-rule may extend instead of committing.
// Returns the cut of the next commit, or nothing when the window reached the endpoint.
std::optional<std::size_t> PwphPipeline::runExtendSegment(std::size_t first_cut,const double eps_lo)
{
    const auto& config=runtime_.config();
    // Once no more commits can come, the window is built straight to the endpoint.
    if(no_more_commits_)
        first_cut=config.eps_breaks.size()-1;
    const double endpoint=config.finalEpsilon();
    const double first=config.eps_breaks.at(first_cut);
    using Clock=std::chrono::steady_clock;
    const bool log=config.extend_worth_log;
    const auto work=[&](const double x)
    {return PackDmtPayback::topWork(runtime_.distanceMatrix(),window_state_,{},x,config.maxdim);};
    // Window times and predicted work, for the log and the rule's calibration.
    const auto note=[&](const char* mode,const double lo,const double hi,const Clock::time_point start,const double prep_s)
    {
        const double seconds=std::chrono::duration<double>(Clock::now()-start).count();
        if(!log && !config.extend_rule)return;
        const double w_lo=std::string(mode)=="extend"?work(lo):0.0,w_hi=work(hi);
        rule_calibration_.add(std::string(mode)=="extend",seconds,w_hi-w_lo);
        if(log)
            std::cout<<std::setprecision(17)<<"[worth-window] mode="<<mode<<" lo="<<lo<<" hi="<<hi<<" core_s="<<seconds
                <<" W_lo="<<w_lo<<" W_hi="<<w_hi<<" active="<<window_state_.activeLabels().size()
                <<" prep_s="<<prep_s<<'\n'<<std::setprecision(6);
    };
    const auto build_start=Clock::now();
    auto edges=window_state_.prepareWindow(runtime_,{.eps_lo=eps_lo,.eps_hi=first,.is_final=first==endpoint});
    // The fixed part of a new window: label distances and edges, before any simplex is enumerated.
    rule_calibration_.prep_s=std::chrono::duration<double>(Clock::now()-build_start).count();
    auto lower=std::make_unique<DimensionFrame>(runtime_,window_state_,std::move(edges),true);
    lower->matchPersistenceForExtension(true);
    std::unique_ptr<DimensionFrame> upper;
    if(config.maxdim==3)
    {
        upper=DimensionFrame::retainedUpperFrame(*lower);
        upper->matchPersistenceForExtension(true);
    }
    note(window_state_.pseudoVertices().empty()?"build":"rebuild",eps_lo,first,build_start,rule_calibration_.prep_s);
    double previous=first;
    for(std::size_t cut=first_cut;cut<config.eps_breaks.size();++cut)
    {
        // Once no more commits can come, the window goes straight to the endpoint.
        if(cut>first_cut && no_more_commits_)
            cut=config.eps_breaks.size()-1;
        const double epsilon=config.eps_breaks[cut];
        const bool final=epsilon==endpoint;
        if(cut>first_cut)
        {
            const auto extend_start=Clock::now();
            window_state_.extendPreparedWindowTo(epsilon,final);
            lower->extendTo(epsilon);
            if(upper)
                upper->extendWith(*lower);
            note("extend",previous,epsilon,extend_start,0.0);
            previous=epsilon;
        }
        if(config.extend_gate)
            checkAgainstDirect(runtime_,window_state_,eps_lo,epsilon,final,*lower,upper.get());
        if(final)
        {
            emitWindow(*lower,eps_lo,endpoint,upper.get());
            upper.reset();
            lower.reset();
            window_state_.invalidateCurrentWindow();
            return std::nullopt;
        }
        // No more commits (the rule stopped, or the commit limit is reached): only extend, without selecting.
        if(no_more_commits_)
            continue;
        const std::size_t live_h1=liveClassCount(*lower);
        const std::size_t live_h2=upper?liveClassCount(*upper):0;
        const auto report=[&](const char* action)
        {
            std::cout<<"[birth-repair] epsilon="<<std::setprecision(17)<<epsilon<<std::setprecision(6)
                     <<" live="<<live_h1;
            if(upper)
                std::cout<<','<<live_h2;
            std::cout<<" action="<<action<<'\n';
        };
        if(live_h1+live_h2>0)
        {
            report("VETO");
            continue;
        }
        const auto batch=(upper?*upper:*lower).snapshotDependencySupportBatch();
        auto preview=previewSelection(runtime_,window_state_,support_processor_,epsilon,batch);
        if(preview.committed_count==0)
        {
            report("ALLOW");
            std::cout<<preview.output;
            if(support_processor_.lastStoppedBeforeFilters())
                no_more_commits_=true;
            continue;
        }
        if(config.extend_rule && !ruleCommits(epsilon,cut))
        {
            // The rule compares committing now with never committing again; once it prefers the latter,
            // the run only extends.
            report("ALLOW");
            no_more_commits_=true;
            continue;
        }
        ++commits_done_;
        if(config.extend_commit_limit>=0 && commits_done_>=static_cast<std::size_t>(config.extend_commit_limit))
            no_more_commits_=true;
        emitWindow(*lower,eps_lo,epsilon,upper.get());
        report("ALLOW");
        upper.reset();
        lower.reset();
        window_state_.invalidateCurrentWindow();
        window_state_=std::move(preview.selected_window);
        std::cout<<preview.output;
        runtime_.ensureBinomialCapacity(window_state_.totalLabelCount());
        return cut;
    }
    throw std::logic_error("EXTEND segment ended without its final endpoint");
}

void PwphPipeline::runExtendAll()
{
    no_more_commits_=runtime_.config().extend_commit_limit==0;
    if(runtime_.config().extend_rule)
        support_processor_.setPreFilterRule(
            [this](const DependencySupportPostProcessor::Worth& w){return ruleDecide(w,"before-filters");});
    std::size_t first_cut=0;
    double eps_lo=0.0;
    for(;;)
    {
        const auto committed=runExtendSegment(first_cut,eps_lo);
        if(!committed)
            return;
        eps_lo=runtime_.config().eps_breaks[*committed];
        first_cut=*committed+1;
    }
}

bool PwphPipeline::ruleCommits(const double epsilon,const std::size_t cut)
{
    const auto& w=support_processor_.lastWorth();
    if(w.epsilon!=epsilon)
        throw std::logic_error("--extend-rule: no predicted work for this cut");
    (void)cut;
    return ruleDecide(w,"after-filters");
}

bool PwphPipeline::ruleDecide(const DependencySupportPostProcessor::Worth& w,const char* stage)
{
    // Rollout of "never commit again": commit now and then build one window to the endpoint (a new window
    // to the next cut, extended from there), against extending the current window to the endpoint. A
    // window built from scale zero costs tb per unit of predicted top-dimension work, an extension te per
    // unit it adds; both are measured on the windows of this run so far. The selection and filter time is
    // already spent when the rule decides, and does not enter.
    // One time per unit of work t for both branches, over every window so far (builds and extensions,
    // weighted by their work, so the large windows dominate); a new window also pays the fixed time of
    // preparing its label distances and edges, measured on the latest build.
    // Stopping cannot be undone, while a commit keeps the choice to commit again: extending to the
    // endpoint must stay cheaper even at the rate the latest window actually ran at, when that was higher
    // than predicted.
    const auto& c=rule_calibration_;
    const double t=c.rate();
    const double margin=std::max(1.0,c.last_ratio);
    const double commit_cost=w.C_s+c.prep_s+t*w.wsel_end;
    const double extend_cost=margin*t*(w.w0_end-w.w0_r);
    const bool commit=commit_cost<extend_cost;
    std::cout<<std::setprecision(17)<<"[rule] epsilon="<<w.epsilon<<" stage="<<stage<<" t="<<t<<" margin="<<margin
             <<" prep_s="<<c.prep_s<<" commit_cost="<<commit_cost<<" extend_cost="<<extend_cost
             <<" action="<<(commit?"COMMIT":"EXTEND")<<'\n'<<std::setprecision(6);
    return commit;
}
