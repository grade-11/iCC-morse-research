#include "Pipelines.hpp"
#include "DimensionFrame.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
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

void emitWindow(const DimensionFrame& frame,const double lo,const double hi)
{
    const auto precision=std::cout.precision();
    std::cout<<"window lower bound, window upper bound\n"
             <<std::setprecision(17)<<lo<<", "<<hi<<'\n'
             <<"dimension, birth weight, death weight\n";
    for(const auto& pair:frame.persistentPairs())
    {
        if(pair.cofacet_weight>=lo || pair.cofacet_weight<0.0)
            std::cout<<frame.dimension()<<", "<<pair.facet_weight<<", "
                     <<pair.cofacet_weight<<'\n';
    }
    std::cout<<'\n';
    std::cout.precision(precision);
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
