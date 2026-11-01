#pragma once
// Enabled-only observations of harness-prescribed branches. No model or action.
#include <algorithm>
#include <array>
#include <bit>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <memory>
#include <iterator>
#include <initializer_list>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <time.h>
#include <unistd.h>
#include <utility>

namespace PolicyCost
{
constexpr size_t N_CAP=1024, EVENT_CAP=64, JSON_CAP=32768, NO_WINDOW=std::numeric_limits<size_t>::max();
static_assert(sizeof(double)==8 && sizeof(uint64_t)==8 && std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53);
inline void need(bool good,const char*why){if(!good)throw std::runtime_error(why);}
inline uint64_t nowNs(){timespec t{};need(clock_gettime(CLOCK_MONOTONIC,&t)==0&&t.tv_sec>=0&&t.tv_nsec>=0&&t.tv_nsec<1000000000,"COST_CLOCK_FAILURE");return uint64_t(t.tv_sec)*1000000000ULL+uint64_t(t.tv_nsec);}
inline uint64_t choose2(uint64_t n){return n<2?0:n*(n-1)/2;}
inline uint64_t choose3(uint64_t n){return n<3?0:n*(n-1)*(n-2)/6;}
inline uint64_t word(double value){return std::bit_cast<uint64_t>(value);}
enum class Script { Close, Defer, Try };
inline Script parseScript(std::string_view s){if(s=="CLOSE_PAID")return Script::Close;if(s=="OPEN_DEFER")return Script::Defer;if(s=="OPEN_TRY_COMMIT")return Script::Try;throw std::runtime_error("COST_SCRIPT");}
inline const char* scriptName(Script s){switch(s){case Script::Close:return "CLOSE_PAID";case Script::Defer:return "OPEN_DEFER";case Script::Try:return "OPEN_TRY_COMMIT";}throw std::runtime_error("COST_SCRIPT_ENUM");}
struct Plan { Script script=Script::Close;size_t n=0;double E=0,q=0;bool forced_defer=false; };
inline void validatePlan(const Plan&p)
{
    need(p.n>=2&&p.n<=N_CAP&&std::isfinite(p.E)&&p.E>0,"COST_PLAN_DOMAIN");
    need(p.forced_defer==(p.script==Script::Defer),"COST_SCRIPT_ACTION");
    if(p.script==Script::Close)need(p.q==0,"COST_CLOSE_Q");
    else need(std::isfinite(p.q)&&p.q>0&&p.q<p.E,"COST_OPEN_Q");
}
struct Metadata {std::string_view run,source,binary,input;uint64_t wrapper=0;};
struct Event {const char*stage=nullptr;uint64_t ns=0,lo=0,hi=0;size_t window=NO_WINDOW,dimension=0;};
struct Opening {bool observed=false,has_q=false;uint64_t q=0,ratio=0,mE=0,mq=0,possible=0,u2sum=0,u2=0,regions=0,reads=0;};
struct Retained
{
    bool observed=false,skipped=false;uint64_t raw=0,protected_count=0,eligible=0,count=0,hash=0,warnings=0,A=0,P=0,L=0;
    uint64_t sum=0,square=0,cube=0,tuples=0,gain=0,proposedA=0,proposedP=0,proposedL=0,pairs=0;
};
struct Actual {bool observed=false;uint64_t committed=0,A=0,P=0,L=0;};
struct Json
{
    std::array<char,JSON_CAP> bytes{};size_t size=0;
    void add(std::string_view s){need(s.size()<=bytes.size()-size,"COST_JSON_CAP");std::copy(s.begin(),s.end(),bytes.begin()+size);size+=s.size();}
    void number(uint64_t x){char b[32];auto r=std::to_chars(b,b+sizeof(b),x);need(r.ec==std::errc{},"COST_NUMBER");add({b,size_t(r.ptr-b)});}
    void boolean(bool x){add(x?"true":"false");}
    void quote(std::string_view s){add("\"");for(unsigned char c:s){need(c>=32&&c<127,"COST_ASCII");if(c=='"'||c=='\\')add("\\");char b=char(c);add({&b,1});}add("\"");}
    void hex(uint64_t x){constexpr char h[]="0123456789abcdef";char b[16];for(size_t i=0;i<16;++i){b[15-i]=h[x&15];x>>=4;}quote({b,16});}
    void nullable(bool yes,uint64_t x){if(yes)number(x);else add("null");}
    void nullableHex(bool yes,uint64_t x){if(yes)hex(x);else add("null");}
};
class Logger;
inline Logger* active=nullptr;
class Logger
{
    int fd_=-1;Plan plan_;uint64_t native_=0,wrapper_=0;
    std::array<char,65> run_{},source_{},binary_{},input_{};
    std::array<Event,EVENT_CAP> events_{};size_t count_=0,window_=NO_WINDOW,windows_=0;uint64_t lo_=0,hi_=0;
    Opening opening_;Retained retained_;Actual actual_;Json json_;
    static bool token(std::string_view s){return !s.empty()&&s.size()<=64&&std::all_of(s.begin(),s.end(),[](unsigned char c){return(c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='.'||c=='_'||c=='-';});}
    static bool sha(std::string_view s){return s.size()==64&&std::all_of(s.begin(),s.end(),[](char c){return(c>='0'&&c<='9')||(c>='a'&&c<='f');});}
    static std::string_view env(const char*n){const char*v=std::getenv(n);need(v&&*v,"COST_REQUIRED_ENV");return v;}
    static void copy(std::array<char,65>&out,std::string_view in){need(in.size()<out.size(),"COST_METADATA_SIZE");std::copy(in.begin(),in.end(),out.begin());}
    static bool known(std::string_view stage)
    {
        constexpr std::string_view names[]={"PIPELINE_BEGIN","GEO_PREPASS_BEGIN","GEO_PREPASS_END","OPENING_READY","EDGE_REBUILD_BEGIN","EDGE_REBUILD_END","FRAME_BUILD_BEGIN","FRAME_BUILD_END","MATCH_BEGIN","MATCH_END","ADVANCE_BEGIN","ADVANCE_END","WINDOW_CORE_END","SELECTOR_BEGIN","FILTER_END","RETAINED_PREVALUATION","VALUATION_BEGIN","VALUATION_END","DECISION_READY","MUTATION_BEGIN","MUTATION_END","DSP_RETURN","BINOMIAL_BEGIN","BINOMIAL_END","BASELINE_RESULT","PIPELINE_SUMMARY_BEGIN","PIPELINE_SUMMARY_END","GEO_FINISH_BEGIN","GEO_FINISH_END"};
        return std::find(std::begin(names),std::end(names),stage)!=std::end(names);
    }
    void record(const char*stage,size_t dimension,uint64_t ns)
    {
        need(known(stage)&&count_<events_.size(),"COST_EVENT_CAP_OR_NAME");
        need(count_==0||ns>=events_[count_-1].ns,"COST_EVENT_CLOCK_ORDER");
        const std::string_view s(stage);
        const bool interface=s.starts_with("FRAME_BUILD_")||s.starts_with("MATCH_")||s.starts_with("ADVANCE_");
        need(dimension==(interface?2:0),"COST_EVENT_DIMENSION");
        events_[count_++]={stage,ns,lo_,hi_,window_,dimension};
    }
    void checkOrder()const
    {
        size_t i=0;
        const auto take=[&](const char*stage){need(i<count_&&std::string_view(events_[i].stage)==stage,"COST_EVENT_ORDER");++i;};
        for(const char*s:{"PIPELINE_BEGIN","GEO_PREPASS_BEGIN","GEO_PREPASS_END","OPENING_READY"})take(s);
        const size_t expected_windows=plan_.script==Script::Close?1:2;
        need(windows_==expected_windows&&opening_.observed,"COST_WINDOW_OR_OPENING_CENSUS");
        for(size_t w=0;w<expected_windows;++w)
        {
            for(const char*s:{"EDGE_REBUILD_BEGIN","EDGE_REBUILD_END","FRAME_BUILD_BEGIN","FRAME_BUILD_END","MATCH_BEGIN","MATCH_END","ADVANCE_BEGIN","ADVANCE_END","WINDOW_CORE_END"})take(s);
            if(expected_windows==2&&w==0)
            {
                for(const char*s:{"SELECTOR_BEGIN","FILTER_END","RETAINED_PREVALUATION","VALUATION_BEGIN","VALUATION_END","DECISION_READY","MUTATION_BEGIN","MUTATION_END","DSP_RETURN"})take(s);
                if(actual_.committed){take("BINOMIAL_BEGIN");take("BINOMIAL_END");}
                take("BASELINE_RESULT");
            }
        }
        for(const char*s:{"PIPELINE_SUMMARY_BEGIN","PIPELINE_SUMMARY_END","GEO_FINISH_BEGIN","GEO_FINISH_END"})take(s);
        need(i==count_,"COST_EXTRA_EVENTS");
        need((plan_.script==Script::Close&&!retained_.observed&&!actual_.observed)||
             (plan_.script!=Script::Close&&retained_.observed&&actual_.observed),"COST_SELECTOR_CENSUS");
    }
    void render(uint64_t finalize)
    {
        auto&j=json_;j.size=0;
        j.add("{\"schema\":\"TWO_GATE_COST_V1\",\"status\":\"COMPLETE\",\"mode\":\"HARNESS_PRESCRIBED_INTERVENTION_ONLY\",\"self_verified\":false,\"native_pid\":");j.number(native_);
        j.add(",\"run_id\":");j.quote(run_.data());j.add(",\"declared_source_sha256\":");j.quote(source_.data());j.add(",\"declared_binary_sha256\":");j.quote(binary_.data());j.add(",\"declared_input_sha256\":");j.quote(input_.data());j.add(",\"wrapper_parent_pid\":");j.number(wrapper_);
        j.add(",\"script\":");j.quote(scriptName(plan_.script));j.add(",\"n\":");j.number(plan_.n);j.add(",\"E_bits\":");j.hex(word(plan_.E));j.add(",\"q_bits\":");j.nullableHex(plan_.script!=Script::Close,word(plan_.q));
        const auto&o=opening_;j.add(",\"opening_features\":{\"observed\":true,\"candidate_q_bits\":");j.nullableHex(o.has_q,o.q);j.add(",\"q_over_E_bits\":");j.nullableHex(o.has_q,o.ratio);
        j.add(",\"positive_edges_E\":");j.number(o.mE);j.add(",\"positive_edges_q\":");j.nullable(o.has_q,o.mq);j.add(",\"suffix_positive_edges\":");j.nullable(o.has_q,o.mE-o.mq);j.add(",\"possible_original_pairs\":");j.number(o.possible);
        j.add(",\"density_q_numerator\":");j.nullable(o.has_q,o.mq);j.add(",\"density_E_numerator\":");j.number(o.mE);j.add(",\"density_denominator\":");j.number(o.possible);j.add(",\"U2_numerator\":");j.nullable(o.has_q,o.u2sum);j.add(",\"U2_divisor\":3,\"U2\":");j.nullable(o.has_q,o.u2);j.add(",\"region_count\":");j.number(o.regions);j.add(",\"geometry_pair_reads\":");j.number(o.reads);j.add(",\"additional_distance_reads\":0}");
        const auto&s=retained_;j.add(",\"retained_features\":{\"observed\":");j.boolean(s.observed);
        const auto field=[&](const char*key,uint64_t value){j.add(",");j.quote(key);j.add(":");j.nullable(s.observed,value);};
        field("raw_support_count",s.raw);field("protected_count",s.protected_count);field("eligible_count",s.eligible);field("retained_count",s.count);j.add(",\"retained_hash_bits\":");j.nullableHex(s.observed,s.hash);j.add(",\"remaining_warnings\":");j.nullable(s.observed&&!s.skipped,s.warnings);j.add(",\"skipped_live\":");if(s.observed)j.boolean(s.skipped);else j.add("null");
        field("state_A",s.A);field("state_P",s.P);field("state_L",s.L);field("member_sum",s.sum);field("member_square_sum",s.square);field("member_cube_sum",s.cube);field("block_choose3_sum",s.tuples);field("proposed_gain",s.gain);field("proposed_A",s.proposedA);field("proposed_P",s.proposedP);field("proposed_L",s.proposedL);field("PV_hash_unique_pairs",s.pairs);j.add("}");
        j.add(",\"actual_result\":{\"observed\":");j.boolean(actual_.observed);for(auto entry:{std::pair{"committed_count",actual_.committed},std::pair{"A",actual_.A},std::pair{"P",actual_.P},std::pair{"L",actual_.L}}){j.add(",");j.quote(entry.first);j.add(":");j.nullable(actual_.observed,entry.second);}j.add("},\"events\":[");
        for(size_t i=0;i<count_;++i){const auto&e=events_[i];if(i)j.add(",");j.add("{\"ordinal\":");j.number(i);j.add(",\"stage\":");j.quote(e.stage);j.add(",\"monotonic_ns\":");j.number(e.ns);j.add(",\"window_index\":");j.nullable(e.window!=NO_WINDOW,e.window);j.add(",\"dimension\":");j.nullable(e.dimension!=0,e.dimension);j.add(",\"eps_lo_bits\":");j.nullableHex(e.window!=NO_WINDOW,e.lo);j.add(",\"eps_hi_bits\":");j.nullableHex(e.window!=NO_WINDOW,e.hi);j.add("}");}
        j.add("],\"clock\":\"CLOCK_MONOTONIC\",\"event_cap\":64,\"event_count\":");j.number(count_);j.add(",\"finalization_begin_monotonic_ns\":");j.number(finalize);j.add(",\"finalization_scope\":\"TIMESTAMPS_EXCLUDE_OWN_FINAL_JSON_WRITE_CLOSE_AND_PROCESS_DESTRUCTOR_TAIL; EXTERNAL_WHOLE_PAYS_ALL\",\"forecast_status\":\"UNKNOWN\",\"runtime_model_loaded\":false,\"native_actions_modified\":false}\n");
    }
public:
    Logger(const char*path,Plan plan,Metadata metadata):plan_(plan)
    {
        validatePlan(plan_);native_=uint64_t(getpid());wrapper_=metadata.wrapper;
        need(path&&*path&&path[0]=='/'&&std::strlen(path)<=4096,"COST_OUTPUT_PATH");
        need(token(metadata.run)&&sha(metadata.source)&&sha(metadata.binary)&&sha(metadata.input)&&wrapper_>0&&wrapper_!=native_,"COST_METADATA");
        copy(run_,metadata.run);copy(source_,metadata.source);copy(binary_,metadata.binary);copy(input_,metadata.input);
        fd_=open(path,O_WRONLY|O_CREAT|O_EXCL|O_CLOEXEC,0600);need(fd_>=0,"COST_EXCLUSIVE_OPEN");
    }
    ~Logger()noexcept{if(active==this)active=nullptr;if(fd_>=0)close(fd_);}
    Logger(const Logger&)=delete;Logger&operator=(const Logger&)=delete;
    template<class Metric,class Config>static std::unique_ptr<Logger> create(const Metric&metric,const Config&cfg)
    {
        const char*path=std::getenv("POLICY_COST_FILE");if(!path)return {};
        const uint64_t begin=nowNs();need(active==nullptr,"COST_NESTED_LOGGER");
        const auto geo=env("GEO_SHADOW_FILE");need(std::string_view(path)!=geo,"COST_GEO_SAME_PATH");
        const Script script=parseScript(env("POLICY_COST_SCRIPT"));
        need(cfg.maxdim==2&&cfg.threads==4&&cfg.h1_risk=="selective-fast"&&cfg.pv_cap_scale==.85&&cfg.pv_min_separation==0&&cfg.pack_mode.empty()&&cfg.g_model.empty()&&cfg.ac_policy.empty()&&!cfg.adaptive&&!cfg.verbose&&!cfg.trace_pwph&&!cfg.h1_realizations&&!cfg.cost_telemetry&&!cfg.calibration_targets&&cfg.stop_after_window==0,"COST_NATIVE_DOMAIN");
        need(cfg.eps_breaks.size()==(script==Script::Close?1:2)&&
             (cfg.quotient_mask.empty()||(script!=Script::Close&&cfg.quotient_mask.size()==1&&cfg.quotient_mask[0]==1)),"COST_SCRIPT_SCHEDULE");
        need(cfg.payback_action==(script==Script::Defer?"defer":"auto"),"COST_SCRIPT_CONFIG_ACTION");
        Plan plan{script,metric.vertexCount(),cfg.finalEpsilon(),script==Script::Close?0:cfg.eps_breaks.front(),script==Script::Defer};
        const auto wrapper=env("GEO_SHADOW_WRAPPER_PARENT_PID");uint64_t pid=0;auto result=std::from_chars(wrapper.data(),wrapper.data()+wrapper.size(),pid);need(result.ec==std::errc{}&&result.ptr==wrapper.data()+wrapper.size(),"COST_WRAPPER_PID");
        Metadata metadata{env("GEO_SHADOW_RUN_ID"),env("GEO_SHADOW_SOURCE_SHA256"),env("GEO_SHADOW_BINARY_SHA256"),env("GEO_SHADOW_INPUT_SHA256"),pid};
        auto logger=std::make_unique<Logger>(path,plan,metadata);active=logger.get();logger->record("PIPELINE_BEGIN",0,begin);return logger;
    }
    void mark(const char*stage,size_t dimension=0){record(stage,dimension,nowNs());}
    void setWindow(size_t index,double lo,double hi)
    {
        const bool close=plan_.script==Script::Close;
        need(index==windows_&&index<(close?1:2)&&word(lo)==word(index==0?0:plan_.q)&&word(hi)==word(close||index==1?plan_.E:plan_.q),"COST_ACTUAL_WINDOW");
        ++windows_;window_=index;lo_=word(lo);hi_=word(hi);
    }
    void global(){window_=NO_WINDOW;lo_=hi_=0;}
    template<class Features>void opening(const Features&f)
    {
        need(!opening_.observed&&f.n==plan_.n&&word(f.E)==word(plan_.E),"COST_OPENING_CONTEXT");
        auto&o=opening_;o.observed=true;o.mE=f.positive;o.possible=choose2(plan_.n);o.regions=f.regions.size();o.reads=f.pair_reads;
        need(o.mE<=o.possible&&o.regions>0&&o.regions<=plan_.n,"COST_OPENING_COUNTS");
        for(const auto&a:f.aliases)if(a.name=="q50"&&a.available)
        {
            need(!o.has_q&&std::isfinite(a.successor)&&a.successor>0&&a.successor<plan_.E,"COST_Q50_ALIAS");
            o.has_q=true;o.q=word(a.successor);o.ratio=word(a.successor/plan_.E);bool found=false;
            for(const auto&c:f.candidates)if(word(c.q)==o.q){need(!found,"COST_Q50_DUPLICATE");found=true;o.mq=c.edges;o.u2sum=c.u2numerator;o.u2=c.u2numerator/3;}
            need(found&&o.mq<=o.mE,"COST_Q50_GRAPH");
        }
        need(plan_.script==Script::Close||(o.has_q&&o.q==word(plan_.q)),"COST_SCRIPT_Q50_IDENTITY");
    }
    template<class Selected>void retained(uint64_t raw,uint64_t protected_count,uint64_t eligible,uint64_t warnings,bool skipped,uint64_t hash,const Selected&selected,uint64_t A,uint64_t P,uint64_t L)
    {
        need(plan_.script!=Script::Close&&!retained_.observed&&window_==0&&A==plan_.n&&P==0&&L==plan_.n,"COST_RETAINED_CONTEXT");
        auto&s=retained_;s.observed=true;s.raw=raw;s.protected_count=protected_count;s.eligible=eligible;s.count=selected.size();s.warnings=warnings;s.skipped=skipped;s.hash=hash;s.A=A;s.P=P;s.L=L;
        need(s.count<=eligible&&eligible<=raw&&protected_count<=plan_.n&&skipped==(protected_count>0),"COST_RETAINED_COUNTS");
        std::array<uint8_t,N_CAP> seen{};uint64_t check=1469598103934665603ULL;
        for(const auto&pv:selected)
        {
            const uint64_t n=pv.representatives.size();need(n>0&&n<=plan_.n,"COST_RETAINED_MEMBER_COUNT");s.sum+=n;s.square+=n*n;s.cube+=n*n*n;s.tuples+=choose3(n);check^=n;check*=1099511628211ULL;
            for(auto id:pv.representatives){need(id<plan_.n&&!seen[id],"COST_RETAINED_DISJOINT");seen[id]=1;check^=id;check*=1099511628211ULL;}
        }
        need(check==hash&&s.sum<=plan_.n,"COST_RETAINED_HASH");s.gain=s.sum-s.count;s.proposedA=plan_.n-s.gain;s.proposedP=s.count;s.proposedL=plan_.n+s.count;need(s.proposedA>=s.count,"COST_PROPOSED_PARTITION");s.pairs=s.count*(s.proposedA-s.count)+choose2(s.count);mark("RETAINED_PREVALUATION");
    }
    void actual(uint64_t committed,uint64_t A,uint64_t P,uint64_t L)
    {
        need(!actual_.observed&&retained_.observed&&window_==0,"COST_ACTUAL_CONTEXT");
        const uint64_t expected=plan_.script==Script::Defer?0:retained_.count;
        need(committed==expected&&P==committed&&L==plan_.n+committed&&A==(committed?retained_.proposedA:plan_.n),"COST_ACTUAL_SCRIPT_RESULT");
        actual_={true,committed,A,P,L};mark("BASELINE_RESULT");
    }
    void finish()
    {
        checkOrder();const uint64_t begin=nowNs();need(count_&&begin>=events_[count_-1].ns,"COST_FINAL_CLOCK");render(begin);
        size_t offset=0;while(offset<json_.size){const ssize_t n=write(fd_,json_.bytes.data()+offset,json_.size-offset);if(n<0&&errno==EINTR)continue;need(n>0,"COST_WRITE_FAILURE");offset+=size_t(n);}
        const int fd=fd_;fd_=-1;need(close(fd)==0,"COST_CLOSE_FAILURE");if(active==this)active=nullptr;
    }
};
static_assert(sizeof(Logger)<=65536,"Enabled logger object must remain bounded; not an allocator/RSS claim.");
}
