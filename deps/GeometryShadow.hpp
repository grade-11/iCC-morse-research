#pragma once
// Opt-in research observer. Never changes schedule, supports, guards or state.
#include "PipelineCommon.hpp"
#include "WindowState.hpp"
#include "PolicyCostLog.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fcntl.h>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace GeoShadow
{
constexpr size_t N_CAP=1024, LINE_CAP=1048576, FILE_CAP=3145728;
static_assert(sizeof(double)==8 && sizeof(uint64_t)==8 && std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53);
inline void need(bool good,const char* reason){if(!good)throw std::runtime_error(reason);}
inline uint64_t nowNs(){timespec t{};need(clock_gettime(CLOCK_MONOTONIC,&t)==0&&t.tv_sec>=0&&t.tv_nsec>=0&&t.tv_nsec<1000000000,"GEO_CLOCK_FAILURE");return uint64_t(t.tv_sec)*1000000000ULL+uint64_t(t.tv_nsec);}
inline uint64_t choose2(uint64_t n){return n<2?0:n*(n-1)/2;}
inline std::string bits(double x)
{
    constexpr char hex[]="0123456789abcdef";uint64_t word=std::bit_cast<uint64_t>(x);
    std::string text(16,'0');for(size_t i=0;i<16;++i){text[15-i]=hex[word&15];word>>=4;}return text;
}
inline std::string hexWord(uint64_t word){constexpr char hex[]="0123456789abcdef";std::string text(16,'0');for(size_t i=0;i<16;++i){text[15-i]=hex[word&15];word>>=4;}return text;}

struct Json
{
    std::string text;
    Json(){text.reserve(LINE_CAP);}
    void add(std::string_view s){need(s.size()<=LINE_CAP-text.size(),"GEO_LINE_CAP");text.append(s);}
    void number(uint64_t x){char b[32];auto r=std::to_chars(b,b+sizeof(b),x);need(r.ec==std::errc{},"GEO_INTEGER_FORMAT");add(std::string_view(b,size_t(r.ptr-b)));}
    void boolean(bool x){add(x?"true":"false");}
    void quote(std::string_view s)
    {
        add("\"");for(unsigned char c:s){need(c>=32 && c<127,"GEO_NON_ASCII_TEXT");if(c=='"'||c=='\\')add("\\");char b=char(c);add(std::string_view(&b,1));}add("\"");
    }
    void key(std::string_view s){quote(s);add(":");}
    void nullableWord(bool present,double x){if(present)quote(bits(x));else add("null");}
};

struct DSU
{
    std::vector<size_t> parent;
    explicit DSU(size_t n):parent(n){for(size_t i=0;i<n;++i)parent[i]=i;}
    size_t find(size_t a){while(parent[a]!=a){parent[a]=parent[parent[a]];a=parent[a];}return a;}
    void unite(size_t a,size_t b){a=find(a);b=find(b);if(a>b)std::swap(a,b);parent[b]=a;}
    std::vector<size_t> labels(){std::vector<size_t> out(parent.size());for(size_t i=0;i<out.size();++i)out[i]=find(i);return out;}
};
struct Component{size_t label=0,size=0;};
inline std::vector<Component> components(const std::vector<size_t>& labels)
{
    std::vector<size_t> sizes(labels.size(),0);for(auto label:labels)++sizes.at(label);
    std::vector<Component> out;for(size_t i=0;i<sizes.size();++i)if(sizes[i])out.push_back({i,sizes[i]});return out;
}
struct Alias
{
    std::string_view name,status="UNAVAILABLE_EMPTY";size_t rank=0;
    bool selected=false,available=false;double distance=0,successor=0;
    size_t multiplicity=0,below=0;
};
struct Profile{size_t label=0,size=0;uint64_t within=0,cross=0;};
struct Candidate
{
    double q=0;std::vector<std::string_view> aliases;std::vector<size_t> degrees,labels;
    std::vector<Component> comps,histogram;std::vector<Profile> regions;
    uint64_t edges=0,within=0,cross=0,within_pairs=0,cross_pairs=0,u2numerator=0;
};
struct Features
{
    size_t n=0,k=0,positive=0,distinct=0;double E=0;
    std::array<Alias,3> aliases{{{"degree32"},{"q50"},{"q70"}}};
    std::vector<Candidate> candidates;std::vector<std::vector<size_t>> knn;
    std::vector<size_t> region_labels;std::vector<Component> regions;
    uint64_t mutual_edges=0,pair_reads=0,rank_ns=0,knn_ns=0,graphs_ns=0,total_ns=0;
};

template<class Metric> Features compute(const Metric& metric,double E)
{
    const uint64_t begin=nowNs();Features f;f.n=metric.vertexCount();f.E=E;
    need(f.n>=2&&f.n<=N_CAP,"N_OUT_OF_DOMAIN");need(std::isfinite(E)&&E>0,"E_OUT_OF_DOMAIN");f.k=std::min<size_t>(8,f.n-1);
    const auto read=[&](size_t i,size_t j){++f.pair_reads;const double d=metric.distance(i,j);need(std::isfinite(d),"NONFINITE_DISTANCE");need(d>=0,"NEGATIVE_DISTANCE");return d;};
    const uint64_t rank_start=nowNs();
    {
        std::vector<double> values;values.reserve(choose2(f.n));
        for(size_t i=1;i<f.n;++i)for(size_t j=0;j<i;++j){const double d=read(i,j);if(d>0&&d<E)values.push_back(d);}
        std::sort(values.begin(),values.end());f.positive=values.size();
        std::vector<double> distinct;distinct.reserve(values.size());for(double x:values)if(distinct.empty()||distinct.back()!=x)distinct.push_back(x);f.distinct=distinct.size();
        f.aliases[0].rank=16*f.n;
        if(!distinct.empty()){f.aliases[1].rank=size_t(std::ceil(.5*double(distinct.size())));f.aliases[2].rank=size_t(std::ceil(.70*double(distinct.size())));}
        for(size_t a=0;a<3;++a)
        {
            auto& alias=f.aliases[a];const auto& ranked=a==0?values:distinct;
            if(alias.rank==0){alias.status="UNAVAILABLE_EMPTY";continue;}
            if(alias.rank>ranked.size()){alias.status="UNAVAILABLE_RANK";continue;}
            alias.selected=true;alias.distance=ranked[alias.rank-1];alias.successor=std::min(E,std::nextafter(alias.distance,std::numeric_limits<double>::infinity()));
            const auto equal=std::equal_range(values.begin(),values.end(),alias.distance);alias.below=size_t(equal.first-values.begin());alias.multiplicity=size_t(equal.second-equal.first);
            alias.available=alias.successor<E;alias.status=alias.available?"AVAILABLE":"UNAVAILABLE_AT_E";
            if(alias.available)
            {
                auto found=std::find_if(f.candidates.begin(),f.candidates.end(),[&](const auto& c){return std::bit_cast<uint64_t>(c.q)==std::bit_cast<uint64_t>(alias.successor);});
                if(found==f.candidates.end()){Candidate c;c.q=alias.successor;f.candidates.push_back(std::move(c));found=f.candidates.end()-1;}
                found->aliases.push_back(alias.name);
            }
        }
        std::sort(f.candidates.begin(),f.candidates.end(),[](const auto& a,const auto& b){return a.q<b.q;});
    }
    f.rank_ns=nowNs()-rank_start;
    const uint64_t knn_start=nowNs();
    struct Neighbor{double d;size_t id;};
    const auto less=[](const Neighbor&a,const Neighbor&b){return a.d<b.d||(a.d==b.d&&a.id<b.id);};
    f.knn.resize(f.n);
    for(size_t i=0;i<f.n;++i)
    {
        std::vector<Neighbor> heap;heap.reserve(f.k);
        for(size_t j=0;j<f.n;++j)if(i!=j)
        {
            Neighbor candidate{read(i,j),j};
            if(heap.size()<f.k){heap.push_back(candidate);std::push_heap(heap.begin(),heap.end(),less);}
            else if(less(candidate,heap.front())){std::pop_heap(heap.begin(),heap.end(),less);heap.back()=candidate;std::push_heap(heap.begin(),heap.end(),less);}
        }
        std::sort_heap(heap.begin(),heap.end(),less);f.knn[i].reserve(f.k);for(const auto& neighbor:heap)f.knn[i].push_back(neighbor.id);
    }
    DSU mutual(f.n);
    for(size_t i=0;i<f.n;++i)for(size_t j:f.knn[i])if(i<j&&std::find(f.knn[j].begin(),f.knn[j].end(),i)!=f.knn[j].end()){mutual.unite(i,j);++f.mutual_edges;}
    f.region_labels=mutual.labels();f.regions=components(f.region_labels);f.knn_ns=nowNs()-knn_start;
    const uint64_t graph_start=nowNs();
    for(auto& c:f.candidates)
    {
        c.degrees.assign(f.n,0);DSU connected(f.n);std::vector<uint64_t> within(f.n,0),cross(f.n,0);
        for(size_t i=1;i<f.n;++i)for(size_t j=0;j<i;++j)
        {
            const double d=read(i,j);if(!(d>0&&d<c.q))continue;
            ++c.edges;++c.degrees[i];++c.degrees[j];connected.unite(i,j);
            const auto a=f.region_labels[i],b=f.region_labels[j];
            if(a==b){++c.within;++within[a];}else{++c.cross;++cross[a];++cross[b];}
        }
        c.labels=connected.labels();c.comps=components(c.labels);std::vector<size_t> histogram(f.n,0);
        for(auto degree:c.degrees){++histogram[degree];c.u2numerator+=choose2(degree);}
        for(size_t d=0;d<f.n;++d)if(histogram[d])c.histogram.push_back({d,histogram[d]});
        for(const auto& region:f.regions){c.regions.push_back({region.label,region.size,within[region.label],cross[region.label]});c.within_pairs+=choose2(region.size);}
        c.cross_pairs=choose2(f.n)-c.within_pairs;
    }
    f.graphs_ns=nowNs()-graph_start;f.total_ns=nowNs()-begin;
    need(f.pair_reads==choose2(f.n)*(3+f.candidates.size()),"GEO_ACCESSOR_COUNT_INVARIANT");return f;
}

inline void numbers(Json& j,const std::vector<size_t>& values){j.add("[");for(size_t i=0;i<values.size();++i){if(i)j.add(",");j.number(values[i]);}j.add("]");}
inline void componentJson(Json& j,const std::vector<Component>& values,bool histogram=false)
{
    j.add("[");for(size_t i=0;i<values.size();++i){if(i)j.add(",");j.add(histogram?"{\"degree\":" : "{\"label\":");j.number(values[i].label);j.add(histogram?",\"count\":" : ",\"size\":");j.number(values[i].size);j.add("}");}j.add("]");
}
inline void density(Json& j,uint64_t edges,uint64_t pairs){j.add("{\"numerator\":");j.number(edges);j.add(",\"denominator\":");j.number(pairs);j.add(",\"status\":");j.quote(pairs?"DEFINED":"NA");j.add("}");}
inline void profile(Json&j,uint64_t edges,uint64_t pairs){j.add("{\"edge_count\":");j.number(edges);j.add(",\"pair_count\":");j.number(pairs);j.add(",\"density\":");density(j,edges,pairs);j.add("}");}
inline void geometryJson(Json&j,const Features&f)
{
    j.add("{\"n\":");j.number(f.n);j.add(",\"E_bits\":");j.quote(bits(f.E));j.add(",\"k\":");j.number(f.k);
    j.add(",\"positive_pair_count\":");j.number(f.positive);j.add(",\"distinct_positive_count\":");j.number(f.distinct);j.add(",\"rank_aliases\":[");
    for(size_t i=0;i<3;++i){const auto&a=f.aliases[i];if(i)j.add(",");j.add("{\"name\":");j.quote(a.name);j.add(",\"status\":");j.quote(a.status);j.add(",\"rank_one_based\":");if(a.rank)j.number(a.rank);else j.add("null");j.add(",\"selected_distance_bits\":");j.nullableWord(a.selected,a.distance);j.add(",\"successor_bits\":");j.nullableWord(a.selected,a.successor);j.add(",\"candidate_q_bits\":");j.nullableWord(a.available,a.successor);j.add(",\"selected_multiplicity\":");if(a.selected)j.number(a.multiplicity);else j.add("null");j.add(",\"below_selected\":");if(a.selected)j.number(a.below);else j.add("null");j.add("}");}
    j.add("],\"candidates\":[");
    for(size_t x=0;x<f.candidates.size();++x)
    {
        const auto&c=f.candidates[x];if(x)j.add(",");j.add("{\"q_bits\":");j.quote(bits(c.q));j.add(",\"aliases\":[");for(size_t a=0;a<c.aliases.size();++a){if(a)j.add(",");j.quote(c.aliases[a]);}j.add("],\"degrees\":");numbers(j,c.degrees);
        j.add(",\"degree_histogram\":");componentJson(j,c.histogram,true);j.add(",\"component_label\":");numbers(j,c.labels);j.add(",\"components\":");componentJson(j,c.comps);j.add(",\"edge_count\":");j.number(c.edges);j.add(",\"density\":");density(j,c.edges,choose2(f.n));
        j.add(",\"within_region\":");profile(j,c.within,c.within_pairs);j.add(",\"cross_region\":");profile(j,c.cross,c.cross_pairs);j.add(",\"per_region\":[");
        for(size_t r=0;r<c.regions.size();++r){const auto&p=c.regions[r];if(r)j.add(",");j.add("{\"label\":");j.number(p.label);j.add(",\"size\":");j.number(p.size);j.add(",\"within_edges\":");j.number(p.within);j.add(",\"within_pairs\":");j.number(choose2(p.size));j.add(",\"cross_edges\":");j.number(p.cross);j.add(",\"cross_pairs\":");j.number(p.size*(f.n-p.size));j.add("}");}
        j.add("],\"U2_numerator\":");j.number(c.u2numerator);j.add(",\"U2_divisor\":3,\"U2\":");j.number(c.u2numerator/3);j.add("}");
    }
    j.add("],\"knn_ids\":[");for(size_t i=0;i<f.knn.size();++i){if(i)j.add(",");numbers(j,f.knn[i]);}j.add("],\"mutual_edge_count\":");j.number(f.mutual_edges);j.add(",\"mutual_component_label\":");numbers(j,f.region_labels);j.add(",\"regions\":");componentJson(j,f.regions);
    j.add(",\"timing_ns\":{\"rank\":");j.number(f.rank_ns);j.add(",\"knn\":");j.number(f.knn_ns);j.add(",\"graphs\":");j.number(f.graphs_ns);j.add(",\"prepass_total\":");j.number(f.total_ns);j.add("}}");
}

class Logger;
inline Logger* active=nullptr;
class Logger
{
    int fd_=-1;size_t ordinal_=0,n_=0,windows_=0;double E_=0,first_q_=0;
    bool first_=false,selector_=false,commit_record_=false,result_record_=false;
    uint64_t bytes_=0,write_ns_=0,pair_reads_=0,prepass_ns_=0;
    size_t A_before_=0,P_before_=0,L_before_=0;
    std::string run_,source_,binary_,input_;uint64_t wrapper_=0;
    Json j_;
    static std::string env(const char*name){const char*v=std::getenv(name);need(v&&*v,"GEO_REQUIRED_METADATA_ENV");return v;}
    static bool token(std::string_view s){return !s.empty()&&s.size()<=64&&std::all_of(s.begin(),s.end(),[](unsigned char c){return (c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='.'||c=='_'||c=='-';});}
    static bool sha(std::string_view s){return s.size()==64&&std::all_of(s.begin(),s.end(),[](char c){return(c>='0'&&c<='9')||(c>='a'&&c<='f');});}
    void begin(std::string_view type)
    {
        j_.text.clear();j_.add("{\"schema\":\"GEO_SHADOW_JSONL_V1\",\"record_type\":");j_.quote(type);j_.add(",\"ordinal\":");j_.number(ordinal_);j_.add(",\"run_id\":");j_.quote(run_);j_.add(",\"native_pid\":");j_.number(uint64_t(getpid()));j_.add(",\"monotonic_ns\":");j_.number(nowNs());j_.add(",\"payload\":");
    }
    void emit()
    {
        j_.add("}\n");need(j_.text.size()<=FILE_CAP-bytes_,"GEO_FILE_CAP");const uint64_t start=nowNs();size_t offset=0;
        while(offset<j_.text.size()){const ssize_t wrote=write(fd_,j_.text.data()+offset,j_.text.size()-offset);if(wrote<0&&errno==EINTR)continue;need(wrote>0,"GEO_WRITE_FAILURE");offset+=size_t(wrote);}write_ns_+=nowNs()-start;bytes_+=j_.text.size();++ordinal_;
    }
    void missingCommit(std::string_view reason,bool q)
    {
        begin("COMMIT_SHADOW");j_.add("{\"observation\":\"NOT_OBSERVED\",\"reason\":");j_.quote(reason);j_.add(",\"q_bits\":");j_.nullableWord(q,first_q_);
        j_.add(",\"raw_support_count\":null,\"protected_count\":null,\"eligible_count\":null,\"retained_count\":null,\"retained_hash_bits\":null,\"remaining_warnings\":null,\"native_filter_result\":\"NOT_OBSERVED\",\"state_before\":null,\"full_S\":null,\"proposal\":null,\"forecast_status\":\"NOT_REACHED\",\"shadow_commit\":\"NOT_REACHED\",\"lb_remaining_DEFER_ns\":null,\"ub_remaining_Q_ns\":null,\"AS_forecast_bytes\":null,\"RSS_forecast_bytes\":null,\"native_action_modified\":false}");emit();commit_record_=true;
    }
public:
    explicit Logger(const char*path)
    {
        need(path&&*path&&std::char_traits<char>::length(path)<=4096,"GEO_OUTPUT_PATH");
        run_=env("GEO_SHADOW_RUN_ID");source_=env("GEO_SHADOW_SOURCE_SHA256");binary_=env("GEO_SHADOW_BINARY_SHA256");input_=env("GEO_SHADOW_INPUT_SHA256");const auto pid=env("GEO_SHADOW_WRAPPER_PARENT_PID");
        need(token(run_)&&sha(source_)&&sha(binary_)&&sha(input_),"GEO_METADATA_FORMAT");auto parsed=std::from_chars(pid.data(),pid.data()+pid.size(),wrapper_);need(parsed.ec==std::errc{}&&parsed.ptr==pid.data()+pid.size()&&wrapper_>0&&wrapper_!=uint64_t(getpid()),"GEO_WRAPPER_PID_FORMAT");
        fd_=open(path,O_WRONLY|O_CREAT|O_EXCL|O_CLOEXEC,0600);need(fd_>=0,"GEO_EXCLUSIVE_OPEN_FAILURE");
    }
    ~Logger()noexcept{if(active==this)active=nullptr;if(fd_>=0)close(fd_);}
    Logger(const Logger&)=delete;Logger&operator=(const Logger&)=delete;
    static std::unique_ptr<Logger> create(const DistanceMatrix&metric,const PipelineConfig&config)
    {
        const char*path=std::getenv("GEO_SHADOW_FILE");if(!path)return {};
        need(active==nullptr,"GEO_NESTED_LOGGER");need(config.maxdim==2&&config.h1_risk=="selective-fast"&&config.pv_cap_scale==.85&&config.pv_min_separation==0&&config.pack_mode.empty(),"GEO_SHADOW_NATIVE_DOMAIN");
        auto log=std::make_unique<Logger>(path);log->initialize(metric,config.finalEpsilon());active=log.get();return log;
    }
    void initialize(const DistanceMatrix&metric,double E)
    {
        n_=metric.vertexCount();E_=E;need(n_>=2&&n_<=N_CAP,"N_OUT_OF_DOMAIN");need(std::isfinite(E_)&&E_>0,"E_OUT_OF_DOMAIN");
        begin("HEADER");j_.add("{\"mode\":\"SHADOW_ONLY\",\"model_status\":\"UNFIT\",\"self_verified\":false,\"declared_source_sha256\":");j_.quote(source_);j_.add(",\"declared_binary_sha256\":");j_.quote(binary_);j_.add(",\"declared_input_sha256\":");j_.quote(input_);j_.add(",\"wrapper_parent_pid\":");j_.number(wrapper_);j_.add(",\"n\":");j_.number(n_);j_.add(",\"E_bits\":");j_.quote(bits(E_));j_.add(",\"max_simplex_dimension\":2,\"k_cap\":8,\"n_cap\":1024,\"line_cap_bytes\":1048576,\"file_cap_bytes\":3145728}");emit();
        const auto f=compute(metric,E_);pair_reads_=f.pair_reads;prepass_ns_=f.total_ns;
        if(auto* cost=PolicyCost::active)cost->opening(f);
        begin("GEOMETRY");geometryJson(j_,f);emit();begin("OPENING");j_.add("{\"model_status\":\"UNFIT\",\"forecast_status\":\"UNKNOWN\",\"shadow_action\":\"UNKNOWN\",\"reference\":\"ONE_WINDOW_NOQ\",\"candidate_q_bits\":[");for(size_t i=0;i<f.candidates.size();++i){if(i)j_.add(",");j_.quote(bits(f.candidates[i].q));}j_.add("],\"chosen_q_bits\":null,\"lb_NOQ_ns\":null,\"ub_GEO_ns\":null,\"score_ns\":null,\"AS_forecast_bytes\":null,\"RSS_forecast_bytes\":null,\"native_action_modified\":false}");emit();
    }
    void noteWindow(){++windows_;}
    void firstDecision(double q,bool selector,size_t A,size_t P,size_t L)
    {
        if(first_)return;first_=true;selector_=selector;first_q_=q;A_before_=A;P_before_=P;L_before_=L;
        if(!selector)missingCommit("BASELINE_BYPASSED_SELECTOR",true);
    }
    bool awaitingCommit()const{return first_&&selector_&&!commit_record_;}
    void retained(double q,size_t raw,size_t protected_count,size_t eligible,size_t warnings,bool skipped_live,uint64_t hash,const std::vector<SelectedPV>&selected)
    {
        if(!awaitingCommit())return;need(bits(q)==bits(first_q_)&&P_before_==0&&A_before_==n_&&L_before_==n_,"GEO_FIRST_COMMIT_CONTEXT");
        size_t gain=0;std::vector<bool> seen(n_,false);for(const auto&pv:selected){need(!pv.representatives.empty()&&pv.representatives.size()<=n_,"GEO_EMPTY_OR_OVERSIZED_RETAINED_SUPPORT");gain+=pv.representatives.size()-1;for(auto id:pv.representatives){need(id<n_&&!seen[id],"GEO_RETAINED_PARTITION_OVERLAP");seen[id]=true;}}
        const size_t P=selected.size(),A=n_-gain,L=n_+P;need(A>=P,"GEO_PROPOSAL_LABEL_CONTEXT");
        begin("COMMIT_SHADOW");j_.add("{\"observation\":\"OBSERVED\",\"reason\":\"FULL_RETAINED_BEFORE_FORCE_DEFER\",\"q_bits\":");j_.quote(bits(q));
        j_.add(",\"raw_support_count\":");j_.number(raw);j_.add(",\"protected_count\":");j_.number(protected_count);j_.add(",\"eligible_count\":");j_.number(eligible);j_.add(",\"retained_count\":");j_.number(P);j_.add(",\"retained_hash_bits\":");j_.quote(hexWord(hash));j_.add(",\"remaining_warnings\":");if(skipped_live)j_.add("null");else j_.number(warnings);j_.add(",\"native_filter_result\":");j_.quote(P?"RETAINED_NONEMPTY":"EMPTY_AFTER_NATIVE_FILTERS");j_.add(",\"state_before\":{\"A\":");j_.number(A_before_);j_.add(",\"P\":");j_.number(P_before_);j_.add(",\"L\":");j_.number(L_before_);j_.add("},\"full_S\":[");
        for(size_t p=0;p<P;++p){if(p)j_.add(",");numbers(j_,selected[p].representatives);}j_.add("],\"proposal\":{\"P_prime\":");j_.number(P);j_.add(",\"gain\":");j_.number(gain);j_.add(",\"A_prime\":");j_.number(A);j_.add(",\"L_prime\":");j_.number(L);j_.add(",\"active_required_reserve\":");j_.number(A_before_+P);j_.add(",\"PV_hash_unique_pairs\":");j_.number(P*(A-P)+choose2(P));j_.add("},\"forecast_status\":\"UNKNOWN\",\"shadow_commit\":\"UNKNOWN\",\"lb_remaining_DEFER_ns\":null,\"ub_remaining_Q_ns\":null,\"AS_forecast_bytes\":null,\"RSS_forecast_bytes\":null,\"native_action_modified\":false}");emit();commit_record_=true;
    }
    void baselineResult(size_t committed,size_t A,size_t P,size_t L)
    {
        if(result_record_)return;need(first_&&commit_record_,"GEO_BASELINE_RESULT_ORDER");begin("BASELINE_RESULT");j_.add("{\"observed\":true,\"reason\":\"FIRST_NONFINAL_RESULT\",\"q_bits\":");j_.quote(bits(first_q_));j_.add(",\"selector_called\":");j_.boolean(selector_);j_.add(",\"native_do_quotient_before_postprocessor\":");j_.boolean(selector_);j_.add(",\"actual_committed_count\":");j_.number(committed);j_.add(",\"A_before\":");j_.number(A_before_);j_.add(",\"P_before\":");j_.number(P_before_);j_.add(",\"L_before\":");j_.number(L_before_);j_.add(",\"A_after\":");j_.number(A);j_.add(",\"P_after\":");j_.number(P);j_.add(",\"L_after\":");j_.number(L);j_.add(",\"native_schedule_modified\":false}");emit();result_record_=true;
    }
    void finish()
    {
        if(!first_){missingCommit("NO_NONFINAL_WINDOW",false);begin("BASELINE_RESULT");j_.add("{\"observed\":false,\"reason\":\"NO_NONFINAL_WINDOW\",\"q_bits\":null,\"selector_called\":null,\"native_do_quotient_before_postprocessor\":null,\"actual_committed_count\":null,\"A_before\":null,\"P_before\":null,\"L_before\":null,\"A_after\":null,\"P_after\":null,\"L_after\":null,\"native_schedule_modified\":false}");emit();result_record_=true;}
        need(commit_record_&&result_record_&&ordinal_==5,"GEO_FOOTER_ORDER");const uint64_t prior_bytes=bytes_,prior_write=write_ns_;begin("FOOTER");j_.add("{\"status\":\"COMPLETE\",\"record_count\":6,\"first_nonfinal_seen\":");j_.boolean(first_);j_.add(",\"commit_snapshot_observed\":");j_.boolean(selector_&&first_);j_.add(",\"baseline_window_count\":");j_.number(windows_);j_.add(",\"geometry_pair_reads\":");j_.number(pair_reads_);j_.add(",\"prepass_total_ns\":");j_.number(prepass_ns_);j_.add(",\"observer_write_ns_before_footer\":");j_.number(prior_write);j_.add(",\"bytes_before_footer\":");j_.number(prior_bytes);j_.add(",\"native_actions_changed\":false,\"model_fitted\":false}");emit();const int fd=fd_;fd_=-1;need(close(fd)==0,"GEO_CLOSE_FAILURE");if(active==this)active=nullptr;
    }
};
}
