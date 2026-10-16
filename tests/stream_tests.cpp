#include "StreamingF2.hpp"
#include "PipelineCommon.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

int failures=0;
#define EXPECT(value) do {if(!(value)){std::cerr<<__LINE__<<": " #value "\n";++failures;}}while(false)

int main()
{
    DistanceMatrix square({{0.,0.},{1.,0.},{1.,1.},{0.,1.}});
    StreamingF2 incremental(square,2.),oneshot(square,2.);
    auto a=incremental.advance(.9);EXPECT(a.total_edges==0);EXPECT(incremental.diagram(0).size()==4);
    auto b=incremental.advance(1.1);EXPECT(b.new_edges==4 && b.new_triangles==0);
    const std::vector<StreamingF2::Pair> live{{1,1.,-1.}};EXPECT(incremental.diagram(1)==live);
    auto c=incremental.advance(1.4);EXPECT(c.new_edges==0 && c.new_triangles==0);EXPECT(incremental.diagram(1)==live);
    auto d=incremental.advance(1.5);EXPECT(d.new_edges==2 && d.new_triangles==4);
    const std::vector<StreamingF2::Pair> closed{{1,1.,std::sqrt(2.)}};EXPECT(incremental.diagram(1)==closed);
    oneshot.advance(1.5);EXPECT(incremental.stateHash()==oneshot.stateHash());
    EXPECT(incremental.diagram(0)==oneshot.diagram(0));EXPECT(incremental.diagram(1)==oneshot.diagram(1));
    bool rejected=false;try{incremental.advance(1.5);}catch(const std::invalid_argument&){rejected=true;}EXPECT(rejected);
    rejected=false;try{StreamingF2 fresh(square,2.);fresh.advance(-.1);}catch(const std::invalid_argument&){rejected=true;}EXPECT(rejected);
    DistanceMatrix duplicates({{0.,0.},{0.,0.},{1.,0.},{0.,1.}});StreamingF2 duplicate_stream(duplicates,2.);
    duplicate_stream.advance(2.);EXPECT(duplicate_stream.diagram(1).empty());
    int zero=0,survivors=0;for(const auto& p:duplicate_stream.diagram(0)){zero+=p.death==0;survivors+=p.death==-1;}
    EXPECT(zero==1 && survivors==1); // Explicit H0 zero bars retained, unlike Ripser's output policy.
    DistanceMatrix boundary({{0.,0.},{2.,0.}});StreamingF2 strict(boundary,2.);strict.advance(2.);
    EXPECT(strict.diagram(0).size()==2);for(const auto& p:strict.diagram(0))EXPECT(p.death==-1);
    PipelineConfig config;config.maxdim=2;config.threads=4;config.pv_cap_scale=.85;config.eps_breaks={1.,2.};config.quotient_mask={0};config.continue_engine="f2-stream";
    config.validate();config.quotient_mask={1};rejected=false;try{config.validate();}catch(const std::invalid_argument&){rejected=true;}EXPECT(rejected);
    config.quotient_mask={0};config.maxdim=3;rejected=false;try{config.validate();}catch(const std::invalid_argument&){rejected=true;}EXPECT(rejected);
    if(!failures)std::cout<<"Streaming F2: old live birth closes, ties/duplicates/strict cutoff, exact state prefix, refusal scope PASS\n";
    return failures?1:0;
}
