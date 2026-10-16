#include "StreamingF2.hpp"
#include "PipelineCommon.hpp"
#include <algorithm>
#include <bit>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>

StreamingF2::StreamingF2(const DistanceMatrix& distances,double cutoff)
    :n_(distances.vertexCount()),cutoff_(cutoff),adjacency_(n_,std::vector<std::uint8_t>(n_)),
     edge_ids_(n_,std::vector<Id>(n_,std::numeric_limits<Id>::max()))
{
    births_.assign(n_,Birth{0,0.,true});
    for(Id b=1;b<n_;++b)for(Id a=0;a<b;++a)
    {
        const auto weight=distances.distance(a,b);
        if(weight<cutoff_)edges_.push_back({a,b,b*(b-1)/2+a,weight});
    }
    std::sort(edges_.begin(),edges_.end(),[](const auto& a,const auto& b)
    {return a.weight<b.weight || (a.weight==b.weight && a.bindex<b.bindex);});
}

void StreamingF2::reduce(std::vector<Id> column,std::size_t dimension,double weight,Id birth_id)
{
    std::sort(column.begin(),column.end());
    while(!column.empty())
    {
        const auto found=pivots_.find(column.back());
        if(found==pivots_.end())break;
        std::vector<Id> next;next.reserve(column.size()+found->second.size());
        std::set_symmetric_difference(column.begin(),column.end(),found->second.begin(),found->second.end(),std::back_inserter(next));
        column=std::move(next);
    }
    if(column.empty())
    {
        if(dimension==1)births_.at(birth_id).alive=true;
        return; // Top-dimensional births would be H2, explicitly unsupported.
    }
    const auto pivot=column.back();auto& birth=births_.at(pivot);
    if(!birth.alive || birth.dimension+1!=dimension)throw std::logic_error("Invalid streaming pivot/liveness");
    birth.alive=false;finite_.push_back({birth.dimension,birth.weight,weight});
    pivots_.emplace(pivot,std::move(column));
}

StreamingF2::Counts StreamingF2::advance(double epsilon)
{
    if(!(epsilon>=0.) || !(epsilon>epsilon_) || epsilon>cutoff_)throw std::invalid_argument("Streaming checkpoint must increase within cutoff");
    const auto old_edges=cursor_,old_triangles=triangle_count_;
    while(cursor_<edges_.size() && edges_[cursor_].weight<epsilon)
    {
        const auto weight=edges_[cursor_].weight;std::vector<Triangle> triangles;
        do
        {
            const auto& e=edges_[cursor_];const Id id=n_+cursor_;
            births_.push_back(Birth{1,weight,false});
            edge_ids_[e.a][e.b]=edge_ids_[e.b][e.a]=id;
            reduce({e.a,e.b},1,weight,id);
            for(Id c=0;c<n_;++c)if(adjacency_[e.a][c] && adjacency_[e.b][c])
            {
                std::array<Id,3> v{e.a,e.b,c};std::sort(v.begin(),v.end());
                const Id index=v[0]+v[1]*(v[1]-1)/2+v[2]*(v[2]-1)*(v[2]-2)/6;
                triangles.push_back({v,index});
            }
            adjacency_[e.a][e.b]=adjacency_[e.b][e.a]=1;++cursor_;
        }while(cursor_<edges_.size() && edges_[cursor_].weight==weight);
        std::sort(triangles.begin(),triangles.end(),[](const auto& a,const auto& b){return a.bindex<b.bindex;});
        for(const auto& triangle:triangles)
        {
            const auto& v=triangle.vertices;
            reduce({edge_ids_[v[0]][v[1]],edge_ids_[v[0]][v[2]],edge_ids_[v[1]][v[2]]},2,weight,0);
            ++triangle_count_;
        }
    }
    epsilon_=epsilon;
    return {cursor_-old_edges,triangle_count_-old_triangles,cursor_,triangle_count_,pivots_.size()};
}

std::vector<StreamingF2::Pair> StreamingF2::diagram(std::size_t dimension) const
{
    if(dimension>1)throw std::invalid_argument("StreamingF2 supports H0/H1 only");
    std::vector<Pair> result;
    for(const auto& p:finite_)if(p.dimension==dimension && (dimension==0 || p.birth!=p.death))result.push_back(p);
    for(const auto& b:births_)if(b.alive && b.dimension==dimension)result.push_back({dimension,b.weight,-1.});
    std::sort(result.begin(),result.end(),[](const auto& a,const auto& b){return a.birth<b.birth || (a.birth==b.birth && a.death<b.death);});
    return result;
}

std::uint64_t StreamingF2::stateHash() const
{
    std::uint64_t hash=1469598103934665603ULL;
    const auto add=[&](std::uint64_t x){hash^=x;hash*=1099511628211ULL;};
    add(cursor_);add(triangle_count_);std::vector<Id> keys;
    for(const auto& [key,column]:pivots_)keys.push_back(key);
    std::sort(keys.begin(),keys.end());
    for(auto key:keys){add(key);add(pivots_.at(key).size());for(auto value:pivots_.at(key))add(value);}
    for(const auto& b:births_){add(b.alive);add(b.dimension);add(std::bit_cast<std::uint64_t>(b.weight));}
    for(const auto& p:finite_){add(p.dimension);add(std::bit_cast<std::uint64_t>(p.birth));add(std::bit_cast<std::uint64_t>(p.death));}
    return hash;
}

void runStreamingF2(const DistanceMatrix& distances,const PipelineConfig& config)
{
    std::cout<<"[stream-engine] kernel=explicit-F2 max_homology=1 quotient=unsupported mode="<<config.continue_engine<<'\n';
    // In f2-rebuild, each checkpoint recomputes the same F2 engine from scratch.
    // In f2-stream, the one engine is retained. Both pay preprocessing and all
    // checkpoint/logging/teardown costs in whole-process time.
    const bool rebuild=config.continue_engine=="f2-rebuild";
    std::unique_ptr<StreamingF2> engine;
    double lo=0.;
    for(std::size_t i=0;i<config.eps_breaks.size();++i)
    {
        const auto hi=config.eps_breaks[i];
        if(!engine || rebuild)engine=std::make_unique<StreamingF2>(distances,config.finalEpsilon());
        const auto counts=engine->advance(hi);
        std::cout<<"window lower bound, window upper bound\n"<<lo<<", "<<hi<<"\ndimension, birth weight, death weight\n";
        for(const auto& p:engine->diagram(1))if(p.death<0 || p.death>=lo)
            std::cout<<"2, "<<std::setprecision(17)<<p.birth<<", "<<p.death<<'\n'<<std::setprecision(6);
        std::cout<<"[stream-state] eps="<<hi<<" mode="<<config.continue_engine
            <<" new_edges="<<counts.new_edges<<" new_triangles="<<counts.new_triangles
            <<" total_edges="<<counts.total_edges<<" total_triangles="<<counts.total_triangles<<" pivots="<<counts.pivots;
        if(config.verbose || config.cost_telemetry)std::cout<<" state_hash="<<engine->stateHash();
        std::cout<<'\n';
        if(config.verbose)
            for(const auto& p:engine->diagram(0))std::cout<<std::setprecision(17)<<"[stream-h0-pair] eps="<<hi<<" birth="<<p.birth<<" death="<<p.death<<'\n'<<std::setprecision(6);
        lo=hi;
    }
}
