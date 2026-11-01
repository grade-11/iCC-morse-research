#pragma once
#include "WindowState.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

// Conservative two-block triangle-connectivity forecast, NOT a PH certificate.
// Reads original distances and immutable candidate blocks. No trial WindowState.
struct H1EdgeRisk
{
    struct Block {std::vector<size_t> vertices;size_t label;bool proposed;};
    struct Event {double w;size_t a,b;};
    struct Warning {
        size_t lhs=0,rhs=0,components=0;
        double epsilon=0,first_weight=0,resolved_at=-1;
        std::vector<std::pair<size_t,size_t>> witness_edges;
    };
    struct Result {
        size_t pairs=0,edges=0,event_batches=0;
        double seconds=0;
        std::vector<Block> blocks;
        std::vector<Warning> warnings;
    };
    static Result inspect(const DistanceMatrix& d,const WindowState& state,
                          const std::vector<SelectedPV>& selected,double q,double endpoint)
    {
        const auto start=std::chrono::steady_clock::now();Result result;
        const size_t n=state.originalVertexCount();std::vector<bool> absorbed(n,false);
        for(const auto& pv:selected)for(auto v:pv.representatives)absorbed.at(v)=true;
        for(size_t label:state.activeLabels()){
            if(label<n){if(!absorbed[label])result.blocks.push_back({{label},label,false});}
            else result.blocks.push_back({state.pseudoVertices().at(label-n).representatives,label,false});
        }
        for(size_t p=0;p<selected.size();++p)
            result.blocks.push_back({selected[p].representatives,state.totalLabelCount()+p,true});
        for(size_t a=0;a<result.blocks.size();++a)for(size_t b=a+1;b<result.blocks.size();++b){
            const auto& A=result.blocks[a];const auto& B=result.blocks[b];
            if(!A.proposed && !B.proposed)continue;
            ++result.pairs;std::vector<Event> events;
            for(size_t i=0;i<A.vertices.size();++i)for(size_t j=0;j<B.vertices.size();++j){
                const double w=d.distance(A.vertices[i],B.vertices[j]);
                if(w>0 && w<endpoint)events.push_back({w,i,A.vertices.size()+j});
            }
            result.edges+=events.size();if(events.size()<2)continue;
            std::sort(events.begin(),events.end(),[](const Event& x,const Event& y){
                if(x.w!=y.w)return x.w<y.w;if(x.a!=y.a)return x.a<y.a;return x.b<y.b;});
            const size_t m=A.vertices.size()+B.vertices.size();
            std::vector<size_t> parent(m),size(m,1);std::vector<bool> active(m,false);
            std::iota(parent.begin(),parent.end(),0);size_t components=0;
            const auto find=[&](size_t v){while(parent[v]!=v){parent[v]=parent[parent[v]];v=parent[v];}return v;};
            const auto add=[&](const Event& e){
                for(auto v:{e.a,e.b})if(!active[v]){active[v]=true;++components;}
                size_t x=find(e.a),y=find(e.b);if(x!=y){if(size[x]<size[y])std::swap(x,y);parent[y]=x;size[x]+=size[y];--components;}
            };
            size_t pos=0;while(pos<events.size() && events[pos].w<q)add(events[pos++]);
            bool warned=false;Warning warning;
            const auto observe=[&](double epsilon,size_t end){
                if(components>1 && !warned){
                    warned=true;warning={a,b,components,epsilon,events.front().w,-1,{}};
                    std::vector<size_t> roots;
                    for(size_t e=0;e<end;++e){
                        const size_t root=find(events[e].a);
                        if(std::find(roots.begin(),roots.end(),root)!=roots.end())continue;
                        roots.push_back(root);warning.witness_edges.push_back({A.vertices[events[e].a],B.vertices[events[e].b-A.vertices.size()]});
                    }
                }else if(warned && components==1 && warning.resolved_at<0)warning.resolved_at=epsilon;
            };
            observe(q,pos);
            while(pos<events.size()){
                const double w=events[pos].w;
                do{add(events[pos++]);}while(pos<events.size() && events[pos].w==w);
                ++result.event_batches;
                const double epsilon=std::nextafter(w,std::numeric_limits<double>::infinity());
                // If epsilon equals E the event is still admitted by the shared E.
                observe(std::min(epsilon,endpoint),pos);
            }
            if(warned)result.warnings.push_back(std::move(warning));
        }
        result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        return result;
    }
    static std::vector<size_t> greedyDrops(const Result& r,const std::vector<SelectedPV>& proposed,size_t base) {
        std::vector<size_t> incident(proposed.size(),0),removed;
        std::vector<std::vector<size_t>> adjacency(proposed.size());
        for(size_t e=0;e<r.warnings.size();++e)for(auto block:{r.warnings[e].lhs,r.warnings[e].rhs})
            if(r.blocks[block].proposed){const size_t id=r.blocks[block].label-base;++incident.at(id);adjacency.at(id).push_back(e);}
        std::vector<bool> edge_live(r.warnings.size(),true),support_live(proposed.size(),true);size_t left=r.warnings.size();
        while(left){
            size_t drop=proposed.size();double best=std::numeric_limits<double>::infinity();
            for(size_t p=0;p<proposed.size();++p)if(support_live[p] && incident[p]){
                const double cost=double(proposed[p].representatives.size()-1)/incident[p];
                if(cost<best){best=cost;drop=p;}
            }
            if(drop==proposed.size())throw std::logic_error("Warning graph has no removable support");
            support_live[drop]=false;removed.push_back(drop);
            for(size_t e:adjacency[drop])if(edge_live[e]){
                edge_live[e]=false;--left;
                for(auto block:{r.warnings[e].lhs,r.warnings[e].rhs})if(r.blocks[block].proposed)
                    --incident.at(r.blocks[block].label-base);
            }
        }return removed;
    }
    static void report(const Result& r,double q,double endpoint,size_t protected_count,
                       size_t eligible,size_t committed,const std::string& mode,bool details)
    {
        const bool live=protected_count>0;
        std::cout<<std::setprecision(17)<<"[h1-risk] epsilon="<<q<<" horizon="<<endpoint
            <<" mode="<<mode<<" eligible="<<eligible<<" committed="<<committed
            <<" live_protected="<<protected_count<<" risky_pairs="<<r.warnings.size()
            <<" pairs_checked="<<r.pairs<<" edges_checked="<<r.edges<<" event_batches="<<r.event_batches
            <<" seconds="<<r.seconds<<" precommit_realization=NA/not-implemented"
            <<" reason="<<(eligible==0?"no-useful-support":live?"live-H1-history":!r.warnings.empty()?"disconnected-edge-fiber":"no-risk-detected")<<'\n';
        if(details)for(const auto& w:r.warnings){
            const auto& a=r.blocks[w.lhs];const auto& b=r.blocks[w.rhs];
            std::cout<<"[h1-risk-pair] quotient_epsilon="<<q<<" epsilon="<<w.epsilon
                <<" lhs="<<a.label<<" rhs="<<b.label<<" components="<<w.components
                <<" first_edge="<<w.first_weight<<" resolved_at="<<w.resolved_at<<" lhs_members=";
            for(auto v:a.vertices)std::cout<<v<<',';
            std::cout<<" rhs_members=";for(auto v:b.vertices)std::cout<<v<<',';
            std::cout<<" witness_edges=";for(auto [u,v]:w.witness_edges)std::cout<<u<<':'<<v<<',';std::cout<<'\n';
        }
        std::cout<<std::setprecision(6);
    }
};
