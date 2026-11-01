#pragma once
#include "PipelineCommon.hpp"
#include "WindowState.hpp"
#include "SimplexUtility.hpp"
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>

// Post-rebuild observation only. Uses the native recorded triangle witness and
// a deterministic minimum-distance original edge witness (ties: original colex).
// A zero count describes these canonical witnesses, never all realizations.
struct H1RealizationDiagnostics
{
    struct EdgeWitness {double weight;std::array<size_t,2> labels,vertices;};
    static void report(const PipelineRuntime& runtime,const WindowState& window,
                       const SimplexList& edges,const SimplexList& triangles,
                       const robin_hood::unordered_map<SimplexBindex,uint64_t>& triangle_witnesses)
    {
        const auto start=std::chrono::steady_clock::now();
        const auto old_precision=std::cout.precision();
        const auto& d=runtime.distanceMatrix();const auto& table=runtime.binomialTable();
        const size_t n=window.originalVertexCount();
        const auto reps=[&](size_t label){return label<n?std::vector<size_t>{label}:window.pseudoVertices().at(label-n).representatives;};
        std::map<SimplexBindex,EdgeWitness> canonical;
        size_t pv_edges=0;
        for(const auto& [index,weight]:edges){
            const auto labels=SimplexUtility::getSimplexVertices(table,index,window.totalLabelCount(),1);
            EdgeWitness e{std::numeric_limits<double>::infinity(),{labels[0],labels[1]},{0,0}};
            std::pair<size_t,size_t> best_order{SIZE_MAX,SIZE_MAX};
            for(size_t a:reps(labels[0]))for(size_t b:reps(labels[1])){
                const double w=d.distance(a,b);const std::pair<size_t,size_t> order{std::max(a,b),std::min(a,b)};
                if(w<e.weight || (w==e.weight && order<best_order)){e.weight=w;e.vertices={a,b};best_order=order;}
            }
            if(e.weight!=weight)throw std::logic_error("H1 canonical edge weight disagrees with native rebuilt edge");
            canonical.emplace(index,e);pv_edges+=labels[0]>=n;
        }
        size_t pv_triangles=0,incidences=0,different=0,multi_different=0,positive_gap=0,overshoot_count=0,emitted=0;
        double max_gap=0,max_overshoot=0;
        for(const auto& [index,weight]:triangles){
            const auto labels=SimplexUtility::getSimplexVertices(table,index,window.totalLabelCount(),2);
            if(labels[0]<n)continue;
            ++pv_triangles;
            const auto it=triangle_witnesses.find(index);
            if(it==triangle_witnesses.end())throw std::logic_error("Missing native H1 triangle witness");
            std::array<size_t,3> vertices;
            for(size_t i=0;i<3;++i){const auto group=reps(labels[i]);const size_t local=(it->second>>(8*i))&255;
                if(local>=group.size())throw std::logic_error("H1 triangle witness local index out of range");vertices[i]=group[local];}
            const double realized=std::max({d.distance(vertices[0],vertices[1]),d.distance(vertices[0],vertices[2]),d.distance(vertices[1],vertices[2])});
            if(realized!=weight)throw std::logic_error("Native H1 triangle realization weight mismatch");
            for(size_t i=0;i<3;++i)for(size_t j=i+1;j<3;++j){
                ++incidences;
                const auto edge_index=SimplexUtility::getEdgeBinomialIndex(table,labels[i],labels[j]);
                const auto& edge=canonical.at(edge_index);
                if(edge.labels[0]!=labels[i] || edge.labels[1]!=labels[j])throw std::logic_error("H1 label alignment mismatch");
                const size_t changed=(edge.vertices[0]!=vertices[i])+(edge.vertices[1]!=vertices[j]);
                different+=changed>0;multi_different+=changed==2;
                const double restricted=d.distance(vertices[i],vertices[j]),gap=restricted-edge.weight;
                if(gap<0)throw std::logic_error("Canonical minimum edge has negative realization gap");
                positive_gap+=gap>0;max_gap=std::max(max_gap,gap);
                const std::array<size_t,4> carrier{vertices[i],vertices[j],edge.vertices[0],edge.vertices[1]};
                double diameter=0;for(size_t a=0;a<4;++a)for(size_t b=a+1;b<4;++b)diameter=std::max(diameter,d.distance(carrier[a],carrier[b]));
                const double overshoot=std::max(0.,diameter-weight);overshoot_count+=overshoot>0;max_overshoot=std::max(max_overshoot,overshoot);
                if((gap>0 || changed || overshoot>0) && emitted<32){
                    ++emitted;std::cout<<std::setprecision(17)<<"[h1-realization-incidence] epsilon="<<window.bounds().eps_hi
                        <<" triangle="<<labels[0]<<','<<labels[1]<<','<<labels[2]<<" triangle_weight="<<weight
                        <<" triangle_vertices="<<vertices[0]<<','<<vertices[1]<<','<<vertices[2]
                        <<" edge="<<labels[i]<<','<<labels[j]<<" edge_weight="<<edge.weight
                        <<" edge_vertices="<<edge.vertices[0]<<','<<edge.vertices[1]<<" restricted_weight="<<restricted
                        <<" gap="<<gap<<" changed_coordinates="<<changed<<" carrier_overshoot="<<overshoot<<'\n';
                }
            }
        }
        std::cout<<std::setprecision(17)<<"[h1-realization] epsilon="<<window.bounds().eps_hi
            <<" status=SUPPORTED_CANONICAL_ONLY precommit=0 pv_edges="<<pv_edges<<" pv_triangles="<<pv_triangles
            <<" incidences="<<incidences<<" different="<<different<<" multi_different="<<multi_different
            <<" positive_gap="<<positive_gap<<" max_gap="<<max_gap<<" carrier_overshoot_count="<<overshoot_count
            <<" max_carrier_overshoot="<<max_overshoot<<" missing=0 emitted_examples="<<emitted
            <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
            <<" all_realizations=NA no_preservation_claim=1\n";
        std::cout.precision(old_precision);
    }
};
