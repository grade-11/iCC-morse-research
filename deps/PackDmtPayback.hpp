#pragma once
#include "WindowState.hpp"
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// DESIGN-LOCK v1: metric-only predictions; this object never creates a PV.
struct PackDmtPayback {
    struct Model {
        double aL=0,aM=0,aU=0,bP=0,bF=0,bC=0;
        void load(const std::string& file) {
            std::ifstream in(file);std::string tag;
            if(!(in>>tag>>aL>>aM>>aU>>bP>>bF>>bC) || tag!="PACK_DMT_PAYBACK_V1")
                throw std::invalid_argument("Invalid packing cost model");
            for(double x:{aL,aM,aU,bP,bF,bC})if(!std::isfinite(x)||x<=0)
                throw std::invalid_argument("Missing/nonpositive packing cost coefficient");
            if(aL>aM || aM>aU)throw std::invalid_argument("Unordered native cost quantiles");
        }
    } model;
    struct Candidate {
        double p=0,q=0,Hq=0,Hend=0,delta=0,X=0,J=0;
        size_t rank=0,blocks=0,groups=0,labels_reduced=0;
        bool admissible=false;
    };
    std::vector<Candidate> candidates;
    double endpoint=0,total_s=0,chosen=0,H0E=0;
    bool full_endpoint=false,duplicates=false;
    std::string reason="NO_CANDIDATE";
    static double completeH(size_t n) {
        const double x=n;return x+x*(x-1)/2+x*(x-1)*(x-2)/2;
    }
    static bool fullEndpoint(const DistanceMatrix& d,double e,bool& duplicate) {
        bool full=true;duplicate=false;
        for(size_t i=0;i<d.vertexCount();++i)for(size_t j=i+1;j<d.vertexCount();++j){
            const double w=d.distance(i,j);if(w==0)duplicate=true;
            if(!(w>0 && std::isfinite(w) && w<e))full=false;
        }return full;
    }
    static double graphH(const DistanceMatrix& d,const std::vector<size_t>& owner,size_t count,double x) {
        std::vector<unsigned char> adjacent(count*count,0);std::vector<size_t> degree(count,0);size_t edges=0;
        for(size_t i=0;i<owner.size();++i)for(size_t j=i+1;j<owner.size();++j){
            const size_t a=owner[i],b=owner[j];if(a==b)continue;
            const double w=d.distance(i,j);
            // d=0 is excluded by the native solver. Duplicates are outside the exact contract.
            if(w>0 && w<x && !adjacent[a*count+b]){
                adjacent[a*count+b]=adjacent[b*count+a]=1;++degree[a];++degree[b];++edges;
            }
        }
        double h=count+edges;for(size_t k:degree)h+=double(k)*double(k-1)/2;return h;
    }
    static double actualH(const DistanceMatrix& d,const WindowState& state,
                          const std::vector<SelectedPV>& selected,double x,bool full) {
        if(full){
            size_t blocks=state.activeLabels().size();
            for(const auto& p:selected)blocks-=p.representatives.size()-1;
            return completeH(blocks);
        }
        const size_t n=state.originalVertexCount();std::vector<size_t> owner(n,n),absorbed(n,0);size_t count=0;
        for(const auto& p:selected)for(auto v:p.representatives)absorbed.at(v)=1;
        for(auto label:state.activeLabels()){
            if(label<n){if(!absorbed[label])owner[label]=count++;}
            else {for(auto v:state.pseudoVertices().at(label-n).representatives)owner.at(v)=count;++count;}
        }
        for(const auto& p:selected){for(auto v:p.representatives)owner.at(v)=count;++count;}
        if(std::find(owner.begin(),owner.end(),n)!=owner.end())throw std::logic_error("Incomplete partition");
        return full?completeH(count):graphH(d,owner,count,x);
    }
    void load(const DistanceMatrix& d,double e,const std::string& mode,const std::string& file) {
        const auto start=std::chrono::steady_clock::now();endpoint=e;chosen=e;
        if(mode=="plan")model.load(file);
        const size_t n=d.vertexCount();full_endpoint=fullEndpoint(d,e,duplicates);
        std::vector<double> values;std::vector<std::vector<size_t>> rows(n);
        for(size_t i=0;i<n;++i){
            rows[i].reserve(n-1);
            for(size_t j=0;j<n;++j)if(i!=j){
                rows[i].push_back(j);const double w=d.distance(i,j);
                if(!std::isfinite(w))throw std::invalid_argument("Nonfinite packing distance");
                if(i<j && w>0 && w<e)values.push_back(w);
            }
            std::sort(rows[i].begin(),rows[i].end(),[&](size_t a,size_t b){
                const double x=d.distance(i,a),y=d.distance(i,b);return x!=y?x<y:a<b;});
        }
        std::sort(values.begin(),values.end());values.erase(std::unique(values.begin(),values.end()),values.end());
        std::vector<size_t> singleton(n);std::iota(singleton.begin(),singleton.end(),0);
        H0E=full_endpoint?completeH(n):graphH(d,singleton,n,e);
        // Integer nearest rank avoids accidental p*m rounding at exact rank boundaries.
        for(size_t tenth=1;tenth<=8 && !values.empty();++tenth){
            const size_t rank=(tenth*values.size()+9)/10;
            const double q=std::nextafter(values[rank-1],std::numeric_limits<double>::infinity());
            if(q>=e || (!candidates.empty() && candidates.back().q==q))continue;
            Candidate c;c.p=double(tenth)/10;c.rank=rank;c.q=q;
            const double radius=.85*q/2;std::vector<size_t> counts(n),centers(n);std::iota(centers.begin(),centers.end(),0);
            for(size_t i=0;i<n;++i)counts[i]=std::upper_bound(rows[i].begin(),rows[i].end(),radius,
                [&](double v,size_t id){return v<d.distance(i,id);})-rows[i].begin();
            std::sort(centers.begin(),centers.end(),[&](size_t a,size_t b){return counts[a]!=counts[b]?counts[a]>counts[b]:a<b;});
            std::vector<size_t> owner(n,n);size_t block=0;
            for(auto center:centers){
                if(owner[center]!=n)continue;std::vector<size_t> members{center};
                for(size_t j:rows[center]){
                    if(d.distance(center,j)>radius)break;
                    if(owner[j]!=n)continue;bool valid=true;
                    for(size_t k:members){double w=d.distance(j,k);if(!(w<=.85*q && w<q)){valid=false;break;}}
                    if(valid)members.push_back(j);if(members.size()==63)break;
                }
                if(members.size()<3)continue;
                for(size_t j:members)owner[j]=block;++block;++c.groups;c.labels_reduced+=members.size()-1;
            }
            for(size_t j=0;j<n;++j)if(owner[j]==n)owner[j]=block++;
            c.blocks=block;c.Hq=graphH(d,singleton,n,q);
            c.Hend=full_endpoint?completeH(block):graphH(d,owner,block,e);c.delta=H0E-c.Hend;candidates.push_back(c);
        }
        double best=-std::numeric_limits<double>::infinity();
        for(auto& c:candidates){
            if(mode=="plan"){
                const double n2=double(n)*n;
                const double P=model.bP*n2*(1+std::log2(n+2)+candidates.size());
                const double F=model.bF*n2*std::log2(n+2),C=model.bC*n2;
                c.X=model.aU*c.Hq+P+F+C;c.J=model.aL*c.delta-1.25*c.X;
                c.admissible=c.X<=.50*model.aM*H0E && c.J>0;
            }else if(mode=="probe"){
                c.J=c.delta-c.Hq;c.admissible=true;
            }
            if(c.admissible && c.J>best){best=c.J;chosen=c.q;}
        }
        reason=candidates.empty()?"NO_CANDIDATE":chosen==e?"ABSTAIN/predicted-no-payback":mode=="probe"?"GRAPH_PROXY_DIAGNOSTIC":"predicted-payback";
        total_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    void report(const std::string& mode,bool details)const {
        std::cout<<std::setprecision(17)<<"[pack] mode="<<mode<<" candidates="<<candidates.size()
            <<" q="<<chosen<<" q_bits="<<std::bit_cast<uint64_t>(chosen)<<" E="<<endpoint
            <<" E_bits="<<std::bit_cast<uint64_t>(endpoint)<<" H0E="<<H0E<<" full_endpoint="<<full_endpoint
            <<" duplicates="<<duplicates<<" prep_s="<<total_s<<" reason="<<reason<<'\n';
        // Small fixed grid is also emitted in compact calibration logs.
        for(const auto& c:candidates)std::cout<<"[pack-candidate] p="<<c.p<<" rank="<<c.rank
            <<" q="<<c.q<<" q_bits="<<std::bit_cast<uint64_t>(c.q)<<" blocks="<<c.blocks
            <<" groups="<<c.groups<<" labels_reduced="<<c.labels_reduced<<" Hq="<<c.Hq<<" Hend="<<c.Hend
            <<" delta="<<c.delta<<" X="<<c.X<<" J="<<c.J<<" admissible="<<c.admissible<<'\n';
        std::cout<<std::setprecision(6);
    }
};
