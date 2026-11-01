#pragma once
#include "DistanceMatrix.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>

// Geometry proposes times only. The full VR metric and DMT remain authoritative.
struct GeometryController
{
    using Clock = std::chrono::steady_clock;
    struct UnionFind {
        std::vector<size_t> p, size;
        explicit UnionFind(size_t n):p(n),size(n,1){std::iota(p.begin(),p.end(),0);}
        size_t find(size_t x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}
        void join(size_t a,size_t b){a=find(a);b=find(b);if(a==b)return;if(size[a]<size[b])std::swap(a,b);p[b]=a;size[a]+=size[b];}
    };
    struct Edge {size_t a,b;double w;};
    struct Cluster {
        std::vector<size_t> vertices;
        double diameter=0,density=0,outside=std::numeric_limits<double>::infinity();
        size_t age=1,growth=0;
    };
    struct Snapshot {double epsilon=0;std::vector<Cluster> clusters;};
    struct Candidate {double epsilon=0,gain=0,stability=0,growth=0,work=0;};
    std::vector<Snapshot> snapshots;
    std::vector<Candidate> candidates;
    std::vector<double> schedule;
    double knn_s=0,cluster_s=0,selection_s=0,total_s=0;
    size_t mutual_edges=0,requested_cuts=0;
    bool fallback=false;

    static std::vector<double> rankedSubset(const std::vector<double>& values,size_t cap)
    {
        if(values.size()<=cap)return values;
        std::vector<double> out;
        for(size_t j=1;j<=cap;++j)out.push_back(values[(j*values.size()+cap-1)/cap-1]);
        return out;
    }
    static double wedges(const DistanceMatrix& d,double epsilon)
    {
        std::vector<size_t> degree(d.vertexCount(),0);
        for(size_t a=1;a<degree.size();++a)for(size_t b=0;b<a;++b)
            if(d.distance(a,b)>0 && d.distance(a,b)<epsilon){++degree[a];++degree[b];}
        double w=0;for(auto v:degree)w+=double(v)*double(v?v-1:0)/2;return w;
    }
    void load(const DistanceMatrix& d,double endpoint,double alpha,size_t cuts)
    {
        const auto start=Clock::now();requested_cuts=cuts;
        const size_t n=d.vertexCount(),k=std::min<size_t>(8,n?n-1:0);
        std::vector<std::vector<size_t>> neighbors(n);
        for(size_t a=0;a<n;++a){
            std::vector<size_t> row;row.reserve(n?n-1:0);
            for(size_t b=0;b<n;++b)if(a!=b)row.push_back(b);
            std::partial_sort(row.begin(),row.begin()+k,row.end(),[&](size_t b,size_t c){
                const double x=d.distance(a,b),y=d.distance(a,c);return x!=y?x<y:b<c;});
            row.resize(k);std::sort(row.begin(),row.end());neighbors[a]=std::move(row);
        }
        std::vector<Edge> edges;std::vector<double> scales;
        for(size_t a=0;a<n;++a)for(size_t b:neighbors[a])
            if(a<b && std::binary_search(neighbors[b].begin(),neighbors[b].end(),a)){
                const double w=d.distance(a,b);
                if(w>0 && w<endpoint){edges.push_back({a,b,w});scales.push_back(w);}
            }
        mutual_edges=edges.size();std::sort(scales.begin(),scales.end());
        scales.erase(std::unique(scales.begin(),scales.end()),scales.end());
        scales=rankedSubset(scales,32);
        for(auto& w:scales)w=std::nextafter(w,std::numeric_limits<double>::infinity());
        scales.erase(std::remove_if(scales.begin(),scales.end(),[&](double x){return x>=endpoint;}),scales.end());
        knn_s=std::chrono::duration<double>(Clock::now()-start).count();
        const auto cluster_start=Clock::now();
        std::map<std::vector<size_t>,size_t> previous;
        std::vector<double> proposed;
        for(double t:scales){
            UnionFind uf(n);for(const auto& e:edges)if(e.w<t)uf.join(e.a,e.b);
            std::map<size_t,std::vector<size_t>> components;
            std::vector<size_t> root_of(n),cluster_of(n);
            for(size_t a=0;a<n;++a){root_of[a]=uf.find(a);components[root_of[a]].push_back(a);}
            Snapshot snap;snap.epsilon=t;std::map<std::vector<size_t>,size_t> ages;
            for(auto& [root,vertices]:components){
                Cluster c;c.vertices=std::move(vertices);c.age=previous.contains(c.vertices)?previous.at(c.vertices)+1:1;
                for(size_t a:c.vertices)cluster_of[a]=snap.clusters.size();
                ages[c.vertices]=c.age;
                snap.clusters.push_back(std::move(c));
            }
            std::vector<size_t> internal_edges(snap.clusters.size(),0);
            for(size_t a=0;a<n;++a)for(size_t b=a+1;b<n;++b){
                const double w=d.distance(a,b);const size_t ca=cluster_of[a],cb=cluster_of[b];
                if(ca==cb){snap.clusters[ca].diameter=std::max(snap.clusters[ca].diameter,w);internal_edges[ca]+=w>0 && w<t;}
                else{snap.clusters[ca].outside=std::min(snap.clusters[ca].outside,w);snap.clusters[cb].outside=std::min(snap.clusters[cb].outside,w);}
            }
            for(size_t i=0;i<snap.clusters.size();++i){
                const auto& c=snap.clusters[i];const size_t m=c.vertices.size();
                snap.clusters[i].density=m>1?2.*internal_edges[i]/(double(m)*(m-1)):0;
                if(m>=3 && m<64){
                    const double q=std::nextafter(std::max(t,c.diameter/alpha),std::numeric_limits<double>::infinity());
                    if(q<endpoint)proposed.push_back(q);
                }
            }
            snapshots.push_back(std::move(snap));previous=std::move(ages);
        }
        for(size_t s=0;s+1<snapshots.size();++s){
            for(auto& c:snapshots[s].clusters){
                std::vector<bool> in(n,false);for(size_t a:c.vertices)in[a]=true;
                for(size_t b=0;b<n;++b)if(!in[b]){
                    bool joins=false;for(size_t a:c.vertices)if(d.distance(a,b)<snapshots[s+1].epsilon){joins=true;break;}
                    c.growth+=joins;
                }
            }
        }
        cluster_s=std::chrono::duration<double>(Clock::now()-cluster_start).count();
        const auto selection_start=Clock::now();
        std::sort(proposed.begin(),proposed.end());proposed.erase(std::unique(proposed.begin(),proposed.end()),proposed.end());
        // Bounded candidates prevent all component appearances from causing O(n^3) scans.
        proposed=rankedSubset(proposed,64);
        std::vector<Edge> work_edges;std::vector<size_t> final_degree(n,0),prefix_degree(n,0);
        for(size_t a=0;a<n;++a)for(size_t b=a+1;b<n;++b){
            const double w=d.distance(a,b);if(w>0 && w<endpoint){work_edges.push_back({a,b,w});++final_degree[a];++final_degree[b];}
        }
        double final_work=0;for(size_t v:final_degree)final_work+=double(v)*double(v?v-1:0)/2.;final_work=std::max(1.,final_work);
        std::sort(work_edges.begin(),work_edges.end(),[](const Edge& a,const Edge& b){return a.w<b.w;});
        size_t work_position=0;double prefix_work=0;
        for(double q:proposed){
            const Snapshot* snap=nullptr;
            for(const auto& s:snapshots){if(s.epsilon>q)break;snap=&s;}
            if(!snap)continue;
            Candidate c;c.epsilon=q;size_t eligible=0;
            for(const auto& cluster:snap->clusters){
                const size_t m=cluster.vertices.size();
                if(m>=3 && m<64 && cluster.diameter<=alpha*q && cluster.outside>=q){
                    c.gain+=m-1;c.stability+=double(cluster.age)/(cluster.age+1);c.growth+=cluster.growth;++eligible;
                }
            }
            if(eligible)c.stability/=eligible;
            while(work_position<work_edges.size() && work_edges[work_position].w<q){
                const auto& edge=work_edges[work_position++];
                prefix_work+=prefix_degree[edge.a]+prefix_degree[edge.b];++prefix_degree[edge.a];++prefix_degree[edge.b];
            }
            c.work=prefix_work/final_work;candidates.push_back(c);
        }
        for(size_t j=1;j<=cuts;++j){
            double best=0;const Candidate* chosen=nullptr;
            for(const auto& c:candidates){
                if(std::find(schedule.begin(),schedule.end(),c.epsilon)!=schedule.end())continue;
                const double gain=c.gain/(1.+c.growth/std::max(1.,c.gain));
                const double score=gain/std::max<size_t>(1,n)*c.stability*(1.-c.work)-.25*std::abs(c.work-double(j)/(cuts+1));
                if(score>best){best=score;chosen=&c;}
            }
            if(chosen)schedule.push_back(chosen->epsilon);
        }
        fallback=schedule.empty(); // Caller invokes the unchanged author's generator, scale 2.
        std::sort(schedule.begin(),schedule.end());
        if(!fallback)schedule.push_back(endpoint);
        selection_s=std::chrono::duration<double>(Clock::now()-selection_start).count();
        total_s=std::chrono::duration<double>(Clock::now()-start).count();
    }
    void report(bool details) const
    {
        std::cout<<std::setprecision(17)<<"[geometry] k=8 max_scales=32 max_candidates=64 requested_cuts="<<requested_cuts
            <<" proposed_cuts="<<(schedule.empty()?0:schedule.size()-1)<<" mutual_edges="<<mutual_edges
            <<" snapshots="<<snapshots.size()<<" candidates="<<candidates.size()<<" fallback="<<fallback
            <<" knn_s="<<knn_s<<" cluster_s="<<cluster_s<<" selection_s="<<selection_s<<" total_s="<<total_s<<'\n';
        if(details){
            for(const auto& s:snapshots)for(const auto& c:s.clusters){
                std::cout<<"[geometry-cluster] epsilon="<<s.epsilon<<" size="<<c.vertices.size()<<" diameter="<<c.diameter
                    <<" density="<<c.density<<" outside="<<c.outside<<" age="<<c.age<<" growth="<<c.growth<<" members=";
                for(auto v:c.vertices)std::cout<<v<<',';std::cout<<'\n';
            }
            for(const auto& c:candidates)std::cout<<"[geometry-candidate] epsilon="<<c.epsilon<<" gain="<<c.gain
                <<" stability="<<c.stability<<" growth="<<c.growth<<" work_fraction="<<c.work<<'\n';
        }
        for(size_t i=0;i<schedule.size();++i)std::cout<<"[geometry-schedule] index="<<i<<" epsilon="<<schedule[i]<<" final="<<(i+1==schedule.size())<<'\n';
        std::cout<<std::setprecision(6);
    }
};
