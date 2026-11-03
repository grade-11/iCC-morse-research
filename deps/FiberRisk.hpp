#pragma once
#include "WindowState.hpp"
#include <algorithm>
#include <array>
#include <iterator>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// Fiber check. The quotient maps the original VR complex onto its image on the labels; by the fiber
// lemma it keeps homology up to degree k when every fiber is acyclic up to degree k. A fiber with one
// pseudo vertex is a cone (its realizing vertex in the PV is adjacent to all of it), and H1EdgeRisk
// covers H1 over pairs. The default check covers the fiber over each triangle spanned by three
// multi-vertex blocks (at least one proposed): the original VR complex on their union must carry no
// H2 at any scale >= q where that triangle exists. With `full` it also checks the fibers over pairs of
// multi-vertex blocks and, for H2, over a pair with one singleton, and it reports H1 as well as H2
// bars (H1 only without `h2`). Development diagnostic, NOT a certificate: fibers over tetrahedra
// (four labels) are not checked.
struct FiberRisk
{
    struct Block {std::vector<size_t> vertices;bool proposed;size_t proposed_id;};
    static constexpr size_t NONE=std::numeric_limits<size_t>::max();
    // core: the vertices left after dropping dominated ones, which carry every bar from `from` on.
    // c is NONE for a pair fiber; a singleton is a block with one vertex.
    struct Warning {size_t a,b,c;double epsilon;std::vector<size_t> core;};
    struct Result {size_t triples=0,pairs=0,pair_vertex=0,vertices=0,undominated=0;double seconds=0;std::vector<Block> blocks;std::vector<Warning> warnings;};

    // Removes, one at a time, every vertex v dominated by another vertex u at all scales >= from: uv is
    // an edge by `from`, and every edge vx of the complex (weight below `endpoint`) has ux by
    // max(d(v,x), from). Then N[v] lies in N[u] at each such scale, so removing v is a strong collapse
    // that commutes with the inclusions: the module on [from, inf) is unchanged, and each H2 bar keeps
    // its death and max(birth, from).
    static std::vector<size_t> undominated(const DistanceMatrix& d,const std::vector<size_t>& vs,double from,double endpoint)
    {
        const size_t m=vs.size();
        std::vector<double> dm(m*m);
        for(size_t i=0;i<m;++i)for(size_t j=0;j<m;++j)dm[i*m+j]=d.distance(vs[i],vs[j]);
        const auto edge=[&](double w){return w>0 && w<endpoint;};
        std::vector<char> alive(m,1);
        for(bool changed=true;changed;){
            changed=false;
            for(size_t v=0;v<m;++v)if(alive[v])for(size_t u=0;u<m;++u){
                const double uv=dm[u*m+v];if(u==v || !alive[u] || !edge(uv) || uv>from)continue;
                bool dominated=true;
                for(size_t x=0;x<m && dominated;++x){
                    if(x==u || x==v || !alive[x] || !edge(dm[v*m+x]))continue;
                    const double ux=dm[u*m+x];dominated=edge(ux) && ux<=std::max(dm[v*m+x],from);
                }
                if(dominated){alive[v]=0;changed=true;break;}
            }
        }
        std::vector<size_t> out;for(size_t i=0;i<m;++i)if(alive[i])out.push_back(vs[i]);
        return out;
    }

    // Birth of the first H2 bar of the flag filtration on `vs` alive at some scale >= from, else -1.
    static double firstVoid(const DistanceMatrix& d,const std::vector<size_t>& vs,double from,double endpoint)
    {
        return firstBar(d,vs,from,endpoint,false,true);
    }

    // Birth of the first H1 (with h1) or H2 (with h2) bar of the flag filtration on `vs` alive at some
    // scale >= from, else -1.
    static double firstBar(const DistanceMatrix& d,const std::vector<size_t>& vs,double from,double endpoint,bool h1,bool h2)
    {
        const size_t m=vs.size();
        struct S {double w;int dim;std::vector<int> faces;};
        std::vector<S> simplices;
        const auto dist=[&](size_t i,size_t j){return d.distance(vs[i],vs[j]);};
        // Index simplices by vertex tuple through small lookup tables. An edge's faces are its vertices.
        std::vector<int> edge_id(m*m,-1),tri_id;
        std::vector<std::array<size_t,3>> tris;
        for(size_t i=0;i<m;++i)for(size_t j=i+1;j<m;++j){
            const double w=dist(i,j);if(!(w>0 && w<endpoint))continue;
            edge_id[i*m+j]=static_cast<int>(simplices.size());simplices.push_back({w,1,{static_cast<int>(i),static_cast<int>(j)}});
        }
        const auto edge=[&](size_t i,size_t j){return edge_id[std::min(i,j)*m+std::max(i,j)];};
        tri_id.assign(m*m*m,-1);
        for(size_t i=0;i<m;++i)for(size_t j=i+1;j<m;++j)for(size_t k=j+1;k<m;++k){
            const int a=edge(i,j),b=edge(i,k),c=edge(j,k);if(a<0||b<0||c<0)continue;
            const double w=std::max({simplices[a].w,simplices[b].w,simplices[c].w});
            tri_id[(i*m+j)*m+k]=static_cast<int>(simplices.size());simplices.push_back({w,2,{a,b,c}});
        }
        const auto tri=[&](size_t i,size_t j,size_t k){return tri_id[(i*m+j)*m+k];};
        if(h2)for(size_t i=0;i<m;++i)for(size_t j=i+1;j<m;++j)for(size_t k=j+1;k<m;++k)for(size_t l=k+1;l<m;++l){
            const int a=tri(i,j,k),b=tri(i,j,l),c=tri(i,k,l),e=tri(j,k,l);if(a<0||b<0||c<0||e<0)continue;
            const double w=std::max({simplices[a].w,simplices[b].w,simplices[c].w,simplices[e].w});
            simplices.push_back({w,3,{a,b,c,e}});
        }
        // Filtration order: weight, then dimension (faces first).
        std::vector<int> order(simplices.size());for(size_t i=0;i<order.size();++i)order[i]=static_cast<int>(i);
        std::sort(order.begin(),order.end(),[&](int x,int y){
            if(simplices[x].w!=simplices[y].w)return simplices[x].w<simplices[y].w;
            if(simplices[x].dim!=simplices[y].dim)return simplices[x].dim<simplices[y].dim;return x<y;});
        std::vector<int> pos(simplices.size());for(size_t r=0;r<order.size();++r)pos[order[r]]=static_cast<int>(r);
        // Reduce the triangle columns first: a triangle with a nonzero reduced column kills an H1 class
        // (negative); the others create voids. Then reduce the tetrahedron columns; each pivot is a
        // positive triangle and gives the H2 bar [w(triangle), w(tetrahedron)).
        std::vector<std::vector<int>> by_pivot(order.size());
        std::vector<char> is_pivot(order.size(),0),negative(order.size(),0);
        std::vector<std::vector<int>> low2(order.size());std::vector<char> pivot1(order.size(),0);
        double first=-1;
        const auto bar=[&](double birth,double death){if(death>birth && death>from && (first<0 || birth<first))first=birth;};
        for(int r:order){
            const auto& s=simplices[r];if(s.dim!=2)continue;
            std::vector<int> col;for(int f:s.faces)col.push_back(pos[f]);std::sort(col.begin(),col.end());
            while(!col.empty()){
                const int p=col.back();if(!pivot1[p])break;
                std::vector<int> out;std::set_symmetric_difference(col.begin(),col.end(),low2[p].begin(),low2[p].end(),std::back_inserter(out));col.swap(out);
            }
            if(!col.empty()){
                pivot1[col.back()]=1;low2[col.back()]=col;negative[pos[r]]=1;
                if(h1)bar(simplices[order[col.back()]].w,s.w);
            }
        }
        if(h1){
            // An edge that closes a cycle and is no triangle's pivot is an H1 class alive at the endpoint.
            std::vector<size_t> parent(m);for(size_t i=0;i<m;++i)parent[i]=i;
            const auto find=[&](size_t v){while(parent[v]!=v){parent[v]=parent[parent[v]];v=parent[v];}return v;};
            for(size_t r=0;r<order.size();++r){
                const auto& s=simplices[order[r]];if(s.dim!=1)continue;
                const size_t x=find(s.faces[0]),y=find(s.faces[1]);
                if(x!=y)parent[x]=y;
                else if(!pivot1[r] && (first<0 || s.w<first))first=s.w;
            }
        }
        if(!h2)return first;
        for(int r:order){
            const auto& s=simplices[r];if(s.dim!=3)continue;
            std::vector<int> col;for(int f:s.faces)col.push_back(pos[f]);std::sort(col.begin(),col.end());
            while(!col.empty()){
                const int p=col.back();if(!is_pivot[p])break;
                std::vector<int> out;std::set_symmetric_difference(col.begin(),col.end(),by_pivot[p].begin(),by_pivot[p].end(),std::back_inserter(out));col.swap(out);
            }
            if(col.empty())continue;
            const int p=col.back();is_pivot[p]=1;by_pivot[p]=col;
            bar(simplices[order[p]].w,s.w);
        }
        // A void still alive at the endpoint also counts.
        for(size_t r=0;r<order.size();++r){
            const auto& s=simplices[order[r]];
            if(s.dim==2 && !negative[r] && !is_pivot[r] && (first<0 || s.w<first))first=s.w;
        }
        return first;
    }

    // focus: only fibers with that proposed support (index into `selected`), else every fiber with one.
    // full, h2: see the header comment.
    static Result inspect(const DistanceMatrix& d,const WindowState& state,const std::vector<SelectedPV>& selected,
                          double q,double endpoint,size_t focus=std::numeric_limits<size_t>::max(),
                          bool full=false,bool h2=true)
    {
        const auto start=std::chrono::steady_clock::now();Result result;
        const size_t n=state.originalVertexCount();
        std::vector<char> absorbed(n,0);
        for(const auto& pv:selected)for(auto v:pv.representatives)absorbed[v]=1;
        std::vector<size_t> singles;
        for(size_t label:state.activeLabels()){
            if(label<n){if(!absorbed[label])singles.push_back(label);continue;}
            // Existing PVs absorbed by a proposed support (--absorb-pv) belong to that block.
            const auto& reps=state.pseudoVertices().at(label-n).representatives;
            if(absorbed[reps.front()])continue;
            if(reps.size()>1)result.blocks.push_back({reps,false,0});else singles.push_back(reps.front());
        }
        for(size_t p=0;p<selected.size();++p)if(selected[p].representatives.size()>1)
            result.blocks.push_back({selected[p].representatives,true,p});
        const size_t B=result.blocks.size();
        const auto is_focus=[&](const Block& x){return x.proposed && x.proposed_id==focus;};
        const auto wanted=[&](std::initializer_list<size_t> ids){
            bool proposed=false,focused=false;
            for(size_t id:ids){proposed|=result.blocks[id].proposed;focused|=is_focus(result.blocks[id]);}
            return proposed && (focus>=selected.size() || focused);
        };
        // far[v][b]: the largest distance from vertex v to block b; near[a][b]: the smallest distance
        // between blocks a and b. A fiber is a cone at every scale >= its radius, the least over its
        // vertices of the largest distance to the others, so only scales below the radius can carry a
        // bar. The radius comes from `far` in one pass, and a fiber whose radius is at most a lower
        // bound of its first scale (q or a crossing distance) is skipped before its triangle scale is
        // computed. This prunes exactly the fibers the radius test below prunes.
        const double inf=std::numeric_limits<double>::infinity();
        std::vector<std::vector<double>> far(n);
        const auto fill=[&](size_t v){
            if(!far[v].empty())return;
            far[v].assign(B,0.0);
            for(size_t b=0;b<B;++b)for(auto u:result.blocks[b].vertices)far[v][b]=std::max(far[v][b],d.distance(v,u));
        };
        for(size_t b=0;b<B;++b)for(auto v:result.blocks[b].vertices)fill(v);
        std::vector<double> near(B*B,inf);
        for(size_t a=0;a<B;++a)for(size_t b=a+1;b<B;++b){
            double w=inf;
            for(auto i:result.blocks[a].vertices)for(auto j:result.blocks[b].vertices)w=std::min(w,d.distance(i,j));
            near[a*B+b]=near[b*B+a]=w;
        }
        const auto radiusOf=[&](std::initializer_list<size_t> ids,size_t x){
            double radius=inf;
            for(size_t id:ids)for(auto v:result.blocks[id].vertices){
                double r=x==NONE?0.0:d.distance(v,x);for(size_t k:ids)r=std::max(r,far[v][k]);
                radius=std::min(radius,r);
            }
            if(x!=NONE){double r=0;for(size_t k:ids)r=std::max(r,far[x][k]);radius=std::min(radius,r);}
            return radius;
        };
        // The fiber over the image simplex on `ids` exists from w_tau; check it from max(q, w_tau) on.
        const auto check=[&](size_t a,size_t b,size_t c,double w_tau,double radius,size_t& count){
            if(!(w_tau<endpoint))return;
            const double from=std::max(q,w_tau);
            if(radius<=from)return;
            std::vector<size_t> vs;
            for(size_t id:{a,b,c})if(id!=NONE)vs.insert(vs.end(),result.blocks[id].vertices.begin(),result.blocks[id].vertices.end());
            ++count;
            const double end=std::min(endpoint,std::nextafter(radius,inf));
            result.vertices+=vs.size();vs=undominated(d,vs,from,end);result.undominated+=vs.size();
            const double v=full?firstBar(d,vs,from,end,true,h2):firstVoid(d,vs,from,end);
            if(v>=0)result.warnings.push_back({a,b,c,std::max(v,from),std::move(vs)});
        };
        for(size_t a=0;a<B;++a)for(size_t b=a+1;b<B;++b)for(size_t c=b+1;c<B;++c){
            if(!wanted({a,b,c}))continue;
            const double radius=radiusOf({a,b,c},NONE);
            if(radius<=std::max({q,near[a*B+b],near[a*B+c],near[b*B+c]}))continue;
            const auto& A=result.blocks[a];const auto& Bk=result.blocks[b];const auto& C=result.blocks[c];
            // Scale at which the triangle (A, B, C) exists in the image: best transversal triangle.
            double w_tau=inf;
            for(auto i:A.vertices)for(auto j:Bk.vertices){
                const double ij=d.distance(i,j);if(ij>=w_tau)continue;
                for(auto k:C.vertices)w_tau=std::min(w_tau,std::max({ij,d.distance(i,k),d.distance(j,k)}));
            }
            check(a,b,c,w_tau,radius,result.triples);
        }
        if(!full){
            result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            return result;
        }
        // Singletons become blocks after the multi-vertex ones; triples above never see them.
        if(h2)for(size_t v:singles){result.blocks.push_back({{v},false,0});fill(v);}
        for(size_t a=0;a<B;++a)for(size_t b=a+1;b<B;++b){
            if(!wanted({a,b}))continue;
            const auto& A=result.blocks[a];const auto& Bk=result.blocks[b];
            const double w_ab=near[a*B+b];
            check(a,b,NONE,w_ab,radiusOf({a,b},NONE),result.pairs);
            if(!(w_ab<endpoint))continue;
            for(size_t s=B;s<result.blocks.size();++s){
                const size_t x=result.blocks[s].vertices.front();
                double near_a=inf,near_b=inf;
                for(auto i:A.vertices)near_a=std::min(near_a,d.distance(i,x));
                for(auto j:Bk.vertices)near_b=std::min(near_b,d.distance(j,x));
                if(!(std::max({w_ab,near_a,near_b})<endpoint))continue;
                const double radius=radiusOf({a,b},x);
                if(radius<=std::max({q,w_ab,near_a,near_b}))continue;
                double w_tau=inf;
                for(auto i:A.vertices){
                    const double ix=d.distance(i,x);if(ix>=w_tau)continue;
                    for(auto j:Bk.vertices)w_tau=std::min(w_tau,std::max({ix,d.distance(i,j),d.distance(j,x)}));
                }
                check(a,b,s,w_tau,radius,result.pair_vertex);
            }
        }
        result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        return result;
    }

    // Proposed supports to drop so that no warning keeps all its proposed blocks.
    static std::vector<size_t> greedyDrops(const Result& r,size_t proposed_count)
    {
        std::vector<char> live(r.warnings.size(),1),dropped(proposed_count,0);std::vector<size_t> out;
        while(true){
            std::vector<size_t> incident(proposed_count,0);bool any=false;
            for(size_t w=0;w<r.warnings.size();++w)if(live[w]){
                any=true;
                for(size_t blk:{r.warnings[w].a,r.warnings[w].b,r.warnings[w].c})
                    if(blk!=NONE && r.blocks[blk].proposed)++incident[r.blocks[blk].proposed_id];
            }
            if(!any)return out;
            const size_t drop=static_cast<size_t>(std::max_element(incident.begin(),incident.end())-incident.begin());
            if(incident[drop]==0)throw std::logic_error("Fiber warning has no proposed support");
            dropped[drop]=1;out.push_back(drop);
            for(size_t w=0;w<r.warnings.size();++w)if(live[w])
                for(size_t blk:{r.warnings[w].a,r.warnings[w].b,r.warnings[w].c})
                    if(blk!=NONE && r.blocks[blk].proposed && r.blocks[blk].proposed_id==drop)live[w]=0;
        }
    }
};
