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

// Three-block fiber check for H2. The quotient maps the original VR complex onto its image on the
// labels; by the fiber lemma it keeps homology when every fiber is acyclic. Fibers with one pseudo
// vertex are cones, and H1EdgeRisk covers pairs. This checks the fiber over each triangle spanned by
// three multi-vertex blocks (at least one proposed): the original VR complex on their union must
// carry no H2 at any scale >= q where that triangle exists. Development diagnostic, NOT a certificate
// for fibers with singletons or over tetrahedra.
struct FiberRisk
{
    struct Block {std::vector<size_t> vertices;bool proposed;size_t proposed_id;};
    // core: the vertices left after dropping dominated ones, which carry every void from `from` on.
    struct Warning {size_t a,b,c;double epsilon;std::vector<size_t> core;};
    struct Result {size_t triples=0,vertices=0,undominated=0;double seconds=0;std::vector<Block> blocks;std::vector<Warning> warnings;};

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
        const size_t m=vs.size();
        struct S {double w;int dim;std::vector<int> faces;};
        std::vector<S> simplices;
        const auto dist=[&](size_t i,size_t j){return d.distance(vs[i],vs[j]);};
        // Index simplices by vertex tuple through small lookup tables.
        std::vector<int> edge_id(m*m,-1),tri_id;
        std::vector<std::array<size_t,3>> tris;
        for(size_t i=0;i<m;++i)for(size_t j=i+1;j<m;++j){
            const double w=dist(i,j);if(!(w>0 && w<endpoint))continue;
            edge_id[i*m+j]=static_cast<int>(simplices.size());simplices.push_back({w,1,{}});
        }
        const auto edge=[&](size_t i,size_t j){return edge_id[std::min(i,j)*m+std::max(i,j)];};
        tri_id.assign(m*m*m,-1);
        for(size_t i=0;i<m;++i)for(size_t j=i+1;j<m;++j)for(size_t k=j+1;k<m;++k){
            const int a=edge(i,j),b=edge(i,k),c=edge(j,k);if(a<0||b<0||c<0)continue;
            const double w=std::max({simplices[a].w,simplices[b].w,simplices[c].w});
            tri_id[(i*m+j)*m+k]=static_cast<int>(simplices.size());simplices.push_back({w,2,{a,b,c}});
        }
        const auto tri=[&](size_t i,size_t j,size_t k){return tri_id[(i*m+j)*m+k];};
        for(size_t i=0;i<m;++i)for(size_t j=i+1;j<m;++j)for(size_t k=j+1;k<m;++k)for(size_t l=k+1;l<m;++l){
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
        for(int r:order){
            const auto& s=simplices[r];if(s.dim!=2)continue;
            std::vector<int> col;for(int f:s.faces)col.push_back(pos[f]);std::sort(col.begin(),col.end());
            while(!col.empty()){
                const int p=col.back();if(!pivot1[p])break;
                std::vector<int> out;std::set_symmetric_difference(col.begin(),col.end(),low2[p].begin(),low2[p].end(),std::back_inserter(out));col.swap(out);
            }
            if(!col.empty()){pivot1[col.back()]=1;low2[col.back()]=col;negative[pos[r]]=1;}
        }
        double first=-1;
        for(int r:order){
            const auto& s=simplices[r];if(s.dim!=3)continue;
            std::vector<int> col;for(int f:s.faces)col.push_back(pos[f]);std::sort(col.begin(),col.end());
            while(!col.empty()){
                const int p=col.back();if(!is_pivot[p])break;
                std::vector<int> out;std::set_symmetric_difference(col.begin(),col.end(),by_pivot[p].begin(),by_pivot[p].end(),std::back_inserter(out));col.swap(out);
            }
            if(col.empty())continue;
            const int p=col.back();is_pivot[p]=1;by_pivot[p]=col;
            const double birth=simplices[order[p]].w,death=s.w;
            if(death>birth && death>from && (first<0 || birth<first))first=birth;
        }
        // A void still alive at the endpoint also counts.
        for(size_t r=0;r<order.size();++r){
            const auto& s=simplices[order[r]];
            if(s.dim==2 && !negative[r] && !is_pivot[r] && (first<0 || s.w<first))first=s.w;
        }
        return first;
    }

    // focus: only triples with that proposed support (index into `selected`), else every triple with one.
    static Result inspect(const DistanceMatrix& d,const WindowState& state,const std::vector<SelectedPV>& selected,
                          double q,double endpoint,size_t focus=std::numeric_limits<size_t>::max())
    {
        const auto start=std::chrono::steady_clock::now();Result result;
        const size_t n=state.originalVertexCount();
        std::vector<char> absorbed(n,0);
        for(const auto& pv:selected)for(auto v:pv.representatives)absorbed[v]=1;
        for(size_t label:state.activeLabels())if(label>=n){
            // Existing PVs absorbed by a proposed support (--absorb-pv) belong to that block.
            const auto& reps=state.pseudoVertices().at(label-n).representatives;
            if(reps.size()>1 && !absorbed[reps.front()])result.blocks.push_back({reps,false,0});
        }
        for(size_t p=0;p<selected.size();++p)if(selected[p].representatives.size()>1)
            result.blocks.push_back({selected[p].representatives,true,p});
        const size_t B=result.blocks.size();
        for(size_t a=0;a<B;++a)for(size_t b=a+1;b<B;++b)for(size_t c=b+1;c<B;++c){
            const auto& A=result.blocks[a];const auto& Bk=result.blocks[b];const auto& C=result.blocks[c];
            if(!A.proposed && !Bk.proposed && !C.proposed)continue;
            const auto is_focus=[&](const Block& x){return x.proposed && x.proposed_id==focus;};
            if(focus<selected.size() && !is_focus(A) && !is_focus(Bk) && !is_focus(C))continue;
            // Scale at which the triangle (A, B, C) exists in the image: best transversal triangle.
            double w_tau=std::numeric_limits<double>::infinity();
            for(auto i:A.vertices)for(auto j:Bk.vertices){
                const double ij=d.distance(i,j);if(ij>=w_tau)continue;
                for(auto k:C.vertices)w_tau=std::min(w_tau,std::max({ij,d.distance(i,k),d.distance(j,k)}));
            }
            if(!(w_tau<endpoint))continue;
            ++result.triples;
            std::vector<size_t> vs=A.vertices;vs.insert(vs.end(),Bk.vertices.begin(),Bk.vertices.end());vs.insert(vs.end(),C.vertices.begin(),C.vertices.end());
            const double from=std::max(q,w_tau);
            // Once one vertex is adjacent to all others the fiber is a cone and stays one: only
            // scales below that radius can carry a void.
            double radius=std::numeric_limits<double>::infinity();
            for(size_t i=0;i<vs.size();++i){
                double far=0;for(size_t j=0;j<vs.size();++j)if(j!=i)far=std::max(far,d.distance(vs[i],vs[j]));
                radius=std::min(radius,far);
            }
            if(radius<=from){--result.triples;continue;}
            const double end=std::min(endpoint,std::nextafter(radius,std::numeric_limits<double>::infinity()));
            result.vertices+=vs.size();vs=undominated(d,vs,from,end);result.undominated+=vs.size();
            const double v=firstVoid(d,vs,from,end);
            if(v>=0)result.warnings.push_back({a,b,c,std::max(v,from),std::move(vs)});
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
                    if(r.blocks[blk].proposed)++incident[r.blocks[blk].proposed_id];
            }
            if(!any)return out;
            const size_t drop=static_cast<size_t>(std::max_element(incident.begin(),incident.end())-incident.begin());
            if(incident[drop]==0)throw std::logic_error("Fiber warning has no proposed support");
            dropped[drop]=1;out.push_back(drop);
            for(size_t w=0;w<r.warnings.size();++w)if(live[w])
                for(size_t blk:{r.warnings[w].a,r.warnings[w].b,r.warnings[w].c})
                    if(r.blocks[blk].proposed && r.blocks[blk].proposed_id==drop)live[w]=0;
        }
    }
};
