#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

// Costs are seconds; memory predictions are KiB, never added to seconds.
struct GCostController
{
    using Features=std::array<double,8>;
    Features rebuild{},residual{},memory{};
    std::array<double,3> select{};
    double max_labels=0,max_cofacets=0,scale=1;
    size_t dimension=0;
    struct Decision
    {
        double horizon=0,cont=0,q=0,r=0,f=0,g=0,threshold=0,mem_continue=0,mem_quotient=0;
        bool quotient=false;
        std::string reason="fallback-no-model";
    };
    static Features features(double a,double c,double pvs,double raw,double eps,double h,size_t d,bool q)
    {
        const double growth=std::log(h/eps), work=c/1e5, action=q?1.:0.;
        return {1.,work,work*growth,std::pow(a/100.,d+1)*growth,
                pvs/std::max(1.,a)*work,action*work,action*raw/100.,action};
    }
    static double dot(const Features& a,const Features& b)
    {
        double sum=0;for(size_t i=0;i<a.size();++i)sum+=a[i]*b[i];return sum;
    }
    void load(const std::string& path,size_t expected_dimension)
    {
        std::ifstream in(path);std::string tag;
        if (!(in>>tag>>dimension>>max_labels>>max_cofacets) || tag!="PWPH_G_V1" || dimension!=expected_dimension)
            throw std::invalid_argument("Invalid G model header/dimension");
        for(auto* values:{&rebuild,&residual,&memory})
            for(auto& x:*values)if(!(in>>x)||!std::isfinite(x))throw std::invalid_argument("Invalid G coefficients");
        for(auto& x:select)if(!(in>>x)||!std::isfinite(x))throw std::invalid_argument("Invalid G select coefficients");
        if (!(max_labels>0 && max_cofacets>0))throw std::invalid_argument("Invalid G domain");
    }
    Decision decide(double a,double c,double pvs,double raw,double entries,double eps,double h,
                    double peak_kib,double budget_kib,double margin) const
    {
        Decision out;out.horizon=h;
        const auto x0=features(a,c,pvs,raw,eps,h,dimension,false),x1=features(a,c,pvs,raw,eps,h,dimension,true);
        out.cont=scale*(std::max(0.,dot(rebuild,x0))+std::max(0.,dot(residual,x0)));
        out.q=scale*std::max(0.,select[0]+select[1]*raw/100.+select[2]*entries/1000.);
        out.r=scale*std::max(0.,dot(rebuild,x1));out.f=scale*std::max(0.,dot(residual,x1));
        out.g=out.cont-out.q-out.r-out.f;
        out.threshold=std::max(.0005,margin*out.cont);
        out.mem_continue=std::max(peak_kib,dot(memory,x0));out.mem_quotient=std::max(peak_kib,dot(memory,x1));
        const bool feasible0=out.mem_continue<=.7*budget_kib,feasible1=out.mem_quotient<=.7*budget_kib;
        if (!feasible0 && !feasible1)throw std::runtime_error("G_MEMORY_LIMIT: neither forecast action fits headroom");
        if (a>1.5*max_labels || c>1.5*max_cofacets)
        {
            if (!feasible0)throw std::runtime_error("G_MEMORY_LIMIT: out-of-domain continuation forecast infeasible");
            out.reason="fallback-out-of-domain";return out;
        }
        if (raw==0)
        {
            if (!feasible0)throw std::runtime_error("G_MEMORY_LIMIT: no supports for feasible quotient");
            out.reason="no-supports";return out;
        }
        if (!feasible0 && feasible1){out.quotient=true;out.reason="memory-constraint";return out;}
        out.quotient=feasible1 && out.g>out.threshold;
        out.reason=!feasible1?"memory-reject-quotient":out.quotient?"positive-G":"G-below-margin";
        return out;
    }
    void update(double observed,double predicted)
    {
        if (predicted>1e-5 && observed>0)
            scale=std::clamp(scale*(.8+.2*std::clamp(observed/predicted,.25,4.)),.25,4.);
    }
};
