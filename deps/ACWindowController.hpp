#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

// A decision made AFTER the current legal quotient and BEFORE next rebuild.
// Each action predicts all remaining work through the SAME final cutoff.
// Costs are seconds, memory is KiB. No PH endpoints are controller inputs.
struct ACWindowController
{
    using Features=std::array<double,8>;
    Features mean{},scale{};
    std::array<Features,4> cost{},memory{};
    struct Decision {
        double next=2,controller_s=0;
        Features x{};
        std::array<double,4> seconds{},kib{};
        int selected=3;
        std::string reason="cost-minimum";
    };
    static Features features(double n,double c,double a,double pv,double raw,double edges)
    {
        return {1.,std::log(n/200.),std::log1p(c/1e5),a/n,pv/n,
            std::log1p(raw/100.),edges/std::max(1.,n*(n-1.)/2.),(n-a+pv)/n};
    }
    void load(const std::string& path)
    {
        std::ifstream in(path);std::string tag;size_t nf=0,na=0;
        if (!(in>>tag>>nf>>na) || tag!="AC_WINDOW_V1" || nf!=8 || na!=4)
            throw std::invalid_argument("Invalid AC model header");
        for(auto* row:{&mean,&scale})for(auto& v:*row)
            if(!(in>>v)||!std::isfinite(v))throw std::invalid_argument("Invalid AC feature scale");
        for(double s:scale)if(s<=0)throw std::invalid_argument("Invalid AC scale");
        for(auto* rows:{&cost,&memory})for(auto& row:*rows)for(auto& v:row)
            if(!(in>>v)||!std::isfinite(v))throw std::invalid_argument("Invalid AC coefficient");
    }
    Decision decide(double n,double c,double a,double pv,double raw,double edges,double peak,double budget) const
    {
        Decision out;out.x=features(n,c,a,pv,raw,edges);
        double best=1e300;
        constexpr std::array<double,4> candidates{.9,1.2,1.5,2.};
        for(size_t j=0;j<4;++j)
        {
            double t=0,m=0;
            for(size_t k=0;k<8;++k){double z=(out.x[k]-mean[k])/scale[k];t+=z*cost[j][k];m+=z*memory[j][k];}
            out.seconds[j]=std::exp(std::clamp(t,-20.,20.));
            // Include current peak and reserve 50% above predicted total peak;
            // this is an empirical guard, NOT a hard memory guarantee.
            out.kib[j]=std::max(peak,1.5*std::exp(std::clamp(m,0.,30.)));
            if(out.kib[j]<budget && out.seconds[j]<best)
                {best=out.seconds[j];out.next=candidates[j];out.selected=j;}
        }
        if(best==1e300)out.reason="memory-fallback-final";
        return out;
    }
};
