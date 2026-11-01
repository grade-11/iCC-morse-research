#pragma once
#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

// An explicit falsifiable cost extrapolation, separate from quotient eligibility.
struct PaybackAnalyticController
{
    double first_fraction=.15,margin=.35,min_gain_s=.002,loading_s=0;
    struct Decision
    {
        std::array<double,3> seconds{};
        std::array<double,2> predicted_active{};
        double next=0,reduction=0;
        size_t selected=0;
        std::string reason="no-nonfinal-window";
    };
    void load(const std::string& file)
    {
        const auto start=std::chrono::steady_clock::now();
        std::ifstream in(file);std::string growth,tag;double late=0,unused=0;
        if(!(in>>growth>>first_fraction>>late>>unused>>tag>>margin>>min_gain_s) || growth!="GROWTH_V1" || tag!="PAYBACK_ANALYTIC_V1" ||
           !(first_fraction>0 && first_fraction<.35 && margin>=0 && margin<1 && min_gain_s>=0) || !std::isfinite(min_gain_s))
            throw std::invalid_argument("Invalid analytic payback model");
        loading_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    static double triangles(double a){return a>=3?std::max(0.,a*(a-1)*(a-2)/6.):0.;}
    static double burden(double c){c=std::max(1.,c);return c*std::log2(2+c);}
    Decision decide(size_t before,size_t after,size_t pv,size_t cofacets,double work_s,double quotient_s,
                    double current,double middle,double late,double cutoff) const
    {
        Decision out;out.next=cutoff;
        if(!cofacets || !(work_s>0)){out.reason="no-reliable-work-rate";return out;}
        out.reduction=std::clamp(1.-static_cast<double>(after)/std::max<size_t>(1,before),0.,1.);
        const double full=triangles(after),price=work_s/burden(cofacets);
        out.seconds[0]=price*burden(full);
        const std::array<double,2> fractions{.35,.70};
        const std::array<double,3> candidates{cutoff,middle,late};
        for(size_t j=0;j<2;++j)
        {
            out.predicted_active[j]=pv+(after-pv)*std::pow(1.-out.reduction,fractions[j]/first_fraction);
            out.seconds[j+1]=price*(burden(fractions[j]*full)+burden(triangles(out.predicted_active[j])))+quotient_s;
        }
        size_t best=0;
        for(size_t j=1;j<3;++j)
            if(candidates[j]>current && candidates[j]<cutoff && out.seconds[j]<out.seconds[best])best=j;
        if(best && out.seconds[0]-out.seconds[best]>std::max(min_gain_s,margin*out.seconds[0]))
        {out.selected=best;out.next=candidates[best];out.reason="predicted-payback-with-margin";}
        else out.reason=out.reduction==0?"no-observed-reduction":"direct-below-payback-margin";
        return out;
    }
};
