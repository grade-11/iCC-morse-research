#pragma once
#include "DistanceMatrix.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>

// Window placement only. Matching, PV selection and witness weights are unchanged.
struct GrowthController
{
    double first_fraction=0,late_fraction=.70,reduction_threshold=.2;
    double first=0,middle=0,late=0,initialization_s=0,controller_s=0,next=0,reduction=0;
    size_t samples=0;
    bool deferred=false;
    std::string reason="unobserved";
    std::string placement_reason="fixed-sampled-quantile";
    double closure_threshold=0,coverage_threshold=0;
    struct Candidate{double fraction=0,epsilon=0,wedges=0,triangles=0,closure=0,coverage=0;bool accepted=false;};
    std::vector<Candidate> candidates;
    void load(const std::string& file,const DistanceMatrix& distances,double cutoff)
    {
        const auto start=std::chrono::steady_clock::now();
        std::ifstream in(file);std::string tag;
        if(!(in>>tag>>first_fraction>>late_fraction>>reduction_threshold) || tag!="GROWTH_V1" ||
           !(first_fraction>0 && first_fraction<late_fraction && late_fraction<1 && reduction_threshold>=0 && reduction_threshold<=1))
            throw std::invalid_argument("Invalid sampled-growth model");
        const size_t n=distances.vertexCount();
        if(n<3){first=middle=late=cutoff;reason="fewer-than-three-points";}
        else
        {
            uint64_t state=0x9e3779b97f4a7c15ULL;
            const auto draw=[&](){state^=state>>12;state^=state<<25;state^=state>>27;return (state*2685821657736338717ULL)%n;};
            std::vector<double> appearance;appearance.reserve(4096);
            for(size_t t=0;t<4096;++t)
            {
                const size_t i=draw();size_t j=draw(),k=draw();
                while(j==i)j=draw();while(k==i || k==j)k=draw();
                appearance.push_back(std::max({distances.distance(i,j),distances.distance(i,k),distances.distance(j,k)}));
            }
            samples=appearance.size();std::sort(appearance.begin(),appearance.end());
            const auto cut=[&](double q){const size_t rank=static_cast<size_t>(std::ceil(q*appearance.size()))-1;
                return std::min(cutoff,std::nextafter(appearance[rank],std::numeric_limits<double>::infinity()));};
            first=cut(first_fraction);middle=cut(.35);late=cut(late_fraction);
            std::string extension;
            if(in>>extension;extension=="GROWTH_GATE_V1")
            {
                if(!(in>>closure_threshold>>coverage_threshold) || !(closure_threshold>=0 && closure_threshold<=1 && coverage_threshold>=0 && coverage_threshold<=1))
                    throw std::invalid_argument("Invalid wedge closure model");
                first=cutoff;first_fraction=1.;placement_reason="skip-no-closed-covered-candidate";
                for(double fraction:std::array<double,4>{.02,.06,.15,.35})
                {
                    Candidate c;c.fraction=fraction;c.epsilon=cut(fraction);
                    std::vector<size_t> degree(n,0);
                    for(size_t i=1;i<n;++i)for(size_t j=0;j<i;++j)
                        if(distances.distance(i,j)<c.epsilon){++degree[i];++degree[j];}
                    size_t covered=0;
                    for(size_t d:degree){c.wedges+=static_cast<double>(d)*(d?d-1:0)/2.;covered+=d>=2;}
                    c.coverage=static_cast<double>(covered)/n;
                    const double observed=std::lower_bound(appearance.begin(),appearance.end(),c.epsilon)-appearance.begin();
                    c.triangles=observed/appearance.size()*(static_cast<double>(n)*(n-1)*(n-2)/6.);
                    c.closure=c.wedges>0?std::clamp(3*c.triangles/c.wedges,0.,1.):0.;
                    c.accepted=c.epsilon<cutoff && c.closure>=closure_threshold && c.coverage>=coverage_threshold;
                    candidates.push_back(c);
                    if(first==cutoff && c.accepted){first=c.epsilon;first_fraction=fraction;placement_reason="first-closed-covered-candidate";}
                }
                std::string ablation;
                if(in>>ablation)
                {
                    if(ablation=="FIXED15" || (ablation=="NO_SKIP" && first==cutoff))
                    {first=cut(.15);first_fraction=.15;placement_reason=ablation=="FIXED15"?"ablation-fixed15-with-features":"ablation-fallback15-instead-of-skip";}
                    else if(ablation!="NO_SKIP")throw std::invalid_argument("Unknown growth gate ablation");
                }
            }
        }
        next=cutoff;
        initialization_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    double decide(const std::string& policy,size_t before,size_t after,bool was_deferred,double cutoff)
    {
        const auto start=std::chrono::steady_clock::now();
        deferred=was_deferred;reduction=1.-static_cast<double>(after)/std::max<size_t>(1,before);
        next=cutoff;reason="direct-final";
        if(policy=="growth-extra"){next=late;reason="fixed-extra-ablation";}
        else if(policy=="growth-rule" && reduction<reduction_threshold)
        {next=late;reason=deferred?"first-deferred":"small-first-reduction";}
        if(!(next>first && next<cutoff)){next=cutoff;reason="skip-nonincreasing-or-terminal";}
        controller_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        return next;
    }
};
