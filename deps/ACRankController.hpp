#pragma once
#include "DistanceMatrix.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct ACRankController
{
    double first_fraction=.2,second_fraction=.5,threshold=.3;
    double first=0,second=0,late=0,initialization_s=0,reduction=0,next=2,controller_s=0;
    size_t first_protected=0;
    bool first_deferred=false;
    std::string reason="unobserved";
    void load(const std::string& file,const DistanceMatrix& distances,double cutoff)
    {
        const auto start=std::chrono::steady_clock::now();
        std::ifstream in(file);std::string tag;
        if(!(in>>tag>>first_fraction>>second_fraction>>threshold) || tag!="AC_RANK_V1" ||
            !(first_fraction>0 && first_fraction<second_fraction && second_fraction<1 && threshold>=0 && threshold<=1))
            throw std::invalid_argument("Invalid rank model");
        std::vector<double> values;const size_t n=distances.vertexCount();values.reserve(n*(n-1)/2);
        for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j)
        {
            const double d=distances.distance(i,j);
            if(!std::isfinite(d))throw std::invalid_argument("Nonfinite rank distance");
            if(d>0 && d<cutoff)values.push_back(d);
        }
        std::sort(values.begin(),values.end());values.erase(std::unique(values.begin(),values.end()),values.end());
        if(values.empty()) {first=cutoff;second=cutoff;late=cutoff;reason="no-positive-distance";}
        else
        {
            const auto rank=[&](double p){const auto k=static_cast<size_t>(std::ceil(p*values.size()))-1;
                return std::min(cutoff,std::nextafter(values[k],std::numeric_limits<double>::infinity()));};
            first=rank(first_fraction);second=rank(second_fraction);late=rank(.70);
        }
        initialization_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    double decide(const std::string& policy,size_t before,size_t after,size_t original,double peak_kib,double budget_kib,double cutoff,bool deferred=false,size_t protected_labels=0)
    {
        reduction=1.-static_cast<double>(after)/std::max<size_t>(1,before);
        first_deferred=deferred;first_protected=protected_labels;
        const double a=after,n=original;
        const double forecast_kib=(128.*std::max(0.,a*(a-1)*(a-2)/6.)+32.*n*n)/1024.;
        // A conservative empirical rebuild+matching proxy, not a proven bound.
        // Reject rather than silently process an infeasible final window.
        if(peak_kib+forecast_kib>budget_kib)throw std::runtime_error("AC rank predicted memory budget exceeded");
        if(policy.starts_with("guard"))
        {
            next=cutoff;reason="enough-reduction";
            if(policy=="guard50") {next=second;reason="fixed-middle-ablation";}
            else if(policy=="guard70") {next=late;reason="fixed-late-ablation";}
            else if(policy=="guard-base")reason="fixed-direct-ablation";
            else if(deferred) {next=late;reason=protected_labels?"deferred-unresolved-birth":"deferred-no-support";}
            else if(reduction<threshold && after>3) {next=second;reason="small-reduction";}
            if(!(next>first && next<cutoff))next=cutoff;
            return next;
        }
        const bool extra=policy=="rank-extra" || (policy=="rank" && reduction<threshold && after>3);
        next=extra && second>first && second<cutoff?second:cutoff;
        reason=next==cutoff?"direct-final":"low-first-reduction";
        if(policy=="rank-extra")reason="fixed-extra";
        if(policy=="rank-base")reason="fixed-direct";
        return next;
    }
};
