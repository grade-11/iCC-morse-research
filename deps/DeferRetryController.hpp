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

// Research-only bounded scheduler. It never changes support eligibility or the
// global/selective guards. The caller supplies ACTUAL committed-PV outcome.
struct DeferRetryController
{
    std::size_t degree=32, early_rank=0, positive_pairs=0, distinct_pairs=0;
    double retry_fraction=.5, first=0, retry=0, next=0;
    double initialization_s=0, controller_s=0;
    std::size_t decisions=0;
    bool early_available=false, retry_available=false, retry_scheduled=false;
    bool last_committed=false;
    std::string reason="uninitialized";

    void load(const std::string& file,const DistanceMatrix& distances,double cutoff)
    {
        const auto start=std::chrono::steady_clock::now();
        std::ifstream in(file);std::string tag,extra;long long requested_degree=0;
        double requested_fraction=0;
        if(!(in>>tag>>requested_degree>>requested_fraction) || tag!="AC_DEFER_RETRY_V1" ||
           requested_degree<=0 || !std::isfinite(requested_fraction) ||
           !(requested_fraction>0 && requested_fraction<1) || (in>>extra) ||
           !std::isfinite(cutoff) || cutoff<=0)
            throw std::invalid_argument("Invalid AC_DEFER_RETRY_V1 model or cutoff");
        degree=static_cast<std::size_t>(requested_degree);retry_fraction=requested_fraction;
        decisions=0;early_rank=0;retry_scheduled=false;last_committed=false;controller_s=0;
        early_available=false;retry_available=false;first=retry=next=cutoff;
        reason="ABSTAIN/degree-unavailable";
        std::vector<double> values;
        const auto n=distances.vertexCount();
        for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j)
        {
            const double value=distances.distance(i,j);
            if(!std::isfinite(value) || value<0)
                throw std::invalid_argument("Nonfinite or negative DEFER_RETRY distance");
            if(value>0 && value<cutoff)values.push_back(value);
        }
        std::sort(values.begin(),values.end());positive_pairs=values.size();
        const long double requested_rank=std::ceil(static_cast<long double>(n)*degree/2.L);
        if(requested_rank>=1 && requested_rank<=static_cast<long double>(values.size()))
        {
            early_rank=static_cast<std::size_t>(requested_rank);
            first=std::min(cutoff,std::nextafter(values[early_rank-1],std::numeric_limits<double>::infinity()));
            early_available=first<cutoff;
            reason=early_available?"early-ready":"ABSTAIN/early-at-final";
        }
        values.erase(std::unique(values.begin(),values.end()),values.end());distinct_pairs=values.size();
        if(!values.empty())
        {
            const auto k=static_cast<std::size_t>(std::ceil(retry_fraction*values.size()))-1;
            retry=std::min(cutoff,std::nextafter(values[k],std::numeric_limits<double>::infinity()));
        }
        retry_available=early_available && retry>first && retry<cutoff;
        next=first;
        initialization_s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }

    // Only one retry is possible. A nonzero support count is not a commit;
    // Pipelines supplies true only after processPwph returns committed_count>0.
    double decide(bool actual_committed,double current,double cutoff)
    {
        last_committed=actual_committed;
        next=cutoff;
        if(current>=cutoff){reason="final";return next;}
        if(decisions++==0)
        {
            if(actual_committed)reason="early-actual-COMMIT/direct-final";
            else if(retry_available && retry>current)
            {next=retry;retry_scheduled=true;reason="early-actual-DEFER/retry-prior50";}
            else reason="early-actual-DEFER/no-later-checkpoint";
        }
        else reason=actual_committed?"retry-actual-COMMIT/direct-final":"retry-actual-DEFER/direct-final";
        return next;
    }
};
