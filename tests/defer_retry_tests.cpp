#include "DeferRetryController.hpp"
#include "PipelineCommon.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace
{
int failures=0;
#define EXPECT(value) do{if(!(value)){std::cerr<<__LINE__<<": " #value "\n";++failures;}}while(false)
struct ModelFile
{
    std::filesystem::path path=std::filesystem::temp_directory_path()/
        ("icc-defer-retry-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".txt");
    ~ModelFile(){std::error_code ignored;std::filesystem::remove(path,ignored);}
    void write(const std::string& content){std::ofstream out(path);out<<content;}
};
bool rejected_model(ModelFile& file,const DistanceMatrix& matrix,const std::string& model,double cutoff=11.)
{
    file.write(model);
    try{DeferRetryController controller;controller.load(file.path.string(),matrix,cutoff);}
    catch(const std::invalid_argument&){return true;}
    return false;
}
bool rejected_config(const PipelineConfig& config)
{
    try{config.validate();}catch(const std::invalid_argument&){return true;}
    return false;
}
}

int main()
{
    ModelFile model;
    const double inf=std::numeric_limits<double>::infinity();
    DistanceMatrix line({{0.},{1.},{4.},{10.}});
    model.write("AC_DEFER_RETRY_V1 1 0.5\n");
    DeferRetryController deferred;deferred.load(model.path.string(),line,11.);
    EXPECT(deferred.early_rank==2 && deferred.positive_pairs==6 && deferred.distinct_pairs==6);
    EXPECT(deferred.first==std::nextafter(3.,inf));
    EXPECT(deferred.retry==std::nextafter(4.,inf));
    EXPECT(deferred.initialization_s>=0 && deferred.retry_available);
    const auto late=deferred.decide(false,deferred.first,11.);
    EXPECT(late==deferred.retry && deferred.retry_scheduled);
    EXPECT(deferred.reason=="early-actual-DEFER/retry-prior50");
    EXPECT(deferred.decide(false,late,11.)==11.);
    EXPECT(deferred.reason=="retry-actual-DEFER/direct-final");
    EXPECT(deferred.decide(false,11.,11.)==11. && deferred.reason=="final");

    DeferRetryController committed;committed.load(model.path.string(),line,11.);
    // Actual commit is the sole feedback, including a hypothetical singleton
    // commit with zero label reduction. No gain threshold substitutes for it.
    EXPECT(committed.decide(true,committed.first,11.)==11.);
    EXPECT(!committed.retry_scheduled && committed.last_committed);
    EXPECT(committed.reason=="early-actual-COMMIT/direct-final");
    DeferRetryController empty;empty.load(model.path.string(),line,11.);
    // raw/eligible-empty and selective-empty both give actual false and retry.
    EXPECT(empty.decide(false,empty.first,11.)==empty.retry);
    EXPECT(empty.decide(true,empty.retry,11.)==11.);
    EXPECT(empty.reason=="retry-actual-COMMIT/direct-final");

    DistanceMatrix square({{0.,0.},{1.,0.},{1.,1.},{0.,1.}});
    DeferRetryController ties;ties.load(model.path.string(),square,2.);
    EXPECT(ties.early_rank==2 && ties.positive_pairs==6 && ties.distinct_pairs==2);
    EXPECT(ties.first==std::nextafter(1.,inf) && ties.retry==ties.first);
    EXPECT(!ties.retry_available && ties.decide(false,ties.first,2.)==2.);
    EXPECT(ties.reason=="early-actual-DEFER/no-later-checkpoint");
    model.write("AC_DEFER_RETRY_V1 2 0.5");
    DeferRetryController backwards;backwards.load(model.path.string(),line,11.);
    EXPECT(backwards.first==std::nextafter(6.,inf) && backwards.retry<backwards.first);
    EXPECT(!backwards.retry_available && backwards.decide(false,backwards.first,11.)==11.);
    model.write("AC_DEFER_RETRY_V1 3 0.5");
    DeferRetryController at_final;at_final.load(model.path.string(),line,std::nextafter(10.,inf));
    EXPECT(!at_final.early_available && at_final.first==std::nextafter(10.,inf));
    EXPECT(at_final.reason=="ABSTAIN/early-at-final");
    model.write("AC_DEFER_RETRY_V1 32 0.5");
    DeferRetryController unavailable;unavailable.load(model.path.string(),line,11.);
    EXPECT(!unavailable.early_available && unavailable.first==11. && unavailable.decisions==0);
    DistanceMatrix duplicates({{0.},{0.},{1.}});
    model.write("AC_DEFER_RETRY_V1 1 0.5");
    DeferRetryController duplicate_rank;duplicate_rank.load(model.path.string(),duplicates,2.);
    EXPECT(duplicate_rank.positive_pairs==2 && duplicate_rank.early_rank==2 && duplicate_rank.first==std::nextafter(1.,inf));
    DeferRetryController strict;strict.load(model.path.string(),duplicates,1.);
    EXPECT(strict.positive_pairs==0 && strict.first==1. && !strict.early_available);

    for(const auto& text:{"WRONG 1 0.5","AC_DEFER_RETRY_V1 0 0.5","AC_DEFER_RETRY_V1 -1 0.5",
                         "AC_DEFER_RETRY_V1 1 0","AC_DEFER_RETRY_V1 1 1","AC_DEFER_RETRY_V1 1 nan",
                         "AC_DEFER_RETRY_V1 1 0.5 trailing","AC_DEFER_RETRY_V1 1.5 0.5",""})
        EXPECT(rejected_model(model,line,text));
    EXPECT(rejected_model(model,line,"AC_DEFER_RETRY_V1 1 0.5",inf));
    EXPECT(rejected_model(model,line,"AC_DEFER_RETRY_V1 1 0.5",0.));

    PipelineConfig config;config.maxdim=2;config.threads=4;config.pv_cap_scale=.85;
    config.pv_min_separation=0;config.eps_breaks={11.};config.ac_policy="guard-defer-retry";
    config.ac_model=model.path.string();config.h1_risk="selective-fast";
    EXPECT(!rejected_config(config));
    auto control=config;control.ac_policy="guard-defer-retry-control";
    EXPECT(!rejected_config(control));
    auto invalid=config;invalid.h1_risk="observe";EXPECT(rejected_config(invalid));
    invalid=config;invalid.h1_risk="";EXPECT(rejected_config(invalid));
    invalid=config;invalid.maxdim=3;EXPECT(rejected_config(invalid));
    invalid=config;invalid.eps_breaks={2.,11.};EXPECT(rejected_config(invalid));
    invalid=config;invalid.pack_mode="plan";EXPECT(rejected_config(invalid));
    invalid=config;invalid.quotient_mask={1};EXPECT(rejected_config(invalid));
    invalid=config;invalid.ac_model="";EXPECT(rejected_config(invalid));
    invalid=config;invalid.pv_cap_scale=.9;EXPECT(rejected_config(invalid));
    invalid=config;invalid.pv_min_separation=.01;EXPECT(rejected_config(invalid));
    if(!failures)std::cout<<"DEFER_RETRY_V1 ranks/ties/duplicates/strict cutoff, actual outcomes, bounded retry and guard validation PASS\n";
    return failures?1:0;
}
