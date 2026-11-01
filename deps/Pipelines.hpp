#pragma once

#include "DependencySupportPostProcessor.hpp"
#include "PipelineCommon.hpp"
#include "WindowState.hpp"

class PersistentHomologyPipeline
{
public:
    PersistentHomologyPipeline(
        DistanceMatrix distance_matrix,
        PipelineConfig config);

    void run();

private:
    PipelineRuntime runtime_;
};

class PwphPipeline
{
public:
    PwphPipeline(
        DistanceMatrix distance_matrix,
        PipelineConfig config);

    void run();

private:
    void runGuardExtendLadder();

    PipelineRuntime runtime_;
    WindowState window_state_;
    DependencySupportPostProcessor support_processor_;
};

class ApparentPairPipeline
{
public:
    ApparentPairPipeline(
        DistanceMatrix distance_matrix,
        PipelineConfig config);

    void run();

private:
    PipelineRuntime runtime_;
};
