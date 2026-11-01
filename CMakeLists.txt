cmake_minimum_required(VERSION 3.22)

project(iCC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
include(FetchContent)

find_package(OpenMP REQUIRED)

find_package(CLI11 CONFIG QUIET)
if(NOT CLI11_FOUND)
    FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG v2.6.2
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(CLI11)
endif()

add_library(icc_core STATIC
    deps/readInput.cpp
    deps/DistanceMatrix.cpp
    deps/PipelineCommon.cpp
    deps/WindowState.cpp
    deps/DimensionFrame.cpp
    deps/ImplicitMorseMatching.cpp
    deps/DependencySupportPostProcessor.cpp
    deps/Pipelines.cpp
    deps/SimplexEnumerator.cpp
)

target_include_directories(icc_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/deps
)

target_link_libraries(icc_core
    PUBLIC
        OpenMP::OpenMP_CXX
)

add_executable(iCC main.cpp)
target_link_libraries(iCC PRIVATE icc_core CLI11::CLI11)

if(BUILD_TESTING)
    add_executable(icc_defer_retry_tests tests/defer_retry_tests.cpp)
    target_link_libraries(icc_defer_retry_tests PRIVATE icc_core)
    add_test(NAME icc_defer_retry_tests COMMAND icc_defer_retry_tests)
    add_executable(icc_pack_filter_tests tests/pack_filter_tests.cpp)
    target_link_libraries(icc_pack_filter_tests PRIVATE icc_core)
    add_test(NAME icc_pack_filter_tests COMMAND icc_pack_filter_tests)
    add_executable(icc_geometry_tests tests/geometry_tests.cpp)
    target_link_libraries(icc_geometry_tests PRIVATE icc_core)
    add_test(NAME icc_geometry_tests COMMAND icc_geometry_tests)
    add_executable(icc_core_tests tests/core_tests.cpp)
    target_link_libraries(icc_core_tests PRIVATE icc_core)
    add_test(NAME icc_core_tests COMMAND icc_core_tests)

    add_test(
        NAME icc_cli_ph_smoke
        COMMAND iCC
            --tool ph
            --file-name ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/irregular5.csv
            --max-dim 3
            --max-eps 3.0
            --threads 1
    )
    set_tests_properties(
        icc_cli_ph_smoke
        PROPERTIES
            PASS_REGULAR_EXPRESSION "\\(1.1401754250991381, 1.2041594578792296\\)"
            FAIL_REGULAR_EXPRESSION "birth facet labels:"
    )

    add_test(
        NAME icc_cli_ph_verbose
        COMMAND iCC
            --tool ph
            --file-name ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/irregular5.csv
            --max-dim 3
            --max-eps 3.0
            --threads 1
            --verbose
    )
    set_tests_properties(
        icc_cli_ph_verbose
        PROPERTIES
            PASS_REGULAR_EXPRESSION "birth facet labels: 3 1"
    )
endif()
