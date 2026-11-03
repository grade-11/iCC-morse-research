#include "DistanceMatrix.hpp"
#include "Pipelines.hpp"
#include "readInput.hpp"
#include "GeometryController.hpp"
#include "PackDmtPayback.hpp"

#include <CLI/CLI.hpp>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr const char* PIECEWISE_PH_TOOL = "piecewise";
constexpr const char* MORSE_PH_TOOL = "ph";
constexpr const char* APPARENT_PAIR_TOOL = "apparent";
constexpr double DEFAULT_EPS_INTERVAL_SCALE = 1.0;
// SimplexEnumerator stores group coverage in a uint64_t, with group indices 0..maxdim.
constexpr size_t MAX_SIMPLEX_DIMENSION = 63;

class HelpSpacingFormatter : public CLI::Formatter
{
public:
    std::string make_option(const CLI::Option* option, bool is_positional) const override
    {
        std::string text = CLI::Formatter::make_option(option, is_positional);
        const auto& long_names = option->get_lnames();

        if (!is_positional &&
            std::find(long_names.begin(), long_names.end(), "help") != long_names.end())
        {
            text += '\n';
        }

        return text;
    }
};

std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";

    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::string normalizeTool(std::string tool)
{
    tool = trim(tool);
    std::transform(tool.begin(), tool.end(), tool.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });

    if (tool == "1" || tool == "piecewise" || tool == "morsepiecewiseph" || tool == "morse-piecewise-ph")
        return PIECEWISE_PH_TOOL;

    if (tool == "2" || tool == "ph" || tool == "morseph" || tool == "morse-ph")
        return MORSE_PH_TOOL;

    if (tool == "3" || tool == "apparent" || tool == "apparent-pairs" || tool == "apparentpairs")
        return APPARENT_PAIR_TOOL;

    throw std::invalid_argument(
        "Unknown tool '" + tool + "'. Use 'piecewise', 'ph', or 'apparent'.");
}

bool isPiecewisePHTool(const std::string& tool)
{
    return normalizeTool(tool) == PIECEWISE_PH_TOOL;
}

const char* displayToolName(const std::string& tool)
{
    const std::string normalized_tool = normalizeTool(tool);
    if (normalized_tool == PIECEWISE_PH_TOOL)
        return "morsePiecewisePH";
    if (normalized_tool == MORSE_PH_TOOL)
        return "morsePH";
    return "apparentPairs";
}

std::string promptRequiredLine(const std::string& prompt)
{
    while (true)
    {
        std::cout << prompt;

        std::string line;
        if (!std::getline(std::cin, line))
            throw std::runtime_error("Input stream closed.");

        line = trim(line);
        if (!line.empty())
            return line;

        std::cout << "  Please enter a value.\n";
    }
}

std::string promptFilename()
{
    while (true)
    {
        const std::string filename = promptRequiredLine("Input CSV filename: ");
        const std::filesystem::path filepath(filename);

        if (!std::filesystem::exists(filepath))
        {
            std::cout << "  File not found. Please enter an existing path.\n";
            continue;
        }

        if (!std::filesystem::is_regular_file(filepath))
        {
            std::cout << "  That path is not a regular file.\n";
            continue;
        }

        std::ifstream file_stream(filepath);
        if (!file_stream.good())
        {
            std::cout << "  Unable to open that file for reading.\n";
            continue;
        }

        return filename;
    }
}

std::string promptTool()
{
    while (true)
    {
        std::cout << "Select tool:\n";
        std::cout << "  1) morsePiecewisePH\n";
        std::cout << "  2) morsePH\n";
        std::cout << "  3) apparentPairs\n";

        const std::string answer = promptRequiredLine("Input tool [1/2/3]: ");
        try
        {
            return normalizeTool(answer);
        }
        catch (const std::invalid_argument&)
        {
            std::cout << "  Please enter 1 for morsePiecewisePH, 2 for morsePH, or 3 for apparentPairs.\n";
        }
    }
}

size_t promptMaxDimension()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input max dimension: ");
        std::stringstream parser(line);
        long long value = 0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || value <= 0 ||
            value > static_cast<long long>(MAX_SIMPLEX_DIMENSION))
        {
            std::cout << "  Please enter an integer from 1 to " << MAX_SIMPLEX_DIMENSION << ".\n";
            continue;
        }

        return static_cast<size_t>(value);
    }
}

std::vector<double> promptEpsilonBreaks()
{
    while (true)
    {
        const std::string raw_line = promptRequiredLine("Input epsilon breaks (space or comma separated): ");
        std::string normalized = raw_line;
        std::replace(normalized.begin(), normalized.end(), ',', ' ');

        std::stringstream parser(normalized);
        std::vector<double> eps_breaks;
        std::string token;
        bool valid = true;

        while (parser >> token)
        {
            size_t parsed_length = 0;
            double value = 0.0;

            try
            {
                value = std::stod(token, &parsed_length);
            }
            catch (const std::exception&)
            {
                valid = false;
                break;
            }

            if (parsed_length != token.size() || value <= 0.0)
            {
                valid = false;
                break;
            }

            eps_breaks.push_back(value);
        }

        if (!valid || eps_breaks.empty())
        {
            std::cout << "  Please enter one or more positive numbers, for example: 1.0 1.4 1.7\n";
            continue;
        }

        bool strictly_increasing = true;
        for (size_t i = 1; i < eps_breaks.size(); ++i)
        {
            if (eps_breaks[i] <= eps_breaks[i - 1])
            {
                strictly_increasing = false;
                break;
            }
        }

        if (!strictly_increasing)
        {
            std::cout << "  Epsilon breaks must be strictly increasing.\n";
            continue;
        }

        return eps_breaks;
    }
}

bool promptUseAutomaticEpsilonBreaks()
{
    while (true)
    {
        std::cout << "Select epsilon input method:\n";
        std::cout << "  1) explicit epsilon breaks\n";
        std::cout << "  2) number of intervals from sorted pairwise distances\n";

        const std::string answer = promptRequiredLine("Input epsilon method [1/2]: ");
        if (answer == "1" || answer == "explicit" || answer == "breaks")
            return false;
        if (answer == "2" || answer == "auto" || answer == "count" || answer == "interval-count" || answer == "break-count" || answer == "windows")
            return true;

        std::cout << "  Please enter 1 for explicit breaks or 2 for automatic intervals.\n";
    }
}

size_t promptEpsilonIntervalCount()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input epsilon interval count: ");
        std::stringstream parser(line);
        long long value = 0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || value <= 0)
        {
            std::cout << "  Please enter a positive integer.\n";
            continue;
        }

        return static_cast<size_t>(value);
    }
}

bool isValidEpsilonIntervalScale(const double value)
{
    return value >= 1.0 && std::isfinite(value);
}

double promptEpsilonIntervalScale()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input epsilon interval scale (>=1.0; 1.0 = linear, larger values make earlier intervals larger): ");
        std::stringstream parser(line);
        double value = 0.0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || !isValidEpsilonIntervalScale(value))
        {
            std::cout << "  Please enter a finite number greater than or equal to 1.0.\n";
            continue;
        }

        return value;
    }
}

int promptThreadNumber()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input thread number: ");
        std::stringstream parser(line);
        int value = 0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || value <= 0)
        {
            std::cout << "  Please enter a positive integer.\n";
            continue;
        }

        return value;
    }
}

double promptMaxEpsilon()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input max epsilon: ");
        std::stringstream parser(line);
        double value = 0.0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || value <= 0.0)
        {
            std::cout << "  Please enter a positive number.\n";
            continue;
        }

        return value;
    }
}

double promptPVCapScale()
{
    while (true)
    {
        const std::string line = promptRequiredLine("Input PV diameter cap scale: ");
        std::stringstream parser(line);
        double value = 0.0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || !std::isfinite(value) || value <= 0.0)
        {
            std::cout << "  Please enter a positive number.\n";
            continue;
        }

        return value;
    }
}

double promptPVMinSeparation()
{
    while (true)
    {
        const std::string line = promptRequiredLine(
            "Input minimum inter-PV separation scale (>=0; 0 disables): ");
        std::stringstream parser(line);
        double value = 0.0;
        char extra = '\0';

        if (!(parser >> value) || (parser >> extra) || !std::isfinite(value) || value < 0.0)
        {
            std::cout << "  Please enter a nonnegative number.\n";
            continue;
        }

        return value;
    }
}

std::vector<double> collectSortedUniquePositiveDistances(const DistanceMatrix& distance_matrix)
{
    const size_t npts = distance_matrix.vertexCount();
    std::vector<double> distances;
    distances.reserve((npts > 1) ? (npts * (npts - 1)) / 2 : 0);

    for (size_t i = 0; i + 1 < npts; ++i)
    {
        for (size_t j = i + 1; j < npts; ++j)
        {
            const double distance = distance_matrix.distance(i, j);
            if (!std::isfinite(distance))
                throw std::runtime_error("Encountered a non-finite pairwise distance while generating epsilon breaks.");
            if (distance > 0.0)
                distances.push_back(distance);
        }
    }

    if (distances.empty())
        throw std::runtime_error("Automatic epsilon breaks require at least one positive pairwise distance.");

    std::sort(distances.begin(), distances.end());
    distances.erase(std::unique(distances.begin(), distances.end()), distances.end());
    return distances;
}

std::vector<double> generateEpsilonBreaksFromIntervalCount(const DistanceMatrix& distance_matrix,
                                                           const size_t interval_count,
                                                           const double interval_scale)
{
    if (interval_count == 0)
        throw std::invalid_argument("Epsilon interval count must be positive.");
    if (!isValidEpsilonIntervalScale(interval_scale))
        throw std::invalid_argument("Epsilon interval scale must be a finite number greater than or equal to 1.0.");

    const auto distances = collectSortedUniquePositiveDistances(distance_matrix);
    const size_t distance_count = distances.size();

    if (interval_count > distance_count)
        throw std::invalid_argument("Epsilon interval count exceeds the number of distinct positive pairwise distances.");

    std::vector<double> eps_breaks;
    eps_breaks.reserve(interval_count);

    size_t previous_index = 0;
    bool has_previous_index = false;

    for (size_t j = 1; j <= interval_count; ++j)
    {
        const double linear_fraction = static_cast<double>(j) / static_cast<double>(interval_count);

        // Scale 1 gives linear distance ranks. Larger values move earlier cut
        // points upward, so earlier intervals contain more distinct distances.
        double scaled_fraction = 1.0 - std::pow(1.0 - linear_fraction, interval_scale);
        scaled_fraction = std::clamp(scaled_fraction, 0.0, 1.0);

        size_t raw_index = 0;
        if (scaled_fraction >= 1.0)
        {
            raw_index = distance_count - 1;
        }
        else if (scaled_fraction > 0.0)
        {
            raw_index = static_cast<size_t>(std::ceil(scaled_fraction * static_cast<double>(distance_count))) - 1;
        }

        const size_t remaining_intervals = interval_count - j;
        const size_t lower_bound_index = has_previous_index ? previous_index + 1 : 0;
        const size_t upper_bound_index = distance_count - 1 - remaining_intervals;

        if (lower_bound_index > upper_bound_index)
            throw std::logic_error("Internal epsilon break index selection error.");

        const size_t selected_index = std::clamp(raw_index, lower_bound_index, upper_bound_index);
        previous_index = selected_index;
        has_previous_index = true;

        // The VR enumerators use strict upper bounds. Moving to the next
        // representable double includes the selected discrete distance.
        eps_breaks.push_back(std::nextafter(distances[selected_index], std::numeric_limits<double>::infinity()));
    }

    return eps_breaks;
}

std::vector<double> resolveEpsilonBreaks(const DistanceMatrix& distance_matrix,
                                         const std::vector<double>& explicit_eps_breaks,
                                         const size_t eps_interval_count,
                                         const double eps_interval_scale)
{
    if (!explicit_eps_breaks.empty())
        return explicit_eps_breaks;

    return generateEpsilonBreaksFromIntervalCount(distance_matrix, eps_interval_count, eps_interval_scale);
}

void printRunConfiguration(const std::string& tool,
                           const std::string& filename,
                           const size_t maxdim,
                           const int threadnumber,
                           const std::vector<double>& eps_breaks,
                           const size_t eps_interval_count,
                           const double eps_interval_scale,
                           const double maxeps,
                           const double pv_cap_scale,
                           const double pv_min_separation)
{
    const std::string normalized_tool = normalizeTool(tool);

    std::cout << "\nRun " << displayToolName(normalized_tool) << " with:\n";
    std::cout << "  input file: " << filename << '\n';
    std::cout << "  max dimension: " << maxdim << '\n';
    std::cout << "  thread number: " << threadnumber << '\n';

    if (normalized_tool == PIECEWISE_PH_TOOL)
    {
        if (!eps_breaks.empty())
        {
            std::cout << "  epsilon breaks:";
            for (const double eps : eps_breaks)
                std::cout << ' ' << eps;
            std::cout << '\n';
        }
        else
        {
            std::cout << "  epsilon interval count: " << eps_interval_count << '\n';
            std::cout << "  epsilon interval scale: " << eps_interval_scale << '\n';
        }
        std::cout << "  pv cap scale: " << pv_cap_scale << '\n';
        std::cout << "  pv minimum separation scale: " << pv_min_separation << '\n';
    }
    else
    {
        std::cout << "  max epsilon: " << maxeps << '\n';
    }
}

bool promptConfirmation(const std::string& tool,
                        const std::string& filename,
                        const size_t maxdim,
                        const int threadnumber,
                        const std::vector<double>& eps_breaks,
                        const size_t eps_interval_count,
                        const double eps_interval_scale,
                        const double maxeps,
                        const double pv_cap_scale,
                        const double pv_min_separation)
{
    while (true)
    {
        printRunConfiguration(tool, filename, maxdim, threadnumber, eps_breaks, eps_interval_count,
                              eps_interval_scale, maxeps, pv_cap_scale, pv_min_separation);

        const std::string answer = promptRequiredLine("Proceed? [y/n]: ");
        if (answer == "y" || answer == "Y" || answer == "yes" || answer == "YES" || answer == "Yes")
            return true;
        if (answer == "n" || answer == "N" || answer == "no" || answer == "NO" || answer == "No")
            return false;

        std::cout << "  Please answer y or n.\n";
    }
}

void validateRunOptions(const std::string& tool,
                        const std::string& filename,
                        const size_t maxdim,
                        const int threadnumber,
                        const std::vector<double>& eps_breaks,
                        const size_t eps_interval_count,
                        const double eps_interval_scale,
                        const double maxeps)
{
    const std::string normalized_tool = normalizeTool(tool);

    if (filename.empty())
        throw std::invalid_argument("Input filename is required.");

    const std::filesystem::path filepath(filename);
    if (!std::filesystem::exists(filepath))
        throw std::invalid_argument("Input file does not exist: " + filename);
    if (!std::filesystem::is_regular_file(filepath))
        throw std::invalid_argument("Input path is not a regular file: " + filename);

    if (maxdim == 0)
        throw std::invalid_argument("Maximum dimension must be a positive integer.");
    if (maxdim > MAX_SIMPLEX_DIMENSION)
        throw std::invalid_argument("Maximum dimension must be less than 64.");
    if (threadnumber <= 0)
        throw std::invalid_argument("Thread number must be a positive integer.");

    if (normalized_tool == PIECEWISE_PH_TOOL)
    {
        const bool has_explicit_breaks = !eps_breaks.empty();
        const bool has_automatic_intervals = (eps_interval_count > 0);

        if (has_explicit_breaks == has_automatic_intervals)
            throw std::invalid_argument("When --tool piecewise is selected, specify exactly one of --eps-breaks or --eps-interval-count.");

        if (has_automatic_intervals && !isValidEpsilonIntervalScale(eps_interval_scale))
            throw std::invalid_argument("--eps-interval-scale must be a finite number greater than or equal to 1.0.");
    }
    else
    {
        if (!std::isfinite(maxeps) || maxeps <= 0.0)
        {
            if (normalized_tool == MORSE_PH_TOOL)
            {
                throw std::invalid_argument(
                    "--max-eps is required and must be finite and positive when --tool ph is selected.");
            }
            throw std::invalid_argument(
                "--max-eps is required and must be finite and positive when --tool apparent is selected.");
        }
    }
}

int runTool(const std::string& tool,
            const std::string& filename,
            const size_t maxdim,
            const int threadnumber,
            const std::vector<double>& eps_breaks,
            const size_t eps_interval_count,
            const double eps_interval_scale,
            const double maxeps,
            const double pv_cap_scale,
            const double pv_min_separation,
            const bool print_configuration,
            const bool verbose,
            const bool adaptive = false,
            const double adaptive_target = 0.35,
            const bool trace_pwph = false,
            const std::string& g_model = "",
            const double g_margin = .05,
            const double g_memory_mib = 4096,
            const bool cost_telemetry = false,
            const std::vector<int>& quotient_mask = {},
            const size_t stop_after_window = 0,
            const std::string& ac_policy = "",
            const std::string& ac_model = "",
            const double ac_auto_cutoff = 0,
            const size_t geometry_cuts = 0,
            const std::string& h1_risk = "",
            const bool h1_realizations = false,
            const std::string& pack_mode = "",
            const std::string& pack_model = "",
            const std::string& payback_action = "auto",
            const bool calibration_targets = false,
            const std::string& extend_dump = "",
            const bool extend_force_chain = false,
            const bool extend_direct_probe = false,
            const bool birth_repair = false,
            const bool live_guard = false,
            const bool live_guard_extend = false,
            const bool extend_gate = false,
            const bool extend_commit_once = false,
            const double extend_min_gain = 0.0,
            const bool extend_cost_rule = false,
            const bool fiber_risk = false,
            const bool absorb_pv = false,
            const bool shrink_risky = false,
            const bool shrink_search = false,
            const bool joint_drops = false,
            const bool fiber_full = false)
{
    const std::string normalized_tool = normalizeTool(tool);
    if(ac_policy=="guard-extend-ladder" || extend_direct_probe || extend_force_chain || !extend_dump.empty())
    {
        if(normalized_tool!=PIECEWISE_PH_TOOL || maxdim!=2 || eps_breaks.empty() ||
           eps_interval_count || geometry_cuts || ac_auto_cutoff>0)
            throw std::invalid_argument("EXTEND and its gate probes require piecewise --max-dim 2 with explicit --eps-breaks");
    }
    // The H1 risk filter composes with the live guard; it only drops supports.
    if((birth_repair || live_guard) && (normalized_tool!=PIECEWISE_PH_TOOL || !ac_policy.empty() || (!h1_risk.empty() && !live_guard) || adaptive ||
        !g_model.empty() || geometry_cuts || !pack_mode.empty() || !quotient_mask.empty() || eps_breaks.empty()))
        throw std::invalid_argument("--birth-repair and --live-guard require plain piecewise with explicit --eps-breaks");
    if((live_guard_extend || extend_gate) && (!live_guard || maxdim<2 || maxdim>3 || eps_interval_count || ac_auto_cutoff>0 ||
        stop_after_window || cost_telemetry))
        throw std::invalid_argument("--extend requires --live-guard, --max-dim 2 or 3 and explicit --eps-breaks only");
    if(extend_gate && !live_guard_extend)
        throw std::invalid_argument("--extend-gate requires --extend");
    if(extend_commit_once && !live_guard_extend)
        throw std::invalid_argument("--extend-commit-once requires --extend");
    if(extend_min_gain!=0.0 && (!live_guard_extend || extend_commit_once || !(extend_min_gain>0.0 && extend_min_gain<1.0)))
        throw std::invalid_argument("--extend-min-gain requires --extend, no --extend-commit-once, and a value in (0, 1)");
    if(fiber_risk && (!live_guard || maxdim!=3))
        throw std::invalid_argument("--fiber-risk requires --live-guard and --max-dim 3");
    if(fiber_full && (!live_guard || maxdim<2 || maxdim>3))
        throw std::invalid_argument("--fiber-full requires --live-guard and --max-dim 2 or 3");
    if((shrink_search || joint_drops) && !shrink_risky)
        throw std::invalid_argument("--shrink-search and --joint-drops require --shrink-risky");
    if(shrink_risky && h1_risk!="selective-fast")
        throw std::invalid_argument("--shrink-risky requires --h1-risk selective-fast");
    if(absorb_pv && (!live_guard || pv_min_separation>0))
        throw std::invalid_argument("--absorb-pv requires --live-guard and --pv-min-separation 0");
    if(extend_cost_rule && (!live_guard_extend || extend_commit_once || extend_min_gain!=0.0))
        throw std::invalid_argument("--extend-cost-rule requires --extend and excludes --extend-commit-once and --extend-min-gain");
    if(h1_realizations && normalized_tool!=PIECEWISE_PH_TOOL)
        throw std::invalid_argument("H1 realization diagnostics require piecewise");
    if((geometry_cuts || !h1_risk.empty()) && (normalized_tool!=PIECEWISE_PH_TOOL ||
        (maxdim!=2 && !(maxdim==3 && live_guard && !geometry_cuts))))
        throw std::invalid_argument("Geometry and H1 risk require piecewise --max-dim 2 (H1 risk: or 3 with --live-guard)");
    if(geometry_cuts && (geometry_cuts>2 || eps_breaks.size()!=1 || eps_interval_count || !ac_policy.empty() || adaptive || !g_model.empty() || !quotient_mask.empty()))
        throw std::invalid_argument("Geometry requires cuts=1 or 2, one explicit final epsilon, no other scheduler/mask");
    if ((!ac_policy.empty() || ac_auto_cutoff>0) && normalized_tool!=PIECEWISE_PH_TOOL)
        throw std::invalid_argument("AC options require piecewise");
    if (adaptive && normalized_tool != PIECEWISE_PH_TOOL)
        throw std::invalid_argument("--adaptive requires --tool piecewise.");
    if (normalized_tool != PIECEWISE_PH_TOOL && (!g_model.empty() || !quotient_mask.empty() || cost_telemetry))
        throw std::invalid_argument("G/quotient/cost telemetry options require piecewise");

    if(!pack_mode.empty() && (normalized_tool!=PIECEWISE_PH_TOOL || maxdim!=2 || geometry_cuts || adaptive || !g_model.empty() ||
        (pack_mode!="plan" && pack_mode!="measure" && pack_mode!="probe" && pack_mode!="replay") ||
        ((pack_mode=="plan" || pack_mode=="probe") && (eps_breaks.size()!=1 || !ac_policy.empty() || eps_interval_count))))
        throw std::invalid_argument("Invalid packing mode or conflicting scheduler");
    if(!pack_mode.empty() && (pv_cap_scale!=.85 || pv_min_separation!=0))throw std::invalid_argument("Packing locks alpha=.85 and delta=0");
    if((pack_mode=="plan" || pack_mode=="replay") && pack_model.empty())throw std::invalid_argument("Packing requires frozen model");
    if(payback_action!="auto" && payback_action!="commit" && payback_action!="defer")throw std::invalid_argument("Invalid payback action");
    if((pack_mode=="plan" || pack_mode=="probe" || pack_mode=="replay") && h1_risk!="selective-fast" && h1_risk!="selective")
        throw std::invalid_argument("Packing actions require locked selective filter");
    validateRunOptions(normalized_tool, filename, maxdim, threadnumber, eps_breaks,
                       eps_interval_count, eps_interval_scale, maxeps);

    if (print_configuration && verbose)
        printRunConfiguration(normalized_tool, filename, maxdim, threadnumber, eps_breaks,
                              eps_interval_count, eps_interval_scale, maxeps, pv_cap_scale,
                              pv_min_separation);

    const auto input_start=std::chrono::steady_clock::now();
    auto point_cloud = readInput::readCSV(filename);
    const auto distance_start=std::chrono::steady_clock::now();
    omp_set_num_threads(threadnumber);
    DistanceMatrix distance_matrix(point_cloud);
    const auto distance_end=std::chrono::steady_clock::now();
    decltype(point_cloud){}.swap(point_cloud);

    PackDmtPayback packing;
    std::vector<double> resolved_eps_breaks;
    if (normalized_tool == PIECEWISE_PH_TOOL)
    {
        resolved_eps_breaks = resolveEpsilonBreaks(
            distance_matrix, eps_breaks, eps_interval_count, eps_interval_scale);
        if(ac_auto_cutoff>0)
        {
            if(eps_interval_count==0 || !std::isfinite(ac_auto_cutoff) || ac_auto_cutoff<resolved_eps_breaks.back())
                throw std::invalid_argument("Auto cutoff must extend native generated final bound");
            resolved_eps_breaks.back()=ac_auto_cutoff;
        }
        if(geometry_cuts)
        {
            const double endpoint=resolved_eps_breaks.back();
            GeometryController geometry;
            geometry.load(distance_matrix,endpoint,pv_cap_scale,geometry_cuts);
            if(geometry.fallback)
            {
                const auto fallback_start=std::chrono::steady_clock::now();
                const size_t count=std::min(geometry_cuts+1,collectSortedUniquePositiveDistances(distance_matrix).size());
                if(count)geometry.schedule=generateEpsilonBreaksFromIntervalCount(distance_matrix,count,2.);
                geometry.schedule.erase(std::remove_if(geometry.schedule.begin(),geometry.schedule.end(),[&](double x){return x>=endpoint;}),geometry.schedule.end());
                geometry.schedule.push_back(endpoint);
                const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-fallback_start).count();
                geometry.selection_s+=elapsed;geometry.total_s+=elapsed;
            }
            resolved_eps_breaks=geometry.schedule;
            geometry.report(verbose);
            std::cout<<std::setprecision(17)<<"[geometry-input] read_s="<<std::chrono::duration<double>(distance_start-input_start).count()
                <<" distance_s="<<std::chrono::duration<double>(distance_end-distance_start).count()<<'\n'<<std::setprecision(6);
        }

        if(!pack_mode.empty() && pack_mode!="replay"){
            packing.load(distance_matrix,resolved_eps_breaks.back(),pack_mode,pack_model);
            if(pack_mode!="measure"){
                resolved_eps_breaks={packing.endpoint};
                if(packing.chosen<packing.endpoint)resolved_eps_breaks.insert(resolved_eps_breaks.begin(),packing.chosen);
            }
            packing.report(pack_mode,verbose);
        }
        if(pack_mode=="replay")packing.model.load(pack_model);
        bool duplicate_domain=false;
        packing.full_endpoint=PackDmtPayback::fullEndpoint(distance_matrix,resolved_eps_breaks.back(),duplicate_domain);
        if (verbose && eps_breaks.empty())
        {
            std::cout << "  generated epsilon breaks:";
            for (const double eps : resolved_eps_breaks)
                std::cout << ' ' << eps;
            std::cout << '\n';
        }
    }

    if (verbose)
        std::cout << "\nstarted running...\n";

    PipelineConfig config{
        .eps_breaks = normalized_tool == PIECEWISE_PH_TOOL
            ? resolved_eps_breaks
            : std::vector<double>{maxeps},
        .maxdim = maxdim,
        .threads = threadnumber,
        .pv_cap_scale = normalized_tool == PIECEWISE_PH_TOOL ? pv_cap_scale : 1.0,
        .pv_min_separation = normalized_tool == PIECEWISE_PH_TOOL ? pv_min_separation : 0.0,
        .verbose = verbose,
        .adaptive = adaptive,
        .adaptive_target = adaptive_target,
        .trace_pwph = trace_pwph,
        .g_model = g_model,
        .g_margin = g_margin,
        .g_memory_mib = g_memory_mib,
        .cost_telemetry = cost_telemetry,
        .quotient_mask = quotient_mask,
        .stop_after_window = stop_after_window,
        .ac_policy = ac_policy,
        .ac_model = ac_model,
        .h1_risk = h1_risk,
        .h1_realizations = h1_realizations,
        .pack_mode = pack_mode, .pack_model = pack_model, .payback_action = payback_action,
        .full_endpoint = packing.full_endpoint, .calibration_targets = calibration_targets,
        .payback_aL = packing.model.aL, .payback_aM = packing.model.aM, .payback_bC = packing.model.bC,
        .extend_dump = extend_dump, .extend_force_chain = extend_force_chain,
        .extend_direct_probe = extend_direct_probe, .birth_repair = birth_repair || live_guard, .live_guard = live_guard,
        .live_guard_extend = live_guard_extend, .extend_gate = extend_gate,
        .extend_commit_once = extend_commit_once, .extend_min_gain = extend_min_gain,
        .extend_cost_rule = extend_cost_rule, .fiber_risk = fiber_risk || fiber_full, .absorb_pv = absorb_pv,
        .shrink_risky = shrink_risky, .shrink_search = shrink_search,
        .joint_drops = joint_drops, .fiber_full = fiber_full};

    const auto st0 = std::chrono::steady_clock::now();
    if (normalized_tool == PIECEWISE_PH_TOOL)
    {
        PwphPipeline pipeline(std::move(distance_matrix), std::move(config));
        pipeline.run();
    }
    else if (normalized_tool == MORSE_PH_TOOL)
    {
        PersistentHomologyPipeline pipeline(
            std::move(distance_matrix), std::move(config));
        pipeline.run();
    }
    else
    {
        ApparentPairPipeline pipeline(
            std::move(distance_matrix), std::move(config));
        pipeline.run();
    }
    const auto st1 = std::chrono::steady_clock::now();

    const auto pt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(st1 - st0);
    std::cout << "run time = " << pt_ms.count() << '\n';

    return 0;
}

int runInteractive()
{
    std::cout << "Morse PH interactive runner\n\n";

    const std::string filename = promptFilename();
    const std::string tool = promptTool();
    const size_t maxdim = promptMaxDimension();
    const int threadnumber = promptThreadNumber();

    std::vector<double> eps_breaks;
    size_t eps_interval_count = 0;
    double eps_interval_scale = DEFAULT_EPS_INTERVAL_SCALE;
    double maxeps = 0.0;
    double pv_cap_scale = 0.0;
    double pv_min_separation = 0.0;

    if (isPiecewisePHTool(tool))
    {
        if (promptUseAutomaticEpsilonBreaks())
        {
            eps_interval_count = promptEpsilonIntervalCount();
            eps_interval_scale = promptEpsilonIntervalScale();
        }
        else
        {
            eps_breaks = promptEpsilonBreaks();
        }

        pv_cap_scale = promptPVCapScale();
        pv_min_separation = promptPVMinSeparation();
    }
    else
    {
        maxeps = promptMaxEpsilon();
    }

    if (!promptConfirmation(tool, filename, maxdim, threadnumber, eps_breaks, eps_interval_count,
                            eps_interval_scale, maxeps, pv_cap_scale, pv_min_separation))
    {
        std::cout << "Run cancelled.\n";
        return 0;
    }

    return runTool(tool, filename, maxdim, threadnumber, eps_breaks, eps_interval_count,
                   eps_interval_scale, maxeps, pv_cap_scale, pv_min_separation, false, true);
}

int runCommandLine(int argc, char** argv)
{
    CLI::App app{"Morse PH command-line runner"};
    app.formatter(std::make_shared<HelpSpacingFormatter>());

    std::string tool;
    std::string filename;
    size_t maxdim = 0;
    int threadnumber = 1;
    std::vector<double> eps_breaks;
    size_t eps_interval_count = 0;
    double eps_interval_scale = DEFAULT_EPS_INTERVAL_SCALE;
    double maxeps = 0.0;
    double pv_cap_scale = 0.0;
    double pv_min_separation = 0.0;
    bool verbose = false;

    app.add_option("-t,--tool", tool, "Tool to run: piecewise, ph, or apparent")
        ->required();

    app.add_option("-f,--file-name", filename, "Input CSV filename")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("-d,--max-dim", maxdim, "Maximum simplex dimension (1-63)")
        ->required()
        ->check(CLI::Range(size_t{1}, MAX_SIMPLEX_DIMENSION));

    app.add_option("-n,--threads", threadnumber, "Number of worker threads")
        ->default_val(1)
        ->check(CLI::PositiveNumber);

    app.add_flag("-v,--verbose", verbose, "Show diagnostic output");
    bool adaptive = false;
    std::string ac_policy,ac_model;
    std::string extend_dump;
    bool extend_force_chain=false,extend_direct_probe=false;
    app.add_option("--extend-dump",extend_dump,"Gate-only JSONL snapshot file; never enable for timing");
    app.add_flag("--extend-force-chain",extend_force_chain,"Gate-only continuation through every explicit cut; requires --extend-dump");
    app.add_flag("--extend-direct-probe",extend_direct_probe,"Gate-only fresh first-window oracle using ordinary matching; requires --extend-dump");
    bool birth_repair=false,live_guard=false;
    app.add_flag("--live-guard",live_guard,"Veto every quotient while any class of any dimension is alive at the cut");
    bool live_guard_extend=false,extend_gate=false;
    app.add_flag("--extend",live_guard_extend,"With --live-guard: retain the first window across vetoed cuts and extend it to the next cut");
    app.add_flag("--extend-gate",extend_gate,"Gate only: with --extend, extend through every cut and compare with a direct window at each");
    bool extend_commit_once=false;
    app.add_flag("--extend-commit-once",extend_commit_once,"With --extend: after the first committed quotient, close with one window up to the endpoint");
    double extend_min_gain=0.0;
    app.add_option("--extend-min-gain",extend_min_gain,"With --extend: keep cutting only while a commit removes at least this share of the next window's predicted top-dimension work");
    bool fiber_risk=false;
    app.add_flag("--fiber-risk",fiber_risk,"With --live-guard --max-dim 3: drop supports whose three-block fiber carries a void");
    bool absorb_pv=false;
    app.add_flag("--absorb-pv",absorb_pv,"With --live-guard: a support may absorb existing PVs; its flattened vertex set forms one new PV");
    bool shrink_risky=false;
    app.add_flag("--shrink-risky",shrink_risky,"With --h1-risk selective-fast: shrink a flagged support to its safe part instead of dropping it");
    bool shrink_search=false;
    app.add_flag("--shrink-search",shrink_search,"With --shrink-risky: first try removing a single unit of the offered part");
    bool joint_drops=false;
    app.add_flag("--joint-drops",joint_drops,"With --shrink-risky: choose drops for both filters together, weighted by support size");
    bool fiber_full=false;
    app.add_flag("--fiber-full",fiber_full,"With --live-guard --max-dim 2 or 3: the fiber check also covers pairs and (H2) pairs with one vertex, for H1 and H2 bars; implies --fiber-risk");
    bool extend_cost_rule=false;
    app.add_flag("--extend-cost-rule",extend_cost_rule,"With --extend: keep cutting only while the next cut is predicted to save more work than its extra window costs");
    app.add_flag("--birth-repair",birth_repair,"Keep the pre-quotient birth of a class alive across a quotient; veto a quotient two classes of one dimension would survive");
    double ac_auto_cutoff=0;
    size_t geometry_cuts=0;
    std::string h1_risk;
    std::string pack_mode,pack_model,payback_action="auto";
    bool calibration_targets=false;
    bool h1_realizations=false;
    app.add_flag("--h1-realizations",h1_realizations,"Post-rebuild H1 canonical edge/triangle witness diagnostics; requires verbose, no action filtering");
    app.add_flag("--calibration-targets",calibration_targets,"Development-only native workload features, excluded from native timer");
    app.add_option("--pack-mode",pack_mode,"Metric-only packing: measure, plan, probe, or same-prefix replay");
    app.add_option("--pack-model",pack_model,"Frozen PACK_DMT_PAYBACK_V1 cost model");
    app.add_option("--payback-action",payback_action,"Development marginal ablation: auto, commit, defer");
    app.add_option("--geometry-cuts",geometry_cuts,"Mutual-kNN geometry schedule: 1 or 2 intermediate quotients; supply E as sole --eps-breaks");
    app.add_option("--h1-risk",h1_risk,"Precommit H1 filter: edge-fiber, selective, or observe (diagnostic without rejection)");
    app.add_option("--ac-policy",ac_policy,"Research window controller: cost");
    app.add_option("--ac-model",ac_model,"Development-fitted AC controller model");
    app.add_option("--ac-auto-cutoff",ac_auto_cutoff,"Extend native auto final bound to common research cutoff");
    bool trace_pwph = false;
    std::string g_model;
    double g_margin=.05,g_memory_mib=4096;
    bool cost_telemetry=false;
    std::vector<int> quotient_mask;
    size_t stop_after_window=0;
    double adaptive_target = 0.35;
    app.add_flag("--adaptive", adaptive, "Pilot online scheduling on the explicit candidate grid");
    app.add_flag("--trace-pwph", trace_pwph, "Diagnostic simplex and PV membership trace (not for benchmarks)");
    app.add_option("--g-model",g_model,"Cost model file; enables online quotient/defer G controller")->check(CLI::ExistingFile);
    app.add_option("--g-margin",g_margin,"Minimum fractional predicted benefit; also requires >=0.5ms")->default_val(.05);
    app.add_option("--g-memory-mib",g_memory_mib,"Forecast memory constraint, MiB, with 30% headroom")->default_val(4096);
    app.add_flag("--cost-telemetry",cost_telemetry,"Cheap nonverbose window cost records for calibration");
    app.add_option("--quotient-mask",quotient_mask,"Forced 0/1 actions at nonfinal checkpoints; diagnostic/offline replay")
        ->delimiter(',')->check(CLI::Range(0,1));
    app.add_option("--stop-after-window",stop_after_window,"Calibration only: stop after this many completed windows");
    app.add_option("--adaptive-target", adaptive_target, "Target fractional active-label reduction in (0,1]")
        ->default_val(0.35);

    // For a std::vector option, CLI11 already accepts one or more values when
    // the option is present. Conditional requirement and resolved-value
    // validation are handled after parsing.
    app.add_option("--eps-breaks", eps_breaks,
                   "For --tool piecewise. Strictly increasing epsilon breaks. Accepts spaces or commas, e.g. '--eps-breaks 1.0 1.4 1.7' or '--eps-breaks=1.0,1.4,1.7'")
        ->delimiter(',')
        ->check(CLI::PositiveNumber);

    app.add_option("--eps-interval-count,--eps-break-count,--eps-window-count", eps_interval_count,
                   "For --tool piecewise. Alternative to --eps-breaks. Builds this many epsilon intervals from sorted distinct pairwise distances")
        ->check(CLI::PositiveNumber);

    app.add_option("--eps-interval-scale,--eps-break-scale,--eps-window-scale", eps_interval_scale,
                   "For --tool piecewise with --eps-interval-count. Must be >= 1.0; 1.0 gives linear distance ranks and larger values make earlier intervals larger")
        ->default_val(DEFAULT_EPS_INTERVAL_SCALE)
        ->check(CLI::PositiveNumber);

    app.add_option("--pv-cap-scale", pv_cap_scale,
                   "Required when --tool piecewise. PV diameter cap scale")
        ->check(CLI::PositiveNumber);

    app.add_option("--pv-min-separation", pv_min_separation,
                   "For --tool piecewise. Minimum inter-PV distance as a fraction of the final epsilon; 0 disables the constraint")
        ->default_val(0.0)
        ->check(CLI::NonNegativeNumber);

    app.add_option("--max-eps", maxeps,
                   "Required when --tool ph or apparent. Maximum epsilon")
        ->check(CLI::PositiveNumber);

    try
    {
        app.parse(argc, argv);
        return runTool(tool, filename, maxdim, threadnumber, eps_breaks, eps_interval_count,
                       eps_interval_scale, maxeps, pv_cap_scale, pv_min_separation, true, verbose,
                       adaptive, adaptive_target, trace_pwph,g_model,g_margin,g_memory_mib,cost_telemetry,quotient_mask,stop_after_window,ac_policy,ac_model,ac_auto_cutoff,geometry_cuts,h1_risk,h1_realizations,pack_mode,pack_model,payback_action,calibration_targets,extend_dump,extend_force_chain,extend_direct_probe,birth_repair,live_guard,live_guard_extend,extend_gate,extend_commit_once,extend_min_gain,extend_cost_rule,fiber_risk,absorb_pv,shrink_risky,shrink_search,joint_drops,fiber_full);
    }
    catch (const CLI::ParseError& err)
    {
        return app.exit(err);
    }
}
}

int main(int argc, char** argv)
{
    try
    {
        // Backward-compatible behavior: running with no arguments keeps the old interactive front end.
        if (argc == 1)
            return runInteractive();

        return runCommandLine(argc, argv);
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
