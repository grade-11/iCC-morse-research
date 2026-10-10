#include "DistanceMatrix.hpp"
#include "Pipelines.hpp"
#include "readInput.hpp"

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
            const bool verbose)
{
    const std::string normalized_tool = normalizeTool(tool);

    validateRunOptions(normalized_tool, filename, maxdim, threadnumber, eps_breaks,
                       eps_interval_count, eps_interval_scale, maxeps);

    if (print_configuration && verbose)
        printRunConfiguration(normalized_tool, filename, maxdim, threadnumber, eps_breaks,
                              eps_interval_count, eps_interval_scale, maxeps, pv_cap_scale,
                              pv_min_separation);

    auto point_cloud = readInput::readCSV(filename);
    omp_set_num_threads(threadnumber);
    DistanceMatrix distance_matrix(point_cloud);
    decltype(point_cloud){}.swap(point_cloud);

    std::vector<double> resolved_eps_breaks;
    if (normalized_tool == PIECEWISE_PH_TOOL)
    {
        resolved_eps_breaks = resolveEpsilonBreaks(
            distance_matrix, eps_breaks, eps_interval_count, eps_interval_scale);

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
        .verbose = verbose};

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
                       eps_interval_scale, maxeps, pv_cap_scale, pv_min_separation, true, verbose);
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
