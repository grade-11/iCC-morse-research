#pragma once

#include <algorithm>
#include <vector>

// Window-level pilot heuristic, not the proposed cost-benefit function G.
// More contraction permits a wider next interval. No future PH or timers enter.
inline double chooseAdaptiveEpsilon(const std::vector<double>& grid,
                                    double current,
                                    double reduction,
                                    double target_reduction)
{
    const double width = std::clamp(0.5 * reduction / target_reduction, 0.2, 0.5);
    auto next = std::upper_bound(grid.begin(), grid.end(), current);
    if (next == grid.end())
        return current;
    double chosen = *next; // Always advance, including a final gap > width.
    for (; next != grid.end() && *next <= current + width + 1e-12; ++next)
        chosen = *next;
    return chosen;
}
