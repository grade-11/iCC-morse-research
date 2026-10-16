# iCC

Incremental Critical Cells (iCC) computes persistent-homology information for
dense Euclidean point clouds. The executable provides three concrete pipelines:

- `ph`: direct implicit Vietoris–Rips persistent homology.
- `piecewise`: piecewise persistent homology (PwPH) with pseudo-vertices (PVs).
- `apparent`: apparent-pair-only reduction, without augmenting-path matching.

This repository contains the C++ implementation used by the HPC lab at the
University of Cincinnati.

## Requirements

- A C++23 compiler
- CMake 3.22 or newer
- OpenMP
- CLI11 (CMake fetches CLI11 v2.6.2 when no installed package is found)

CGAL and Eigen are not required.

## Building

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The executable is written to `build/iCC`.

To build and run the core regression tests:

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Running

Running without arguments opens the interactive runner:

```bash
./build/iCC
```

Direct persistent homology:

```bash
./build/iCC --tool ph --file-name <input.csv> --max-dim <maxDim> --max-eps <epsilon> [-n <thread_count>] [-v]
```

PwPH with explicit epsilon breaks:

```bash
./build/iCC --tool piecewise --file-name <input.csv> --max-dim <maxDim> --eps-breaks <eps1> <eps2> ... --pv-cap-scale <scale> [--pv-min-separation <delta>] [-n <thread_count>] [-v]
```

PwPH with automatically generated epsilon breaks:

```bash
./build/iCC --tool piecewise --file-name <input.csv> --max-dim <maxDim> --eps-interval-count <count> [--eps-interval-scale <scale>] --pv-cap-scale <scale> [--pv-min-separation <delta>] [-n <thread_count>] [-v]
```

Apparent-pair-only reduction:

```bash
./build/iCC --tool apparent --file-name <input.csv> --max-dim <maxDim> --max-eps <epsilon> [-n <thread_count>] [-v]
```

The apparent pipeline prints deterministic `dimension, unmatched weight` rows.
It includes unmatched top-dimensional cofacets and does not print H0 cells.

The numeric tool aliases `1`, `2`, and `3` select `piecewise`, `ph`, and
`apparent`, respectively. Existing PH aliases remain accepted; the apparent
pipeline also accepts `apparent-pairs` and `apparentpairs`.

For every option and alias:

```bash
./build/iCC --help
```

## Epsilon schedules

`--eps-breaks` and `--eps-interval-count` are mutually exclusive. The aliases
`--eps-break-count` and `--eps-window-count` are accepted for the interval
count. `--eps-break-scale` and `--eps-window-scale` are accepted for the
interval scale.

Automatic scheduling sorts the distinct positive pairwise distances and chooses
the requested number of upper bounds from those ranks. The interval scale must
be finite and at least `1.0`; `1.0` gives linear distance ranks, while larger
values place more distinct distances in earlier windows. Selected distances are
advanced by one representable floating-point value because filtration admission
uses a strict upper bound.

`--pv-cap-scale` limits a new PV's diameter to `scale * eps_hi`, where `eps_hi`
is its formation-window upper bound. `--pv-min-separation` requires every pair
of PVs to have minimum cross-distance at least `delta * eps_max`; its default is
`0`, which disables this constraint.

The executable accepts any finite positive PV cap scale for experiments. The
quotient and single-coordinate safety arguments assume the cap does not exceed
the formation scale. Because admission is strict, use a scale below `1.0` to
avoid the equality boundary.

The cap alone does not guarantee preservation of the full PH diagram after a
quotient. A deterministic eight-point unit-S2 example at cap scale `0.85` loses
a future H1 interval with breaks `0.5 2.0`, including on untouched upstream
commit `734be3d`. The [provenance study](experiments/pwph2026-provenance-continue/README.md)
links the input, direct version comparisons and per-configuration results.

## Output and timing

Direct PH retains the existing per-dimension persistent-pair report. Nonverbose
PwPH writes a CSV-like window header followed by
`dimension, birth weight, death weight` rows. In verbose PH and PwPH output,
every persistence pair includes decoded birth-facet and death-cofacet labels.
Verbose PwPH additionally reports pipeline counts, PV statistics, and available
FFI diagnostics.

For PwPH windows containing PVs in dimensions 3 through 7, verbose mode reports
canonical-witness false-facet-identification statistics for the top interface.
Each report describes the complex rebuilt at that window's upper scale; it is
not restricted to cofacets whose weights lie strictly inside the recording
interval.

Every command-line run prints `run time = ...` in milliseconds. Input reading,
distance-matrix construction, and automatic epsilon-schedule generation happen
before this timer starts.

## Experimental no-quotient continuation

`--continue-engine f2-stream` selects an explicit F2 boundary-reduction engine
that retains reduced columns and live births between checkpoints. The control
mode `f2-rebuild` starts the same F2 engine afresh at every checkpoint. Both
support `--max-dim 2` (H0/H1), require an explicit all-zero `--quotient-mask`
for multiple windows, and reject quotient, G, adaptive v1, H2, trace and early
stop options. The default remains `native-rebuild` with the Morse pipeline.

```bash
./build/iCC --tool piecewise --file-name tests/data/irregular5.csv \
 --max-dim 2 --threads 4 --eps-breaks 1.0 1.5 2.0 --pv-cap-scale .85 \
 --quotient-mask 0 0 --continue-engine f2-stream --verbose
```

The F2 engine sorts input edges once per engine, generates triangles when their
last edge arrives, and retains pivot columns by stable edge IDs. Its reduction
is currently serial; `--threads` still bounds the process's configured worker
count. H0 prefix output (including zero bars) is available with `--verbose`.
Every startup, enumeration, checkpoint, reduction and logging cost is included
in whole-process measurements. This experimental kernel does not implement
continuation of native Morse interface workspaces or migration through quotient.
