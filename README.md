# Artifact: A Run-Time Window Controller for Piecewise Persistent Homology

This branch holds the protocols, locks, raw measurements and analysis scripts behind the paper *A Run-Time Window Controller for Piecewise Persistent Homology* (FJCAI 2026 submission). The source code is on another branch of this repository (see Code).

- Every number in the paper can be traced to a file here.
- Nothing in `e10/` or `dev-stop-rule/` was edited after the runs; the files are copied unchanged.
- Protocol and deviation notes are written in Vietnamese.

## Code

| What | Where |
|---|---|
| Kernel used for every PwPH run | Commit [`f4733ef5e`](https://github.com/grade-11/iCC-morse-research/commit/f4733ef5eb5184d13c9a826c8625ed44fea31cf2) (full: `f4733ef5eb5184d13c9a826c8625ed44fea31cf2`) on branch [`work/pwph-h2-birth-repair-20261102`](https://github.com/grade-11/iCC-morse-research/tree/work/pwph-h2-birth-repair-20261102). The next commit, `3a39ba4c1`, changes only a code comment in `deps/H1EdgeRisk.hpp`. |
| Original PwPH implementation | [wilseypa/iCC](https://github.com/wilseypa/iCC) (X. Yang and P. A. Wilsey); this repository is built on it |
| Build | CMake Release, g++ 13.3.0, `-O3 -DNDEBUG` |
| Binary | `e10/bin/B_X31`, sha256 `9a0f76aa…e9b311d`, built from `f4733ef5e` |

All 450 PwPH runs of E10 used `B_X31`; the run arguments are in `e10/RUNS.jsonl`. `e10/bin/B_X22` (commit `0746d0815`) is pinned in the lock but was not used by any E10 run.

The reference implementation is Ripser 0.6.15 (Python package, one thread), run through `e10/ripser_run.py`.

## Names used in the paper and in the data

| Paper | Arm code in `e10/` |
|---|---|
| CTRL (the controller) | `UNI` |
| ALL | `CTRL` |
| EXTEND | `PAPER` |
| NoQ | `NOQ` |
| authors' schedule, H2 with filters | `AUTHOR-F` |
| authors' schedule, H2 without filters / H1 | `AUTHOR` |
| Ripser | `RIPSER` |

## Layout

```
e10/            evaluation cohort (Section IV): protocol, locks, cohort, raw runs, analyses
e9-protocol/    PROTOCOL.txt of E9, which e10/PROTOCOL.txt refers to ("như E9" = "as in E9")
dev-stop-rule/  development data behind the stop-rule sensitivity sentence (Section III-D)
paper/          scripts added for the paper: Table II memory/CPU columns, Fig. 3 data, per-cloud sensitivity
```

## From the paper to the files

**Protocol and locks (Section IV-A)**
- Protocol: `e10/PROTOCOL.txt`.
- Hash lock of the scripts, binaries and parameters, written before the cohort was generated: `e10/LOCK.json`.
- Cohort hashes: `e10/COHORT-LOCK.json`. Clouds: `e10/cohort/*.csv`. Cuts per cloud: `e10/COHORT.json`.
- Generators: `e10/e10_generate.py`, `e10/new_shapes.py`, and `e10/authors-generator/`, the unmodified sphere generator of the original authors.
- Deviations from the protocol: `e10/DEVIATIONS.txt`.

**Raw measurements**
- One line per run, with wall-clock seconds, user plus system CPU seconds and peak RSS: `e10/RUNS.jsonl`.
- Full stdout and stderr of every run: `e10/RAW/<run>/`.
- Mains-power log during timing: `e10/power.log`.

**Table II and Table III**
- Corrected analysis (`log t ~ arm + block`), which the paper reports: `e10/e10_corrected.py` → `e10/SUMMARY-CORRECTED.json`.
- Locked analysis, kept as run: `e10/e10_analyze.py` → `e10/SUMMARY.json`.
- Print the table rows: `python3 e10/e10_tables.py`, which uses the corrected analysis by default (pass `SUMMARY.json` for the locked one).
- Why the locked model is rank-deficient on two H2 clouds: `e10/e10_identifiability.py` → `e10/IDENTIFIABILITY.json`.
- Other sensitivity analyses (locked model on the identified clouds, per-cloud medians): `e10/e10_sensitivity.py` → `e10/SENSITIVITY.json`.
- Mem. and CPU columns: `python3 paper/artifact/cpu_memory_table.py`, output in `paper/artifact/CPU-MEMORY-TABLE.txt`. CPU time is measured, never derived from wall time.

**Correctness (Section IV-B)**
- Agreement with Ripser at float32, H1 and H2: `e10/e10_check.py` → `e10/CHECKS.json`.
- Equality with the uncut run of the same kernel: `e10/e10_uncut_equality.py` → `e10/UNCUT-EQUALITY.json`.
- H0 is not checked.

**Fig. 3:** `paper/figures/make_pd_data.py`. It reruns Ripser and asserts equality before writing the `.dat` files.

**Stop-rule sensitivity (Section III-D)**
- `dev-stop-rule/sensitivity_beta.py` → `dev-stop-rule/SENSITIVITY-BETA.json`.
- Per cloud: `paper/artifact/stop_rule_sensitivity_per_cloud.py`.
- How the replay works (only the stop is replayed; the times are measured): `paper/artifact/README.md`.
- The development runs were made with binary B_X30 (commit `2028f8e56`), whose run script is `dev-stop-rule/extend_oracle_v5.py`. B_X31 adds the margin, which these files evaluate offline.

## Running the scripts

You need Python 3 with numpy, plus `ripser==0.6.15` for `make_pd_data.py` and `e10_check.py`.

The scripts in `paper/` find `e10/` and `dev-stop-rule/` relative to themselves. To point them elsewhere, set `PWPH_E10` or `PWPH_DEV`.

The sealed scripts keep the absolute paths of the machine they ran on, so that their hashes still match `e10/LOCK.json`. Edit these lines before running them:

| File | Line with the absolute path |
|---|---|
| `e10/e10_common.py` | `PYTHON` |
| `e10/e10_generate.py` | `E9` |
| `e10/run_e10.sh` | `P` and the `cd` line |
| `dev-stop-rule/extend_oracle.py` | the `sys.path` entry for E8 |

The analysis scripts (`e10_analyze.py`, `e10_corrected.py`, `e10_tables.py`, `e10_sensitivity.py`, `e10_identifiability.py`, `e10_uncut_equality.py`) read only files in `e10/` and need no change.

Rerunning the timing (`e10_time.py`) needs the binaries; on the machine of Section IV-A the 596 runs took about two hours of wall-clock time. Times will not repeat exactly.
