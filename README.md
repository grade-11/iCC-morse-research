# Artifact: A Run-Time Window Controller for Piecewise Persistent Homology

This branch holds the protocols, locks, raw measurements and analysis scripts behind the paper *A Run-Time Window Controller for Piecewise Persistent Homology* (FJCAI 2026 submission). The source code is on another branch of this repository (see Code).

- Every number in the paper can be traced to a file here.
- Nothing in `evaluation/` or `development-stop-rule/` was edited after the runs; the files are copied unchanged.
- Protocol and deviation notes are written in Vietnamese.

**Internal code names.** Inside the sealed folders, the evaluation keeps the code name it had when it was run, `E10`. It appears:
- in file names such as `e10_analyze.py`;
- in cloud names such as `E10-H2-author-s2-n200`;
- in the protocol.

These names are hashed by `evaluation/LOCK.json`, so they are left as they are. Earlier code names (E8, E9) occur only where the protocol refers to the framework it builds on, which is in `evaluation-base-protocol/`.

## Code

| What | Where |
|---|---|
| Kernel used for every PwPH run | Commit [`f4733ef5e`](https://github.com/grade-11/iCC-morse-research/commit/f4733ef5eb5184d13c9a826c8625ed44fea31cf2) (full: `f4733ef5eb5184d13c9a826c8625ed44fea31cf2`) on branch [`work/pwph-h2-birth-repair-20261102`](https://github.com/grade-11/iCC-morse-research/tree/work/pwph-h2-birth-repair-20261102). The next commit, `3a39ba4c1`, changes only a code comment in `deps/H1EdgeRisk.hpp`. |
| Original PwPH implementation | [wilseypa/iCC](https://github.com/wilseypa/iCC) (X. Yang and P. A. Wilsey); this repository is built on it |
| Build | CMake Release, g++ 13.3.0, `-O3 -DNDEBUG` |
| Binary | `evaluation/bin/B_X31`, sha256 `9a0f76aa…e9b311d`, built from `f4733ef5e` |

All 450 PwPH runs of the evaluation used `B_X31`; the run arguments are in `evaluation/RUNS.jsonl`. `evaluation/bin/B_X22` (commit `0746d0815`) is pinned in the lock but was not used by any run.

The reference implementation is Ripser 0.6.15 (Python package, one thread), run through `evaluation/ripser_run.py`.

## Names used in the paper and in the data

| Paper | Arm code in the data |
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
evaluation/                evaluation cohort (Section IV): protocol, locks, cohort, raw runs, analyses
evaluation-base-protocol/  the protocol that evaluation/PROTOCOL.txt builds on ("như E9" = "as in E9")
development-stop-rule/     development data behind the stop-rule sensitivity sentence (Section III-D)
paper/                     scripts written for the paper: Table II memory/CPU columns, Fig. 3 data, per-cloud sensitivity
```

## From the paper to the files

**Protocol and locks (Section IV-A)**
- Protocol: `evaluation/PROTOCOL.txt`.
- Hash lock of the scripts, binaries and parameters, written before the cohort was generated: `evaluation/LOCK.json`.
- Cohort hashes: `evaluation/COHORT-LOCK.json`. Clouds: `evaluation/cohort/*.csv`. Cuts per cloud: `evaluation/COHORT.json`.
- Generators: `evaluation/e10_generate.py`, `evaluation/new_shapes.py`, and `evaluation/authors-generator/`, the unmodified sphere generator of the original authors.
- Deviations from the protocol: `evaluation/DEVIATIONS.txt`.

**Raw measurements**
- One line per run, with wall-clock seconds, user plus system CPU seconds and peak RSS: `evaluation/RUNS.jsonl`.
- Full stdout and stderr of every run: `evaluation/RAW/<run>/`.
- Mains-power log during timing: `evaluation/power.log`.

**Table II and Table III**
- Corrected analysis (`log t ~ arm + block`), which the paper reports: `evaluation/e10_corrected.py` → `evaluation/SUMMARY-CORRECTED.json`.
- Locked analysis, kept as run: `evaluation/e10_analyze.py` → `evaluation/SUMMARY.json`.
- Print the table rows: `python3 evaluation/e10_tables.py`, which uses the corrected analysis by default (pass `SUMMARY.json` for the locked one).
- Why the locked model is rank-deficient on two H2 clouds: `evaluation/e10_identifiability.py` → `evaluation/IDENTIFIABILITY.json`.
- Other sensitivity analyses (locked model on the identified clouds, per-cloud medians): `evaluation/e10_sensitivity.py` → `evaluation/SENSITIVITY.json`.
- Mem. and CPU columns: `python3 paper/artifact/cpu_memory_table.py`, output in `paper/artifact/CPU-MEMORY-TABLE.txt`. CPU time is measured, never derived from wall time.

**Correctness (Section IV-B)**
- Agreement with Ripser at float32, H1 and H2: `evaluation/e10_check.py` → `evaluation/CHECKS.json`.
- Equality with the uncut run of the same kernel: `evaluation/e10_uncut_equality.py` → `evaluation/UNCUT-EQUALITY.json`.
- H0 is not checked.

**Fig. 3:** `paper/figures/make_pd_data.py`. It reruns Ripser and asserts equality before writing the `.dat` files.

**Stop-rule sensitivity (Section III-D)**
- `development-stop-rule/sensitivity_beta.py` → `development-stop-rule/SENSITIVITY-BETA.json`.
- Per cloud: `paper/artifact/stop_rule_sensitivity_per_cloud.py`.
- How the replay works (only the stop is replayed; the times are measured): `paper/artifact/README.md`.
- The development runs were made with binary B_X30 (commit `2028f8e56`), whose run script is `development-stop-rule/extend_oracle_v5.py`. B_X31 adds the margin, which these files evaluate offline.

## Running the scripts

You need Python 3 with numpy, plus `ripser==0.6.15` for `make_pd_data.py` and `e10_check.py`.

The scripts in `paper/` find `evaluation/` and `development-stop-rule/` relative to themselves. To point them elsewhere, set `PWPH_EVAL` or `PWPH_DEV`.

The sealed scripts keep the absolute paths of the machine they ran on, so that their hashes still match `evaluation/LOCK.json`. Edit these lines before running them:

| File | Line with the absolute path |
|---|---|
| `evaluation/e10_common.py` | `PYTHON` |
| `evaluation/e10_generate.py` | `E9`, the earlier cohort checked for overlapping seeds |
| `evaluation/run_e10.sh` | `P` and the `cd` line |
| `development-stop-rule/extend_oracle.py` | the `sys.path` entry |

The analysis scripts (`e10_analyze.py`, `e10_corrected.py`, `e10_tables.py`, `e10_sensitivity.py`, `e10_identifiability.py`, `e10_uncut_equality.py`) read only files in `evaluation/` and need no change.

Rerunning the timing (`e10_time.py`) needs the binaries; on the machine of Section IV-A the 596 runs took about two hours of wall-clock time. Times will not repeat exactly.
