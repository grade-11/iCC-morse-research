# Scripts added for the paper

All scripts here only read existing data. They write nothing outside this folder and start no PwPH or Ripser run.

## `cpu_memory_table.py` → `CPU-MEMORY-TABLE.txt`

This script recomputes the Mem. and CPU columns of Table II from the evaluation's `SUMMARY-CORRECTED.json`.

**Per cloud and arm**
- It takes the median, over the MEASURED and COMPLETE runs, of:
  - peak RSS;
  - CPU time, which is user plus system time from `wait4` on the whole process.
- Most arms have three runs per cloud, one per block. NoQ has a single run, in block 1.

**Per row**
- Each row is the GM over clouds of (median of CTRL) / (median of the arm).
- Nothing is derived from wall time or thread count.
- PwPH ran with 4 OpenMP threads; Ripser runs on one.

**Name mapping**

| Paper | Arm code in the data |
|---|---|
| CTRL | `UNI` |
| ALL | `CTRL` |
| EXTEND | `PAPER` |
| authors (H2) | `AUTHOR-F` |
| authors (H1) | `AUTHOR` |

The file also prints:
- the CPU ratio to Ripser for each cloud size;
- the capability runs (n = 400).

## `stop_rule_sensitivity_per_cloud.py` → `STOP-RULE-SENSITIVITY-PER-CLOUD.json`

This script gives a per-cloud view of the sensitivity sentence in Section III-D:

> fixing β = 1 or dropping the margin changed the GM time by at most 1.4%, although single clouds moved by up to 31%.

It reuses `sensitivity_beta.py` and `extend_oracle_table.py` in `development-stop-rule/` unchanged. It uses the data of that experiment's step 3: `RUNS-extend-oracle-v5.jsonl` and `extend-oracle-v5-runs/`, with binary `B_X30` on the ladder q0.1 … q0.9.

### What was run (development, 2026-11-05; no new runs here)

For each development cloud (96 for H1, 17 for H2), these runs were made with the same binary and ladder:

| Run | Flags | What it does |
|---|---|---|
| `LOG` | `--extend-all --extend-worth-log` | Keeps the frame and commits at every admitted cut. At each commit it logs the inputs of the stop rule: the work counts, per-window times, preparation time and collapse time. |
| `K<k>`, for k = 0 … K−1 | `--extend-all --extend-commit-limit k` | Commits at the first k admitted cuts as LOG does, then only extends the frame to E. K is the number of commits of LOG. |
| `ALL` | `--extend-all` | The same as LOG without logging. Its time stands for k = K. |

### How the "time" of a variant is obtained

The time is **measured, not estimated**. Only the stopping point is replayed.

1. **The decision is replayed.** The rule variant is evaluated offline on the logged inputs of LOG at each of its commits:
   - chosen rule: β clamped to [1, 3], with margin;
   - variant 1: β = 1;
   - variant 2: no margin.

   The variant stops at the first commit where the predicted saving is not positive. If that is the (j+1)-th commit, the variant's path is "commit j times, then extend to E".
2. **The time is measured.** That path is exactly what run `K<j>` executed, so the variant's time is the measured wall-clock time of `K<j>`. If the variant never stops, its time is that of `ALL`.

### Why no unlogged state is needed

A stop rule can only:
- accept the same commit that LOG made, or
- stop for good.

Up to its stop, a variant therefore follows LOG's path: same windows, same supports, same commits. After the stop, it follows the path of a measured `K<j>` run. Every path a variant can take is one of the measured runs, so no state outside the logs is needed.

### Limits of this replay

- The decisions use the timings logged by LOG, a different process from the `K<j>` run whose time is charged. In a real run, the rule would see that run's own timings. Timing noise can therefore move a decision; in step 3b, 24 of 96 H1 clouds changed their number of commits over three repetitions.
- Each `K<j>` is a single measured run, so the per-cloud ratios include run-to-run noise. This is one reason single clouds can move by up to 31% while the GM moves by at most 1.4%.
- The binary is `B_X30` (commit 2028f8e56). The paper's binary `B_X31` (commit f4733ef5e) adds the margin m, which is evaluated offline here, and computes only the work the rule uses.
- One H2 cloud (`H2-author-s2-n400-d0`) has no `K0` run. A variant that would stop at the first commit there is reported as missing. With these three variants this does not occur.

### Results (per-cloud ratio = variant time / chosen-rule time)

| Variant | Dim | Clouds | Decision changed | GM | Per-cloud range | Largest change |
|---|---|---|---|---|---|---|
| β = 1 | H1 | 96 | 3 | 0.998 | [0.918, 1.000] | TEST-P2b-rp2-veronese-r5-n384-d1, 0.918 |
| β = 1 | H2 | 17 | 2 | 1.011 | [1.000, 1.190] | H2-author-s2-n400-d0, 1.190 |
| no margin | H1 | 96 | 19 | 0.987 | [0.694, 1.237] | TEST-P2b-hopf-link-r3-n256-d1, 0.694 |
| no margin | H2 | 17 | 8 | 1.013 | [0.917, 1.201] | H2-author-s2-n200-d0, 1.201 |

The aggregates of `SENSITIVITY-BETA.json` are reproduced exactly: GM of T(stop)/base 0.9366 (H1) and 0.9344 (H2) for the chosen rule. The largest GM change is 1.32%.
