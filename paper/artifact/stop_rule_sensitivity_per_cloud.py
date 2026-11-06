"""Per-cloud view of the stop-rule sensitivity reported in the paper (written 2026-11-06, offline, no new timing).
Reuses sensitivity_beta.py of the development experiment pwph-h2-extend-dev-20261102 unchanged: for each variant
(beta fixed to 1; no margin) and each development cloud, the measured time T(stop) of the K series at the decision
where that variant stops, divided by T(stop) of the chosen rule (beta clamped to [1,3], with margin). T(stop) is the MEASURED wall-clock time of the
development run that commits that many times and then only extends (K<j>, or ALL if the rule never stops); only the
stopping decision is replayed, from the logged rule inputs of the LOG run. See README.md. Also the GM of
T(stop)/base as in SENSITIVITY-BETA.json, to check that the aggregate is reproduced. Reads the development
experiment only; prints a table and writes JSON to the path given (default: stdout only)."""
import collections
import json
import math
import os
import sys
from pathlib import Path

# Development directory: $PWPH_DEV, else ../../dev-stop-rule (layout of the artifact branch), else the original.
_HERE = Path(__file__).resolve().parent
DEV = Path(os.environ.get("PWPH_DEV") or next(
    (p for p in (_HERE.parents[1] / "dev-stop-rule",) if p.is_dir()),
    "/home/alex_1709/projects/iCC-team/experiments/pwph-h2-extend-dev-20261102"))
sys.path.insert(0, str(DEV))
out_path = sys.argv[1] if len(sys.argv) > 1 else None
import sensitivity_beta as sb  # noqa: E402  (sets argv to v5 for extend_oracle_table)

t = sb.t


def stop_times(beta_fn, rule):
    t.local_beta = beta_fn
    runs = collections.defaultdict(dict)
    for line in open(DEV / "RUNS-extend-oracle-v5.jsonl"):
        r = json.loads(line)
        runs[(r["exp"], r["case"])][(r["ladder"], r["arm"])] = r
    res = {}
    for (exp, case), a in sorted(runs.items()):
        log = a.get(("P9", "LOG"))
        if not log or log["status"] != "COMPLETE":
            continue
        K = log["commits"]
        T = [a[("P9", f"K{k}")]["seconds"] if ("P9", f"K{k}") in a and a[("P9", f"K{k}")]["status"] == "COMPLETE"
             else math.inf for k in range(K)] + [a[("P9", "ALL")]["seconds"]]
        base = a[("ONCE", "PAPER-B22")] if exp == "H1" else a[("EVEN9", "B22")]
        d = t.decisions(DEV / "extend-oracle-v5-runs" / exp / case / "P9" / "LOG" / "stdout.txt")
        stop = next((j for j, x in enumerate(d[:K]) if x[rule] <= 0), K)
        res[(exp, case)] = (stop, T[stop], base["seconds"])
    return res


def main():
    chosen = stop_times(*sb.VARIANTS["beta [1,3] (chosen)"])
    out = {}
    for name in ("beta = 1", "beta [1,3], no margin"):
        var = stop_times(*sb.VARIANTS[name])
        out[name] = {}
        for exp in ("H1", "H2"):
            rows = []
            for (e, case), (s0, t0, base) in chosen.items():
                if e != exp or (e, case) not in var:
                    continue
                s1, t1, _ = var[(e, case)]
                if not (math.isfinite(t0) and math.isfinite(t1)):
                    continue
                rows.append(dict(case=case, stop_chosen=s0 + 1, stop_variant=s1 + 1, ratio=t1 / t0,
                                 chosen_over_base=t0 / base, variant_over_base=t1 / base))
            gm = lambda xs: math.exp(sum(map(math.log, xs)) / len(xs))
            changed = [r for r in rows if r["stop_chosen"] != r["stop_variant"]]
            worst = max(rows, key=lambda r: abs(math.log(r["ratio"])))
            out[name][exp] = dict(clouds=len(rows), decision_changed=len(changed),
                                  gm_ratio_to_chosen=gm([r["ratio"] for r in rows]),
                                  gm_chosen_over_base=gm([r["chosen_over_base"] for r in rows]),
                                  gm_variant_over_base=gm([r["variant_over_base"] for r in rows]),
                                  min_ratio=min(r["ratio"] for r in rows), max_ratio=max(r["ratio"] for r in rows),
                                  largest_change=worst, per_cloud=rows)
            o = out[name][exp]
            print(f"{name:22s} {exp}: {o['clouds']} clouds, decision changed on {o['decision_changed']}; "
                  f"GM variant/chosen {o['gm_ratio_to_chosen']:.4f}; per-cloud range [{o['min_ratio']:.3f}, "
                  f"{o['max_ratio']:.3f}]; largest {worst['case']} {worst['ratio']:.3f} "
                  f"(stop {worst['stop_chosen']} -> {worst['stop_variant']}); GM/base chosen "
                  f"{o['gm_chosen_over_base']:.4f}, variant {o['gm_variant_over_base']:.4f}")
    if out_path:
        Path(out_path).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
