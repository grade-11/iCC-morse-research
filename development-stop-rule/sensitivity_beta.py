"""Review response (2026-11-05), offline only, no new timing. Two questions on the stop rule of B_X31 (localr):

1. Sensitivity to the clamp of beta and to the margin. From the step 3 K series (RUNS-extend-oracle-v5.jsonl,
   P9 ladder, the same evaluation as variant_choice.py): T(stop) / base, GM over H1 (base: PAPER-B22) and H2
   (base: B22-EVEN9), for beta clamped to [1, 2], [1, 3] (chosen), [1, 4], [1, inf), beta = 1, and without the
   margin. A stop at decision 1 on the n400 cloud (no K0 run) is reported as missing.
2. Stability of the decision across repetitions: in step 3b (RUNS-step3b.jsonl), the number of commits of RULE-P9
   per cloud over its rounds.
Output: SENSITIVITY-BETA.json."""
import collections
import json
import math
import sys
from pathlib import Path

sys.argv = sys.argv[:1] + ["v5"]
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R))
import extend_oracle_table as t  # noqa: E402

ORIG_BETA = t.local_beta


def clamped(lo, hi):
    def f(obs):
        pts = []
        for i, (a, b, sec, build) in enumerate(obs):
            if i and build:
                break
            if b > a and sec > 0:
                pts.append((math.log((a + b) / 2), math.log(sec / (b - a))))
        if len(pts) < 3:
            return 1.0
        mx = sum(x for x, _ in pts) / len(pts)
        my = sum(y for _, y in pts) / len(pts)
        sxx = sum((x - mx) ** 2 for x, _ in pts)
        slope = sum((x - mx) * (y - my) for x, y in pts) / sxx if sxx else 0.0
        return min(hi, max(lo, 1.0 + slope))
    return f


VARIANTS = {
    "beta [1,3] (chosen)": (clamped(1.0, 3.0), "localr"),
    "beta [1,2]": (clamped(1.0, 2.0), "localr"),
    "beta [1,4]": (clamped(1.0, 4.0), "localr"),
    "beta [1,inf)": (clamped(1.0, math.inf), "localr"),
    "beta = 1": (lambda obs: 1.0, "localr"),
    "beta [1,3], no margin": (clamped(1.0, 3.0), "local"),
}


def gm(xs):
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float("nan")


def main():
    runs = collections.defaultdict(dict)
    for line in open(R / "RUNS-extend-oracle-v5.jsonl"):
        r = json.loads(line)
        runs[(r["exp"], r["case"])][(r["ladder"], r["arm"])] = r
    out = {"sensitivity": {}, "stability": {}}
    for name, (beta_fn, rule) in VARIANTS.items():
        t.local_beta = beta_fn
        ratios, missing, stops = {"H1": [], "H2": []}, [], collections.Counter()
        for (exp, case), a in sorted(runs.items()):
            log = a.get(("P9", "LOG"))
            if not log or log["status"] != "COMPLETE":
                continue
            K = log["commits"]
            T = [a[("P9", f"K{k}")]["seconds"] if ("P9", f"K{k}") in a and a[("P9", f"K{k}")]["status"] == "COMPLETE"
                 else math.inf for k in range(K)] + [a[("P9", "ALL")]["seconds"]]
            base = a[("ONCE", "PAPER-B22")] if exp == "H1" else a[("EVEN9", "B22")]
            d = t.decisions(R / "extend-oracle-v5-runs" / exp / case / "P9" / "LOG" / "stdout.txt")
            stop = next((j for j, x in enumerate(d[:K]) if x[rule] <= 0), K)
            stops[f"{exp} stop at decision {stop + 1}" if stop < K else f"{exp} never stops"] += 1
            if not math.isfinite(T[stop]):
                missing.append(case)
                continue
            ratios[exp].append(T[stop] / base["seconds"])
        out["sensitivity"][name] = {"H1": gm(ratios["H1"]), "n_H1": len(ratios["H1"]), "H2": gm(ratios["H2"]),
                                    "n_H2": len(ratios["H2"]), "missing": missing}
        print(f"{name:24s} H1 {gm(ratios['H1']):.3f} (n={len(ratios['H1'])})  H2 {gm(ratios['H2']):.3f} "
              f"(n={len(ratios['H2'])})  missing {missing}")
    t.local_beta = ORIG_BETA

    commits = collections.defaultdict(list)
    for line in open(R / "RUNS-step3b.jsonl"):
        r = json.loads(line)
        if r["arm"] == "RULE-P9" and r["status"] == "COMPLETE":
            commits[(r["exp"], r["case"])].append(r["commits"])
    for exp in ("H1", "H2"):
        cases = {c: v for (e, c), v in commits.items() if e == exp}
        varying = {c: v for c, v in cases.items() if len(set(v)) > 1}
        out["stability"][exp] = {"clouds": len(cases), "rounds": sorted({len(v) for v in cases.values()}),
                                 "clouds_with_varying_commits": len(varying), "examples": dict(list(varying.items())[:8])}
        print(f"stability {exp}: {len(varying)}/{len(cases)} clouds changed their number of commits across rounds "
              f"{sorted({len(v) for v in cases.values()})}; e.g. {dict(list(varying.items())[:5])}")
    (R / "SENSITIVITY-BETA.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
