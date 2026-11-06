"""E10 exploratory sensitivity (written after the locked analysis, not part of it). In E10-H2-genus2-r3-n300 the
shuffled arm order happened to be the same in all three blocks, so arm and position are collinear for UNI, CTRL
and RIPSER, and the position slope is set by the single NOQ run; the locked OLS then gives CTRL/UNI 0.27 although
the medians differ by 6%. Here every contrast is recomputed (a) by OLS without the position term and (b) from
per-cloud medians, with the locked bootstrap seed and reps. Output: SENSITIVITY.json."""
import json
import math
import statistics as st
import sys
from collections import defaultdict

import numpy as np

import e10_analyze as A
from e10_common import R, load


def fit_nopos(rows):
    present = sorted({r["arm"] for r in rows} - {"UNI"})
    blocks = sorted({r["block"] for r in rows})
    X, y = [], []
    for r in rows:
        X.append([1.0] + [1.0 if r["arm"] == a else 0.0 for a in present]
                 + [1.0 if r["block"] == b else 0.0 for b in blocks[1:]])
        y.append(math.log(r["seconds"]))
    beta, *_ = np.linalg.lstsq(np.array(X), np.array(y), rcond=None)
    return {"UNI": 0.0, **{a: float(beta[1 + k]) for k, a in enumerate(present)}}


def fit_median(rows):
    med = {a: st.median(r["seconds"] for r in rows if r["arm"] == a) for a in {r["arm"] for r in rows}}
    return {a: math.log(v / med["UNI"]) for a, v in med.items()}


def main():
    plan = load(R / "LOCK.json")["analysis"]
    cohort = {c["case"]: c for c in load(R / "COHORT.json")["rows"]}
    runs = [json.loads(l) for l in (R / "RUNS.jsonl").read_text().splitlines()]
    family_of = {c: v["family"] for c, v in cohort.items()}
    tests = {"H2": [("H2_N", "NOQ"), ("H2_C", "CTRL"), ("H2_S", "AUTHOR-F")],
             "H1": [("H1_N", "NOQ"), ("H1_C", "CTRL"), ("H1_E", "PAPER"), ("H1_S", "AUTHOR")]}
    out = {}
    for method, fitter in (("locked_ols", A.fit), ("ols_without_position", fit_nopos), ("medians", fit_median)):
        rng = np.random.default_rng(plan["bootstrap_seed"])
        out[method] = {}
        for e, hs in tests.items():
            clouds = sorted(c for c, v in cohort.items() if v["experiment"] == e and v["role"] == "MAIN")
            per = {c: fitter([r for r in runs if r["case"] == c and r["stage"] == "MEASURED"
                              and r["status"] == "COMPLETE"]) for c in clouds}
            for name, den in hs:
                logs = {c: -per[c][den] for c in clouds if den in per[c]}
                res = A.summarize(f"UNI/{den}", logs, family_of, rng, plan["bootstrap_reps"])
                res["slow"] = sorted(f for f, v in res["families"].items() if v["gm"] > plan["slow_threshold"])
                out[method][name] = res
                print(method, name, res["name"], round(res["gm"], 4), [round(x, 3) for x in res["cloud_ci"]],
                      [round(x, 3) for x in res["family_ci"]], res["slow"], flush=True)
    (R / "SENSITIVITY.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
