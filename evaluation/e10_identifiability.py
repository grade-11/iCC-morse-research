"""E10 exploratory identifiability check (written 06/11/2026 after an external review, not part of the locked
analysis). For every cloud: rank of the locked design matrix (arm + position + block) and whether each arm
contrast against UNI is estimable (the contrast vector lies in the row space of X). Then, with the locked
bootstrap seed and reps:
  (a) the locked OLS on the clouds where the contrast is estimable;
  (b) UNI/NOQ from block 1 alone, the paired ratio of the two runs measured in that block, without any model.
Output: IDENTIFIABILITY.json (written once)."""
import json
import math
import sys

import numpy as np

import e10_analyze as A
from e10_common import R, load, write_new


def design(rows):
    present = sorted({r["arm"] for r in rows} - {"UNI"})
    blocks = sorted({r["block"] for r in rows})
    X = [[1.0] + [1.0 if r["arm"] == a else 0.0 for a in present] + [float(r["position"])]
         + [1.0 if r["block"] == b else 0.0 for b in blocks[1:]] for r in rows]
    return np.array(X), present


def estimable(X, k):
    c = np.zeros(X.shape[1])
    c[1 + k] = 1.0
    return bool(np.allclose(c @ np.linalg.pinv(X) @ X, c, atol=1e-8))


def main():
    plan = load(R / "LOCK.json")["analysis"]
    cohort = {c["case"]: c for c in load(R / "COHORT.json")["rows"]}
    runs = [json.loads(l) for l in (R / "RUNS.jsonl").read_text().splitlines()]
    family_of = {c: v["family"] for c, v in cohort.items()}
    tests = {"H2": [("H2_N", "NOQ"), ("H2_C", "CTRL"), ("H2_S", "AUTHOR-F")],
             "H1": [("H1_N", "NOQ"), ("H1_C", "CTRL"), ("H1_E", "PAPER"), ("H1_S", "AUTHOR")]}
    out = {"design": {}, "locked_ols_estimable_only": {}, "noq_block1_paired": {}}
    rng = np.random.default_rng(plan["bootstrap_seed"])
    for e, hs in tests.items():
        clouds = sorted(c for c, v in cohort.items() if v["experiment"] == e and v["role"] == "MAIN")
        per, est = {}, {}
        for c in clouds:
            rows = [r for r in runs if r["case"] == c and r["stage"] == "MEASURED" and r["status"] == "COMPLETE"]
            X, present = design(rows)
            s = np.linalg.svd(X, compute_uv=False)
            est[c] = {a: estimable(X, k) for k, a in enumerate(present)}
            out["design"][c] = dict(columns=X.shape[1], rank=int(np.linalg.matrix_rank(X)),
                                    condition=float(s[0] / s[-1]) if s[-1] > 0 else None, estimable=est[c])
            per[c] = A.fit(rows)
        for name, den in hs:
            logs = {c: -per[c][den] for c in clouds if den in per[c] and est[c].get(den)}
            res = A.summarize(f"UNI/{den}", logs, family_of, rng, plan["bootstrap_reps"])
            res["dropped"] = sorted(c for c in clouds if den in per[c] and not est[c].get(den))
            res["slow"] = sorted(f for f, v in res["families"].items() if v["gm"] > plan["slow_threshold"])
            out["locked_ols_estimable_only"][name] = res
            print("estimable-only", name, round(res["gm"], 4), [round(x, 3) for x in res["cloud_ci"]],
                  [round(x, 3) for x in res["family_ci"]], "slow", res["slow"], "dropped", res["dropped"], flush=True)
        logs = {}
        for c in clouds:
            b1 = {r["arm"]: r["seconds"] for r in runs if r["case"] == c and r["stage"] == "MEASURED"
                  and r["status"] == "COMPLETE" and r["block"] == 1}
            if "UNI" in b1 and "NOQ" in b1:
                logs[c] = math.log(b1["UNI"] / b1["NOQ"])
        res = A.summarize("UNI/NOQ block 1", logs, family_of, rng, plan["bootstrap_reps"])
        res["slow"] = sorted(f for f, v in res["families"].items() if v["gm"] > plan["slow_threshold"])
        out["noq_block1_paired"][f"{e}_N"] = res
        print("block1", e, round(res["gm"], 4), [round(x, 3) for x in res["cloud_ci"]],
              [round(x, 3) for x in res["family_ci"]], "slow", res["slow"], "clouds", res["clouds"], flush=True)
    deficient = sorted(c for c, d in out["design"].items() if d["rank"] < d["columns"])
    out["rank_deficient"] = deficient
    print("rank-deficient", deficient, flush=True)
    write_new(R / "IDENTIFIABILITY.json", out)


if __name__ == "__main__":
    sys.exit(main())
