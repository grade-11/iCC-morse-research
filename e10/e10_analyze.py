"""E10 locked analysis (PROTOCOL.txt). Per cloud: OLS log t ~ arm + position + block over the complete measured
runs, reference UNI; equal-cloud geometric means; cloud and family-cluster bootstrap. Exploratory additions go
in a separate file."""
import json
import math
import sys
from collections import Counter, defaultdict

import numpy as np

from e10_common import ARMS, R, load, require, utc, verify_lock, write_new

REFERENCE = "UNI"


def fit(rows):
    """{arm: log time relative to UNI} from one cloud's complete measured runs."""
    present = sorted({r["arm"] for r in rows} - {REFERENCE})
    blocks = sorted({r["block"] for r in rows})
    X, y = [], []
    for r in rows:
        x = [1.0] + [1.0 if r["arm"] == a else 0.0 for a in present] + [float(r["position"])]
        x += [1.0 if r["block"] == b else 0.0 for b in blocks[1:]]
        X.append(x)
        y.append(math.log(r["seconds"]))
    beta, *_ = np.linalg.lstsq(np.array(X), np.array(y), rcond=None)
    out = {REFERENCE: 0.0}
    out.update({a: float(beta[1 + k]) for k, a in enumerate(present)})
    return out


def gm(logs):
    return math.exp(sum(logs) / len(logs))


def interval(samples):
    s = np.sort(np.array(samples))
    return [float(s[int(0.025 * len(s))]), float(s[int(math.ceil(0.975 * len(s))) - 1])]


def cloud_boot(logs, rng, reps):
    a = np.array(logs)
    return interval([math.exp(a[rng.integers(0, len(a), len(a))].mean()) for _ in range(reps)])


def family_boot(by_family, rng, reps):
    fams = sorted(by_family)
    out = []
    for _ in range(reps):
        pick = rng.integers(0, len(fams), len(fams))
        pooled = [v for k in pick for v in by_family[fams[k]]]
        out.append(math.exp(sum(pooled) / len(pooled)))
    return interval(out)


def summarize(name, logs, family_of, rng, reps):
    by_family = defaultdict(list)
    for c, v in logs.items():
        by_family[family_of[c]].append(v)
    fam = {f: dict(gm=gm(v), clouds=len(v)) for f, v in sorted(by_family.items())}
    return dict(name=name, clouds=len(logs), gm=gm(list(logs.values())), cloud_ci=cloud_boot(list(logs.values()), rng, reps),
                family_ci=family_boot(by_family, rng, reps), families=fam,
                faster_clouds=sum(v < 0 for v in logs.values()))


def main():
    lock = verify_lock()
    plan = lock["analysis"]
    require((R / "TIMING-END.json").exists(), "timing has not closed")
    require((R / "CHECKS.json").exists(), "correctness check has not run")
    cohort = {c["case"]: c for c in load(R / "COHORT.json")["rows"]}
    runs = [json.loads(l) for l in (R / "RUNS.jsonl").read_text().splitlines()]
    checks = load(R / "CHECKS.json")["runs"]
    mismatch = {(x["case"], x["arm"]) for x in checks if x["check"] != "PASS" and x["status"] == "COMPLETE"}
    family_of = {c: v["family"] for c, v in cohort.items()}
    rng = np.random.default_rng(plan["bootstrap_seed"])
    reps = plan["bootstrap_reps"]

    def decide(result):
        """Secondary: UNI against the per-dimension controller (and Extend in H1): FASTER, else NONINFERIOR (both upper bounds below 1.05), else
        NOT_CONFIRMED."""
        slow = sorted(f for f, v in result["families"].items() if v["gm"] > plan["slow_threshold"])
        result["slow_families_over_5pct"] = slow
        if result["cloud_ci"][1] < 1 and result["family_ci"][1] < 1 and len(slow) <= plan["max_slow_families"]:
            return "FASTER"
        bound = plan["noninferiority_bound"]
        if result["cloud_ci"][1] < bound and result["family_ci"][1] < bound:
            return "NONINFERIOR"
        return "NOT_CONFIRMED"

    def decide_noq(result):
        """Primary: UNI against the uncut computation. CONFIRMED iff both upper bounds are below 1 and at most
        max_slow_families families have GM above slow_threshold."""
        slow = sorted(f for f, v in result["families"].items() if v["gm"] > plan["slow_threshold"])
        result["slow_families_over_5pct"] = slow
        ok = result["cloud_ci"][1] < 1 and result["family_ci"][1] < 1 and len(slow) <= plan["max_slow_families"]
        return "CONFIRMED" if ok else "NOT_CONFIRMED"

    def decide_sphere(result):
        slow = sorted(f for f, v in result["families"].items() if v["gm"] > plan["slow_threshold"])
        result["slow_families_over_5pct"] = slow
        if result["cloud_ci"][1] >= 1:
            return "NOT_CONFIRMED"
        return "CONFIRMED" if not slow else "PARTIAL"

    summary = dict(utc=utc(), verdicts={}, hypotheses={}, descriptive={}, by_size={}, memory_gm_ratio={},
                   cpu_gm_ratio={}, per_cloud={}, correctness=load(R / "CHECKS.json")["counts"],
                   run_status={f"{k[0]}:{k[1]}:{k[2]}": v for k, v in
                               Counter((r["experiment"], r["stage"], r["status"]) for r in runs).items()})
    for experiment in ("H2", "H1"):
        clouds = sorted(c for c, v in cohort.items() if v["experiment"] == experiment and v["role"] == "MAIN")
        per_cloud, med_rss, med_cpu, med_s = {}, {}, {}, {}
        for case in clouds:
            rows = [r for r in runs if r["case"] == case and r["stage"] == "MEASURED" and r["status"] == "COMPLETE"]
            per_cloud[case] = fit(rows)
            arms = sorted({r["arm"] for r in rows})
            med_s[case] = {a: float(np.median([r["seconds"] for r in rows if r["arm"] == a])) for a in arms}
            med_rss[case] = {a: float(np.median([r["maxrss_kib"] for r in rows if r["arm"] == a])) / 1024 for a in arms}
            med_cpu[case] = {a: float(np.median([r["cpu_seconds"] for r in rows if r["arm"] == a])) for a in arms}

        def contrast(num, den, cs, clean=False):
            return {c: per_cloud[c][num] - per_cloud[c][den] for c in cs if num in per_cloud[c] and den in per_cloud[c]
                    and not (clean and ((c, num) in mismatch or (c, den) in mismatch))}

        spheres = [c for c in clouds if "AUTHOR" in cohort[c]["cuts"]]
        if experiment == "H2":
            main = ("H2_N", "UNI", "NOQ", clouds)
            secondary = [("H2_C", "UNI", "CTRL", clouds)]
            extra = [("H2_S", "UNI", "AUTHOR-F", spheres)]
            desc = [("CTRL/NOQ", "CTRL", "NOQ", clouds),
                    ("UNI/RIPSER", "UNI", "RIPSER", clouds), ("CTRL/RIPSER", "CTRL", "RIPSER", clouds),
                    ("CTRL/AUTHOR-F spheres", "CTRL", "AUTHOR-F", spheres),
                    ("AUTHOR-F/AUTHOR spheres", "AUTHOR-F", "AUTHOR", spheres),
                    ("UNI/AUTHOR spheres", "UNI", "AUTHOR", spheres),
                    ("AUTHOR-F/RIPSER spheres", "AUTHOR-F", "RIPSER", spheres)]
        else:
            main = ("H1_N", "UNI", "NOQ", clouds)
            secondary = [("H1_C", "UNI", "CTRL", clouds), ("H1_E", "UNI", "PAPER", clouds)]
            extra = [("H1_S", "UNI", "AUTHOR", spheres)]
            desc = [("UNI/PAPER", "UNI", "PAPER", clouds), ("CTRL/PAPER", "CTRL", "PAPER", clouds),
                    ("CTRL/NOQ", "CTRL", "NOQ", clouds), ("UNI/RIPSER", "UNI", "RIPSER", clouds),
                    ("CTRL/RIPSER", "CTRL", "RIPSER", clouds), ("PAPER/RIPSER", "PAPER", "RIPSER", clouds),
                    ("CTRL/AUTHOR spheres", "CTRL", "AUTHOR", spheres), ("PAPER/AUTHOR spheres", "PAPER", "AUTHOR", spheres)]
        name, num, den, cs = main
        res = summarize(f"{num}/{den}", contrast(num, den, cs), family_of, rng, reps)
        summary["verdicts"][name] = decide_noq(res)
        res["mismatch_free"] = summarize(f"{num}/{den} mismatch-free", contrast(num, den, cs, True), family_of, rng, reps)
        summary["hypotheses"][name] = res
        for name, num, den, cs in secondary:
            res = summarize(f"{num}/{den}", contrast(num, den, cs), family_of, rng, reps)
            summary["verdicts"][name] = decide(res)
            summary["hypotheses"][name] = res
        for name, num, den, cs in extra:
            res = summarize(f"{num}/{den}", contrast(num, den, cs), family_of, rng, reps)
            summary["verdicts"][name] = decide_sphere(res)
            summary["hypotheses"][name] = res
        summary["descriptive"][experiment] = [summarize(n, contrast(a, b, cs), family_of, rng, reps)
                                              for n, a, b, cs in desc if contrast(a, b, cs)]
        summary["by_size"][experiment] = {}
        for n in sorted({cohort[c]["n"] for c in clouds}):
            same = [c for c in clouds if cohort[c]["n"] == n]
            summary["by_size"][experiment][str(n)] = {f"{a}/{b}": gm(list(contrast(a, b, same).values()))
                                                      for _, a, b, _ in [main] + desc[:4] if contrast(a, b, same)}
        for store, table in ((summary["memory_gm_ratio"], med_rss), (summary["cpu_gm_ratio"], med_cpu)):
            store[experiment] = {f"{a}/{b}": gm([math.log(table[c][a] / table[c][b]) for c in cs
                                                 if a in table[c] and b in table[c]])
                                 for _, a, b, cs in [main] + desc if any(a in table[c] and b in table[c] for c in cs)}
        summary["per_cloud"][experiment] = [dict(case=c, family=family_of[c], n=cohort[c]["n"],
                                                 ratio_to_UNI={a: math.exp(v) for a, v in per_cloud[c].items()},
                                                 median_seconds=med_s[c], median_rss_mib=med_rss[c],
                                                 median_cpu_seconds=med_cpu[c]) for c in clouds]
    summary["capability"] = [dict(case=r["case"], family=r["family"], n=r["n"], arm=r["arm"], status=r["status"],
                                  seconds=r["seconds"], cpu_seconds=r["cpu_seconds"], maxrss_mib=r["maxrss_kib"] / 1024,
                                  check=next((x["check"] for x in checks if x["run_id"] == r["run_id"]), None))
                             for r in runs if r["stage"] == "CAPABILITY"]
    commits = defaultdict(Counter)
    for x in checks:
        if x["check"] in ("PASS", "PH_MISMATCH") and x["stage"] == "MEASURED":
            commits[f"{x['experiment']}:{x['arm']}"][x["commits"]] += 1
    summary["commits"] = {k: dict(sorted(v.items())) for k, v in commits.items()}
    write_new(R / "SUMMARY.json", summary)
    for k, res in summary["hypotheses"].items():
        print("E10", k, summary["verdicts"][k], res["name"], round(res["gm"], 4), res["cloud_ci"], res["family_ci"],
              "slow", res["slow_families_over_5pct"], flush=True)


if __name__ == "__main__":
    sys.exit(main())
