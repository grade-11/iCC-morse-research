"""Tables for extend_oracle.py version 4 (RUNS-extend-oracle-v4.jsonl, extend-oracle-v4-runs/).

T(k): time of "commit k times, then only extend" (K<k>; k = K is ALL), per cloud and ladder. Oracle: min over k
(biased low by noise: a minimum over K + 1 noisy times).
A stop rule commits until its first "extend"; if that is at the (j+1)-th commit decision, its cost is T(j), which
the K series measures for any rule computed from the LOG run alone. Rules (all with the LOG run's [worth] line at
the cut and the windows before it):
  B_X26   (also B_X27) one rate t = sum time / sum work over all windows; commit iff C + prep + t*Wsel(end) < t*(W0(end) - W0(r))
          (prep: the latest build's label-distance and edge time). This is --extend-rule in B_X26.
  noprep  the same without prep.
  B_X24   build rate tb for a new window, extension rate te: C + tb*Wsel(next) + te*(Wsel(end) - Wsel(next)).
  never   never commit (K0); always: commit at every admitted cut (ALL).
Usage: python3 extend_oracle_table.py [H1|H2]"""
import collections
import json
import math
import statistics as st
import sys
from pathlib import Path

R = Path(__file__).resolve().parent
RULES = ("B_X26", "noprep", "B_X24", "power", "chain", "margin", "local", "localm", "localr")
V5 = "v5" in sys.argv[1:]


def local_beta(obs):
    """beta from how the per-unit rate of the first window grows along its own extensions: rate_i = sec_i / dW_i at
    the midpoint of step i; least squares of log rate on log W gives beta - 1. Clamped to [1, 3]; 1 with fewer than
    three steps. Uses only what the program has seen before the decision."""
    pts = []
    for i, (lo, hi, sec, build) in enumerate(obs):
        if i and build:
            break
        if hi > lo and sec > 0:
            pts.append((math.log((lo + hi) / 2), math.log(sec / (hi - lo))))
    if len(pts) < 3:
        return 1.0
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    slope = sum((x - mx) * (y - my) for x, y in pts) / sxx if sxx else 0.0
    return min(3.0, max(1.0, 1.0 + slope))


def local_cost(obs, beta, w_from, w_to):
    """Cost of covering work w_from..w_to when the rate grows as W**(beta-1), scaled to the latest window: an
    extension step [lo, hi] has its average rate at its midpoint; a build [0, hi] has average rate r(hi)/beta."""
    lo, hi, sec, build = obs[-1]
    if build or lo == 0:
        r_hi = beta * sec / hi
    else:
        r_hi = sec / (hi - lo) * (hi / ((lo + hi) / 2)) ** (beta - 1)
    return r_hi * hi / beta * ((w_to / hi) ** beta - (w_from / hi) ** beta)


def chain_beta(obs):
    """beta from the first window's own growth: cumulative time against W over its build and extensions
    (least squares in log-log), clamped to [1, 2]; 1 with fewer than three points."""
    pts, total = [], 0.0
    for lo, hi, sec, build in obs:
        if build and pts:
            break
        total += sec
        pts.append((math.log(hi), math.log(total)))
    if len(pts) < 3:
        return 1.0
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    beta = sum((x - mx) * (y - my) for x, y in pts) / sxx if sxx else 1.0
    return min(2.0, max(1.0, beta))


def fit_power(obs):
    """obs: (W_lo, W_hi, seconds, build?) per window so far. A window built to W costs t*W**beta and an
    extension from W_lo to W_hi costs t*(W_hi**beta - W_lo**beta); beta on a grid in [1, 2], t closed form,
    least squares on relative error."""
    best = (float("inf"), 1.0, 0.0)
    for k in range(21):
        beta = 1.0 + 0.05 * k
        g = [(hi ** beta - (0.0 if build else lo ** beta), sec) for lo, hi, sec, build in obs if hi > lo]
        if not g:
            return 1.0, 0.0
        # minimise sum((t*g/s - 1)^2): t = sum(g/s) / sum((g/s)^2)
        num = sum(x / y for x, y in g)
        den = sum((x / y) ** 2 for x, y in g)
        t = num / den if den else 0.0
        err = sum((t * x / y - 1) ** 2 for x, y in g)
        if err < best[0]:
            best = (err, beta, t)
    return best[1], best[2]


def gmean(xs):
    xs = [x for x in xs if x > 0 and math.isfinite(x)]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float("nan")


def fields(line):
    return {k: v for k, v in (t.split("=", 1) for t in line.split()[1:] if "=" in t)}


def decisions(path):
    """Per commit of the LOG run, the predicted saving (extend cost - commit cost) under each rule."""
    bs = bw = es = ew = 0.0
    prep = 0.0
    obs = []
    out = []
    last_ratio = 1.0  # the latest window's time over its prediction from the pooled rate before it
    local_ratio = 1.0  # the latest window's time over its prediction from the growing-rate model before it
    for line in path.read_text().splitlines():
        if line.startswith("[worth-window] "):
            f = fields(line)
            s, dw = float(f["core_s"]), float(f["W_hi"]) - float(f["W_lo"])
            new = (float(f["W_lo"]), float(f["W_hi"]), s, f["mode"] != "extend")
            if obs and dw > 0:
                pred = local_cost(obs, local_beta(obs), 0.0 if new[3] else new[0], new[1])
                if pred > 0:
                    local_ratio = s / pred
            obs.append(new)
            if bw + ew > 0 and dw > 0:
                last_ratio = s / ((bs + es) / (bw + ew) * dw)
            if f["mode"] != "extend":
                prep = float(f.get("prep_s", 0.0))
            if dw > 0:
                if f["mode"] == "extend":
                    es, ew = es + s, ew + dw
                else:
                    bs, bw = bs + s, bw + dw
        elif line.startswith("[worth] "):
            f = {k: float(v) for k, v in fields(line).items()}
            if f["committed"] == 0:
                continue
            t = (bs + es) / (bw + ew) if bw + ew else 0.0
            tb = bs / bw if bw else 0.0
            te = es / ew if ew else tb
            extend_pooled = t * (f["W0_end"] - f["W0_r"])
            beta, tp = fit_power(obs)
            bc = chain_beta(obs)
            g = [(hi ** bc - (0.0 if build else lo ** bc), sec) for lo, hi, sec, build in obs if hi > lo]
            tc = sum(sec for _, sec in g) / sum(x for x, _ in g)
            lb = local_beta(obs)
            ext_local = local_cost(obs, lb, f["W0_r"], f["W0_end"])
            commit_local = f["C_s"] + prep + local_cost(obs, lb, 0.0, f["Wsel_end"])
            out.append({
                "local": ext_local - commit_local,
                "localm": max(1.0, last_ratio) * ext_local - commit_local,
                "localr": max(1.0, local_ratio) * ext_local - commit_local,
                "margin": max(1.0, last_ratio) * extend_pooled - (f["C_s"] + prep + t * f["Wsel_end"]),
                "chain": tc * (f["W0_end"] ** bc - f["W0_r"] ** bc) - (f["C_s"] + prep + tc * f["Wsel_end"] ** bc),
                "power": tp * (f["W0_end"] ** beta - f["W0_r"] ** beta) - (f["C_s"] + prep + tp * f["Wsel_end"] ** beta),
                "B_X26": extend_pooled - (f["C_s"] + prep + t * f["Wsel_end"]),
                "noprep": extend_pooled - (f["C_s"] + t * f["Wsel_end"]),
                "B_X24": te * (f["W0_end"] - f["W0_r"])
                - (f["C_s"] + tb * f["Wsel_next"] + te * (f["Wsel_end"] - f["Wsel_next"])),
            })
    return out


def main():
    which = next((a for a in sys.argv[1:] if a in ("H1", "H2")), None)
    if V5:  # B_X30, one session (extend_oracle_v5.py); PAPER-B22 is under ladder ONCE; n400 has no K0
        rows = [json.loads(l) for l in open(R / "RUNS-extend-oracle-v5.jsonl")]
        rows = [dict(r, ladder="PAPER") if r["arm"] == "PAPER-B22" else r for r in rows if r["arm"] not in ("AUTHOR", "AUTHOR-F")]
        base = lambda case: R / "extend-oracle-v5-runs"  # noqa: E731
    else:
        # The n400 cloud was rerun in one session after WSL shut down (oracle_n400_rerun.py); its own files.
        rows = [json.loads(l) for l in open(R / "RUNS-extend-oracle-v4.jsonl")]
        extra = R / "RUNS-extend-oracle-v4-n400.jsonl"
        if extra.exists():
            rows += [json.loads(l) for l in open(extra)]
        base = lambda case: R / ("extend-oracle-v4-n400-runs" if case == "H2-author-s2-n400-d0" else "extend-oracle-v4-runs")  # noqa: E731
    runs = collections.defaultdict(dict)
    for r in rows:
        runs[(r["exp"], r["case"], r["ladder"])][r["arm"]] = r
    sec = lambda r: r["seconds"] if r and r["status"] == "COMPLETE" else float("inf")  # noqa: E731
    for exp in ("H1", "H2"):
        if which and exp != which:
            continue
        bad = sorted({(c, l, a) for (e, c, l), arms in runs.items() if e == exp for a, r in arms.items()
                      if r["status"] != "COMPLETE" or not r.get("exact")})
        print(f"===== {exp}; runs not complete or not exact: {bad[:10]}")
        for ladder in ("EVEN9", "P9"):
            keys = [k for k in runs if k[0] == exp and k[2] == ladder and "ALL" in runs[k]]
            if not keys:
                continue
            ratio = collections.defaultdict(list)
            counts = collections.defaultdict(collections.Counter)
            truth = collections.defaultdict(list)
            for key in sorted(keys):
                a = runs[key]
                paper = runs.get((exp, key[1], "PAPER"), {}).get("PAPER-B22")
                ratio["ALL / B22 (extend every window vs E8 controller)"].append(sec(a["ALL"]) / sec(a["B22"]))
                ratio["RULE (B_X26, run) / ALL"].append(sec(a["RULE"]) / sec(a["ALL"]))
                ratio["RULE / B22"].append(sec(a["RULE"]) / sec(a["B22"]))
                counts["ALL"][a["ALL"].get("commits")] += 1
                counts["RULE"][a["RULE"].get("commits")] += 1
                if paper:
                    ratio["B22 / PAPER-B22"].append(sec(a["B22"]) / sec(paper))
                    ratio["ALL / PAPER-B22"].append(sec(a["ALL"]) / sec(paper))
                    ratio["RULE / PAPER-B22"].append(sec(a["RULE"]) / sec(paper))
                if "LOG" not in a:
                    continue
                K = a["LOG"]["commits"]
                T = [sec(a.get(f"K{k}")) for k in range(K)] + [sec(a["ALL"])]
                best = min(range(K + 1), key=lambda k: T[k])
                counts["oracle"][best] += 1
                ratio["oracle / ALL"].append(T[best] / T[K])
                ratio["never (K0) / ALL"].append(T[0] / T[K])
                if K:
                    ratio["commit once (K1) / ALL"].append(T[1] / T[K])
                if paper:
                    ratio["oracle / PAPER-B22"].append(T[best] / sec(paper))
                d = decisions(base(key[1]) / exp / key[1] / ladder / "LOG" / "stdout.txt")
                for rule in RULES:
                    stop = next((j for j, x in enumerate(d[:K]) if x[rule] <= 0), K)
                    counts[rule][stop] += 1
                    ratio[f"rule {rule} (offline, T at its stop) / ALL"].append(T[stop] / T[K])
                    ratio[f"rule {rule} (offline) / oracle"].append(T[stop] / T[best])
                    if paper:
                        ratio[f"rule {rule} (offline) / PAPER-B22"].append(T[stop] / sec(paper))
                for j, x in enumerate(d[:K]):
                    if math.isfinite(T[j]) and math.isfinite(T[j + 1]):
                        truth[j].append(((T[j] - T[j + 1]) / T[K], x))
            print(f"--- {exp} {ladder}: {len(keys)} clouds")
            for k, v in ratio.items():
                print(f"  {k:58s} {gmean(v):.3f}  (n={len(v)}, median {st.median(v):.3f})")
            for k, c in counts.items():
                print(f"  commits {k:8s}", dict(sorted(c.items(), key=lambda kv: (kv[0] is None, kv[0] or 0))))
            for j, v in sorted(truth.items()):
                line = f"    decision {j + 1}: committing pays (truth) {sum(t > 0 for t, _ in v)}/{len(v)}, " \
                       f"median true saving {st.median(t for t, _ in v):+.3f} of ALL; rule says commit: "
                line += ", ".join(f"{r} {sum(x[r] > 0 for _, x in v)} (agree {sum((t > 0) == (x[r] > 0) for t, x in v)})"
                                  for r in RULES)
                print(line)


if __name__ == "__main__":
    main()
