"""CPU and memory columns of Table II (written 2026-11-06). Reads only SUMMARY-CORRECTED.json of the evaluation.
For each dimension and arm: GM over clouds of (per-cloud median of CTRL) / (per-cloud median of the arm), where the
medians are over the MEASURED, COMPLETE runs of that cloud (three blocks; NoQ: block 1 only). CPU time is user plus
system time from wait4 of the whole process, as recorded; nothing is derived from wall time or thread count.
Paper names: CTRL = arm UNI, ALL = arm CTRL, EXTEND = arm PAPER; authors = AUTHOR-F in H2, AUTHOR in H1.
Also the CPU ratio to Ripser by cloud size, and the capability runs. Usage: python3 cpu_memory_table.py"""
import json
import math
import os
from pathlib import Path

# Evaluation directory: $PWPH_EVAL, else ../../evaluation (layout of the artifact branch), else the original location.
_HERE = Path(__file__).resolve().parent
E10 = Path(os.environ.get("PWPH_EVAL") or next(
    (p for p in (_HERE.parents[1] / "evaluation",) if p.is_dir()),
    "/home/alex_1709/projects/iCC-team/experiments/pwph-e10-even-ladder-confirm-20261105"))
ROWS = {"H2": [("NoQ", "NOQ"), ("authors", "AUTHOR-F"), ("All", "CTRL"), ("Ripser", "RIPSER")],
        "H1": [("NoQ", "NOQ"), ("authors", "AUTHOR"), ("All", "CTRL"), ("Extend", "PAPER"), ("Ripser", "RIPSER")]}


def gm_ratio(clouds, key, arm):
    logs = [math.log(c[key]["UNI"] / c[key][arm]) for c in clouds if arm in c[key]]
    return math.exp(sum(logs) / len(logs)), len(logs)


def main():
    s = json.loads((E10 / "SUMMARY-CORRECTED.json").read_text())
    print("dim  vs.      clouds  mem   cpu   (wall, ratio of medians; Table II uses the regression)")
    for e, rows in ROWS.items():
        pc = s["per_cloud"][e]
        for label, arm in rows:
            mem, n = gm_ratio(pc, "median_rss_mib", arm)
            cpu, _ = gm_ratio(pc, "median_cpu_seconds", arm)
            wall, _ = gm_ratio(pc, "median_seconds", arm)
            print(f"{e}   {label:8s} {n:6d}  {mem:.3f} {cpu:.3f} ({wall:.3f})")
        for n in sorted({c["n"] for c in pc}):
            same = [c for c in pc if c["n"] == n]
            cpu, k = gm_ratio(same, "median_cpu_seconds", "RIPSER")
            wall, _ = gm_ratio(same, "median_seconds", "RIPSER")
            print(f"{e}   Ripser at n={n}: {k} clouds, cpu {cpu:.3f}, wall (medians) {wall:.3f}")
    print("capability runs (n=400, one run each):")
    for r in s["capability"]:
        print(f"  {r['case']:36s} {r['arm']:9s} {r['seconds']:8.1f} s  cpu {r['cpu_seconds']:8.1f} s  "
              f"{r['maxrss_mib'] / 1024:5.1f} GiB  {r['check']}")


if __name__ == "__main__":
    main()
