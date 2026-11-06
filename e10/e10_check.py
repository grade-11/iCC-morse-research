"""E10 correctness, after timing closed: every complete PwPH run (measured and capability) against Ripser 0.6.15
in the same cloud, H1 and H2 for H2 runs and H1 for H1 runs, after projecting PwPH endpoints to float32. H0 is
not checked."""
import collections
import json
import sys

import numpy as np

from e10_common import MAXDIM, R, diff, load, parse, project, require, utc, verify_lock, verify_pins, write_new


def main():
    import ripser
    verify_pins()
    verify_lock()
    require(ripser.__version__ == "0.6.15", "Ripser version differs")
    require((R / "TIMING-END.json").exists(), "timing has not closed")
    require(not (R / "CHECKS.json").exists(), "E10 check already ran")
    runs = [json.loads(l) for l in (R / "RUNS.jsonl").read_text().splitlines()]
    out = []
    for cloud in load(R / "COHORT.json")["rows"]:
        top = MAXDIM[cloud["experiment"]] - 1
        points = np.loadtxt(R / cloud["path"], delimiter=",")
        dgms = ripser.ripser(points, maxdim=top, thresh=float(cloud["E_decimal17"]), coeff=2,
                             metric="euclidean", do_cocycles=False, n_perm=None)["dgms"]
        ref = {h: collections.Counter((float(b), float(d)) for b, d in dgms[h].tolist()) for h in range(1, top + 1)}
        for r in runs:
            if r["case"] != cloud["case"] or r["arm"] == "RIPSER":
                continue
            row = {k: r[k] for k in ("run_id", "case", "experiment", "family", "n", "arm", "stage", "block", "status")}
            if r["status"] != "COMPLETE":
                row["check"] = "NOT_CHECKED"
            else:
                try:
                    p = parse(R / r["capture_directory"] / "stdout.txt")
                    res = {h: diff(ref[h], project(p["bars"][h])) for h in ref}
                    row.update({f"H{h}": res[h] for h in res})
                    row["check"] = "PASS" if all(v["status"] == "PASS" for v in res.values()) else "PH_MISMATCH"
                    row["commits"] = sum(x["committed"] > 0 for x in p["posts"])
                    row["windows"] = len(p["windows"])
                except Exception as exc:  # parse failures are reported, never dropped
                    row.update(check="CHECK_ERROR", reason=str(exc))
            out.append(row)
        print("CHECKED", cloud["case"], flush=True)
    counts = collections.Counter((r["experiment"], r["arm"], r["check"]) for r in out)
    write_new(R / "CHECKS.json", dict(utc=utc(), H0="NOT_CHECKED",
                                      counts={f"{e}:{a}:{c}": v for (e, a, c), v in sorted(counts.items())},
                                      runs=out))
    for k, v in sorted(counts.items()):
        print("E10 CHECK", *k, v, flush=True)


if __name__ == "__main__":
    sys.exit(main())
