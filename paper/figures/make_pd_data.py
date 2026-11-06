"""Data of Fig. 3 (persistence diagrams of the controller and Ripser), written 2026-11-06.
The controller's bars come from the block-1 UNI runs of E10 (run 114: E10-H2-author-s2-n200; run 333:
E10-H1-clifford-torus-r4-n512), read with e10_common.parse and projected to float32 with e10_common.project, exactly
as e10_check.py compares them. Ripser 0.6.15 is run on the same cloud with the arguments of e10_check.py. Both sides
are printed with %g, a class alive at E on the line death = E. Reads E10 only; writes into the directory given
(default: this directory). Usage: <venv python> make_pd_data.py [out_dir]"""
import collections
import json
import os
import sys
from pathlib import Path

import numpy as np
import ripser

# E10 directory: $PWPH_E10, else ../../e10 (layout of the artifact branch), else the original location.
_HERE = Path(__file__).resolve().parent
E10 = Path(os.environ.get("PWPH_E10") or next(
    (p for p in (_HERE.parents[1] / "e10",) if p.is_dir()),
    "/home/alex_1709/projects/iCC-team/experiments/pwph-e10-even-ladder-confirm-20261105"))
sys.path.insert(0, str(E10))
from e10_common import parse, project  # noqa: E402

FIGS = [("s2", "E10-H2-author-s2-n200", 2, (1, 2)), ("cliff", "E10-H1-clifford-torus-r4-n512", 1, (1,))]


def write(path, counter, E):
    rows = sorted(counter.elements())
    path.write_text("b d\n" + "".join("%g %g\n" % (b, E if d == float("inf") else d) for b, d in rows))


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    runs = [json.loads(l) for l in (E10 / "RUNS.jsonl").read_text().splitlines()]
    cohort = {c["case"]: c for c in json.loads((E10 / "COHORT.json").read_text())["rows"]}
    for tag, case, top, dims in FIGS:
        c = cohort[case]
        E = float(c["E_decimal17"])
        run = next(r for r in runs if r["case"] == case and r["arm"] == "UNI" and r["block"] == 1)
        bars = parse(E10 / run["capture_directory"] / "stdout.txt")["bars"]
        points = np.loadtxt(E10 / c["path"], delimiter=",")
        dgms = ripser.ripser(points, maxdim=top, thresh=E, coeff=2, metric="euclidean", do_cocycles=False,
                             n_perm=None)["dgms"]
        for h in dims:
            ctrl = project(bars[h])
            rips = collections.Counter((float(b), float(d)) for b, d in dgms[h].tolist())
            assert ctrl == rips, f"{case} H{h}: controller and Ripser differ"
            write(out / f"{tag}-h{h}-ctrl.dat", ctrl, float(np.float32(E)))
            write(out / f"{tag}-h{h}-ripser.dat", rips, float(np.float32(E)))
            print(case, f"run {run['run_id']}", f"H{h}", sum(ctrl.values()), "bars, identical at float32")


if __name__ == "__main__":
    main()
