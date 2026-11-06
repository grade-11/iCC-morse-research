"""Print the rows of the paper's Table II (results) and Table III (families) from an E10 summary.
Usage: e10_tables.py [SUMMARY-CORRECTED.json | SUMMARY.json]   (default: the corrected analysis).
Reads only; writes nothing. Memory ratios are GMs over clouds of per-cloud median peak RSS."""
import json
import math
import sys
from pathlib import Path

R = Path(__file__).resolve().parent
LETTER = {"CONFIRMED": "C", "FASTER": "F", "NONINFERIOR": "NI", "NOT_CONFIRMED": "NC", "PARTIAL": "P"}
FAMILIES = ["author-s2", "author-s3", "corrugated-torus", "concentric-spheres", "linked-tori", "klein-bottle-r4",
            "cube-surface", "clifford-torus-r4", "ellipsoid-r3", "three-spheres-r3", "genus2-r3", "thick-shell-r3"]


def memory(s, e, arm):
    logs = [math.log(c["median_rss_mib"]["UNI"] / c["median_rss_mib"][arm]) for c in s["per_cloud"][e]
            if arm in c["median_rss_mib"]]
    return math.exp(sum(logs) / len(logs))


def main():
    s = json.loads((R / (sys.argv[1] if len(sys.argv) > 1 else "SUMMARY-CORRECTED.json")).read_text())
    print("Table II:", s.get("specification", "locked: log t ~ arm + position + block"))
    for e, rows in (("H2", [("NoQ", "H2_N", "NOQ"), ("authors", "H2_S", "AUTHOR-F"), ("All", "H2_C", "CTRL"),
                            ("Ripser", None, "RIPSER")]),
                    ("H1", [("NoQ", "H1_N", "NOQ"), ("authors", "H1_S", "AUTHOR"), ("All", "H1_C", "CTRL"),
                            ("Ext", "H1_E", "PAPER"), ("Ripser", None, "RIPSER")])):
        for label, h, arm in rows:
            if h:
                r, v = s["hypotheses"][h], LETTER[s["verdicts"][h]]
            else:
                r, v = next(d for d in s["descriptive"][e] if d["name"] == "UNI/RIPSER"), "--"
            ci = lambda x: "[%.3f, %.3f]" % tuple(x)
            print(f"{e} {label:8s} {r['gm']:.3f} {ci(r['cloud_ci'])} {ci(r['family_ci'])} {v:3s} "
                  f"mem {memory(s, e, arm):.2f}")
    print("Table III: family, H2 vs NoQ, H2 vs All, H1 vs NoQ, H1 vs All")
    for f in FAMILIES:
        g = [s["hypotheses"][h]["families"][f]["gm"] for h in ("H2_N", "H2_C", "H1_N", "H1_C")]
        print(f"{f:20s} " + " ".join(f"{x:.3f}" for x in g))


if __name__ == "__main__":
    sys.exit(main())
