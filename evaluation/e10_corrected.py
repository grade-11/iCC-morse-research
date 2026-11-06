"""E10 corrected analysis (written 06/11/2026, after the locked analysis; DEVIATIONS.txt item 3). The locked model
log t ~ arm + position + block is rank-deficient on two H2 clouds, where no arm contrast is identified
(IDENTIFIABILITY.json). The corrected specification drops the position term, log t ~ arm + block, on all 48
clouds of both dimensions. Everything else is the locked analysis unchanged: e10_analyze.main() with its
contrasts, decision rules, bootstrap seed and reps, descriptive statistics and memory/CPU ratios. The locked
output SUMMARY.json is kept; this writes SUMMARY-CORRECTED.json once."""
import math
import sys

import numpy as np

import e10_analyze as A
from e10_common import R, write_new


def fit_without_position(rows):
    """{arm: log time relative to UNI} by OLS of log t on arm and block."""
    present = sorted({r["arm"] for r in rows} - {A.REFERENCE})
    blocks = sorted({r["block"] for r in rows})
    X, y = [], []
    for r in rows:
        X.append([1.0] + [1.0 if r["arm"] == a else 0.0 for a in present]
                 + [1.0 if r["block"] == b else 0.0 for b in blocks[1:]])
        y.append(math.log(r["seconds"]))
    X = np.array(X)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise RuntimeError("corrected design is rank-deficient")
    beta, *_ = np.linalg.lstsq(X, np.array(y), rcond=None)
    return {A.REFERENCE: 0.0, **{a: float(beta[1 + k]) for k, a in enumerate(present)}}


def main():
    A.fit = fit_without_position
    A.write_new = lambda path, value: write_new(R / "SUMMARY-CORRECTED.json",
                                                dict(value, specification="log t ~ arm + block (no position)"))
    A.main()


if __name__ == "__main__":
    sys.exit(main())
