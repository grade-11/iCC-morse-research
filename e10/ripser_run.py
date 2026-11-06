"""Whole-process Ripser 0.6.15 run: H0..H<maxdim> up to the threshold E, coefficients F2."""
import sys

import numpy as np
import ripser

points = np.loadtxt(sys.argv[1], delimiter=",")
dgms = ripser.ripser(points, maxdim=int(sys.argv[3]), thresh=float(sys.argv[2]), coeff=2, metric="euclidean",
                     do_cocycles=False, n_perm=None)["dgms"]
print("ripser", ripser.__version__, " ".join(f"H{k} {len(d)}" for k, d in enumerate(dgms)))
