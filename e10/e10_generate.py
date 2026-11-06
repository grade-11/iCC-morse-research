"""E10 cohort generation: CSV and SHA only. No metric, q, E, PH or timing here.
H2 main: FAMILIES x H2_SIZES; H1 main: FAMILIES x H1_SIZES; capability: CAPABILITY (H2). One draw each."""
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

from e10_common import (AMBIENT, CAPABILITY, FAMILIES, H1_SIZES, H2_SIZES, R, SPHERE_GEN, load, require, sha,
                       utc, verify_lock, verify_pins, write_new)
from new_shapes import DRAW

E9 = Path("/home/alex_1709/projects/iCC-team/experiments/pwph-e9-unified-confirm-20261105")
SEED_BASE = 20261107300000


def authors_module():
    spec = importlib.util.spec_from_file_location("dSphere_gen", SPHERE_GEN)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def draw(family, n, seed):
    if family.startswith("author-s"):
        k = int(family[len("author-s"):])
        np.random.seed(seed % 2 ** 32)  # the authors' generator draws from numpy's global state
        a = np.array(authors_module().dSphere_uniformRandom(dimension=k + 1, numPoints=n), dtype=float)
    else:
        rng = np.random.Generator(np.random.PCG64(seed))
        if family in DRAW:
            a = DRAW[family](rng, n)
        elif family == "corrugated-torus":  # as in E3-E7
            t = np.where(rng.random(n) < .55, rng.normal(2, .25, n), rng.uniform(0, 2 * np.pi, n))
            p = rng.uniform(0, 2 * np.pi, n)
            tube = .18 + .09 * (1 + np.cos(3 * t))
            center = 1.3 + .15 * np.sin(2 * t)
            a = np.column_stack(((center + tube * np.cos(p)) * np.cos(t),
                                 .8 * (center + tube * np.cos(p)) * np.sin(t), tube * np.sin(p)))
            a += rng.normal(size=(n, 3)) * .006
        elif family == "concentric-spheres":  # as in the H2 development set
            label = rng.integers(0, 2, n)
            u = rng.uniform(-1, 1, n)
            t = rng.uniform(0, 2 * np.pi, n)
            radius = np.array([1., 1.8])[label]
            s = np.sqrt(1 - u * u)
            a = radius[:, None] * np.column_stack((s * np.cos(t), s * np.sin(t), u)) + rng.normal(0, .012, (n, 3))
        elif family == "linked-tori":
            side = rng.integers(0, 2, n)
            u = rng.uniform(0, 2 * np.pi, n)
            v = rng.uniform(0, 2 * np.pi, n)
            ring = 1 + 0.25 * np.cos(v)
            first = np.column_stack((ring * np.cos(u), ring * np.sin(u), 0.25 * np.sin(v)))
            second = np.column_stack((1 + ring * np.cos(u), 0.25 * np.sin(v), ring * np.sin(u)))
            a = np.where(side[:, None] == 0, first, second) + rng.normal(0, 0.005, (n, 3))
        elif family == "klein-bottle-r4":
            u = rng.uniform(0, 2 * np.pi, n)
            v = rng.uniform(0, 2 * np.pi, n)
            a = np.column_stack(((2 + np.cos(v)) * np.cos(u), (2 + np.cos(v)) * np.sin(u),
                                 np.sin(v) * np.cos(u / 2), np.sin(v) * np.sin(u / 2))) + rng.normal(0, .006, (n, 4))
        elif family == "cube-surface":
            face = rng.integers(0, 6, n)
            p = rng.uniform(-1, 1, (n, 3))
            axis = face // 2
            p[np.arange(n), axis] = np.where(face % 2 == 0, -1.0, 1.0)
            a = p + rng.normal(0, 0.01, (n, 3))
        else:
            raise ValueError(family)
    require(a.shape == (n, AMBIENT.get(family, 3)), "shape differs for " + family)
    # Same constructive distinctness transform as E3-E7 (first coordinate only).
    for i in range(n):
        a[i, 0] = math.ldexp(math.floor(float(np.clip(a[i, 0], -16., 16.)) * 2 ** 24) * 512 + i + 1, -33)
    return a


def history():
    """Seeds and CSV hashes of E9 and of everything E9 excluded (E3-E8, P1), which E10 must not reuse."""
    seeds, hashes = set(), set()
    for c in load(E9 / "COHORT-LOCK.json")["rows"]:
        seeds.add(c["seed"])
        hashes.add(c["sha256"])
    sys.path.insert(0, str(E9))
    import e9_generate  # noqa: E402
    old_seeds, old_hashes = e9_generate.history()
    return seeds | old_seeds, hashes | old_hashes


def main():
    verify_pins()
    verify_lock()
    require(not (R / "COHORT-LOCK.json").exists(), "E10 cohort already generated")
    old_seeds, old_hashes = history()
    plan = [("H2", "MAIN", f, n) for f in FAMILIES for n in H2_SIZES]
    plan += [("H1", "MAIN", f, n) for f in FAMILIES for n in H1_SIZES]
    plan += [("H2", "CAPABILITY", f, n) for f, n in CAPABILITY]
    (R / "cohort").mkdir()
    rows, new_hashes = [], set()
    for experiment, role, family, n in plan:
        group = {"H2": 0, "H1": 1}[experiment] if role == "MAIN" else 2
        seed0 = SEED_BASE + 100000 * group + 1000 * FAMILIES.index(family) + 10 * (n // 100)
        case = f"E10-{experiment}-{family}-n{n}" + ("-cap" if role == "CAPABILITY" else "")
        rejected = []
        for seed in range(seed0, seed0 + 10):
            require(seed not in old_seeds, "Seed reused: " + str(seed))
            path = R / "cohort" / (case + ("" if not rejected else f"-seed{seed}") + ".csv")
            with path.open("x") as f:
                np.savetxt(f, draw(family, n, seed), delimiter=",", fmt="%.17g")
            digest = sha(path)
            if digest not in old_hashes | new_hashes:
                break
            rejected.append(dict(seed=seed, sha256=digest, reason="HASH_COLLISION"))
        else:
            raise RuntimeError("No predeclared seed available for " + case)
        new_hashes.add(digest)
        rows.append(dict(case=case, experiment=experiment, role=role, family=family, n=n,
                         ambient_dim=AMBIENT.get(family, 3), seed=seed, path=str(path.relative_to(R)),
                         sha256=digest, rejected=rejected))
        print("GENERATED", case, flush=True)
    write_new(R / "COHORT-LOCK.json", dict(
        status="SEALED_CSV_AND_SHA_ONLY", utc=utc(), generator_sha256=sha(__file__),
        families=FAMILIES, clouds=len(rows), rows=rows, metric_q_E_PH="NOT_COMPUTED_AT_GENERATION"))
    (R / "cohort/SHA256SUMS").write_text(
        "".join(r["sha256"] + "  " + Path(r["path"]).name + "\n" for r in rows))
    print("E10 COHORT LOCKED", len(rows), flush=True)


if __name__ == "__main__":
    sys.exit(main())
