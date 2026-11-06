"""Generators of the E8 families that were never used before (clifford-torus-r4, ellipsoid-r3, three-spheres-r3,
genus2-r3, thick-shell-r3). Each takes a numpy Generator and n and returns an (n, ambient) array."""
import numpy as np


def unit(rng, n, k):
    v = rng.normal(size=(n, k))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def clifford_torus(rng, n):
    """Flat torus in R^4: H1 = 2, H2 = 1."""
    u = rng.uniform(0, 2 * np.pi, n)
    v = rng.uniform(0, 2 * np.pi, n)
    return np.column_stack((np.cos(u), np.sin(u), np.cos(v), np.sin(v))) / np.sqrt(2) + rng.normal(0, .01, (n, 4))


def ellipsoid(rng, n):
    """Ellipsoid with axes 1.5, 1.0, 0.6 (directions uniform on the sphere): H2 = 1."""
    return unit(rng, n, 3) * np.array([1.5, 1.0, 0.6]) + rng.normal(0, .01, (n, 3))


def three_spheres(rng, n):
    """Three disjoint spheres of radii 1.0, 0.75, 0.5, points in proportion to area: H2 = 3."""
    radius = np.array([1.0, 0.75, 0.5])
    centre = np.array([[0, 0, 0], [2.4, 0, 0], [0.6, 2.3, 0]])
    k = rng.choice(3, n, p=radius ** 2 / (radius ** 2).sum())
    return centre[k] + radius[k, None] * unit(rng, n, 3) + rng.normal(0, .01, (n, 3))


def genus2(rng, n):
    """Boundary of a thickened figure eight, (x^4 - x^2 + y^2)^2 + z^2 = 0.1, by Newton projection of uniform
    points in a box: H1 = 4, H2 = 1."""
    out = []
    while len(out) < n:
        p = rng.uniform([-1.3, -0.8, -0.4], [1.3, 0.8, 0.4], (4 * n, 3))
        for _ in range(30):
            x, y, z = p.T
            q = x ** 4 - x ** 2 + y ** 2
            f = q ** 2 + z ** 2 - 0.1
            g = np.column_stack((2 * q * (4 * x ** 3 - 2 * x), 4 * q * y, 2 * z))
            p = p - (f / np.maximum((g ** 2).sum(1), 1e-12))[:, None] * g
        x, y, z = p.T
        ok = np.abs((x ** 4 - x ** 2 + y ** 2) ** 2 + z ** 2 - 0.1) < 1e-9
        out.extend(p[ok][: n - len(out)])
    return np.array(out) + rng.normal(0, .005, (n, 3))


def thick_shell(rng, n):
    """Uniform in the shell 0.8 <= r <= 1: H2 = 1 at larger scales, much small-scale noise."""
    r = (0.8 ** 3 + rng.uniform(0, 1, n) * (1 - 0.8 ** 3)) ** (1 / 3)
    return unit(rng, n, 3) * r[:, None]


DRAW = {"clifford-torus-r4": clifford_torus, "ellipsoid-r3": ellipsoid, "three-spheres-r3": three_spheres,
        "genus2-r3": genus2, "thick-shell-r3": thick_shell}
