"""E10 shared helpers (every PwPH arm on the binary of the unified controller): paths, pins, families, arms, hashing, guarded writes, native capture, output parsing."""
import collections
import datetime
import hashlib
import json
import math
import os
import re
import resource
import struct
import subprocess
import time
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parent
PYTHON = "/home/alex_1709/.venvs/icc-ripser-0.6.15-VnIso6/bin/python"
# Two kernels of one branch (work/pwph-h2-birth-repair-20261102, local commits only): the draft's controller
# (E8) and the unified controller (one ladder plus the stop rule with the growing-rate estimate).
BINARIES = {"B_X22": (R / "bin/B_X22", "2a94c159cd5f487828a87edf3c0650e0400adcf75116007ee6428b4ba54e1123", "0746d0815"),
            "B_X31": (R / "bin/B_X31", "9a0f76aa6a40ac74f519e135d6b067f082dae059b085df02d43fb1fa7e9b311d", "f4733ef5e")}
METRIC = R / "bin/native-metric-export"
METRIC_SHA = "4624800161db468ebd50e1061a0f8f5440be469f8b1652275b50e43f818672f0"
# The authors' d-sphere generator, unmodified (as in E7).
SPHERE_GEN = R / "authors-generator/dSphere_gen.py"
SPHERE_GEN_SHA = "c3d95a94fc9482e2e4dadb000af01b38fc9054a6ed9c76ce1c5fec4ecf5f01d0"
SPHERE_UTIL = R / "authors-generator/utilities.py"
SPHERE_UTIL_SHA = "37c4cfbc1d53335e1cfbb3ea2abf18db367fd39c02742311946da18c1d474c43"
RIPSER_RUN = R / "ripser_run.py"

# Twelve families: the authors' S^2 and S^3, five families of the H2 development set (fresh draws), and five
# shapes never used before (new_shapes.py).
FAMILIES = ["author-s2", "author-s3", "corrugated-torus", "concentric-spheres", "linked-tori", "klein-bottle-r4",
            "cube-surface", "clifford-torus-r4", "ellipsoid-r3", "three-spheres-r3", "genus2-r3", "thick-shell-r3"]
AMBIENT = {"author-s3": 4, "klein-bottle-r4": 4, "clifford-torus-r4": 4}
AUTHOR_WINDOWS = {"author-s2": ["1.2", "1.7", "2.0"], "author-s3": ["1.3", "1.6", "1.8", "2.0"]}
H2_SIZES = [200, 300]
H1_SIZES = [512, 768]
CAPABILITY = [("author-s2", 400), ("corrugated-torus", 400)]
EVEN9 = [10, 20, 30, 40, 50, 60, 70, 80, 90]
PAPER = [20, 25, 30, 40, 50]
P9 = [20, 25, 30, 40, 50, 60, 70, 80, 90]
PCTS = sorted(set(EVEN9 + PAPER + P9))

BASE = ["--tool", "piecewise", "--threads", "4", "--pv-cap-scale", "0.85", "--pv-min-separation", "0"]
GUARD = ["--live-guard"]
FILTERS = ["--h1-risk", "selective-fast", "--fiber-full"]
# The controller: Extend before the first commit, then a commit at every admitted cut, absorbing old pseudo
# vertices; both filters; a flagged support is shrunk to its safe part instead of dropped.
CONTROLLER = GUARD + ["--extend", "--absorb-pv", "--shrink-risky"] + FILTERS
# The unified controller: the same flags, every window retained and extended (--extend-all), and the stop rule
# (--extend-rule) on one ladder for H1 and H2.
UNIFIED = CONTROLLER + ["--extend-all", "--extend-rule"]
UNIFIED_LADDER = "EVEN9"  # q10..q90: ten equal parts of the distinct distances (PROTOCOL.txt)
H1_CONTROLLER_LADDER = "PAPER"  # E8
# arm: (experiment, ladder, flags, binary). Ladders: EVEN9, P9, PAPER, AUTHOR, NOQ (just E).
ARMS = {
    "H2": {"UNI": (UNIFIED_LADDER, UNIFIED, "B_X31"),
           "CTRL": ("EVEN9", CONTROLLER, "B_X31"),
           "AUTHOR-F": ("AUTHOR", GUARD + FILTERS, "B_X31"),
           "AUTHOR": ("AUTHOR", GUARD, "B_X31"),
           "NOQ": ("NOQ", [], "B_X31"),
           "RIPSER": (None, [], None)},
    "H1": {"UNI": (UNIFIED_LADDER, UNIFIED, "B_X31"),
           "CTRL": (H1_CONTROLLER_LADDER, CONTROLLER, "B_X31"),
           # The paper's Extend: its ladder, commit once, pair filter only.
           "PAPER": ("PAPER", GUARD + ["--extend", "--extend-commit-once", "--h1-risk", "selective-fast"], "B_X31"),
           "AUTHOR": ("AUTHOR", GUARD + ["--h1-risk", "selective-fast"], "B_X31"),
           "NOQ": ("NOQ", [], "B_X31"),
           "RIPSER": (None, [], None)},
}
MAXDIM = {"H2": 3, "H1": 2}
# Arms timed once per cloud (first block) instead of in every block.
ONCE_PER_CLOUD = {"H2": {"NOQ"}, "H1": {"NOQ"}}
BLOCKS = 3
CAP = {"H2": 18 * 1024 ** 3, "H1": 9 * 1024 ** 3}
TIMEOUT = {"H2": 2400, "H1": 600, "CAPABILITY": 3600}
ENV = {"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C", "OMP_NUM_THREADS": "4",
       "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "4", "OPENBLAS_NUM_THREADS": "1",
       "MKL_NUM_THREADS": "1"}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def word(value):
    return struct.pack(">d", float(value)).hex()


def d17(value):
    return format(float(value), ".17g")


def load(path):
    return json.loads(Path(path).read_text())


def write_new(path, value):
    """Sealed artifacts are written once; reruns must not overwrite them."""
    with Path(path).open("x") as f:
        json.dump(value, f, indent=1, ensure_ascii=False, allow_nan=False)
        f.write("\n")


def verify_pins():
    require(UNIFIED_LADDER in ("EVEN9", "P9"), "unified ladder not set")
    for name, (path, digest, _) in BINARIES.items():
        require(sha(path) == digest, name + " pin differs")
    require(sha(METRIC) == METRIC_SHA, "metric helper pin differs")
    require(sha(SPHERE_GEN) == SPHERE_GEN_SHA, "authors' generator pin differs")
    require(sha(SPHERE_UTIL) == SPHERE_UTIL_SHA, "authors' utilities pin differs")


def verify_lock():
    lock = load(R / "LOCK.json")
    for name, digest in lock["script_sha256"].items():
        require(sha(R / name) == digest, "Locked file changed: " + name)
    return lock


def argv_for(cloud, arm):
    experiment = cloud["experiment"]
    if arm == "RIPSER":
        return [PYTHON, "-I", RIPSER_RUN, R / cloud["path"], cloud["E_decimal17"], str(MAXDIM[experiment] - 1)]
    ladder, flags, binary = ARMS[experiment][arm]
    return ([BINARIES[binary][0]] + BASE + ["--max-dim", str(MAXDIM[experiment]), "--file-name", R / cloud["path"]]
            + flags + ["--eps-breaks"] + cloud["cuts"][ladder])


def capture(argv, directory, timeout, cap):
    """Whole-process wall clock: perf_counter around Popen..wait4; RLIMIT_AS = cap bytes. Also user and
    system CPU seconds from wait4."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)

    def limit():
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    with (directory / "stdout.txt").open("w") as out, (directory / "stderr.txt").open("w") as err:
        start = time.perf_counter()
        proc = subprocess.Popen([str(a) for a in argv], stdout=out, stderr=err, env=ENV, cwd=directory,
                                preexec_fn=limit)
        status = "COMPLETE"
        deadline = start + timeout
        while True:
            pid, code, usage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if time.perf_counter() > deadline:
                proc.kill()
                pid, code, usage = os.wait4(proc.pid, 0)
                status = "TIMEOUT"
                break
            time.sleep(0.005)
        seconds = time.perf_counter() - start
    rc = os.waitstatus_to_exitcode(code)
    if status == "COMPLETE" and rc != 0:
        status = "FAILED"
    record = dict(argv=[str(a) for a in argv], status=status, returncode=rc, seconds=seconds,
                  cpu_seconds=usage.ru_utime + usage.ru_stime, maxrss_kib=usage.ru_maxrss, cap_bytes=cap, utc=utc())
    (directory / "status.json").write_text(json.dumps(record, indent=1) + "\n")
    return record


def parse(path):
    """Emitted H1/H2 (as in the H2 survey): a finite death is owned by its window [lo, hi); -1 in a non-final
    window is a temporary survivor and is dropped; -1 in the final window is a class alive at the endpoint.
    Also the [postcommit] lines."""
    lines = Path(path).read_text().splitlines()
    require(lines and lines[-1].startswith("run time = "), "incomplete output: " + str(path))
    wins, rows, posts = [], [], []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "window lower bound, window upper bound":
            lo, hi = map(float, lines[i + 1].split(","))
            wins.append((lo, hi))
            i += 3
            continue
        m = re.fullmatch(r"([23]),\s*([^,]+),\s*([^,]+)", line)
        if m:
            rows.append((len(wins) - 1, int(m[1]) - 1, float(m[2]), float(m[3])))
        elif line.startswith("[postcommit] "):
            fields = dict(t.split("=", 1) for t in line.split()[1:])
            posts.append(dict(epsilon=float(fields["epsilon"]), committed=int(fields["committed"]),
                              eligible=int(fields["eligible"]), active=int(fields["active"])))
        i += 1
    bars = {1: collections.Counter(), 2: collections.Counter()}
    for w, h, b, d in rows:
        lo, hi = wins[w]
        if d == -1 and w != len(wins) - 1:
            continue
        require(d == -1 or lo <= d < hi, "death outside its window: " + str(path))
        bars[h][(b, d)] += 1
    return dict(bars=bars, windows=wins, posts=posts)


def project(counter):
    """PwPH bars at float32 resolution, as Ripser reports them; -1 (alive at the endpoint) becomes inf."""
    out = collections.Counter()
    for (b, d), c in counter.items():
        out[(float(np.float32(b)), math.inf if d == -1 else float(np.float32(d)))] += c
    return out


def diff(reference, candidate):
    missing, extra = reference - candidate, candidate - reference
    return dict(status="PASS" if not missing and not extra else "MISMATCH",
                reference_count=sum(reference.values()),
                missing=[[b, d, c] for (b, d), c in sorted(missing.items())][:20],
                extra=[[b, d, c] for (b, d), c in sorted(extra.items())][:20])
