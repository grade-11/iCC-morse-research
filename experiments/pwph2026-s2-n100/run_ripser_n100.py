#!/usr/bin/env python3
"""Create a new Ripser reference for the fixed PwPH S2/n100 CSV.

Run from the iCC repository root with ripser==0.6.15 installed.
This script does not run iCC or compare diagrams. Exit 0 means a reference
was saved successfully, not that iCC/Ripser agreement has passed.
Outputs go to a new experiments/pwph2026-s2-n100/runs/ripser-reference-*.
"""

import argparse
import contextlib
import csv
import datetime
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import tarfile
import tempfile
import traceback

DATA = Path("experiments/pwph2026-s2-n100/data/sphere_S2_R3_n100_seed20261007.csv")
DATA_SHA = "e39f71c26b3f71c96860c04316b9eae7ec192891ebdd39c80bbebb2029f76b48"
RUNS = Path("experiments/pwph2026-s2-n100/runs")
PARAMETERS = dict(maxdim=2, thresh=2.0, coeff=2, metric="euclidean",
                  distance_matrix=False, do_cocycles=False, n_perm=None)


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8")


def timestamp():
    return datetime.datetime.now().astimezone().isoformat()


def command_output(args):
    env = dict(os.environ, GIT_OPTIONAL_LOCKS="0")
    return subprocess.check_output(args, text=True, stderr=subprocess.STDOUT,
                                   env=env).strip()


def validate_diagram(np, diagram, dimension):
    if diagram.ndim != 2 or diagram.shape[1] != 2:
        raise ValueError(f"H{dimension}: expected shape (n, 2)")
    births, deaths = diagram[:, 0], diagram[:, 1]
    if not np.all(np.isfinite(births) & (births >= 0) & (births <= 2.0)):
        raise ValueError(f"H{dimension}: invalid birth endpoint")
    finite = np.isfinite(deaths)
    if not np.all(finite | np.isposinf(deaths)):
        raise ValueError(f"H{dimension}: NaN or negative infinity death")
    if not np.all((deaths[finite] >= births[finite]) & (deaths[finite] <= 2.0)):
        raise ValueError(f"H{dimension}: invalid finite interval")
    return dict(total=int(len(diagram)), finite=int(finite.sum()),
                surviving_at_cutoff=int((~finite).sum()),
                zero_length=int(np.sum(finite & (births == deaths))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    run_dir = None
    try:
        root = Path.cwd().resolve()
        if Path(command_output(["git", "rev-parse", "--show-toplevel"])).resolve() != root:
            raise ValueError("Run this script from the iCC repository root")
        source_bytes = DATA.read_bytes()
        if sha_bytes(source_bytes) != DATA_SHA:
            raise ValueError("CSV checksum differs from the fixed S2/n100 dataset")
        if importlib.metadata.version("ripser") != "0.6.15":
            raise ValueError("This reference requires ripser==0.6.15")
        import numpy as np
        from ripser import ripser

        points = np.loadtxt(io.BytesIO(source_bytes), delimiter=",", dtype=np.float64,
                            ndmin=2)
        if points.shape != (100, 3) or not np.all(np.isfinite(points)):
            raise ValueError("CSV must contain exactly 100 finite rows of 3 coordinates")

        packages = {name: importlib.metadata.version(name) for name in
                    ("ripser", "numpy", "scipy", "scikit-learn", "persim")}
        freeze = command_output([sys.executable, "-m", "pip", "freeze", "--all"])
        script_bytes = Path(__file__).read_bytes()
        metadata = dict(
            schema_version=1, status="RUNNING", started_at=timestamp(),
            purpose="New PH reference; no iCC run, diagram comparison, or performance benchmark",
            command=shlex.join([sys.executable, *sys.argv]), cwd=str(root),
            head=command_output(["git", "rev-parse", "HEAD"]),
            branch=command_output(["git", "branch", "--show-current"]),
            git_status_before=command_output(["git", "status", "--short"]),
            python_executable=sys.executable, python_version=sys.version,
            python_prefix=sys.prefix, python_base_prefix=sys.base_prefix,
            platform=platform.platform(), packages=packages,
            dataset_path=str(DATA.resolve()), dataset_sha256=DATA_SHA,
            input_shape=list(points.shape), input_dtype=str(points.dtype),
            script_sha256=sha_bytes(script_bytes), parameters=PARAMETERS,
            notes=["All intervals and multiplicities are retained",
                   "inf is recorded as surviving at cutoff; no claim about a longer filtration",
                   "This is a new WSL reference, not recovery of the historical Codespaces run",
                   "No iCC agreement has been checked, including H0"],
        )
        RUNS.mkdir(parents=True, exist_ok=True)
        run_dir = Path(tempfile.mkdtemp(prefix="ripser-reference-", dir=RUNS))
        write_json(run_dir / "context.json", metadata)
        (run_dir / "input.csv").write_bytes(source_bytes)
        (run_dir / "runner.py").write_bytes(script_bytes)
        (run_dir / "requirements-freeze.txt").write_text(freeze + "\n", encoding="utf-8")
        summary = None
        error = None
        with (run_dir / "ripser.log").open("w", encoding="utf-8") as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                try:
                    print("Python-level run log. Purpose: PH reference only.", flush=True)
                    print(json.dumps(PARAMETERS, sort_keys=True), flush=True)
                    result = ripser(points, **PARAMETERS)
                    diagrams = [np.asarray(d) for d in result["dgms"]]
                    # Save raw diagrams before validation, including repeated intervals.
                    np.savez_compressed(run_dir / "diagrams.npz",
                                        **{f"H{d}": x for d, x in enumerate(diagrams)})
                    if len(diagrams) != 3:
                        raise ValueError("Expected H0, H1 and H2 from maxdim=2")
                    counts = {f"H{d}": validate_diagram(np, x, d)
                              for d, x in enumerate(diagrams)}
                    distances = np.asarray(result["dperm2all"])
                    if distances.shape != (100, 100) or not np.all(np.isfinite(distances)):
                        raise ValueError("Unexpected returned distance matrix")
                    # Keep the API-returned matrix without changing its precision.
                    np.save(run_dir / "dperm2all.npy", distances, allow_pickle=False)
                    with (run_dir / "intervals.csv").open("w", newline="", encoding="utf-8") as f:
                        writer = csv.writer(f)
                        writer.writerow(["homology_dimension", "birth", "death", "death_status"])
                        for d, diagram in enumerate(diagrams):
                            for birth, death in diagram:
                                alive = bool(np.isposinf(death))
                                writer.writerow([d, format(float(birth), ".17g"),
                                                 "inf" if alive else format(float(death), ".17g"),
                                                 "surviving_at_cutoff" if alive else "finite"])
                    after_hash = sha_bytes(DATA.read_bytes())
                    metadata["dataset_sha256_after"] = after_hash
                    if after_hash != DATA_SHA:
                        raise ValueError("Source CSV changed while the reference was running")
                    summary = dict(status="REFERENCE_CREATED", comparison_status="NOT_RUN",
                                   counts=counts, num_edges=int(result["num_edges"]),
                                   cutoff=2.0, coefficient_field=2,
                                   H2=[[format(float(v), ".17g") for v in row]
                                       for row in diagrams[2]])
                    metadata.update(status="REFERENCE_CREATED", exit_status=0,
                                    diagram_dtypes=[str(x.dtype) for x in diagrams],
                                    returned_distance_dtype=str(distances.dtype))
                    write_json(run_dir / "summary.json", summary)
                    print(json.dumps(summary, indent=2), flush=True)
                except Exception as exc:
                    error = f"{type(exc).__name__}: {exc}"
                    metadata.update(status="ERROR", exit_status=2, error=error)
                    traceback.print_exc()
        metadata["finished_at"] = timestamp()
        write_json(run_dir / "context.json", metadata)
        # The manifest covers evidence files; it deliberately does not hash itself.
        manifest = {p.name: sha_bytes(p.read_bytes()) for p in sorted(run_dir.iterdir())
                    if p.is_file()}
        write_json(run_dir / "checksums.json", manifest)
        archive = run_dir.with_suffix(".tar.gz")
        with tarfile.open(archive, "x:gz") as tar:
            tar.add(run_dir, arcname=run_dir.name)
        print(json.dumps(dict(status=metadata["status"], output_directory=str(run_dir.resolve()),
                              archive=str(archive.resolve()), summary=summary, error=error), indent=2))
        return 2 if error else 0
    except Exception as exc:
        print(json.dumps(dict(status="ERROR", message=f"{type(exc).__name__}: {exc}",
                              output_directory=str(run_dir) if run_dir else None), indent=2),
              file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
