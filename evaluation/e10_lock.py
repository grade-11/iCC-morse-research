"""Writes LOCK.json before the cohort is generated: hashes of every script and of PROTOCOL.txt, the pins, and the
timing and analysis parameters."""
import sys

from e10_common import (ARMS, BINARIES, BLOCKS, CAP, FAMILIES, H1_CONTROLLER_LADDER, H1_SIZES, H2_SIZES, METRIC_SHA,
                       ONCE_PER_CLOUD, R, SPHERE_GEN_SHA, UNIFIED_LADDER, SPHERE_UTIL_SHA, TIMEOUT, sha, utc,
                       verify_pins, write_new)

FILES = ["e10_common.py", "new_shapes.py", "e10_generate.py", "e10_prepare.py", "e10_time.py", "e10_check.py",
         "e10_analyze.py", "ripser_run.py", "e10_lock.py", "run_e10.sh", "PROTOCOL.txt"]


def main():
    verify_pins()
    write_new(R / "LOCK.json", dict(
        utc=utc(), status="LOCKED_BEFORE_COHORT_GENERATION", binaries={k: dict(sha256=d, commit=c) for k, (_, d, c) in BINARIES.items()},
        unified_ladder=UNIFIED_LADDER, metric_helper_sha256=METRIC_SHA,
        authors_generator_sha256={"dSphere_gen.py": SPHERE_GEN_SHA, "utilities.py": SPHERE_UTIL_SHA},
        script_sha256={name: sha(R / name) for name in FILES},
        families=FAMILIES, sizes=dict(H2=H2_SIZES, H1=H1_SIZES), h1_controller_ladder=H1_CONTROLLER_LADDER,
        arms={e: {a: dict(ladder=l, flags=f, binary=b) for a, (l, f, b) in v.items()} for e, v in ARMS.items()},
        timing=dict(blocks=BLOCKS, once_per_cloud={k: sorted(v) for k, v in ONCE_PER_CLOUD.items()},
                    order_seed=2026110691, block_seed=2026110692, capability_seed=2026110693,
                    cap_bytes=CAP, timeout_s=TIMEOUT),
        analysis=dict(bootstrap_seed=2026110694, bootstrap_reps=5000, slow_threshold=1.05, max_slow_families=1,
                      noninferiority_bound=1.05)))
    print("E10 LOCKED", flush=True)


if __name__ == "__main__":
    sys.exit(main())
