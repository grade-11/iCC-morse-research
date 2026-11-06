"""E10 timing, sequential, one process at a time. H2 main cohort, then H1 main cohort: clouds in seeded order,
BLOCKS blocks per cloud, every arm once per block in seeded order (arms in ONCE_PER_CLOUD only in block 1).
An arm that FAILED or TIMED OUT is not run again for that cloud. Then the capability set: every H2 arm except
NOQ once, in seeded order. No automatic rerun."""
import json
import random
import sys

from e10_common import (ARMS, BLOCKS, CAP, ONCE_PER_CLOUD, R, TIMEOUT, argv_for, capture, load, require, sha, utc,
                       verify_lock, verify_pins, write_new)


def arms_for(cloud):
    return [a for a, (ladder, _, _) in ARMS[cloud["experiment"]].items() if ladder is None or ladder in cloud["cuts"]]


def main():
    verify_pins()
    lock = verify_lock()
    cohort = load(R / "COHORT.json")["rows"]
    for c in cohort:
        require(sha(R / c["path"]) == c["sha256"], "cohort CSV changed: " + c["case"])
    require(not (R / "RUNS.jsonl").exists(), "E10 timing already started; no automatic rerun")
    plan = lock["timing"]
    rng = random.Random(plan["order_seed"])
    order = {}
    for experiment in ("H2", "H1"):
        rows = [c for c in cohort if c["experiment"] == experiment and c["role"] == "MAIN"]
        rng.shuffle(rows)
        order[experiment] = rows
    write_new(R / "TIMING-START.json", dict(utc=utc(), cloud_order={k: [c["case"] for c in v] for k, v in order.items()}))
    run_id = 0
    with (R / "RUNS.jsonl").open("x") as log:
        def run(cloud, arm, stage, block, position, timeout, cap):
            nonlocal run_id
            run_id += 1
            name = f"{run_id:05d}-{stage}-{cloud['case']}-b{block}-{arm}"
            rec = capture(argv_for(cloud, arm), R / "RAW" / name, timeout, cap)
            row = dict(run_id=run_id, case=cloud["case"], experiment=cloud["experiment"], role=cloud["role"],
                       family=cloud["family"], n=cloud["n"], arm=arm, stage=stage, block=block, position=position,
                       capture_directory=str((R / "RAW" / name).relative_to(R)), **rec)
            log.write(json.dumps(row) + "\n")
            log.flush()
            return row

        for experiment in ("H2", "H1"):
            for k, cloud in enumerate(order[experiment]):
                local = random.Random(plan["block_seed"] * 1000 + 100 * ("H2", "H1").index(experiment) + k)
                failed = set()
                for b in range(1, BLOCKS + 1):
                    block = [a for a in arms_for(cloud) if a not in failed
                             and (b == 1 or a not in ONCE_PER_CLOUD[experiment])]
                    local.shuffle(block)
                    for p, arm in enumerate(block):
                        row = run(cloud, arm, "MEASURED", b, p, TIMEOUT[experiment], CAP[experiment])
                        if row["status"] != "COMPLETE":
                            failed.add(arm)
                print("TIMED", cloud["case"], utc(), flush=True)
        crng = random.Random(plan["capability_seed"])
        for cloud in [c for c in cohort if c["role"] == "CAPABILITY"]:
            arms = [a for a in arms_for(cloud) if a != "NOQ"]
            crng.shuffle(arms)
            for p, arm in enumerate(arms):
                run(cloud, arm, "CAPABILITY", 0, p, TIMEOUT["CAPABILITY"], CAP["H2"])
            print("CAPABILITY", cloud["case"], utc(), flush=True)
    write_new(R / "TIMING-END.json", dict(utc=utc(), runs=run_id, runs_sha256=sha(R / "RUNS.jsonl")))


if __name__ == "__main__":
    sys.exit(main())
