"""Step 3 of the unified rule (PROTOCOL-step3.txt): oracle remeasured with B_X30 and development timing against
the draft (B_X22) and the authors' schedule. Clouds and ladders as extend_oracle.py (version 4).
Per cloud and ladder (H1: EVEN9 and P9; H2: P9 with the K series, EVEN9 with ALL, RULE, B22), LOG first, then
in a fixed shuffled order (seed 20261105):
  LOG    B_X30 --extend-all --extend-worth-log
  K<k>   B_X30 --extend-all --extend-commit-limit k, k = 0 .. K-1 (K = commits of LOG; n400 without K0)
  ALL    B_X30 --extend-all
  RULE   B_X30 --extend-all --extend-rule   (the unified controller U)
  B22    B_X22 on this ladder
Once per cloud, in a fixed shuffled order: H1 PAPER-B22 (B_X22, q20 q25 q30 q40 q50); sphere clouds: AUTHOR
(H1, H2) and AUTHOR-F (H2), B_X22 on 1.2, 1.7, endpoint.
Writes extend-oracle-v5-runs/, RUNS-extend-oracle-v5.jsonl. Resumable per cloud (DONE marker)."""
import json
import os
import random
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.argv = sys.argv[:1]
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R))
import extend_oracle as eo  # noqa: E402  clouds(), cuts, with_cuts, count, committed, parse, ENV

B30, B22 = str(R / "bin/B_X30"), str(R / "bin/B_X22")
OUT = R / "extend-oracle-v5-runs"
LOGFILE = R / "RUNS-extend-oracle-v5.jsonl"
TIMEOUT = {"H1": 600, "H2": 1800}
CONTROLLER_ONLY = {"--extend", "--absorb-pv", "--shrink-risky"}
NO_K0 = {"H2-author-s2-n400-d0"}


def run(argv, directory, timeout):
    directory.mkdir(parents=True, exist_ok=False)
    killed = []
    with (directory / "stdout.txt").open("w") as o, (directory / "stderr.txt").open("w") as e:
        start = time.perf_counter()
        proc = subprocess.Popen([str(a) for a in argv], stdout=o, stderr=e, env=eo.ENV, cwd=directory)
        timer = threading.Timer(timeout, lambda: (killed.append(1), proc.kill()))
        timer.start()
        _, code, usage = os.wait4(proc.pid, 0)
        seconds = time.perf_counter() - start
        timer.cancel()
    rc = os.waitstatus_to_exitcode(code)
    status = "TIMEOUT" if killed else "FAILED" if rc else "COMPLETE"
    (directory / "status.json").write_text(json.dumps(dict(
        argv=[str(a) for a in argv], status=status, seconds=seconds, maxrss_kib=usage.ru_maxrss,
        cpu_seconds=usage.ru_utime + usage.ru_stime), indent=1) + "\n")
    return status, seconds, usage.ru_maxrss


def author_argv(base, exp, filtered):
    """B_X22 argv on the authors' cuts: the controller-only flags removed; filters as in E8."""
    end = eo.cuts(base)[-1]
    assert float(end) > 1.7, end
    argv = eo.with_cuts(base, ["1.2", "1.7", end])
    flags_at = argv.index("--live-guard")
    keep = [a for a in argv[flags_at:] if a not in CONTROLLER_ONLY]
    if exp == "H1" or not filtered:
        keep = [a for a in keep if a not in ("--fiber-risk", "--fiber-full")]
    if not filtered:
        keep = [a for a in keep if a == "--live-guard"]
    return [B22] + argv[1:flags_at] + keep


def main():
    rng = random.Random(20261105)
    OUT.mkdir(exist_ok=True)
    with open(LOGFILE, "a") as log:
        for exp, case, base, ladders, extra in eo.clouds():
            top = OUT / exp / case
            if (top / "DONE").exists():
                continue
            if top.exists():
                raise SystemExit(f"partial cloud {top}: move it aside by hand")
            rows, ref = [], None

            def record(ladder, arm, d, status, seconds, rss):
                nonlocal ref
                row = dict(exp=exp, case=case, ladder=ladder, arm=arm, status=status, seconds=seconds,
                           gib=rss / 1048576)
                if status == "COMPLETE":
                    out = d / "stdout.txt"
                    row.update(commits=eo.count(out, "[postcommit]", eo.committed),
                               windows=eo.count(out, "window lower bound"))
                    bars = eo.parse(out)["bars"]
                    ref = ref or bars
                    row["exact"] = bars == ref
                rows.append(row)
                return row

            for ladder, full in (("EVEN9", exp == "H1"), ("P9", True)):
                argv = eo.with_cuts(base, ladders[ladder])
                d = top / ladder
                jobs = []
                if full:
                    row = record(ladder, "LOG", d / "LOG", *run([B30] + argv[1:] + ["--extend-all", "--extend-worth-log"],
                                                               d / "LOG", TIMEOUT[exp]))
                    if row["status"] == "COMPLETE":
                        K = row["commits"]
                        jobs += [(f"K{k}", B30, ["--extend-all", "--extend-commit-limit", str(k)])
                                 for k in range(K) if not (k == 0 and case in NO_K0)]
                jobs += [("ALL", B30, ["--extend-all"]), ("RULE", B30, ["--extend-all", "--extend-rule"]),
                         ("B22", B22, [])]
                rng.shuffle(jobs)
                for arm, binary, flags in jobs:
                    record(ladder, arm, d / arm, *run([binary] + argv[1:] + flags, d / arm, TIMEOUT[exp]))
            once = []
            if exp == "H1":
                once.append(("PAPER-B22", [B22] + extra["PAPER-B22"][1][1:]))
            sphere = "sphere-s2" in case if exp == "H1" else "author-s2" in case
            if sphere:
                once.append(("AUTHOR", author_argv(base, exp, exp == "H1")))
                if exp == "H2":
                    once.append(("AUTHOR-F", author_argv(base, exp, True)))
            rng.shuffle(once)
            for arm, argv in once:
                record("ONCE", arm, top / arm, *run(argv, top / arm, TIMEOUT[exp]))
            for r in rows:
                log.write(json.dumps(r) + "\n")
            log.flush()
            (top / "DONE").write_text("")
            bad = [(r["ladder"], r["arm"]) for r in rows if r["status"] != "COMPLETE" or not r.get("exact")]
            short = {f'{r["ladder"][:2]}-{r["arm"]}': round(r["seconds"], 2) for r in rows
                     if not r["arm"].startswith("K")}
            print(exp, case, short, "NOT EXACT/COMPLETE:" if bad else "", bad, flush=True)


if __name__ == "__main__":
    main()
