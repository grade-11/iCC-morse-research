"""Step 2 of the unified rule: what does each policy actually cost? (development data only)

Version 4 (B_X27: old simplices skipped before the witness search; straight to the endpoint once no more
commits can come). Version 3 (B_X26, extend-oracle-v3-runs/, stopped after 7 H1 clouds) re-enumerated every
simplex at each extension of a window with pseudo vertices. Version 2 (extend-oracle-v2-runs/, stopped after 3 clouds) timed runs with
subprocess.run(timeout=...), which polls about every 50 ms: too coarse for the small H1 clouds. Version 3 waits
with a blocking wait4 and also records the program's own "run time". Version 1 (B_X25, RUNS-extend-oracle.jsonl, extend-oracle-runs/) was stopped after 46 H1
clouds: there the commit limit still selected supports after the limit, and the rule took its rate for a new
window from the small first build. Kept for reference.

Every window is retained and extended through the cuts where it does not commit (--extend-all). "Commit k
times, then only extend" (--extend-commit-limit k) gives, for each k, the true cost T(k) of the policy that stops
committing after k commits; T(j) - T(j+1) is the true saving of committing at the (j+1)-th admitted cut and
never again, against never again, which is what the rule (--extend-rule, a stop rule) predicts. The cheapest k is
an oracle the rule cannot know in advance.

Two ladders, both nine cuts: EVEN9 = q10, q20, ..., q90; P9 = q20, q25, q30, q40, q50, q60, q70, q80, q90 (the
paper's H1 cuts, extended; in H2 the B_X7 development measurement found it as fast as EVEN9).
Per cloud and ladder, one run at a time, LOG first, then the others in a fixed shuffled order:
  LOG    --extend-all --extend-worth-log (commits at every admitted cut; the rule's inputs at each cut)
  K<k>   --extend-all --extend-commit-limit k, k = 0 .. K-1 (K = commits of LOG)
  ALL    --extend-all
  RULE   --extend-all --extend-rule
  B22    B_X22, the controller of E8 (rebuilds the window at a vetoed cut after the first commit)
H1: the 96 development clouds, both ladders with all arms; once per cloud also PAPER-B22 (the paper's five cuts,
    B_X22: the E8 H1 controller). Ladders from the argv of shrink-timing-h1-v3-runs (*-ALL-AS9F, *-ALL-ASF).
H2: the 17 development clouds, smallest first: P9 with all arms; EVEN9 with ALL, RULE and B22 only. Ladders from
    shrink-timing-h2-v3-runs/*-ALL-ASF (EVEN9) and commit-once-runs/*-H1UP-CONT (P9).
Every run of a cloud must give the same H1 and H2 bars as its first LOG run. Writes RUNS-extend-oracle-v2.jsonl.
A development measurement, not a locked timing experiment; nothing else running. Run with the venv python."""
import json
import os
import random
import threading
import re
import subprocess
import sys
import time
from pathlib import Path

R = Path(__file__).resolve().parent
sys.path.insert(0, "/home/alex_1709/projects/iCC-team/experiments/pwph-e8-controller-confirm-20261104")
from e8_common import parse  # noqa: E402  (reads outputs only)

B22, B26 = str(R / "bin/B_X22"), str(R / "bin/B_X27")
OUT = R / "extend-oracle-v4-runs"
LOGFILE = R / "RUNS-extend-oracle-v4.jsonl"
ENV = {"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C", "OMP_NUM_THREADS": "4", "OMP_DYNAMIC": "FALSE",
       "OMP_THREAD_LIMIT": "4", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
TIMEOUT = {"H1": 600, "H2": 1800}


def argv_of(path):
    return json.loads((Path(path) / "status.json").read_text())["argv"]


def cuts(argv):
    i = argv.index("--eps-breaks")
    out = []
    for a in argv[i + 1:]:
        if a.startswith("--"):
            break
        out.append(a)
    return out


def with_cuts(argv, new):
    i = argv.index("--eps-breaks")
    j = i + 1 + len(cuts(argv))
    return argv[: i + 1] + new + argv[j:]


def clouds():
    """(exp, case, base argv with EVEN9 cuts, {ladder: cuts}, extra arms)"""
    out = []
    for d in sorted((R / "shrink-timing-h1-v3-runs").glob("*-ALL-AS9F")):
        case = d.name[: -len("-ALL-AS9F")]
        even = argv_of(d)
        paper = argv_of(R / "shrink-timing-h1-v3-runs" / (case + "-ALL-ASF"))
        e, p = cuts(even), cuts(paper)
        assert float(e[4]) == float(p[4]) and float(e[-1]) == float(p[-1]), case
        out.append(("H1", case, even, {"EVEN9": e, "P9": p[:5] + e[5:]}, {"PAPER-B22": (B22, paper)}))
    h2 = sorted((R / "shrink-timing-h2-v3-runs").glob("*-ALL-ASF"),
                key=lambda d: (int(re.search(r"-n(\d+)", d.name)[1]), d.name))
    for d in h2:
        case = d.name[: -len("-ALL-ASF")]
        even = argv_of(d)
        p9 = cuts(argv_of(R / "commit-once-runs" / (case + "-H1UP-CONT")))
        assert float(p9[-1]) == float(cuts(even)[-1]), case
        out.append(("H2", case, even, {"EVEN9": cuts(even), "P9": p9}, {}))
    return out


def run(argv, directory, timeout):
    """Whole-process wall clock: perf_counter around Popen and a blocking wait4; a timer kills on timeout."""
    directory.mkdir(parents=True, exist_ok=False)
    killed = []
    with (directory / "stdout.txt").open("w") as o, (directory / "stderr.txt").open("w") as e:
        start = time.perf_counter()
        proc = subprocess.Popen([str(a) for a in argv], stdout=o, stderr=e, env=ENV, cwd=directory)
        timer = threading.Timer(timeout, lambda: (killed.append(1), proc.kill()))
        timer.start()
        _, code, usage = os.wait4(proc.pid, 0)
        seconds = time.perf_counter() - start
        timer.cancel()
        proc.returncode = os.waitstatus_to_exitcode(code)
    status = "TIMEOUT" if killed else "FAILED" if proc.returncode else "COMPLETE"
    lines = (directory / "stdout.txt").read_text().splitlines()
    own_ms = float(lines[-1].split("=")[1].split()[0]) if lines and lines[-1].startswith("run time = ") else None
    (directory / "status.json").write_text(json.dumps(dict(argv=[str(a) for a in argv], status=status,
                                                           seconds=seconds, own_ms=own_ms,
                                                           cpu_seconds=usage.ru_utime + usage.ru_stime), indent=1) + "\n")
    return status, seconds


def count(path, prefix, cond=lambda l: True):
    return sum(l.startswith(prefix) and cond(l) for l in path.read_text().splitlines())


def committed(line):
    return " committed=0 " not in line


def main():
    rng = random.Random(20261104)
    with open(LOGFILE, "a") as log:
        for exp, case, base, ladders, extra in clouds():
            top = OUT / exp / case
            if (top / "DONE").exists():
                continue
            rows, ref = [], None
            plans = [("EVEN9", exp == "H1"), ("P9", True)]
            for ladder, full in plans:
                argv = with_cuts(base, ladders[ladder])
                d = top / ladder
                jobs = []
                if full:
                    status, seconds = run([B26] + argv[1:] + ["--extend-all", "--extend-worth-log"], d / "LOG",
                                          TIMEOUT[exp])
                    row = dict(exp=exp, case=case, ladder=ladder, arm="LOG", status=status, seconds=seconds)
                    rows.append(row)
                    if status != "COMPLETE":
                        continue
                    out = d / "LOG" / "stdout.txt"
                    K = count(out, "[postcommit]", committed)
                    row.update(commits=K, windows=count(out, "window lower bound"))
                    bars = parse(out)["bars"]
                    ref = ref or bars
                    row["exact"] = bars == ref
                    jobs += [(f"K{k}", B26, ["--extend-all", "--extend-commit-limit", str(k)]) for k in range(K)]
                jobs += [("ALL", B26, ["--extend-all"]), ("RULE", B26, ["--extend-all", "--extend-rule"]),
                         ("B22", B22, [])]
                rng.shuffle(jobs)
                for arm, binary, flags in jobs:
                    status, seconds = run([binary] + argv[1:] + flags, d / arm, TIMEOUT[exp])
                    row = dict(exp=exp, case=case, ladder=ladder, arm=arm, status=status, seconds=seconds)
                    if status == "COMPLETE":
                        out = d / arm / "stdout.txt"
                        row.update(commits=count(out, "[postcommit]", committed),
                                   windows=count(out, "window lower bound"))
                        bars = parse(out)["bars"]
                        ref = ref or bars
                        row["exact"] = bars == ref
                    rows.append(row)
            for arm, (binary, argv) in extra.items():
                status, seconds = run([binary] + argv[1:], top / arm, TIMEOUT[exp])
                row = dict(exp=exp, case=case, ladder="PAPER", arm=arm, status=status, seconds=seconds)
                if status == "COMPLETE":
                    out = top / arm / "stdout.txt"
                    row.update(commits=count(out, "[postcommit]", committed), windows=count(out, "window lower bound"),
                               exact=parse(out)["bars"] == ref)
                rows.append(row)
            for r in rows:
                log.write(json.dumps(r) + "\n")
            log.flush()
            (top / "DONE").write_text("")
            bad = [(r["ladder"], r["arm"]) for r in rows if r["status"] != "COMPLETE" or not r.get("exact")]
            short = {f'{r["ladder"][:2]}-{r["arm"]}': round(r["seconds"], 2) for r in rows}
            print(exp, case, short, "NOT EXACT/COMPLETE:" if bad else "", bad, flush=True)


if __name__ == "__main__":
    main()
