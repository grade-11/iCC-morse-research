"""E10 preparation: native binary64 metric -> E and q_p (P1 rank rule); ladders per cloud. No PH, no timing."""
import json
import math
import struct
import sys

from e10_common import (AUTHOR_WINDOWS, CAP, EVEN9, METRIC, P9, PAPER, PCTS, R, capture, d17, load, require, sha, utc,
                       verify_lock, verify_pins, word, write_new)


def metric_values(path, n):
    values, footer = [], None
    with path.open() as f:
        header = json.loads(next(f))
        require(header["schema"] == "ICC_ORIGINAL_METRIC_BITS_V1" and header["parsed_n"] == n,
                "metric header differs")
        for line in f:
            row = json.loads(line)
            if row["type"] == "complete":
                footer = row
                continue
            values.append(struct.unpack(">d", bytes.fromhex(row["bits_hex16"]))[0])
    require(footer is not None and footer["actual_records"] == n * (n - 1) // 2 and footer["all_finite"],
            "metric footer differs")
    return values, footer["zero_pairs"]


def main():
    verify_pins()
    verify_lock()
    cohort_lock = load(R / "COHORT-LOCK.json")
    require(not (R / "COHORT.json").exists(), "E10 preparation already ran")
    out = []
    (R / "prep").mkdir(exist_ok=True)
    for row in cohort_lock["rows"]:
        require(sha(R / row["path"]) == row["sha256"], "cohort CSV changed: " + row["case"])
        metric = R / "prep" / (row["case"] + "-metric.jsonl")
        rec = capture([METRIC, "--file-name", R / row["path"], "--threads", "4",
                       "--expect-n", str(row["n"]), "--expect-dim", str(row["ambient_dim"]),
                       "--output", metric], R / "RAW" / ("PREP-" + row["case"]), 600, CAP["H1"])
        require(rec["status"] == "COMPLETE", "metric export failed: " + row["case"])
        values, zero_pairs = metric_values(metric, row["n"])
        metric_sha = sha(metric)
        metric.unlink()  # large; its hash is kept and the helper is deterministic
        endpoint = math.nextafter(max(values), math.inf)
        ds = sorted({v for v in values if 0 < v < endpoint})
        q = {}
        for pc in PCTS:
            rank = (pc * len(ds) + 99) // 100
            value = math.nextafter(ds[rank - 1], math.inf)
            require(0 < value < endpoint, "q outside (0,E)")
            q[str(pc)] = dict(rank_1based=rank, decimal17=d17(value), bits=word(value))
        c = dict(row, E_decimal17=d17(endpoint), E_bits=word(endpoint), q=q, zero_pairs=zero_pairs,
                 distinct_positive=len(ds), metric_sha256=metric_sha)
        cuts = dict(NOQ=[c["E_decimal17"]],
                    EVEN9=[q[str(p)]["decimal17"] for p in EVEN9] + [c["E_decimal17"]],
                    PAPER=[q[str(p)]["decimal17"] for p in PAPER] + [c["E_decimal17"]],
                    P9=[q[str(p)]["decimal17"] for p in P9] + [c["E_decimal17"]])
        if row["family"] in AUTHOR_WINDOWS:
            w = list(AUTHOR_WINDOWS[row["family"]])
            if endpoint <= 2.0:
                w[-1] = c["E_decimal17"]
            else:
                w.append(c["E_decimal17"])
            require(all(float(a) < float(b) for a, b in zip(w, w[1:])), "author windows not increasing")
            cuts["AUTHOR"] = w
        c["cuts"] = cuts
        out.append(c)
        print("PREP", row["case"], c["E_decimal17"], "zero_pairs", zero_pairs, flush=True)
    write_new(R / "COHORT.json", dict(utc=utc(), script_sha256=sha(__file__), rows=out,
                                      rule="distinct positive d<E, rank=ceil(p*m/100), "
                                           "q=nextafter(d_rank,+inf); E=nextafter(max d,+inf)"))


if __name__ == "__main__":
    sys.exit(main())
