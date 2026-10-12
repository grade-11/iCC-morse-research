#!/usr/bin/env python3
"""Compare current iCC verbose ph/piecewise logs, H1..H(maxdim-1).

Standard library only. Exact Decimal endpoint multisets; no tolerance or
Ripser/float32 conversion. Context checks use recorded checksums, not a reread
of the original dataset. This is not a proof of algorithmic correctness.
Exit: 0 match, 1 PH mismatch, 2 invalid/incompatible inputs.
"""
import argparse
from collections import Counter
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import re
import shlex
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(text):
    value = Decimal(text)
    require(value.is_finite(), "Non-finite endpoint or parameter")
    return value


def one(pattern, text, name):
    matches = re.findall(pattern, text, re.M)
    require(len(matches) == 1, f"Expected exactly one {name}")
    return matches[0]


def context(path, tool):
    text = path.read_text()
    command = shlex.split(one(r"^Command: (.+)$", text, "Command"))
    require(command and command[0] == "./build/iCC", "Unexpected executable")
    options = {}
    i = 1
    while i < len(command):
        key = command[i]
        require(key.startswith("--") and key not in options,
                "Malformed or duplicate command option")
        i += 1
        values = []
        while i < len(command) and not command[i].startswith("--"):
            values.append(command[i])
            i += 1
        options[key] = values

    def scalar(key):
        require(len(options.get(key, [])) == 1, f"Missing/invalid {key}")
        return options[key][0]

    require(scalar("--tool") == tool, "Wrong tool in context")
    require(options.get("--verbose") == [], "Verbose log required")
    allowed = {"--tool", "--file-name", "--max-dim", "--threads", "--verbose"}
    allowed |= ({"--max-eps"} if tool == "ph" else
                {"--eps-breaks", "--pv-cap-scale"})
    require(set(options) <= allowed, "Unsupported options; review parser first")
    dataset = scalar("--file-name")
    checksums = re.findall(r"^([0-9a-f]{64})  (.+)$", text, re.M)

    def checksum(name):
        found = [sha for sha, filename in checksums if filename == name]
        require(len(found) == 1, f"Missing/ambiguous checksum: {name}")
        return found[0]

    dim = int(scalar("--max-dim"))
    require(dim >= 2, "No supported H1+ output for max-dim < 2")
    threads = int(scalar("--threads"))
    require(threads > 0, "Invalid thread count")
    breaks = ([number(scalar("--max-eps"))] if tool == "ph" else
              [number(v) for v in options.get("--eps-breaks", [])])
    require(breaks and all(a < b for a, b in zip([Decimal(0)] + breaks, breaks)),
            "Epsilon bounds must be positive and increasing")
    alpha = number(scalar("--pv-cap-scale")) if tool == "piecewise" else None
    return dict(dataset=dataset, dataset_sha=checksum(dataset),
                binary_sha=checksum("./build/iCC"), dim=dim, threads=threads,
                breaks=breaks, alpha=alpha)


def parse_log(path, meta, tool):
    text = path.read_text()
    require(re.findall(r"^(?:iCC )?exit status: (\d+)\s*$", text, re.M) == ["0"],
            f"{path}: missing, repeated, or unsuccessful exit status")
    require(re.fullmatch(r"(?:iCC )?exit status: 0", text.strip().splitlines()[-1]),
            "Exit status must be the final nonblank line")
    one(r"^run time = (\d+)\s*$", text, "runtime completion marker")
    expected_name = "morsePH" if tool == "ph" else "morsePiecewisePH"
    require(one(r"^Run (\w+) with:$", text, "tool header") == expected_name,
            "Context/log tool mismatch")
    require(one(r"^  input file: (.+)$", text, "input header") == meta["dataset"],
            "Context/log dataset path mismatch")
    require(int(one(r"^  max dimension: (\d+)$", text, "dimension")) == meta["dim"],
            "Context/log dimension mismatch")
    require(int(one(r"^  thread number: (\d+)$", text, "threads")) == meta["threads"],
            "Context/log thread mismatch")
    if tool == "ph":
        bounds = [number(one(r"^  max epsilon: (\S+)$", text, "cutoff"))]
    else:
        bounds = [number(v) for v in one(
            r"^  epsilon breaks: (.+)$", text, "breaks").split()]
        require(number(one(r"^  pv cap scale: (\S+)$", text, "cap scale"))
                == meta["alpha"], "Context/log alpha mismatch")
        require(number(one(r"^  pv minimum separation scale: (\S+)$", text,
                           "separation")) == 0, "Only separation=0 supported")
    require(bounds == meta["breaks"], "Context/log epsilon mismatch")
    windows = list(zip([Decimal(0)] + bounds[:-1], bounds))
    expected = [(lo, hi, d) for lo, hi in windows
                for d in range(2, meta["dim"] + 1)]
    blocks = []
    diagrams = {d: Counter() for d in range(1, meta["dim"])}
    temporary = 0
    current = None
    records = 0
    empty = False
    in_pairs = False
    pair_header_seen = False

    def complete_block():
        if current is not None:
            require(pair_header_seen, "Missing persistent-pairs header")
            require((records > 0) != empty,
                    "Each section must have pairs OR an explicit empty marker")

    for line in text.splitlines():
        base = re.fullmatch(r"Processing dim (\d+), cofacet number: \d+, facet number: \d+", line)
        piece = re.fullmatch(r"in eps range (\S+)\s+(\S+)\s+dimension = (\d+)\s+cofacet num = \d+\s+facet num = \d+", line)
        if base or piece:
            complete_block()
            require(bool(base) == (tool == "ph"), "Unexpected section type")
            current = ((Decimal(0), bounds[-1], int(base[1])) if base else
                       (number(piece[1]), number(piece[2]), int(piece[3])))
            blocks.append(current)
            require(blocks == expected[:len(blocks)], "Missing/repeated/out-of-order section")
            records, empty = 0, False
            in_pairs, pair_header_seen = False, False
            continue
        stripped = line.strip()
        if stripped in {"persistent pairs:", "dimensional persistent pairs:"}:
            require(current is not None and not pair_header_seen,
                    "Unexpected persistent-pairs header")
            in_pairs, pair_header_seen = True, True
            continue
        if not in_pairs:
            continue
        if not stripped or stripped.startswith(("birth facet labels:", "death cofacet labels:")):
            continue
        if not stripped.startswith("("):
            in_pairs = False
            continue
        require(current is not None, "Pair outside PH section")
        if stripped in {"(empty)", "(no new interval or surviving interval from previous eps range)"}:
            require(not empty and records == 0, "Conflicting empty marker")
            empty = True
            continue
        pair = re.fullmatch(r"\(([^,]+),\s*([^,]+)\)", stripped)
        require(pair is not None and not empty, f"Malformed pair: {stripped}")
        birth, death = number(pair[1]), number(pair[2])
        lo, hi, dim = current
        require(0 <= birth <= hi and (death == -1 or birth <= death <= hi),
                "Invalid/out-of-range interval")
        records += 1
        if death == -1 and tool == "piecewise" and hi != bounds[-1]:
            temporary += 1
            continue
        require(death == -1 or tool == "ph" or death >= lo,
                "Finite death predates piecewise window")
        diagrams[dim - 1][(birth, death)] += 1
    complete_block()
    require(blocks == expected, "Incomplete PH sections")
    return diagrams, temporary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--baseline-context", required=True, type=Path)
    parser.add_argument("--piecewise", required=True, type=Path)
    parser.add_argument("--piecewise-context", required=True, type=Path)
    args = parser.parse_args()
    try:
        baseline = context(args.baseline_context, "ph")
        piecewise = context(args.piecewise_context, "piecewise")
        for key in ["dataset_sha", "binary_sha", "dim", "threads"]:
            require(baseline[key] == piecewise[key], f"Incompatible {key}")
        require(baseline["breaks"][-1] == piecewise["breaks"][-1], "Different cutoffs")
        a, _ = parse_log(args.baseline, baseline, "ph")
        b, temporary = parse_log(args.piecewise, piecewise, "piecewise")
        passed = a == b

        def entries(counter):
            return [dict(birth=str(birth), death=str(death), multiplicity=count)
                    for (birth, death), count in sorted(counter.items())]

        report = dict(status="PASS" if passed else "FAIL", comparison="exact Decimal multisets",
                      dataset_sha256=baseline["dataset_sha"], binary_sha256=baseline["binary_sha"],
                      cutoff=str(baseline["breaks"][-1]),
                      ignored_intermediate_survivor_records=temporary,
                      notes=["H0 not checked", "Ripser not checked",
                             "death=-1 means surviving at cutoff, not proven infinite",
                             "Checksums are validated from supplied contexts"],
                      input_sha256={k: hashlib.sha256(v.read_bytes()).hexdigest()
                                    for k, v in vars(args).items()},
                      dimensions={f"H{d}": dict(baseline_count=sum(a[d].values()),
                                                piecewise_count=sum(b[d].values()),
                                                missing_in_piecewise=entries(a[d] - b[d]),
                                                extra_in_piecewise=entries(b[d] - a[d])) for d in a})
        print(json.dumps(report, indent=2))
        return 0 if passed else 1
    except (ValueError, OSError, InvalidOperation) as error:
        print(json.dumps(dict(status="ERROR", message=str(error)), indent=2))
        return 2


if __name__ == "__main__":
    sys.exit(main())
