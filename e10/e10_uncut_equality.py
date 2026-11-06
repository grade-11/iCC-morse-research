"""E10 descriptive (after the locked analysis): every complete PwPH run on a main cloud against the NOQ run of the
same cloud, emitted H1 (and H2 for H2 runs) compared exactly in binary64. Output: UNCUT-EQUALITY.json."""
import collections
import json

from e10_common import R, parse

runs = [json.loads(l) for l in (R / "RUNS.jsonl").read_text().splitlines()]
noq = {r["case"]: parse(R / r["capture_directory"] / "stdout.txt")["bars"] for r in runs if r["arm"] == "NOQ"}
count = collections.Counter()
for r in runs:
    if r["arm"] in ("NOQ", "RIPSER") or r["status"] != "COMPLETE" or r["case"] not in noq:
        continue
    bars = parse(R / r["capture_directory"] / "stdout.txt")["bars"]
    dims = (1, 2) if r["experiment"] == "H2" else (1,)
    count[f"{r['experiment']}:{r['arm']}:{'EQUAL' if all(bars[h] == noq[r['case']][h] for h in dims) else 'DIFFERENT'}"] += 1
(R / "UNCUT-EQUALITY.json").write_text(json.dumps(dict(sorted(count.items())), indent=1) + "\n")
print(dict(sorted(count.items())))
