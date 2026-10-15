#!/usr/bin/env python3
"""Regression checks for archive comparison/capture; no PH computation."""
import argparse
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
from evidence_run import capture, write


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', required=True, type=Path)
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    base = Path(__file__).resolve().parent
    old = base/'runs/whole-run-monotonic-ma523yck.tar.gz'
    ref = base/'runs/ripser-reference-noa_bg6w.tar.gz'
    parser = base/'benchmark_icc.py'
    comparator = base/'compare_ripser_n100.py'
    rows = []
    with tarfile.open(old) as t:
        entries = [(m.name, t.extractfile(m).read()) for m in t if m.isfile()]
    target = next(n for n, _ in entries if n.endswith('round-01/ph/program.log'))
    for name, expected in [('success', 0), ('endpoint-mismatch', 1),
                           ('multiplicity-mismatch', 1), ('missing-log', 2),
                           ('truncated-log', 2), ('parser-pin', 2)]:
        archive = a.output/(name+'.tar.gz')
        with tarfile.open(archive, 'w:gz') as t:
            for n, data in entries:
                if n == target:
                    if name == 'missing-log':
                        continue
                    text = data.decode()
                    # Perturb one complete valid record, keeping config/provenance.
                    line = next(x for x in text.splitlines() if x.strip().startswith('(') and ',' in x)
                    if name == 'endpoint-mismatch':
                        import re
                        birth = re.search(r'\(([^,]+),', line)[1]
                        from decimal import Decimal
                        text = text.replace(line, line.replace(birth, str(Decimal(birth)/2)), 1)
                    elif name == 'multiplicity-mismatch':
                        text = text.replace(line, line+'\n'+line, 1)
                    elif name == 'truncated-log':
                        text = text.split('run time =')[0]
                    data = text.encode()
                m = tarfile.TarInfo(n); m.size = len(data); t.addfile(m, io.BytesIO(data))
        pinned = parser
        if name == 'parser-pin':
            pinned = a.output/'unreviewed-parser.py'
            pinned.write_bytes(parser.read_bytes()+b'\n# test change\n')
        cmd = [sys.executable, str(comparator), '--icc-archive', str(archive),
               '--ripser-archive', str(ref), '--benchmark-script', str(pinned)]
        rc = capture(a.output/name, cmd, (archive, ref, pinned, comparator), 60)
        rows.append(dict(test=name, expected_exit=expected, actual_exit=rc, pass_=rc==expected))
    for name, cmd, expected, timeout in [
            ('capture-json', [sys.executable, '-c', 'print("{\\"status\\":\\"PASS\\"}")'], 0, 5),
            ('capture-nonzero', [sys.executable, '-c', 'import sys; print("failure", file=sys.stderr); sys.exit(7)'], 7, 5),
            ('capture-timeout', [sys.executable, '-c', 'import time; time.sleep(5)'], 124, 0.1)]:
        rc = capture(a.output/name, cmd, (), timeout)
        rows.append(dict(test=name, expected_exit=expected, actual_exit=rc, pass_=rc==expected))
    try:
        capture(a.output/'capture-json', [sys.executable, '-c', 'print(1)'])
        preserved = False
    except FileExistsError:
        preserved = True
    rows.append(dict(test='no-overwrite', pass_=preserved))
    rc = capture(a.output/'capture-missing-input', [sys.executable, '-c', 'print(1)'],
                 [a.output/'absent-input.csv'])
    rows.append(dict(test='capture-missing-input', expected_exit=125, actual_exit=rc,
                     pass_=rc==125 and (a.output/'capture-missing-input/stderr.log').is_file()
                     and (a.output/'capture-missing-input/stdout.log').is_file()))
    report = dict(status='PASS' if all(r['pass_'] for r in rows) else 'FAIL', tests=rows)
    write(a.output/'report.json', report)
    print(json.dumps(report, indent=2))
    return 0 if report['status']=='PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
