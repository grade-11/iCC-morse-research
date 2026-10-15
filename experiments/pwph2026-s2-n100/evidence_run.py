#!/usr/bin/env python3
"""Capture a command in a NEW evidence directory; preserve failures and timeouts.
Example: python evidence_run.py --output runs/check --input data.csv -- python check.py
JSON stdout is also saved as report.json. Non-JSON commands get a capture summary.
The wrapper exits with the child status (124 timeout, 125 launch/capture error).
"""
import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shlex
import signal
import subprocess
import sys
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')


def capture(output, argv, inputs=(), timeout=600):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    context = dict(argv=argv, command=shlex.join(argv), cwd=str(Path.cwd()),
                   started_at=now(), timeout_s=timeout, python=sys.version,
                   python_executable=sys.executable, platform=platform.platform(),
                   wrapper_sha256=sha(__file__), input_sha256={},
                   environment={k: os.environ.get(k) for k in
                       ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OMP_DYNAMIC', 'LC_ALL')},
                   packages={})
    for name in ('numpy', 'ripser', 'scipy', 'scikit-learn'):
        try:
            context['packages'][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    for key, cmd in {'head': ['git', 'rev-parse', 'HEAD'],
                     'branch': ['git', 'branch', '--show-current']}.items():
        r = subprocess.run(cmd, capture_output=True, text=True)
        context[key] = dict(exit_status=r.returncode, stdout=r.stdout, stderr=r.stderr)
    (output/'command.txt').write_text(shlex.join(argv)+'\n')
    (output/'stdout.log').touch()
    (output/'stderr.log').touch()
    write(output/'context.json', context)
    start = time.perf_counter()
    rc, error = 125, None
    try:
        context['input_sha256'] = {str(Path(p).resolve()): sha(p) for p in inputs}
        write(output/'context.json', context)
        with (output/'stdout.log').open('wb') as out, (output/'stderr.log').open('wb') as err:
            child = subprocess.Popen(argv, stdout=out, stderr=err, start_new_session=True)
            try:
                rc = child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                rc, error = 124, 'Command timed out; process group killed'
        context['input_sha256_after'] = {str(Path(p).resolve()): sha(p) for p in inputs}
        if context['input_sha256_after'] != context['input_sha256']:
            rc, error = 125, 'Input changed during command'
    except Exception as exc:
        rc, error = 125, str(exc)
    if error:
        with (output/'stderr.log').open('a') as err:
            err.write('\nEvidence capture: '+error+'\n')
    context.update(ended_at=now(), elapsed_monotonic_s=time.perf_counter()-start,
                   exit_status=rc, error=error)
    write(output/'context.json', context)
    (output/'exit-status.txt').write_text(str(rc)+'\n')
    try:
        report = json.loads((output/'stdout.log').read_text())
    except (OSError, ValueError):
        report = dict(status='CAPTURED' if rc == 0 else 'ERROR', exit_status=rc,
                      error=error, note='Command has no JSON stdout; inspect stdout.log/stderr.log')
    write(output/'report.json', report)
    return rc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--input', action='append', default=[], type=Path)
    parser.add_argument('--timeout', type=float, default=600)
    parser.add_argument('argv', nargs=argparse.REMAINDER)
    a = parser.parse_args()
    cmd = a.argv[1:] if a.argv[:1] == ['--'] else a.argv
    if not cmd or a.timeout <= 0:
        parser.error('A command and positive timeout are required')
    rc = capture(a.output, cmd, a.input, a.timeout)
    print(json.dumps(dict(output=str(a.output.resolve()), exit_status=rc)))
    return rc if rc >= 0 else 128-rc


if __name__ == '__main__':
    sys.exit(main())
