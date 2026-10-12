#!/usr/bin/env python3
"""Fixed S2/n100 experiment: run or analyze nonverbose iCC batches.

Python standard library only. Run from the iCC repository root.
analyze never runs iCC and never edits the supplied batch.
Exit codes: 0 = PH comparison passed; 1 = PH mismatch; 2 = invalid input.
H0 and Ripser are NOT checked. Decimal multisets compare printed endpoints.
"""
import argparse
from collections import Counter
import datetime
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import statistics
import subprocess
import sys
import tempfile
import time

DATA = 'experiments/pwph2026-s2-n100/data/sphere_S2_R3_n100_seed20261007.csv'
DATA_SHA = 'e39f71c26b3f71c96860c04316b9eae7ec192891ebdd39c80bbebb2029f76b48'
RUNS = Path('experiments/pwph2026-s2-n100/runs')
BOUNDS = [(Decimal('0'), Decimal('1.2')),
          (Decimal('1.2'), Decimal('1.7')), (Decimal('1.7'), Decimal('2'))]


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def command(mode):
    common = ['./build/iCC', '--tool', mode, '--file-name', DATA,
              '--max-dim', '3', '--threads', '4']
    return common + (['--max-eps', '2.0'] if mode == 'ph' else
                     ['--eps-breaks', '1.2', '1.7', '2.0', '--pv-cap-scale', '0.85'])


def decimal(text):
    x = Decimal(text)
    need(x.is_finite(), 'Non-finite endpoint')
    return x


def parse_log(path, mode):
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    need(len(lines) >= 3 and lines[-1] == 'iCC exit status: 0',
         f'{path}: missing successful completion')
    runtime = re.fullmatch(r'run time = (\d+)', lines[-2])
    need(runtime is not None, f'{path}: missing runtime marker')
    lines = lines[:-2]
    diagrams = {1: Counter(), 2: Counter()}
    ignored = 0

    def add(dim, birth, death, lo, hi):
        nonlocal ignored
        need(dim in (2, 3), 'Only H1/H2 supported')
        b, d = decimal(birth), decimal(death)
        need(0 <= b <= hi and (d == -1 or b <= d <= hi), 'Invalid interval')
        need(d == -1 or d >= lo, 'Death predates current window')
        if d == -1 and mode == 'piecewise' and hi != 2:
            ignored += 1
        else:
            diagrams[dim - 1][(b, d)] += 1

    if mode == 'ph':
        need(lines[0] == 'total point number: 100', 'Unexpected point count')
        pos = 1
        for dim in (2, 3):
            need(pos < len(lines) and re.fullmatch(
                rf'Processing dim {dim}, cofacet number: \d+, facet number: \d+', lines[pos]),
                'Missing or out-of-order baseline section')
            pos += 1
            need(pos < len(lines) and lines[pos] == 'dimensional persistent pairs:',
                 'Missing baseline pair header')
            pos += 1
            count = 0
            if pos < len(lines) and lines[pos] == '(empty)':
                pos += 1
                count = 1
            else:
                while pos < len(lines):
                    pair = re.fullmatch(r'\(([^,]+),\s*([^,]+)\)', lines[pos])
                    if not pair:
                        break
                    add(dim, *pair.groups(), Decimal(0), Decimal(2))
                    pos += 1
                    count += 1
            need(count > 0, 'Missing pairs or empty marker')
        need(pos == len(lines), 'Unexpected baseline content')
    else:
        pos = 0
        for lo, hi in BOUNDS:
            need(pos + 2 < len(lines) and lines[pos] == 'window lower bound, window upper bound',
                 'Missing piecewise window header')
            parts = lines[pos + 1].split(',')
            need(len(parts) == 2 and tuple(map(decimal, parts)) == (lo, hi),
                 'Unexpected window bounds')
            need(lines[pos + 2] == 'dimension, birth weight, death weight',
                 'Missing piecewise pair header')
            pos += 3
            last_dim = 2
            while pos < len(lines) and lines[pos] != 'window lower bound, window upper bound':
                pair = re.fullmatch(r'([23]),\s*([^,]+),\s*([^,]+)', lines[pos])
                need(pair is not None, 'Unexpected piecewise content')
                dim = int(pair[1])
                need(dim >= last_dim, 'Out-of-order pair dimension')
                last_dim = dim
                add(dim, pair[2], pair[3], lo, hi)
                pos += 1
        need(pos == len(lines), 'Extra piecewise window/content')
    return diagrams, ignored, int(runtime[1])


def numeric(value, name, positive=False):
    need(not isinstance(value, bool), f'Invalid {name}')
    x = float(value)
    need(math.isfinite(x) and (x > 0 if positive else x >= 0), f'Invalid {name}')
    return x


def stats(values):
    return dict(median=statistics.median(values), minimum=min(values), maximum=max(values))


def analyze(root):
    meta = read_json(root / 'context.json')
    clock = meta['clock']
    need(clock['monotonic'] is True and clock['adjustable'] is False,
         'Monotonic clock metadata required')
    checks = meta['checksums']
    need(checks == read_json(root / 'checksums-after.json'), 'Checksums changed during batch')
    need(checks.get(DATA) == DATA_SHA, 'Wrong fixed dataset checksum')
    need(re.fullmatch('[0-9a-f]{64}', checks.get('./build/iCC', '')), 'Invalid binary checksum')
    folders = sorted(root.glob('round-*'))
    need(len(folders) >= 2 and [p.name for p in folders] ==
         [f'round-{i:02d}' for i in range(len(folders))], 'Incomplete round sequence')
    if 'measured_rounds' in meta:
        need(len(folders) == meta['measured_rounds'] + 1, 'Missing rounds')
    else:
        need(len(folders) == 11, 'Legacy batch must contain 11 rounds')
    rows, mismatches, warnings = [], [], []
    reference = None
    evidence = {'context.json': sha(root / 'context.json'),
                'checksums-after.json': sha(root / 'checksums-after.json')}
    for i, folder in enumerate(folders):
        row = {'round': i, 'warmup': i == 0}
        for mode in ('ph', 'piecewise'):
            sub = folder / mode
            for name in ('command.txt', 'program.log', 'measurement.json', 'resources.txt'):
                evidence[str((sub / name).relative_to(root))] = sha(sub / name)
            need(shlex.split((sub / 'command.txt').read_text()) == command(mode),
                 f'{sub}: command differs from fixed configuration')
            diagram, ignored, runtime = parse_log(sub / 'program.log', mode)
            if reference is None:
                reference = diagram
            for dim in (1, 2):
                if diagram[dim] != reference[dim]:
                    def pairs(counter):
                        return [{'birth': str(b), 'death': str(d), 'multiplicity': n}
                                for (b, d), n in sorted(counter.items())]
                    mismatches.append(dict(round=i, mode=mode, dimension=dim,
                        missing=pairs(reference[dim] - diagram[dim]),
                        extra=pairs(diagram[dim] - reference[dim])))
            measurement = read_json(sub / 'measurement.json')
            items = [line.split('=', 1) for line in (sub / 'resources.txt').read_text().splitlines()]
            need(all(len(x) == 2 for x in items), 'Malformed resources')
            resources = dict(items)
            need(len(resources) == len(items) and set(resources) ==
                 {'elapsed_s', 'user_cpu_s', 'system_cpu_s', 'max_rss_kib', 'exit_status'},
                 'Missing/duplicate resource fields')
            need(measurement['exit_status'] == 0 and resources['exit_status'] == '0',
                 'Unsuccessful run')
            elapsed = numeric(measurement['elapsed_monotonic_s'], 'monotonic duration', True)
            real = float(measurement['elapsed_realtime_s'])
            gnu = float(resources['elapsed_s'])
            need(math.isfinite(real) and math.isfinite(gnu), 'Nonfinite clock readings')
            # 50 ms is a diagnostic discrepancy threshold, not an exclusion rule.
            if abs(real - elapsed) > 0.05 or abs(gnu - elapsed) > 0.05:
                warnings.append(dict(round=i, mode=mode, monotonic_s=elapsed,
                                     realtime_s=real, gnu_elapsed_s=gnu))
            row[mode] = dict(elapsed_monotonic_s=elapsed,
                max_rss_kib=numeric(resources['max_rss_kib'], 'peak RSS', True),
                user_cpu_s=numeric(resources['user_cpu_s'], 'user CPU'),
                system_cpu_s=numeric(resources['system_cpu_s'], 'system CPU'),
                icc_runtime_ms=runtime, ignored_intermediate_survivors=ignored,
                pair_counts={f'H{d}': sum(diagram[d].values()) for d in (1, 2)})
        rows.append(row)
    summary = {mode: {key: stats([r[mode][key] for r in rows[1:]]) for key in
        ('elapsed_monotonic_s', 'max_rss_kib', 'user_cpu_s', 'system_cpu_s', 'icc_runtime_ms')}
        for mode in ('ph', 'piecewise')}
    passed = not mismatches
    ratios = None
    if passed:
        ratios = dict(ratio_of_time_medians=
            summary['ph']['elapsed_monotonic_s']['median'] / summary['piecewise']['elapsed_monotonic_s']['median'],
            paired_time_ratio=stats([r['ph']['elapsed_monotonic_s'] / r['piecewise']['elapsed_monotonic_s']
                                    for r in rows[1:]]),
            peak_rss_median_reduction_percent=100 * (1 -
                summary['piecewise']['max_rss_kib']['median'] / summary['ph']['max_rss_kib']['median']))
    return dict(status='PASS' if passed else 'FAIL', measured_rounds=len(rows)-1,
        reference='round-00/ph; exact Decimal endpoint multisets', mismatches=mismatches,
        notes=['H0 not checked; Ripser not checked',
               'Agreement is not independent proof of correctness',
               'Intermediate death=-1 excluded; final death=-1 means surviving at cutoff',
               'Warmup round-00 excluded; no other rounds removed',
               'Elapsed includes subprocess/GNU Time launch and completion overhead',
               'RSS is kernel-reported iCC peak; GNU/realtime elapsed not used in summary',
               'Checksums are recorded evidence; analyze does not reread the original dataset'],
        clock_warnings=warnings, summary=summary, ratios=ratios, rows=rows,
        input_sha256=evidence, checksums=checks, head=meta['head'])


def run(rounds):
    need(1 <= rounds <= 99, 'Use 1..99 measured rounds')
    need(Path('/usr/bin/time').is_file(), 'Missing /usr/bin/time')
    need(sha(DATA) == DATA_SHA, 'Wrong dataset')
    clock = time.get_clock_info('perf_counter')
    need(clock.monotonic and not clock.adjustable, 'Unsupported clock')
    RUNS.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix='whole-run-monotonic-', dir=RUNS))
    print(f'Results: {root}', file=sys.stderr, flush=True)
    def capture(args):
        return subprocess.check_output(args, text=True)
    def checks():
        return {name: sha(name) for name in ('./build/iCC', DATA)}
    meta = dict(date=datetime.datetime.now().astimezone().isoformat(), python=sys.version,
        clock=vars(clock), measured_rounds=rounds, script_sha256=sha(__file__),
        head=capture(['git', 'rev-parse', 'HEAD']).strip(),
        branch=capture(['git', 'branch', '--show-current']).strip(),
        status=capture(['git', 'status', '--short']),
        diff=capture(['git', 'diff', 'HEAD', '--', 'deps']), checksums=checks(),
        uname=capture(['uname', '-a']), cpu=capture(['lscpu']), memory=capture(['free', '-h']),
        gnu_time=capture(['/usr/bin/time', '--version']))
    write_json(root / 'context.json', meta)
    env = dict(os.environ, LC_ALL='C')
    for i in range(rounds + 1):
        for mode in (('piecewise', 'ph') if i % 2 == 0 else ('ph', 'piecewise')):
            folder = root / f'round-{i:02d}' / mode
            folder.mkdir(parents=True)
            cmd = command(mode)
            (folder / 'command.txt').write_text(shlex.join(cmd) + '\n')
            measured = ['/usr/bin/time', '-f',
                'elapsed_s=%e\nuser_cpu_s=%U\nsystem_cpu_s=%S\nmax_rss_kib=%M\nexit_status=%x',
                '-o', str(folder / 'resources.txt'), *cmd]
            with (folder / 'program.log').open('w') as log:
                real_start = time.time_ns()
                start = time.perf_counter_ns()
                result = subprocess.run(measured, stdout=log, stderr=subprocess.STDOUT, env=env)
                end = time.perf_counter_ns()
                real_end = time.time_ns()
                log.write(f'\niCC exit status: {result.returncode}\n')
            write_json(folder / 'measurement.json', dict(
                elapsed_monotonic_s=(end-start)/1e9, elapsed_realtime_s=(real_end-real_start)/1e9,
                exit_status=result.returncode, scope='subprocess launch through GNU time completion'))
            need(result.returncode == 0, f'iCC failed: {folder}')
            print(f'round-{i:02d} {mode}: exit 0', file=sys.stderr, flush=True)
    write_json(root / 'checksums-after.json', checks())
    report = analyze(root)
    write_json(root / 'summary.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest='action', required=True)
    a = subs.add_parser('analyze', help='Read existing monotonic batch; does not run iCC')
    a.add_argument('batch', type=Path)
    r = subs.add_parser('run', help='Run fixed configuration and save a fresh batch')
    r.add_argument('--rounds', type=int, default=10)
    args = parser.parse_args()
    try:
        report = analyze(args.batch) if args.action == 'analyze' else run(args.rounds)
        print(json.dumps(report, indent=2))
        return 0 if report['status'] == 'PASS' else 1
    except (ValueError, ArithmeticError, OSError, KeyError, TypeError, IndexError,
            subprocess.SubprocessError) as error:
        print(json.dumps(dict(status='ERROR', message=str(error)), indent=2))
        return 2


if __name__ == '__main__':
    sys.exit(main())
