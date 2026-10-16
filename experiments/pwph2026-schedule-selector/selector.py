#!/usr/bin/env python3
"""Choose a complete schedule before PH, then replace this process with Morse/iCC.

Standard library only. No Ripser, observed future PVs, or PH outputs are inputs.
The enclosing timer includes interpreter startup, features, logging, exec and PH.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import time

SCHEDULES = {
    'noQ': None,
    'q08': [8, 20], 'q10': [10, 20], 'q12': [12, 20],
    'q14': [14, 20], 'q16': [16, 20], 'author': [12, 17, 20],
    'q11_16': [11, 16, 20], 'q13_18': [13, 18, 20],
    'q08_12_16': [8, 12, 16, 20],
}
FEATURE_NAMES = ['intercept', 'log_n_over_100', 'edge_density_08',
                 'quotient_count', 'has_quotient', 'first_quotient_over_2',
                 'sum_edge_densities_at_quotients', 'hasQ_log_n',
                 'hasQ_density_08', 'boundary_distance_proxy']

def raw_features(path):
    raw = Path(path).read_bytes()
    points = [tuple(map(float, row)) for row in csv.reader(io.StringIO(raw.decode())) if row]
    if len(points) < 2 or any(len(p) != 3 or not all(map(math.isfinite, p)) for p in points):
        raise ValueError('Expected at least two finite points in R3')
    counts = {k: 0 for k in range(8, 19)}
    pairs = len(points) * (len(points)-1) // 2
    mean = 0.0
    minimum = math.inf
    for i, a in enumerate(points):
        for b in points[:i]:
            squared = sum((x-y)**2 for x,y in zip(a,b))
            d = math.sqrt(squared)
            mean += d
            minimum = min(minimum,d)
            for k in counts:
                if squared < (k/10)**2:
                    counts[k] += 1
    return dict(n=len(points), density={str(k): v/pairs for k,v in counts.items()},
                mean_distance=mean/pairs, minimum_distance=minimum,
                max_unit_norm_error=max(abs(sum(x*x for x in p)-1) for p in points),
                input_sha256=hashlib.sha256(raw).hexdigest())

def vector(x, schedule):
    q = (SCHEDULES[schedule] or [20])[:-1]
    d = x['density']; logn = math.log(x['n']/100)
    return [1., logn, d['8'], len(q), int(bool(q)), q[0]/20 if q else 0.,
            sum(d[str(k)] for k in q), int(bool(q))*logn,
            int(bool(q))*d['8'], sum(d[str(k)]*(1-d[str(k)])*(2-k/10) for k in q)]

def choose(x, model):
    domain = model['domain']
    inside = (domain['n'][0] <= x['n'] <= domain['n'][1] and
              domain['density08'][0] <= x['density']['8'] <= domain['density08'][1] and
              x['minimum_distance'] > 0 and x['max_unit_norm_error'] < 1e-12)
    predictions = {s: math.exp(max(-15., min(10., sum(a*b for a,b in zip(vector(x,s),model['coefficients'])))))
                   for s in SCHEDULES}
    # Eligibility is an empirical development/selection filter, not a PH certificate.
    chosen = min(model['eligible_schedules'], key=lambda s: (predictions[s], list(SCHEDULES).index(s))) if inside else model['fallback']
    return dict(schedule=chosen, ticks=SCHEDULES[chosen], in_domain=inside,
                reason='minimum_predicted_seconds' if inside else 'predeclared_fallback',
                predictions_s=predictions,
                G_seconds={s:predictions['noQ']-v for s,v in predictions.items()})

def command(binary, data, schedule, verbose=False, threads=4):
    ticks = SCHEDULES[schedule] if isinstance(schedule,str) else schedule
    cmd = [str(binary), '--tool', 'ph' if ticks is None else 'piecewise',
           '--file-name',str(data),'--max-dim','2','--threads',str(threads)]
    if ticks is None: cmd += ['--max-eps','2']
    else: cmd += ['--eps-breaks',*[format(k/10,'.17g') for k in ticks], '--pv-cap-scale','.85']
    if verbose: cmd += ['--verbose']
    return cmd

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--binary',required=True,type=Path)
    p.add_argument('--data',required=True,type=Path)
    p.add_argument('--model',required=True,type=Path)
    p.add_argument('--choice',required=True,type=Path)
    a = p.parse_args(); start = time.perf_counter()
    model = json.loads(a.model.read_text()); x = raw_features(a.data)
    result = choose(x,model)
    result.update(features=x, pre_exec_s=time.perf_counter()-start,
                  model_sha256=hashlib.sha256(a.model.read_bytes()).hexdigest(),
                  selector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    a.choice.parent.mkdir(parents=True,exist_ok=True)
    with a.choice.open('x') as f: json.dump(result,f,indent=2);f.write('\n')
    cmd = command(a.binary,a.data,result['schedule'])
    os.execv(str(a.binary),cmd)

if __name__ == '__main__': main()
