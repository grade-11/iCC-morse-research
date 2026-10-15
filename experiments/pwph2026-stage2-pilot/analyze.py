#!/usr/bin/env python3
"""Read evidence only. Require five successful measured rounds for each row."""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import statistics as st
from pilot import sha, write, SCHEDULES, DEV, now

def csvwrite(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def analyze(run,binary,lock=False):
    binary_hash=sha(binary);rows=[];groups={}
    for p in sorted((run/'cases').glob('*/context.json')):
        c=json.loads(p.read_text());m=re.fullmatch(r'measure-n(100|200)-(\d+)-(\w+)-round(\d+)',c['case'])
        if not m or binary_hash not in c['inputs'].values():continue
        check=p.parent/'ph-check.json'
        if not check.exists():continue
        ph=json.loads(check.read_text());n,seed,method,r=int(m[1]),int(m[2]),m[3],int(m[4])
        row=dict(n=n,seed=seed,split='dev' if seed in DEV else 'held',method=method,round=r,
                 status=ph['status'],exit_status=c['exit_status'],elapsed_s=c['elapsed_monotonic_s'],
                 peak_rss_kib=int(c.get('resources',{}).get('max_rss_kib',-1)),path=str(p.parent))
        rows.append(row)
        if r>0:groups.setdefault((n,seed,method),[]).append(row)
    summary=[]
    for (n,seed,method),g in sorted(groups.items()):
        valid=len(g)==5 and {x['round'] for x in g}==set(range(1,6)) and all(x['status']=='PASS' and x['peak_rss_kib']>0 for x in g)
        t=[x['elapsed_s'] for x in g];mem=[x['peak_rss_kib'] for x in g]
        summary.append(dict(n=n,seed=seed,split='dev' if seed in DEV else 'held',method=method,
                            status='PASS' if valid else 'INCOMPLETE',rounds=len(g),
                            median_s=st.median(t),min_s=min(t),max_s=max(t),
                            median_rss_kib=st.median(mem),min_rss_kib=min(mem),max_rss_kib=max(mem)))
    csvwrite(run/'measurements.csv',rows);csvwrite(run/'summary.csv',summary)
    pairs=[]
    for (n,seed,method),g in sorted(groups.items()):
        if method!='adaptive' or len(g)!=5:continue
        byround={x['round']:x for x in g}
        for fixed in SCHEDULES:
            f=groups.get((n,seed,fixed),[])
            if len(f)!=5 or any(x['status']!='PASS' for x in f+g):continue
            t=[x['elapsed_s']/byround[x['round']]['elapsed_s'] for x in f]
            mem=[x['peak_rss_kib']/byround[x['round']]['peak_rss_kib'] for x in f]
            pairs.append(dict(n=n,seed=seed,fixed=fixed,split='dev' if seed in DEV else 'held',
                              ratio='fixed/adaptive; >1 favors adaptive',median_time_ratio=st.median(t),min_time_ratio=min(t),max_time_ratio=max(t),
                              median_rss_ratio=st.median(mem),min_rss_ratio=min(mem),max_rss_ratio=max(mem)))
    csvwrite(run/'paired-ratios.csv',pairs)
    if lock:
        need=run/'LOCK.json'
        if need.exists():raise RuntimeError('Lock already exists; never retune after held-out')
        scores={}
        for method in SCHEDULES:
            found=[x for x in summary if x['n']==100 and x['seed'] in DEV and x['method']==method and x['status']=='PASS']
            if len(found)==3:scores[method]=math.exp(st.mean(math.log(x['median_s']) for x in found))
        if not scores:raise RuntimeError('No fixed schedule eligible across development seeds')
        best=min(scores,key=scores.get)
        write(need,dict(locked_at=now(),best_fixed=best,schedule=SCHEDULES[best],
                       criterion='Smallest geometric mean of per-seed median whole-run time, PH PASS all three development seeds',
                       scores_s=scores,adaptive_target=.35,controller_sha256=sha(Path('deps/AdaptiveSchedule.hpp')),
                       pilot_sha256=sha(Path(__file__).with_name('pilot.py')),binary_sha256=binary_hash,
                       note='Heuristic unchanged after initial rule. No held-out read before this lock.'))
    print(json.dumps(dict(binary_sha256=binary_hash,summary_rows=len(summary),pairs=len(pairs),lock=json.loads((run/'LOCK.json').read_text()) if (run/'LOCK.json').exists() else None),indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--lock',action='store_true')
    a=ap.parse_args();analyze(a.run,a.binary,a.lock)
