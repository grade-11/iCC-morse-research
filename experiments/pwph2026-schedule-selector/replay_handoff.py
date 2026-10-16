#!/usr/bin/env python3
"""Rebuild bundled Morse source and rerun locked final methods in a fresh output.

No training or retuning. --budget-from keeps the current task's original budget
when validating this handoff; an independent later reproduction gets its own clock.
"""
import argparse
from collections import Counter
from datetime import datetime,timezone,timedelta
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time

def main():
    p=argparse.ArgumentParser();p.add_argument('--bundle',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--dataset')
    p.add_argument('--budget-from',type=Path)
    a=p.parse_args();b=a.bundle.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    helpers=b/'helper-source';sys.path.insert(0,str(helpers/'experiments/pwph2026-stage2-pilot/v2'));sys.path.insert(0,str(b))
    spec=importlib.util.spec_from_file_location('bounded_experiment',b/'offline.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    m.BASE=b;m.REPO=helpers
    plan=json.loads((b/'PLAN.json').read_text())
    if a.dataset:
        final=[x for x in plan['datasets'] if x['split']=='final' and x['id']==a.dataset]
        if len(final)!=1:raise ValueError('Expected one exact reserved final dataset id')
        plan['datasets']=final
    m.write(out/'PLAN.json',plan)
    for name in ['model.json','LOCK.json']:shutil.copy2(b/name,out/name)
    if a.budget_from:budget=json.loads(a.budget_from.read_text())
    else:
        t=datetime.now(timezone.utc);mono=time.monotonic()
        budget=dict(started_at=t.isoformat(),deadline=(t+timedelta(hours=2)).isoformat(),experiments_stop=(t+timedelta(minutes=100)).isoformat(),
                    started_monotonic=mono,deadline_monotonic=mono+7200,experiments_stop_monotonic=mono+6000,boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                    note='Independent reproduction budget; original evidence and original budget remain unchanged')
    m.write(out/'budget.json',budget)
    available=next(int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:'))
    m.write(out/'initial-context.json',dict(head='bundled frozen Morse source 259d4bf; selector hashes in LOCK.json',memory_limit_bytes=min(int(.65*available),5*1024**3)))
    (out/'PROGRESS.md').write_text('# Independent reproduction of locked methods\n')
    shutil.copytree(b/'engine-source',out/'engine-source');shutil.copytree(b/'third-party',out/'third-party')
    (out/'data').mkdir()
    for x in plan['datasets']:
        if x['split']=='final':
            for suffix in ['.csv','.json']:shutil.copy2(b/'data'/(x['id']+suffix),out/'data'/(x['id']+suffix))
    m.build(out);m.final(out)
    rows=json.loads((out/'final-timings.json').read_text())
    print(json.dumps(dict(output=str(out),reproduction_only=True,rows=len(rows),statuses=dict(Counter(x['status'] for x in rows))),indent=2))

if __name__=='__main__':main()
