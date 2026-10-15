#!/usr/bin/env python3
"""Audit saved measurements/identity without rerunning any PH experiment."""
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
from pilot import sha,write,need

def audit(run):
    lock=json.loads((run/'LOCK.json').read_text())
    need(sha(run/'build/iCC')==lock['binary_sha256'],'Frozen binary changed')
    need(sha(Path(__file__).with_name('pilot.py'))==lock['pilot_sha256'],'Frozen harness changed')
    states=Counter();cases=[]
    for p in sorted((run/'cases').glob('*/context.json')):
        c=json.loads(p.read_text());need('exit_status' in c,'Interrupted case: '+str(p))
        need((p.parent/'exit-status.txt').read_text().strip()==str(c['exit_status']),'Exit status file mismatch')
        need(c['inputs']==c['inputs_after'],'Input changed: '+str(p))
        need(all(sha(name)==value for name,value in c['inputs'].items()),'Evidence input changed after case: '+str(p))
        states[c['status']]+=1;cases.append((c['case'],c['key']))
    need(len(set(cases))==len(cases),'Duplicate case identity')
    rows=list(csv.DictReader((run/'summary.csv').open()))
    need(len(rows)==52,'Expected 52 eligible n/seed/method rows after six known n100 mismatches')
    need(all(r['status']=='PASS' and r['rounds']=='5' for r in rows),'Incomplete measurement group')
    ms=list(csv.DictReader((run/'measurements.csv').open()))
    need(len(ms)==52*6 and all(x['status']=='PASS' for x in ms),'Expected warmup plus five rounds per group')
    for n in (100,200):
        for seed in range(20261007,20261012):
            need(any(int(r['n'])==n and int(r['seed'])==seed and r['method']=='adaptive' for r in rows),'Missing adaptive seed')
    names=['experiments/pwph2026-s2-n100/data/sphere_S2_R3_n100_seed20261007.csv',
           'experiments/pwph2026-s2-n100/runs/whole-run-monotonic-ma523yck.tar.gz',
           'experiments/pwph2026-s2-n100/runs/ripser-reference-noa_bg6w.tar.gz']
    # Archive pins are loaded from the reviewed README rather than duplicated by hand.
    import re
    text=Path('experiments/pwph2026-s2-n100/README-n100.md').read_text()
    expected={}
    for name in names:
        m=re.search(r"'([0-9a-f]{64})  "+re.escape(name)+r"'",text)
        need(m is not None,'Missing reviewed evidence pin: '+name)
        expected[name]=m[1];need(sha(name)==m[1],'Reviewed evidence hash changed: '+name)
    start=datetime.fromisoformat(json.loads((run/'initial-context.json').read_text())['started_at'])
    result=dict(status='PASS',processes=len(cases),capture_states=dict(states),measurement_groups=len(rows),
                measured_processes_with_warmups=len(ms),binary_sha256=lock['binary_sha256'],
                elapsed_budget_s=(datetime.now(timezone.utc)-start).total_seconds(),reviewed_evidence=expected)
    write(run/'audit.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();audit(a.run)
