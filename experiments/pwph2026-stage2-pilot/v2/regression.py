#!/usr/bin/env python3
"""Retain known failures, exercise all-DEFER, and verify real G actions."""
import argparse
from lab import *

def run(run,binary):
    prior=BASE.parent/'runs/stage2-20261015T190630Z';results=[]
    schedules={'f15':[1.2,1.5,2.],'f4':[1.2,1.4,1.7,2.],'f5':[1.2,1.4,1.6,1.8,2.]}
    for seed in (20261007,20261008):
        data=prior/'data'/f'S2_R3_n100_seed{seed}.csv'
        for method,grid in schedules.items():
            r=check(run,f'known-full-{seed}-{method}',binary,data,grid)
            old=json.loads((prior/f'validate-S2_R3_n100_seed{seed}-{method}-fc613d5ca2a8.json').read_text())
            need(r['status']=='PH_MISMATCH' and r['diagrams']==old['diagrams'],'Known failure changed; investigate')
            results.append(r)
    models=run/'unit-models';models.mkdir(exist_ok=True)
    for mode,qcost in [('positive',.001),('negative',.1)]:
        model=models/(mode+'.txt');model.write_text('PWPH_G_V1 3 200 1000000\n.01 0 0 0 0 0 0 0\n.02 0 0 0 0 0 0 -.015\n5000 0 0 0 0 0 0 0\n'+str(qcost)+' 0 0\n')
        r=check(run,'unit-G-action-'+mode,binary,run/'data/20261007-f15.csv',[1.2,1.5,2.],extra=['--g-model',str(model)])
        states=parse_states((Path(r['path'])/'stdout.log').read_text())
        need(any(s['action']=='QUOTIENT' for s in states) if mode=='positive' else all(s['action']!='QUOTIENT' for s in states),'G must change action')
        need(r['status']==('PH_MISMATCH' if mode=='positive' else 'PASS'),'G unit fixture PH changed')
        results.append(dict(**r,unit_fixture=True))
    write(run/'regression-v2.json',results)
    print(json.dumps(dict(status='PASS_EXPECTED_RESULTS',known_failures=6,G_actual_action_fixtures=2)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);a=ap.parse_args();guard=runner_lock(a.run);run(a.run,a.binary.resolve())
