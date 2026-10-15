#!/usr/bin/env python3
"""At most two 20-minute point-deletion investigations; not global minima."""
import argparse
from pathlib import Path
import time
import numpy as np
from pilot import execute, parse_compact, command, write, index, serial, sha

def shrink(run,binary,seed,method):
    target=run/'shrunk'/f'{seed}-{method}';target.mkdir(parents=True,exist_ok=True)
    if (target/'result.json').exists():return
    source=run/'data'/f'S2_R3_n100_seed{seed}.csv'
    original=np.loadtxt(source,delimiter=',');ids=list(range(len(original)))
    start=time.monotonic();attempt=0;last_diff=None
    def test(candidate):
        nonlocal attempt,last_diff
        p=target/f'candidate-{attempt:04d}.csv';attempt+=1
        np.savetxt(p,original[candidate],delimiter=',',fmt='%.17g')
        diagrams=[];states=[]
        for m in ('ph',method):
            out,c=execute(run,f'shrink-{seed}-{method}-{attempt:04d}-{m}',command(binary,p,m),[binary,p,Path(__file__)])
            states.append(c['status'])
            if c['status']!='OK':return False
            diagrams.append(parse_compact((out/'stdout.log').read_text(),'ph' if m=='ph' else 'piecewise')[0])
        mismatch=diagrams[0]!=diagrams[1]
        diff={str(d):dict(missing=serial(diagrams[0][d]-diagrams[1][d]),extra=serial(diagrams[1][d]-diagrams[0][d])) for d in (1,2)}
        write(p.with_suffix('.json'),dict(point_ids=candidate,n=len(candidate),mismatch=mismatch,difference=diff))
        if mismatch:last_diff=diff
        return mismatch
    reproduced=test(ids)
    if reproduced:
        granularity=2
        while len(ids)>4 and time.monotonic()-start<20*60:
            size=(len(ids)+granularity-1)//granularity;changed=False
            for offset in range(0,len(ids),size):
                candidate=ids[:offset]+ids[offset+size:]
                if len(candidate)<4:continue
                if time.monotonic()-start>=20*60:break
                if test(candidate):
                    ids=candidate;granularity=max(2,granularity-1);changed=True;break
            if not changed:
                if granularity>=len(ids):break
                granularity=min(len(ids),granularity*2)
    final=target/'retained.csv';np.savetxt(final,original[ids],delimiter=',',fmt='%.17g')
    result=dict(case=f'shrink-{seed}-{method}',status='PH_MISMATCH' if reproduced else 'NOT_REPRODUCED',
                exit_status=0,reason='Explicit reproduction and bounded point deletion of new schedule discrepancy',
                seed=seed,method=method,original_n=100,retained_n=len(ids),point_ids=ids,
                attempts=attempt,elapsed_s=time.monotonic()-start,data_sha256=sha(final),
                difference=last_diff,claim='Reduced case; no claim of global minimality')
    write(target/'result.json',result);index(run,result)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True)
    a=ap.parse_args()
    for seed,method in [(20261007,'f15'),(20261008,'f4')]:shrink(a.run,a.binary.resolve(),seed,method)
