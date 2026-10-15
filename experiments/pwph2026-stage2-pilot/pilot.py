#!/usr/bin/env python3
"""Sequential, resumable stage-2 pilot. Each subprocess has immutable evidence.
Resume identity includes argv, binary/data/script hashes. No failed case retries.
"""
import argparse
from collections import Counter
from decimal import Decimal
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import shlex
import signal
import statistics
import subprocess
import sys
import threading
import time

BASE = Path(__file__).resolve().parent
REPO = BASE.parent.parent
sys.path.insert(0, str(REPO/'experiments/pwph2026-s2-n100'))
from evidence_run import sha, now
from validate_n100 import msf, projected, serial, need
import numpy as np
import compare_icc_ph as verbose_parser

SCHEDULES = {'f2':[1.2,2.0], 'f15':[1.2,1.5,2.0], 'n100':[1.2,1.7,2.0],
             'f4':[1.2,1.4,1.7,2.0], 'f5':[1.2,1.4,1.6,1.8,2.0]}
DEV = [20261007,20261008,20261009]
HELD = [20261010,20261011]

def write(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
    tmp.replace(path)

def index(run, row):
    p=run/'index.jsonl'
    old=p.read_text() if p.exists() else ''
    tmp=p.with_suffix('.tmp'); tmp.write_text(old+json.dumps(row,ensure_ascii=False)+'\n'); tmp.replace(p)
    with (run/'PROGRESS.md').open('a') as f:
        f.write(f"\n- {now()} `{row.get('case')}`: {row.get('status')} (exit {row.get('exit_status')}).\n")

def execute(run, name, argv, inputs=(), timeout=300):
    context={'argv':list(map(str,argv)), 'inputs':{str(p):sha(p) for p in inputs},
             'harness_sha256':sha(__file__)}
    key=hashlib.sha256(json.dumps(context,sort_keys=True).encode()).hexdigest()
    out=run/'cases'/(name+'-'+key[:12])
    if (out/'context.json').exists():
        c=json.loads((out/'context.json').read_text())
        if 'exit_status' in c: return out,c
        raise RuntimeError(f'Interrupted evidence {out}; inspect and retain before explicit retry')
    out.mkdir(parents=True,exist_ok=False)
    context.update(case=name,key=key,started_at=now(),cwd=str(Path.cwd()),
                   timeout_s=timeout,platform=platform.platform(),python=sys.version)
    write(out/'context.json',context)
    (out/'command.txt').write_text(shlex.join(context['argv'])+'\n')
    limit=json.loads((run/'initial-context.json').read_text())['memory_limit_bytes']
    start=time.perf_counter(); stop=threading.Event(); reason=[]; max_tree=[0]
    env=dict(os.environ,OMP_NUM_THREADS='4',OMP_THREAD_LIMIT='4',OPENBLAS_NUM_THREADS='4',
             MKL_NUM_THREADS='4',LC_ALL='C',MPLCONFIGDIR='/tmp/icc-stage2-mpl')
    context['environment']={k:env[k] for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','LC_ALL')}
    command=['/usr/bin/time','-f','max_rss_kib=%M\nuser_s=%U\nsystem_s=%S\nexit_status=%x',
             '-o',str(out/'resources.txt')]+context['argv']
    rc=125
    try:
        with (out/'stdout.log').open('wb') as stdout, (out/'stderr.log').open('wb') as stderr:
            proc=subprocess.Popen(command,stdout=stdout,stderr=stderr,env=env,start_new_session=True)
            def monitor():
                while not stop.is_set():
                    rss=0
                    for p in Path('/proc').glob('[0-9]*'):
                        try:
                            fields=(p/'stat').read_text().rsplit(')',1)[1].split()
                            if int(fields[2])==proc.pid:
                                rss+=int((p/'statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
                        except (OSError,ValueError,IndexError): pass
                    max_tree[0]=max(max_tree[0],rss)
                    if rss>limit or time.perf_counter()-start>timeout:
                        reason.append('MEMORY_LIMIT' if rss>limit else 'TIMEOUT')
                        try: os.killpg(proc.pid,signal.SIGKILL)
                        except ProcessLookupError: pass
                        return
                    stop.wait(.05)
            watcher=threading.Thread(target=monitor);watcher.start()
            rc=proc.wait(); elapsed=time.perf_counter()-start
            stop.set();watcher.join()
    except Exception as e:
        elapsed=time.perf_counter()-start
        (out/'stderr.log').write_text(str(e));reason.append('ERROR')
    context.update(ended_at=now(),elapsed_monotonic_s=elapsed,exit_status=rc,
                   status=reason[0] if reason else ('OK' if rc==0 else 'ERROR'),
                   tree_rss_sampled_bytes=max_tree[0])
    context['inputs_after']={str(p):sha(p) for p in inputs}
    if context['inputs_after']!=context['inputs']:context['status']='ERROR'
    if (out/'resources.txt').exists():
        context['resources']=dict(x.split('=',1) for x in (out/'resources.txt').read_text().splitlines() if '=' in x)
    write(out/'context.json',context);(out/'exit-status.txt').write_text(str(rc)+'\n')
    index(run,context)
    return out,context

def parse_compact(text, mode, cutoff=Decimal('2')):
    """Generalize benchmark_icc.parse_log: derive all actual windows from output.
    Finite records belong to death window; only FINAL-window survivors are bars.
    No deduplication by birth: Counter retains endpoint multiplicities.
    """
    need(re.search(r'run time = \d+\s*$',text) is not None,'Missing completion')
    diagrams={1:Counter(),2:Counter()}; ignored=0; windows=[];dim=None;lo=Decimal(0);hi=cutoff
    lines=text.splitlines();i=0;sections=[]
    while i<len(lines):
        line=lines[i].strip()
        if line=='window lower bound, window upper bound':
            i+=1;lo,hi=map(Decimal,lines[i].split(','));windows.append((lo,hi))
        m=re.match(r'Processing dim ([23]),',line)
        if m: dim=int(m[1]);sections.append(dim)
        pair=(re.fullmatch(r'\(([^,]+),\s*([^,]+)\)',line) if mode=='ph' else
              re.fullmatch(r'([23]),\s*([^,]+),\s*([^,]+)',line))
        if pair:
            if mode=='ph': b,d=map(Decimal,pair.groups())
            else: dim=int(pair[1]);b,d=map(Decimal,pair.groups()[1:])
            need(dim in (2,3) and b.is_finite() and d.is_finite() and 0<=b<=hi and
                 (d==-1 or b<=d<=hi) and (d==-1 or d>=lo),'Invalid interval')
            if d==-1 and hi!=cutoff: ignored+=1
            else:diagrams[dim-1][(b,d)]+=1
        elif line and not (m or line in {'window lower bound, window upper bound',
                'dimension, birth weight, death weight','dimensional persistent pairs:','(empty)'} or
                line.startswith(('total point number:','run time =','[adaptive]'))):
            raise ValueError('Unexpected PH output: '+line)
        i+=1
    if mode!='ph':
        need(windows and windows[0][0]==0 and windows[-1][1]==cutoff and
             all(a<b for a,b in windows) and all(windows[i-1][1]==windows[i][0] for i in range(1,len(windows))), 'Window chain incomplete')
    else: need(sections==[2,3],'Missing/repeated/out-of-order baseline section')
    return diagrams,ignored,[[str(a),str(b)] for a,b in windows]

def command(binary,data,method,verbose=False):
    cmd=[str(binary),'--tool','ph' if method=='ph' else 'piecewise','--file-name',str(data),
         '--max-dim','3','--threads','4']
    if method=='ph':cmd+=['--max-eps','2']
    elif method=='auto':cmd+=['--eps-interval-count','3','--pv-cap-scale','.85']
    else:
        cmd+=['--eps-breaks',*map(str,([1.2,1.4,1.5,1.6,1.7,1.8,2.] if method=='adaptive' else SCHEDULES[method])), '--pv-cap-scale','.85']
        if method=='adaptive':cmd+=['--adaptive','--adaptive-target','.35']
    if verbose:cmd+=['--verbose']
    return cmd

def generate(run,n,seeds):
    data_dir=run/'data';data_dir.mkdir(exist_ok=True)
    sys.path.insert(0,str(REPO/'Sphere gen'))
    import dSphere_gen
    rows=[]
    for seed in seeds:
        p=data_dir/f'S2_R3_n{n}_seed{seed}.csv'
        if not p.exists():
            if n==100 and seed==20261007:
                p.write_bytes((REPO/'experiments/pwph2026-s2-n100/data/sphere_S2_R3_n100_seed20261007.csv').read_bytes())
            else:
                np.random.seed(seed)
                pts=dSphere_gen.dSphere_uniformRandom(dimension=3,numPoints=n,r1=1.,r2=1.,trials=0,selection='min')
                with p.open('w',newline='') as f:csv.writer(f).writerows(pts)
        rows.append(dict(n=n,seed=seed,path=str(p),sha256=sha(p),
                         generator='dSphere_gen.dSphere_uniformRandom',generator_sha256=sha(REPO/'Sphere gen/dSphere_gen.py'),
                         utilities_sha256=sha(REPO/'Sphere gen/utilities.py'),trials=0,selection='min',
                         note='Existing reviewed CSV' if n==100 and seed==20261007 else 'New equivalent sample, not author data'))
    write(run/f'datasets-n{n}-{seeds[0]}.json',rows)

def reference(run,data):
    code="""import sys,numpy as np;from ripser import ripser
p=np.loadtxt(sys.argv[1],delimiter=',')
r=ripser(p,maxdim=2,thresh=2.,coeff=2,metric='euclidean',distance_matrix=False,do_cocycles=False,n_perm=None)
np.savez_compressed(sys.argv[2],**{'H'+str(i):a for i,a in enumerate(r['dgms'])})
print({i:len(a) for i,a in enumerate(r['dgms'])})
"""
    # Stable destination is an output of the captured process, never used for identity.
    dest=run/'references'/(data.stem+'.npz');dest.parent.mkdir(exist_ok=True)
    out,c=execute(run,'ripser-'+data.stem,[sys.executable,'-c',code,str(data),str(dest)],[data])
    need(c['status']=='OK','Ripser failed')
    with np.load(dest,allow_pickle=False) as z:return {d:Counter(map(tuple,z[f'H{d}'].tolist())) for d in range(3)}

def validate(run,binary,data,method,base,refs):
    name=f'validate-{data.stem}-{method}'
    reportfile=run/(name+'-'+sha(binary)[:12]+'.json')
    if reportfile.exists():return json.loads(reportfile.read_text())
    mode='ph' if method=='ph' else 'piecewise'
    out,c=execute(run,'baseline-'+data.stem if method=='ph' else name,command(binary,data,method),[binary,data])
    result=dict(case=name,method=method,data=str(data),binary_sha256=sha(binary),status=c['status'],exit_status=c['exit_status'])
    if c['status']=='OK':
        diag,ignored,windows=parse_compact((out/'stdout.log').read_text(),mode)
        result.update(H1_H2_exact=all(diag[d]==base[d] for d in (1,2)),
                      H1_H2_ripser_float32=all(projected(diag[d])==refs[d] for d in (1,2)),
                      diagrams={str(d):serial(diag[d]) for d in (1,2)},ignored_checkpoint_survivors=ignored,windows=windows)
        vo,vc=execute(run,'baseline-msf-'+data.stem if method=='ph' else name+'-diagnostic',command(binary,data,method,True),[binary,data])
        result['diagnostic_exit_status']=vc['exit_status']
        if vc['status']!='OK':result['status']=vc['status']
        else:
            vtext=(vo/'stdout.log').read_text()
            actual_breaks=[Decimal(w[1]) for w in windows] if mode=='piecewise' else [Decimal('2')]
            # Derived parser input declares ACTUAL schedule; original CLI grid stays in raw log/context.
            if mode=='piecewise':
                vtext=re.sub(r'^  epsilon breaks: .+$','  epsilon breaks: '+' '.join(map(str,actual_breaks)),vtext,flags=re.M)
            derived=vo/'parser-input.log';derived.write_text(vtext+'\niCC exit status: 0\n')
            vd,vi=verbose_parser.parse_log(derived,dict(dataset=str(data),dim=3,threads=4,
                breaks=actual_breaks,alpha=Decimal('.85')),mode)
            need(vd==diag and vi==ignored,'Verbose/nonverbose PH or actual schedule differs')
            frames=msf(vtext);need(frames,'No actual MSF diagnostics')
            first=frames[0];h0=first['pairs']
            h0ok=projected(h0)==refs[0] and (method=='ph' or (first['metadata']['components']=='1' and h0==base[0]))
            result.update(H0='PASS' if h0ok else 'H0_NOT_VERIFIED',H0_intervals=serial(h0))
            result['difference']={str(d):dict(missing=serial(base[d]-diag[d]),extra=serial(diag[d]-base[d])) for d in (1,2)}
            result['status']=('PH_MISMATCH' if not result['H1_H2_exact'] or not result['H1_H2_ripser_float32'] else ('PASS' if h0ok else 'H0_NOT_VERIFIED'))
    write(reportfile,result);index(run,result);return result

def batch(run,binary,n,seeds,methods,rounds):
    generate(run,n,seeds)
    validations={}
    for seed in seeds:
        data=run/'data'/f'S2_R3_n{n}_seed{seed}.csv';refs=reference(run,data)
        bo,bc=execute(run,'baseline-'+data.stem,command(binary,data,'ph'),[binary,data]);need(bc['status']=='OK','Baseline failed')
        base,_,_=parse_compact((bo/'stdout.log').read_text(),'ph')
        vo,vc=execute(run,'baseline-msf-'+data.stem,command(binary,data,'ph',True),[binary,data]);need(vc['status']=='OK','Baseline diagnostics failed')
        base[0]=msf((vo/'stdout.log').read_text())[0]['pairs']
        validations[seed]={m:validate(run,binary,data,m,base,refs) for m in methods}
        eligible=[m for m in methods if validations[seed][m]['status']=='PASS']
        for r in range(rounds+1):
            order=eligible[r%len(eligible):]+eligible[:r%len(eligible)] if eligible else []
            if r%2:order.reverse()
            for m in order:
                out,c=execute(run,f'measure-n{n}-{seed}-{m}-round{r}',command(binary,data,m),[binary,data])
                p=out/'ph-check.json'
                if not p.exists():
                    result=dict(case=c['case'],status=c['status'],exit_status=c['exit_status'],round=r,method=m,n=n,seed=seed)
                    if c['status']=='OK':
                        diag,ignored,windows=parse_compact((out/'stdout.log').read_text(),'ph' if m=='ph' else 'piecewise')
                        result.update(status='PASS' if all(diag[d]==base[d] and projected(diag[d])==refs[d] for d in (1,2)) else 'PH_MISMATCH',windows=windows)
                    write(p,result);index(run,result)
                if json.loads(p.read_text())['status']!='PASS':break
    write(run/(f'validation-summary-n{n}-{seeds[0]}-'+sha(binary)[:12]+'.json'),validations)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True)
    ap.add_argument('--n',type=int,choices=[100,200],default=100);ap.add_argument('--split',choices=['dev','held'],default='dev')
    ap.add_argument('--methods',nargs='+',default=list(SCHEDULES));ap.add_argument('--rounds',type=int,default=5)
    a=ap.parse_args();batch(a.run,a.binary.resolve(),a.n,DEV if a.split=='dev' else HELD,a.methods,a.rounds)
if __name__=='__main__':main()
