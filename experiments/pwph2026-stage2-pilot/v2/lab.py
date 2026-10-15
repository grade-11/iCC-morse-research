#!/usr/bin/env python3
"""Evidence/PH helpers for the four-hour follow-up; preserves the v1 workflow."""
import sys,os
from pathlib import Path
BASE=Path(__file__).resolve().parent
REPO=BASE.parents[2]
sys.path.insert(0,str(BASE.parent))
import pilot as old
from evidence_run import sha,now
from validate_n100 import msf,projected,serial,need
import compare_icc_ph as verbose_parser
from collections import Counter
from decimal import Decimal
import json,time,hashlib,re,fcntl,importlib.metadata
import numpy as np

def write(p,obj):old.write(p,obj)
def run_path():return Path((BASE/'CURRENT-RUN.txt').read_text().strip())
def budget(run):
    b=json.loads((run/'budget.json').read_text())
    need(Path('/proc/sys/kernel/random/boot_id').read_text().strip()==b['boot_id'],'Boot changed: monotonic budget requires review')
    return b
def remaining(run,pack=False):return budget(run)['deadline_monotonic' if pack else 'experiments_stop_monotonic']-time.monotonic()
def runner_lock(run):
    fd=(run/'runner.lock').open('w');fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB);fd.write(str(os.getpid()));fd.flush();return fd
def execute(run,name,argv,inputs=(),timeout=300):
    left=remaining(run)
    need(left>0,'Experiment budget ended; reserve final 20 minutes for packaging')
    model_inputs=[Path(argv[argv.index('--g-model')+1])] if '--g-model' in argv else []
    # Explicit compatibility whitelist for a diagnostic-only harness fix.
    # Child argv and every non-harness input must remain byte-identical.
    rules=run/'COMPATIBLE-HARNESS.json'
    expected={str(p):sha(p) for p in [*inputs,*model_inputs]}
    if rules.exists():
        allowed=json.loads(rules.read_text())['allow_source_hashes']
        for path in (run/'cases').glob(name+'-*/context.json'):
            c=json.loads(path.read_text());previous=c.get('inputs',{});labkey=str(Path(__file__))
            if (c.get('status')=='OK' and c.get('argv')==list(map(str,argv)) and previous.get(labkey) in allowed and
                {k:v for k,v in previous.items() if k!=labkey}==expected):
                return path.parent,c
    return old.execute(run,name,argv,[*inputs,*model_inputs,Path(__file__)],min(timeout,left))
def note(run,name,**fields):
    old.index(run,dict(case=name,at=now(),**fields))

def compact(text,mode,maxdim=3,cutoff=2):
    # Adapt the reviewed parser to dimensions 2/3 and an arbitrary cutoff.
    hi=Decimal(str(cutoff));lo=Decimal(0);final=hi;d=None;windows=[];ignored=0
    bars={i:Counter() for i in range(1,maxdim)};sections=[];lines=text.splitlines();i=0
    need(re.search(r'run time = \d+\s*$',text) is not None,'Missing completion')
    while i<len(lines):
        s=lines[i].strip()
        if s=='window lower bound, window upper bound':
            i+=1;lo,hi=map(Decimal,lines[i].split(','));windows.append((lo,hi))
        m=re.match(r'Processing dim (\d+),',s)
        if m:d=int(m[1]);sections.append(d)
        p=re.fullmatch(r'\(([^,]+),\s*([^,]+)\)',s) if mode=='ph' else re.fullmatch(r'([23]),\s*([^,]+),\s*([^,]+)',s)
        if p:
            if mode=='ph':b,death=map(Decimal,p.groups())
            else:d=int(p[1]);b,death=map(Decimal,p.groups()[1:])
            need(d-1 in bars and b.is_finite() and death.is_finite() and 0<=b<=hi and (death==-1 or b<=death<=hi) and (death==-1 or death>=lo),'Invalid endpoints')
            if death==-1 and hi!=final:ignored+=1
            else:bars[d-1][b,death]+=1
        elif s and not (m or s in {'window lower bound, window upper bound','dimension, birth weight, death weight','dimensional persistent pairs:','(empty)'} or s.startswith(('total point number:','run time =','[adaptive]','[g]','[state]','[simplex]','[pv-members]','[support-trace]'))):raise ValueError('Unexpected output: '+s)
        i+=1
    if mode=='ph':need(sections==list(range(2,maxdim+1)),'Incomplete baseline dimensions')
    else:need(windows and windows[0][0]==0 and windows[-1][1]==final and all(a<b for a,b in windows) and all(windows[i-1][1]==windows[i][0] for i in range(1,len(windows))),'Incomplete window chain')
    return bars,ignored,windows

def command(binary,data,breaks=None,threads=4,maxdim=3,cutoff=2,verbose=False,extra=()):
    cmd=[str(binary),'--tool','ph' if breaks is None else 'piecewise','--file-name',str(data),'--max-dim',str(maxdim),'--threads',str(threads)]
    cmd+=['--max-eps',str(cutoff)] if breaks is None else ['--eps-breaks',*map(str,breaks),'--pv-cap-scale','.85']
    if verbose:cmd+=['--verbose']
    return cmd+list(extra)

def reference(run,data,maxdim=3,cutoff=2):
    config=dict(data_sha256=sha(data),maxdim=maxdim-1,thresh=cutoff,coeff=2,metric='euclidean',distance_matrix=False,do_cocycles=False,n_perm=None,
                ripser=importlib.metadata.version('ripser'),numpy=np.__version__)
    key=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest();dest=run/'references'/(key+'.npz');dest.parent.mkdir(exist_ok=True)
    meta=dest.with_suffix('.json')
    if meta.exists():need(json.loads(meta.read_text())['sha256']==sha(dest),'Reference cache changed')
    else:
        code="import numpy as np,sys;from ripser import ripser;d=ripser(np.loadtxt(sys.argv[1],delimiter=','),maxdim=int(sys.argv[3]),thresh=float(sys.argv[4]),coeff=2,metric='euclidean',distance_matrix=False,do_cocycles=False,n_perm=None)['dgms'];np.savez_compressed(sys.argv[2],**{'H'+str(i):a for i,a in enumerate(d)})"
        o,c=execute(run,'reference-'+key[:16],[sys.executable,'-c',code,str(data),str(dest),str(maxdim-1),str(cutoff)],[data])
        need(c['status']=='OK','Ripser failed');write(meta,dict(config=config,sha256=sha(dest),case=str(o)))
    with np.load(dest) as z:return {d:Counter(map(tuple,z['H'+str(d)].tolist())) for d in range(maxdim)}

def check(run,name,binary,data,breaks=None,maxdim=3,cutoff=2,threads=4,extra=(),diagnostic=True):
    refs=reference(run,data,maxdim,cutoff)
    bo,bc=execute(run,'base-'+sha(data)[:12]+f'-d{maxdim}-c{cutoff}-t{threads}',command(binary,data,maxdim=maxdim,cutoff=cutoff,threads=threads),[binary,data])
    need(bc['status']=='OK','Baseline failed')
    base,_,_=compact((bo/'stdout.log').read_text(),'ph',maxdim,cutoff)
    validation_extra=list(extra)
    if breaks is not None and '--cost-telemetry' not in validation_extra:validation_extra+=['--cost-telemetry']
    o,c=execute(run,name,command(binary,data,breaks,threads,maxdim,cutoff,extra=validation_extra),[binary,data])
    result=dict(case=name,binary_sha256=sha(binary),dataset_sha256=sha(data),dataset=str(data),breaks=breaks,maxdim=maxdim,cutoff=cutoff,threads=threads,
                status=c['status'],exit_status=c['exit_status'],path=str(o))
    if c['status']=='OK':
        bars,ignored,windows=compact((o/'stdout.log').read_text(),'ph' if breaks is None else 'piecewise',maxdim,cutoff)
        result.update(exact={str(d):bars[d]==base[d] for d in bars},ripser={str(d):projected(bars[d])==refs[d] for d in bars},
                      baseline_ripser={str(d):projected(base[d])==refs[d] for d in base},diagrams={str(d):serial(bars[d]) for d in bars},
                      difference={str(d):dict(missing=serial(base[d]-bars[d]),extra=serial(bars[d]-base[d])) for d in bars},ignored_survivors=ignored,windows=[[str(a),str(b)] for a,b in windows])
        result['status']='PASS' if all(result['exact'].values()) and all(result['ripser'].values()) and all(result['baseline_ripser'].values()) else 'PH_MISMATCH'
        result['H0']='H0_NOT_VERIFIED'
        if diagnostic:
            states=parse_states((o/'stdout.log').read_text())
            actions=[int(s['action']=='QUOTIENT') for s in states if s['action'] not in ('FINAL','STOP')]
            replay_breaks=[float(b) for a,b in result['windows']] if breaks is not None else None
            replay_extra=(['--cost-telemetry']+(['--quotient-mask',*map(str,actions)] if actions else [])) if breaks is not None else []
            vo,vc=execute(run,name+'-diagnostic-replay',command(binary,data,replay_breaks,threads,maxdim,cutoff,True,replay_extra),[binary,data])
            if vc['status']!='OK':result['status']=vc['status']
            else:
                text=(vo/'stdout.log').read_text();frames=msf(text);need(frames,'No actual MSF')
                actual=[Decimal(b) for a,b in result['windows']] if breaks is not None else [Decimal(str(cutoff))]
                if breaks is not None:text=re.sub(r'^  epsilon breaks: .+$','  epsilon breaks: '+' '.join(map(str,actual)),text,flags=re.M)
                pp=vo/'parser-input.log';pp.write_text(text+'\niCC exit status: 0\n')
                vb,vi=verbose_parser.parse_log(pp,dict(dataset=str(data),dim=maxdim,threads=threads,breaks=actual,alpha=Decimal('.85')),'ph' if breaks is None else 'piecewise')
                need(vb==bars and vi==ignored,'Verbose/compact mismatch')
                if breaks is not None:
                    replay_states=parse_states((vo/'stdout.log').read_text())
                    need([s['state_hash'] for s in states]==[s['state_hash'] for s in replay_states],'Diagnostic replay state mismatch')
                # Baseline edge-only run calls the same actual MSF clearing path,
                # avoiding an unnecessary repeated high-dimensional baseline.
                ho,hc=execute(run,'base-msf-'+sha(data)[:12]+f'-c{cutoff}-t{threads}',command(binary,data,maxdim=1,cutoff=cutoff,threads=threads,verbose=True),[binary,data])
                need(hc['status']=='OK','Baseline MSF failed');base0=msf((ho/'stdout.log').read_text())[0]['pairs']
                from h0_ledger import ledger
                h0=ledger(frames,result['windows'] if breaks is not None else [['0',str(cutoff)]],base0,refs[0])
                h0ok=h0['status']=='PASS'
                result.update(H0='PASS' if h0ok else 'H0_NOT_VERIFIED',h0_ledger=h0,
                              h0_verifier_sha256=sha(BASE/'h0_ledger.py'),diagnostic=str(vo),actions=actions)
                if result['status']=='PASS' and not h0ok:result['status']='H0_NOT_VERIFIED'
    write(o/'ph-check.json',result);note(run,name,**{k:v for k,v in result.items() if k!='case'})
    return result

def parse_states(text):
    return [dict(x.split('=') for x in line.split()[1:]) for line in text.splitlines() if line.startswith('[state] ')]
