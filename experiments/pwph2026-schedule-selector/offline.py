#!/usr/bin/env python3
"""Bounded sequential experiment; old evidence is read-only, no online Ripser."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import statistics
import subprocess
import sys
import tarfile
import time
import traceback

BASE = Path(__file__).resolve().parent
REPO = BASE.parents[1]
sys.path.insert(0,str(REPO/'experiments/pwph2026-stage2-pilot/v2'))
import lab
from h0_ledger import ledger
import numpy as np
import selector

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now(): return datetime.now(timezone.utc).isoformat()
def write(p,x): lab.write(p,x)
def read(p): return json.loads(Path(p).read_text())
def table(p,rows):
    rows=list(rows)
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def note(run,name,status='OK',**kw): lab.note(run,name,status=status,exit_status=kw.pop('exit_status',0),**kw)
def execute(run,name,argv,inputs=(),timeout=300):
    return lab.execute(run,name,list(map(str,argv)),[*inputs,BASE/'offline.py',BASE/'selector.py'],timeout)
def binary(run): return run/'build/iCC'
def grid(s): return None if selector.SCHEDULES[s] is None else [k/10 for k in selector.SCHEDULES[s]]

def exact_bounds(text,ticks):
    if ticks is None:return text
    bounds=[0]+[k/10 for k in ticks];formatted=[format(x,'.17g') for x in bounds]
    lines=text.splitlines();i=0
    for j,line in enumerate(lines):
        if line=='window lower bound, window upper bound':
            found=list(map(float,lines[j+1].split(',')))
            assert all(a in (b,float(format(b,'.6g'))) for a,b in zip(found,bounds[i:i+2]))
            lines[j+1]=formatted[i]+', '+formatted[i+1];i+=1
        elif line.startswith('  epsilon breaks: '):lines[j]='  epsilon breaks: '+' '.join(formatted[1:])
        elif line.startswith('in eps range '):
            m=re.match(r'in eps range (\S+)\s+(\S+)(\s+dimension = .*)',line)
            if m:
                k=next(k for k in range(len(ticks)) if float(m[1]) in (bounds[k],float(format(bounds[k],'.6g'))) and float(m[2]) in (bounds[k+1],float(format(bounds[k+1],'.6g'))))
                lines[j]='in eps range '+formatted[k]+'  '+formatted[k+1]+m[3]
    assert i==len(ticks) or i==0
    return '\n'.join(lines)+'\n'

def parse(out,s,verbose=False):
    ticks=selector.SCHEDULES[s] if isinstance(s,str) else s
    text=exact_bounds((out/'stdout.log').read_text(),ticks)
    windows=[(Decimal(0),Decimal(2))] if ticks is None else [(Decimal(format(a/10,'.17g')),Decimal(format(b/10,'.17g'))) for a,b in zip([0]+ticks[:-1],ticks)]
    if not verbose:
        bars,ignored,win=lab.compact(text,'ph' if ticks is None else 'piecewise',2,2)
    else:
        p=out/'parser-input.log';p.write_text(text+'\niCC exit status: 0\n')
        argv=read(out/'context.json')['argv']
        meta=dict(dataset=argv[argv.index('--file-name')+1],dim=2,threads=4,breaks=[b for a,b in windows],alpha=Decimal('.85') if ticks is not None else None)
        bars,ignored=lab.verbose_parser.parse_log(p,meta,'ph' if ticks is None else 'piecewise')
    return bars[1],ignored,windows

def prepare(run):
    if (run/'prepare-phase.json').exists():return
    data=[]
    for split,repetitions,offset in [('dev',2,100),('selection',1,200),('final',1,300)]:
        index=0
        for rep in range(repetitions):
            for n in [64,100,200]:
                for geometry in ['uniform','clustered']:
                    seed=2026101600+offset+index;index+=1
                    data.append(dict(id=f'{split}-{geometry}-n{n}-s{seed}',split=split,geometry=geometry,n=n,seed=seed))
    plan=dict(created_at=now(),datasets=data,schedules=selector.SCHEDULES,
        scope=dict(simplex_dimension_max=2,PH='H0-H1',field=2,cutoff=2.,alpha=.85,threads=4),
        survey=dict(warmups=1,measured=3),final=dict(warmups=1,measured=5),
        fallback='noQ',ridge_candidates=[.01,1.],model_features=selector.FEATURE_NAMES,
        model_target='log whole-process seconds; predictions exponentiated to seconds',
        model_selection='lowest selection MAE seconds; dev fit only; empirical eligibility PASS on every dev+selection dataset',
        best_fixed='minimum geometric mean median whole-process runtime on dev+selection among empirically eligible schedules',
        final_policy='Generate final only after LOCK; evaluate A/B/C/D/E, no final oracle or schedule change',
        clustered_generator='six equally likely signed coordinate-axis centers + isotropic N(0,0.18²), normalize each vector; not coordinate scaling',
        uniform_generator='independent isotropic Gaussian triples normalized, uniform surface-area S2',
        timing='whole GNU-time process, G execv preserves one process across prediction and PH; verbose false',
        correctness='Exact C++ Decimal endpoint multisets; Ripser float32 projection without tolerance; actual per-window MSF H0 ledger; diagnostics separate',
        scope_excluded=['CONTINUE development','f2-stream','f2-rebuild','PV rule changes','matching fixes','n300','H2'],
        engine_commit='259d4bf',engine_reason='Frozen pre-f2 Morse implementation, all new methods share this source and binary')
    if not (run/'PLAN.json').exists():write(run/'PLAN.json',plan)
    if not (run/'engine-source').exists():
        (run/'engine-source').mkdir()
        with (run/'engine-source.tar').open('wb') as f:
            subprocess.run(['git','archive','259d4bf','CMakeLists.txt','main.cpp','deps','tests','LICENSE'],cwd=REPO,stdout=f,check=True)
        with tarfile.open(run/'engine-source.tar') as t:t.extractall(run/'engine-source',filter='data')
        (run/'engine-source.tar').unlink()
    cli=REPO/'build/_deps/cli11-src'
    if not (run/'third-party/CLI11').exists():shutil.copytree(cli,run/'third-party/CLI11',ignore=shutil.ignore_patterns('.git','build','__pycache__'))
    for name in ['offline.py','selector.py','README.md']:shutil.copy2(BASE/name,run/name)
    docs=[REPO/'experiments/pwph2026-provenance-continue/runs/provenance-20261016T025050Z/REPORT.md',
          REPO/'experiments/pwph2026-stage2-pilot/v2/CONTROLLER-v2.md',REPO/'experiments/pwph2026-stage2-pilot/v2/PH-ROOT-CAUSE.md',REPO/'experiments/pwph2026-stage2-pilot/v2/H0-LEDGER.md',
          REPO/'experiments/pwph2026-provenance-continue/runs/provenance-20261016T025050Z/case-verdicts.csv',
          REPO/'experiments/pwph2026-s2-n100/runs/n100-dossier-20261015T175443Z/REPORT.md',
          REPO/'experiments/pwph2026-provenance-continue/runs/provenance-20261016T025050Z/n200-dossier.json']
    dest=run/'context-read';dest.mkdir(exist_ok=True)
    for i,p in enumerate(docs):shutil.copy2(p,dest/f'{i:02}-{p.name}')
    write(run/'context-read.json',[dict(path=str(p),sha256=sha(p)) for p in docs])
    # Keep the complete reviewed helper dependency closure in the handoff.
    helpers=[REPO/'experiments/pwph2026-stage2-pilot/v2'/p for p in ['lab.py','h0_ledger.py']]
    helpers += [REPO/'experiments/pwph2026-stage2-pilot'/'pilot.py']
    helpers += [REPO/'experiments/pwph2026-s2-n100'/p for p in ['validate_n100.py','compare_icc_ph.py','evidence_run.py','benchmark_icc.py']]
    for p in helpers:
        q=run/'helper-source'/p.relative_to(REPO);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
    write(run/'helper-source.json',[dict(path=str(p),sha256=sha(p)) for p in helpers])
    note(run,'prepare',datasets=len(data),final_reserved=6)

def generate(run,item):
    dest=run/'data'/(item['id']+'.csv');dest.parent.mkdir(exist_ok=True)
    if dest.exists():return dest
    if item['split']=='final':assert (run/'LOCK.json').exists(),'Final cannot be generated before lock'
    rng=random.Random(item['seed']);points=[]
    centers=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    for _ in range(item['n']):
        if item['geometry']=='uniform':v=[rng.gauss(0,1) for j in range(3)]
        else:
            c=rng.choice(centers);v=[c[j]+rng.gauss(0,.18) for j in range(3)]
        norm=math.sqrt(sum(x*x for x in v));points.append([x/norm for x in v])
    with dest.open('w',newline='') as f:csv.writer(f).writerows([[format(x,'.17g') for x in p] for p in points])
    x=selector.raw_features(dest)
    write(dest.with_suffix('.json'),dict(**item,path=str(dest),sha256=sha(dest),features=x,generator_sha256=sha(BASE/'offline.py'),created_at=now()))
    return dest

def build(run):
    for name,argv in [
        ('configure',['cmake','-S',run/'engine-source','-B',run/'build','-DCMAKE_BUILD_TYPE=Release','-DBUILD_TESTING=ON','-DCMAKE_DISABLE_FIND_PACKAGE_CLI11=ON','-DFETCHCONTENT_FULLY_DISCONNECTED=ON','-DFETCHCONTENT_SOURCE_DIR_CLI11='+str(run/'third-party/CLI11')]),
        ('build',['cmake','--build',run/'build','-j4']),
        ('ctest',['ctest','--test-dir',run/'build','--output-on-failure'])]:
        o,c=execute(run,name,argv,timeout=300)
        assert c['status']=='OK',str(o)
    env={}
    for name,argv in [('compiler',['c++','--version']),('cmake',['cmake','--version']),('cpu',['lscpu']),('python',[sys.executable,'-m','pip','freeze'])]:
        o,c=execute(run,'environment-'+name,argv);env[name]=dict(path=str(o),exit_status=c['exit_status'])
    write(run/'build.json',dict(binary=str(binary(run)),sha256=sha(binary(run)),source_commit='259d4bf',
       source_files={str(p.relative_to(run/'engine-source')):sha(p) for p in (run/'engine-source').rglob('*') if p.is_file()},
       flags=(run/'build/CMakeFiles/iCC.dir/flags.make').read_text(),threads=4,environment=env,
       start_HEAD=read(run/'initial-context.json')['head'],engine_unchanged=True))

def validate(run,item,s,refs,baseline,base0):
    data=run/'data'/(item['id']+'.csv');b=binary(run)
    out,c=execute(run,'diagnostic-'+item['id']+'-'+s,selector.command(b,data,s,True),[b,data])
    r=dict(dataset=item['id'],split=item['split'],schedule=s,status=c['status'],exit_status=c['exit_status'],diagnostic=str(out))
    if c['status']!='OK':return r
    bars,ignored,windows=parse(out,s,True)
    h0=ledger(lab.msf((out/'stdout.log').read_text()),windows,base0,refs[0])
    exact=bars==baseline;rp=lab.projected(bars)==refs[1]
    r.update(H0=h0['status'],H1_exact_CPP=exact,H1_Ripser_float32=rp,H0_ledger=h0,
             bars=lab.serial(bars),ignored_temporary_survivors=ignored,
             missing=lab.serial(baseline-bars),extra=lab.serial(bars-baseline),
             windows=[[str(a),str(b)] for a,b in windows])
    r['status']='PH_MISMATCH' if not exact or not rp else ('PASS' if h0['status']=='PASS' else 'UNVERIFIED')
    write(out/'ph-check.json',r);return r

def baseline_reference(run,item):
    data=generate(run,item);b=binary(run);refs=lab.reference(run,data,2,2)
    out,c=execute(run,'diagnostic-'+item['id']+'-noQ',selector.command(b,data,'noQ',True),[b,data])
    assert c['status']=='OK'
    bars,_,_=parse(out,'noQ',True);frames=lab.msf((out/'stdout.log').read_text());base0=frames[0]['pairs']
    return refs,bars,base0

def survey(run):
    plan=read(run/'PLAN.json');rows=[];checks=[]
    for item in [x for x in plan['datasets'] if x['split']!='final']:
        refs,base,base0=baseline_reference(run,item);data=generate(run,item)
        checkmap={}
        for s in selector.SCHEDULES:
            r=validate(run,item,s,refs,base,base0);checks.append(r);checkmap[s]=r
            write(run/'survey-correctness.json',checks)
        for rep in range(4):
            order=list(selector.SCHEDULES);shift=(item['seed']+rep)%len(order);order=order[shift:]+order[:shift]
            if rep%2:order.reverse()
            for s in order:
                o,c=execute(run,f'survey-{item["id"]}-{s}-r{rep}',selector.command(binary(run),data,s),[binary(run),data])
                status=c['status']
                if status=='OK':
                    bars,_,_=parse(o,s)
                    status=checkmap[s]['status'] if bars==Counter({(Decimal(x['birth']),Decimal(x['death'])):x['multiplicity'] for x in checkmap[s].get('bars',[])}) else 'PH_MISMATCH'
                rows.append(dict(dataset=item['id'],split=item['split'],schedule=s,round=rep,warmup=rep==0,
                    elapsed_s=c['elapsed_monotonic_s'],rss_kib=int(c.get('resources',{}).get('max_rss_kib',-1)),
                    sampled_tree_rss_bytes=c['tree_rss_sampled_bytes'],status=status,exit_status=c['exit_status'],path=str(o)))
            write(run/'survey-timings.json',rows)
        note(run,'survey-complete-'+item['id'],statuses=dict(Counter(checkmap[s]['status'] for s in checkmap)))
    table(run/'survey-timings.csv',rows);table(run/'survey-correctness.csv',[{k:v for k,v in r.items() if k not in ['bars','H0_ledger','missing','extra']} for r in checks])
    # The old 8-point failing schedule is a sentinel, outside the candidate pool.
    old=REPO/'experiments/pwph2026-provenance-continue/runs/provenance-20261016T025050Z/inputs/two-squares-s2-8.csv'
    item=dict(id='sentinel-n8',split='sentinel',geometry='known-counterexample',n=8,seed=None)
    dest=run/'data/sentinel-n8.csv';shutil.copy2(old,dest)
    refs,base,base0=baseline_reference(run,item)
    results=[]
    for s in ['noQ','author']:
        results.append(validate(run,item,s,refs,base,base0))
    out,c=execute(run,'sentinel-n8-q05',selector.command(binary(run),dest,[5,20],True),[binary(run),dest])
    bars,_,windows=parse(out,[5,20],True)
    h0=ledger(lab.msf((out/'stdout.log').read_text()),windows,base0,refs[0])
    results.append(dict(dataset='sentinel-n8',schedule='[0.5,2]',status='PASS' if bars==base and lab.projected(bars)==refs[1] and h0['status']=='PASS' else 'PH_MISMATCH',H0=h0['status'],missing=lab.serial(base-bars),extra=lab.serial(bars-base),path=str(out),exit_status=c['exit_status']))
    write(run/'known-n8-sentinel.json',results)

def fit(run):
    checks=read(run/'survey-correctness.json');rows=read(run/'survey-timings.json');plan=read(run/'PLAN.json')
    items={x['id']:x for x in plan['datasets'] if x['split']!='final'}
    targets={}
    for id in items:
        for s in selector.SCHEDULES:
            values=[r['elapsed_s'] for r in rows if r['dataset']==id and r['schedule']==s and not r['warmup'] and r['exit_status']==0]
            if values:targets[id,s]=statistics.median(values)
    eligible=[s for s in selector.SCHEDULES if all(any(r['dataset']==id and r['schedule']==s and r['status']=='PASS' for r in checks) for id in items)]
    assert 'noQ' in eligible,'Baseline needs correctness investigation, not silent fallback'
    features={id:selector.raw_features(run/'data'/(id+'.csv')) for id in items}
    dev=[(id,s) for id,s in targets if items[id]['split']=='dev'];sel=[(id,s) for id,s in targets if items[id]['split']=='selection']
    X=np.asarray([selector.vector(features[id],s) for id,s in dev]);y=np.log([targets[key] for key in dev])
    mu=X.mean(0);scale=X.std(0);mu[0]=0;scale[0]=1;scale[scale<1e-12]=1;Z=(X-mu)/scale
    candidates=[]
    for ridge in plan['ridge_candidates']:
        penalty=np.eye(Z.shape[1])*ridge;penalty[0,0]=0
        w=np.linalg.solve(Z.T@Z+penalty,Z.T@y);coeff=w/scale;coeff[0]-=np.sum(w*mu/scale)
        pred=[math.exp(float(np.dot(selector.vector(features[id],s),coeff))) for id,s in sel]
        mae=statistics.mean(abs(p-targets[key]) for p,key in zip(pred,sel))
        candidates.append(dict(ridge=ridge,coefficients=coeff.tolist(),selection_MAE_s=mae))
    winner=min(candidates,key=lambda r:r['selection_MAE_s'])
    fixed_scores={s:math.exp(statistics.mean(math.log(targets[id,s]) for id in items)) for s in eligible}
    fixed=min(fixed_scores,key=fixed_scores.get)
    densities=[x['density']['8'] for x in features.values()]
    model=dict(created_at=now(),**winner,features=selector.FEATURE_NAMES,eligible_schedules=eligible,
               empirical_eligibility='PASS every dev+selection input, not a proof for new inputs',
               fallback='noQ',domain=dict(n=[64,200],density08=[max(0,min(densities)-.03),min(1,max(densities)+.03)]),
               reference='noQ',output_unit='seconds',scope=plan['scope'],schedules=selector.SCHEDULES)
    write(run/'model.json',model)
    report=dict(candidates=candidates,fixed_schedule=fixed,fixed_scores_s=fixed_scores,eligible_schedules=eligible,
                training_datasets=[id for id in items if items[id]['split']=='dev'],selection_datasets=[id for id in items if items[id]['split']=='selection'],
                features_by_dataset=features,targets=[dict(dataset=id,schedule=s,seconds=t) for (id,s),t in targets.items()])
    write(run/'fit.json',report)
    predictions=[]
    for (id,s),actual in targets.items():
        decision=selector.choose(features[id],model)
        predictions.append(dict(dataset=id,split=items[id]['split'],schedule=s,actual_s=actual,predicted_s=decision['predictions_s'][s],G_s=decision['G_seconds'][s],eligible=s in eligible,selected=decision['schedule']==s))
    table(run/'development-predictions.csv',predictions)
    lock=dict(locked_at=now(),model_sha256=sha(run/'model.json'),code={p.name:sha(p) for p in [BASE/'selector.py',BASE/'offline.py']},
              engine=read(run/'build.json'),plan_sha256=sha(run/'PLAN.json'),schedules=selector.SCHEDULES,
              fixed_schedule=fixed,fallback='noQ',final_datasets=[x for x in plan['datasets'] if x['split']=='final'],
              no_final_data_exists=not any((run/'data'/(x['id']+'.csv')).exists() for x in plan['datasets'] if x['split']=='final'),
              rule='Choose minimum predicted seconds among empirical eligible schedules; outside domain noQ; never alter after Ripser')
    assert lock['no_final_data_exists'];write(run/'LOCK.json',lock)
    note(run,'model-locked',fixed=fixed,eligible=eligible,ridge=winner['ridge'],selection_MAE_s=winner['selection_MAE_s'])

def final(run):
    lock=read(run/'LOCK.json');assert sha(run/'model.json')==lock['model_sha256']
    for p,h in lock['code'].items():assert sha(BASE/p)==h,'Locked experiment code changed'
    assert sha(binary(run))==lock['engine']['sha256']
    plan=read(run/'PLAN.json');model=read(run/'model.json');rows=[];checks=[];choices=[]
    for item in [x for x in plan['datasets'] if x['split']=='final']:
        data=generate(run,item);features=selector.raw_features(data);choice=selector.choose(features,model)
        # This is preflight only. No PH/reference has yet been computed on this input.
        choice.update(dataset=item['id'],input_sha256=sha(data));choices.append(choice);write(run/'choices.json',choices)
        methods=dict(A='noQ',B='author',C=lock['fixed_schedule'],D=choice['schedule'],E=choice['schedule'])
        first_h1={};method_paths={}
        for rep in range(6):
            order=list(methods);shift=(item['seed']+rep)%5;order=order[shift:]+order[:shift]
            if rep%2:order.reverse()
            for method in order:
                s=methods[method];choicepath=run/'choices'/(item['id']+f'-r{rep}.json')
                argv=([sys.executable,BASE/'selector.py','--binary',binary(run),'--data',data,'--model',run/'model.json','--choice',choicepath]
                      if method=='D' else selector.command(binary(run),data,s))
                o,c=execute(run,f'final-{item["id"]}-{method}-r{rep}',argv,[binary(run),data,run/'model.json'])
                status=c['status'];chosen_ok=True
                if method=='D' and status=='OK':chosen_ok=read(choicepath)['schedule']==choice['schedule']
                if status=='OK':
                    bars,_,_=parse(o,s)
                    if method not in first_h1:first_h1[method]=bars
                    if not chosen_ok or bars!=first_h1[method]:status='PH_MISMATCH'
                rows.append(dict(dataset=item['id'],method=method,schedule=s,round=rep,warmup=rep==0,
                    elapsed_s=c['elapsed_monotonic_s'],rss_kib=int(c.get('resources',{}).get('max_rss_kib',-1)),
                    sampled_tree_rss_bytes=c['tree_rss_sampled_bytes'],status=status,exit_status=c['exit_status'],path=str(o),
                    controller_pre_exec_s=read(choicepath)['pre_exec_s'] if method=='D' and choicepath.exists() else ''))
                method_paths[method]=str(o)
            write(run/'final-timings.json',rows)
        # Independent checks happen only after all locked performance runs for this input.
        refs,base,base0=baseline_reference(run,item)
        for s in dict.fromkeys(methods.values()):
            r=validate(run,item,s,refs,base,base0);checks.append(r)
        for method,s in methods.items():
            check=next(r for r in checks if r['dataset']==item['id'] and r['schedule']==s)
            expected=Counter({(Decimal(x['birth']),Decimal(x['death'])):x['multiplicity'] for x in check.get('bars',[])})
            for r in rows:
                if r['dataset']==item['id'] and r['method']==method and r['status']=='OK':
                    r['status']=check['status'] if first_h1[method]==expected else 'PH_MISMATCH'
        write(run/'final-timings.json',rows);write(run/'final-correctness.json',checks)
        note(run,'final-complete-'+item['id'],choice=choice['schedule'],statuses={m:next(r['status'] for r in rows if r['dataset']==item['id'] and r['method']==m) for m in methods})
    table(run/'final-timings.csv',rows)
    table(run/'choices.csv',[dict(dataset=x['dataset'],schedule=x['schedule'],ticks=x['ticks'],in_domain=x['in_domain'],reason=x['reason'],predicted_s=x['predictions_s'][x['schedule']],G_s=x['G_seconds'][x['schedule']]) for x in choices])
    table(run/'final-predictions.csv',[dict(dataset=x['dataset'],schedule=s,predicted_s=t,G_s=x['G_seconds'][s],chosen=s==x['schedule']) for x in choices for s,t in x['predictions_s'].items()])

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('phase',choices=['prepare','build','survey','fit','final'])
    a=p.parse_args();run=a.run.resolve();start=time.monotonic()
    if a.phase=='prepare':prepare(run)
    elif a.phase=='build':build(run)
    elif a.phase=='survey':survey(run)
    elif a.phase=='fit':fit(run)
    elif a.phase=='final':final(run)
    write(run/(a.phase+'-phase.json'),dict(completed_at=now(),elapsed_s=time.monotonic()-start,exit_status=0))

if __name__=='__main__':
    try:main()
    except Exception:
        args=sys.argv
        if '--run' in args:
            run=Path(args[args.index('--run')+1]);name='tooling-error-'+str(time.time_ns())
            write(run/(name+'.json'),dict(at=now(),argv=args,exit_status=1,traceback=traceback.format_exc()))
            note(run,name,'ERROR',exit_status=1)
        raise
