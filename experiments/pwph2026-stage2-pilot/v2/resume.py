#!/usr/bin/env python3
"""Continue missing stages, never reset the monotonic/boot budget."""
import argparse,fcntl
from lab import *
from evidence_run import capture

def resume(run):
    guard=(run/'workflow.lock').open('a');fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if remaining(run,pack=True)<=0:
        print('Original four-hour deadline expired; no experiments or report rewrites. Existing handoff retained.');return
    binary=run/'build/iCC';need(binary.exists(),'Build binary missing; archived handoff needs a fresh run/build identity')
    lockpath=run/'LOCK-v2.json'
    if lockpath.exists():
        lock=json.loads(lockpath.read_text());need(sha(binary)==lock['binary_sha256'],'Locked binary changed')
        for p,h in {**lock['source_hashes'],**lock['model_hashes']}.items():need(sha(p)==h,'Locked input changed: '+p)
    if (run/'completion.json').exists():
        print('Run already completed with documented PH limitations; no pending stages, no reruns.');return
    stages=[]
    for split in ('dev','validation'):
        stages.append((f'calibration-{split}.json','study.py',['calibrate','--split',split]))
    stages += [('selected-models.json','study.py',['fit']),('active-models.json','ablate.py',[])]
    for split in ('dev','validation'):stages.append((f'screen-{split}.json','study.py',['screen','--split',split]))
    stages.append(('screen-h0-enhanced-validation.json','enhance_h0.py',[]))
    for split in ('dev','validation'):stages.append((f'calibration-global-{split}.json','global_horizon.py',['collect','--split',split]))
    stages.append(('global-models.json','global_horizon.py',['fit']))
    for split in ('dev','validation'):stages.append((f'screen-global-{split}.json','global_horizon.py',['screen','--split',split]))
    stages += [('new-shrink-summary.json','shrink_new.py',[]),('stability.json','stability.py',[]),('ram-study.json','ram_probe.py',[]),
               ('LOCK-v2.json','study.py',['lock']),('final-validation.json','study.py',['final'])]
    for marker,script,args in stages:
        if (run/marker).exists():continue
        if remaining(run)<=0:
            note(run,'resume-budget-stop',status='BUDGET_STOP',exit_status=0,missing_stage=marker);break
        need(not lockpath.exists() or marker=='final-validation.json','Locked run lacks pre-final stage; review rather than retrain after final')
        cmd=[sys.executable,str(BASE/script),*args,'--run',str(run),'--binary',str(binary)]
        if script in ('ablate.py','enhance_h0.py'):cmd=cmd[:-2]
        out=run/('resume-'+script.replace('.py','')+'-'+str(time.monotonic_ns()))
        c=capture(out,cmd,[BASE/script,BASE/'lab.py'],remaining(run));need(c==0,'Stage failed: '+str(out))
    for script in ('summarize.py','report.py'):
        if script=='report.py' and not (run/'final-validation.json').exists():continue
        c=capture(run/('resume-'+script[:-3]+'-'+str(time.monotonic_ns())),[sys.executable,str(BASE/script),'--run',str(run)],[BASE/script],max(1,remaining(run,pack=True)))
        need(c==0,'Summary failed')
    print('Resume complete; existing evidence retained, original deadline unchanged. Packaging is explicit: package.py --run RUN.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,default=run_path());ap.add_argument('--worker',action='store_true',help=argparse.SUPPRESS);a=ap.parse_args();run=a.run.resolve()
    if a.worker:resume(run)
    else:
        out=run/('resume-invocation-'+str(time.monotonic_ns()))
        rc=capture(out,[sys.executable,str(Path(__file__).resolve()),'--run',str(run),'--worker'],
                   [Path(__file__),BASE/'lab.py'],max(30,remaining(run,pack=True)))
        print(f'Resume exit {rc}; evidence: {out}')
        raise SystemExit(rc)
