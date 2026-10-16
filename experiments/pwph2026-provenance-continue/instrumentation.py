#!/usr/bin/env python3
import shutil,difflib
from matrix import *

def variants():
    versions=json.loads((RUN/'versions.json').read_text());original=Path(versions['D']['source'])
    for label in ['D_NO_CLOCK','D_NO_INTERNAL_RAM']:
        if label in versions and Path(versions[label]['binary']).exists():continue
        dst=RUN/'versions'/label;shutil.copytree(original,dst,dirs_exist_ok=True);patch=[]
        for p in [dst/'main.cpp',*list((dst/'deps').glob('*.cpp'))]:
            before=p.read_text();after=before
            if label=='D_NO_CLOCK':
                after=after.replace('Clock::now()','Clock::time_point{}').replace('std::chrono::steady_clock::now()','std::chrono::steady_clock::time_point{}').replace('std::chrono::high_resolution_clock::now()','std::chrono::high_resolution_clock::time_point{}')
            elif p.name=='Pipelines.cpp':
                after=after.replace(': enabled_(enabled), eps_lo_',': enabled_(false), eps_lo_')
                after=after.replace('const auto peak_kib=[]()\n    {','const auto peak_kib=[]()\n    {\n        return 0.; // Controlled fixed-action diagnostic: no internal RSS observation\n')
            if after!=before:p.write_text(after);patch+=list(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile=str(p.relative_to(dst)),tofile=label+'/'+str(p.relative_to(dst))))
        (RUN/'versions'/(label+'.patch')).write_text(''.join(patch));build=RUN/'builds'/label
        for name,argv in [('configure',['cmake','-S',str(dst),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DBUILD_TESTING=OFF','-DFETCHCONTENT_SOURCE_DIR_CLI11='+str(REPO/'build/_deps/cli11-src'),'-DFETCHCONTENT_FULLY_DISCONNECTED=ON']),('build',['cmake','--build',str(build),'-j4'])]:
            out,c=execute(name+'-'+label,argv,[dst/'deps/Pipelines.cpp',dst/'deps/DimensionFrame.cpp']);print(label,name,c['status'],flush=True);assert c['status']=='OK'
        versions[label]=dict(commit=versions['D']['commit'],source=str(dst),binary=str(build/'iCC'),binary_sha256=sha(build/'iCC'),changes='Controlled fixed-action diagnostic only: '+label,source_hashes={str(p.relative_to(dst)):sha(p) for p in dst.rglob('*') if p.is_file() and ('deps' in p.parts or p.name in ('main.cpp','CMakeLists.txt'))});write(RUN/'versions.json',versions)

def direct():
    variants();items=[x for x in setup()['cases'] if x['small']];rows=[]
    for item in items:
        for d in (2,3):
            for t in (1,4):
                ref=one('D',item,d,t);expected=counters(ref)
                settings=[('default','D',[],False),('telemetry','D',['--cost-telemetry'],False),('verbose','D',[],True),('trace','D',['--trace-pwph','--cost-telemetry'],True),('clock-off','D_NO_CLOCK',['--cost-telemetry'],False),('internal-ram-off','D_NO_INTERNAL_RAM',['--cost-telemetry'],True)]
                for repeat in (1,2):
                    for name,v,extra,verbose in settings:
                        mask=['--quotient-mask',*['1']*(len(item['breaks'])-1)]
                        r=one(v,item,d,t,extra=extra+mask,verbose=verbose,tag=f'-instrument-{name}-r{repeat}')
                        equal=r['status']=='OK' and counters(r)==expected
                        r.update(setting=name,repeat=repeat,PH_equal_with_locked_actions=equal,status='PASS' if equal else 'PH_MISMATCH',scope='Instrumentation isolation, not PH correctness: original known failures remain failures versus Ripser')
                        rows.append(r);write(RUN/'instrumentation.json',rows);assert equal,'Instrumentation changed PH with fixed action schedule'
    table(RUN/'instrumentation.csv',[{k:r.get(k,'') for k in ['case','dimension','threads','setting','repeat','status','PH_equal_with_locked_actions','path']} for r in rows]);progress('instrumentation-fixed','PASS','Clock, internalRAM, verbosity/FFI, telemetry/trace toggles preserve printed PH under fixed actions on small cases')

def g_replay():
    old=REPO/'experiments/pwph2026-stage2-pilot/v2/runs/v2-20261015T195640Z';model=old/'models/g-d3-ridge0.01.txt';rows=[]
    # Existing model is used ONLY to isolate G/replay; no claim of new-grid calibration.
    for item in [x for x in setup()['cases'] if x['small'] or x['id']=='standard100']:
        probe=dict(item,breaks=[1.2,1.5,1.7,2.])
        for t in ((1,4) if item['small'] else (4,)):
            for repeat in (1,2):
                for verbose in (False,True):
                    r=one('D',probe,3,t,extra=['--g-model',str(model),'--g-margin','.05'],verbose=verbose,tag=f'-G-r{repeat}-v{int(verbose)}')
                    states=reviewed.parse_states((Path(r['path'])/'stdout.log').read_text());mask=[int(s['action']=='QUOTIENT') for s in states[:-1]]
                    replay=one('D',probe,3,t,extra=['--cost-telemetry','--quotient-mask',*map(str,mask)],verbose=False,tag=f'-G-replay-r{repeat}-v{int(verbose)}')
                    rs=reviewed.parse_states((Path(replay['path'])/'stdout.log').read_text());same=counters(r)==counters(replay);state_equal=[s['state_hash'] for s in states]==[s['state_hash'] for s in rs]
                    assert same and state_equal,'G/replay differs at fixed actual trajectory'
                    rows.append(dict(case=item['id'],threads=t,repeat=repeat,verbose=verbose,mask=mask,PH_equal=same,state_equal=state_equal,G_path=r['path'],replay_path=replay['path'],model_sha256=sha(model),scope='Existing-grid diagnostic only, not fresh controller evaluation'))
                    write(RUN/'g-replay.json',rows)
    differences=[]
    for item in [x for x in setup()['cases'] if x['small'] or x['id']=='standard100']:
        for t in (1,4):
            rr=[r for r in rows if r['case']==item['id'] and r['threads']==t]
            if len({tuple(r['mask']) for r in rr})>1:differences.append(dict(case=item['id'],threads=t,trajectories=[dict(repeat=r['repeat'],verbose=r['verbose'],mask=r['mask']) for r in rr]))
    write(RUN/'g-trajectory-differences.json',differences);table(RUN/'g-replay.csv',[{**r,'mask':''.join(map(str,r['mask']))} for r in rows]);print('Observed trajectory differences',len(differences))

if __name__=='__main__':
    direct();g_replay()
