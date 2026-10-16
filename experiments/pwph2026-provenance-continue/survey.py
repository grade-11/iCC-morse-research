#!/usr/bin/env python3
"""Checkpoint-cost survey, no new controller/model; correctness gates each cell."""
from stream_checks import *
from h0_ledger import ledger
import statistics

def main():
    binary=Path(json.loads((RUN/'versions.json').read_text())['S']['binary'])
    items=[x for x in setup()['cases'] if x['id'] in ('pass32','standard100')]
    items.append(dict(id='two-squares-s2-8',file=str(RUN/'inputs/two-squares-s2-8.csv')))
    for n in (200,300):
        src=REPO/f'experiments/pwph2026-s2-n{n}/data/sphere_S2_R3_n{n}_seed20261007.csv'
        dst=RUN/'inputs'/f'standard{n}.csv';shutil.copy2(src,dst)
        items.append(dict(id=f'standard{n}',file=str(dst),source=str(src),sha256=sha(dst)))
    write(RUN/'survey-plan.json',dict(items=items,rounds=5,warmups=1,grid01=[k/10 for k in range(1,21)],scope='Descriptive checkpoint survey, previously inspected data, no G training/evaluation; F2 kernels serial; native threads4',binary_sha256=sha(binary)))
    validations=[];runs=[];summaries=[]
    for item in items:
      for d in ((2,3) if item['id']=='standard300' else (2,)):
        refs=reviewed.reference(RUN,Path(item['file']),d)
        # One direct bounded native run is both the resource probe and exact baseline.
        base=one('S',item,d,4,'ph',verbose=True,tag='-survey-probe')
        basebars=counters(base) if base['status']=='OK' else None
        basemsf=reviewed.msf((Path(base['path'])/'stdout.log').read_text())[0]['pairs'] if basebars else None
        if basebars and not all(reviewed.projected(v)==refs[k] for k,v in basebars.items()):
            base['status']='PH_MISMATCH';basebars=None
        definitions=[('native-one',None,[]),('native-author',[1.2,1.7,2.],[]),('native-prior-fixed',[1.2,2.],[])]
        if d==2:
            definitions += [('native-defer-grid01',[k/10 for k in range(1,21)],['--quotient-mask',*['0']*19]),('native-Q-grid01',[k/10 for k in range(1,21)],[])]
            definitions += [(engine+'-'+name,grid,['--continue-engine',engine]+(['--quotient-mask',*['0']*(len(grid)-1)] if len(grid)>1 else [])) for engine,name,grid in [('f2-stream','one',[2.]),('f2-stream','sparse',[1.2,1.7,2.]),('f2-stream','grid01',[k/10 for k in range(1,21)]),('f2-rebuild','grid01',[k/10 for k in range(1,21)])]]
        def invoke(definition,verbose,tag):
            method,grid,extra=definition
            if method.startswith('f2-'):
                result=stream(binary,item,grid,extra[1],4,verbose,tag)
                result['diagrams']={'1':reviewed.serial(result['bars'])} if result['status']=='OK' else {}
            else:result=one('S',dict(item,breaks=grid),d,4,'ph' if grid is None else 'quotient',extra=extra,verbose=verbose,tag=tag+'-'+method)
            return result
        eligible=[]
        for definition in definitions:
            method,grid,extra=definition
            if method=='native-one':r=base
            elif basebars is None:
                r=dict(status='NOT_RUN_NO_VALID_NATIVE_BASELINE',path=base['path'])
            else:r=invoke(definition,True,'survey-validation')
            if r['status']=='OK':
                bars=counters(r);exact=bars==basebars;rp=all(reviewed.projected(v)==refs[k] for k,v in bars.items())
                if method.startswith('f2-'):h0=dict(status='PASS' if r['h0']['2']==basemsf and reviewed.projected(r['h0']['2'])==refs[0] else 'H0_NOT_VERIFIED')
                else:
                    frames=reviewed.msf((Path(r['path'])/'stdout.log').read_text());h0=ledger(frames,r.get('windows') or [['0','2']],basemsf,refs[0])
                r.update(H0=h0['status'],H0_ledger=h0,exact_CPP=exact,Ripser_float32=rp)
                r['status']='PASS' if exact and rp and h0['status']=='PASS' else 'PH_MISMATCH'
            clean={k:v for k,v in r.items() if k not in ('bars','h0')}
            clean.update(dataset=item['id'],dimension=d,method=method);validations.append(clean)
            write(RUN/'survey-validation.json',validations)
            if r['status']=='PASS':eligible.append(definition)
        # Rotate order and reverse odd rounds; include one whole-process warmup.
        for rep in range(6):
            order=eligible[rep%len(eligible):]+eligible[:rep%len(eligible)] if eligible else []
            if rep%2:order=list(reversed(order))
            for definition in order:
                r=invoke(definition,False,f'survey-round{rep}')
                if r['status']=='OK':r['status']='PASS' if counters(r)==basebars else 'PH_MISMATCH'
                record=dict(dataset=item['id'],dimension=d,method=definition[0],phase='warmup' if rep==0 else 'measured',round=rep,status=r['status'],elapsed_s=r.get('elapsed_s'),rss_kib=r.get('rss_kib'),path=r['path'],binary_sha256=sha(binary),input_sha256=sha(item['file']))
                runs.append(record);write(RUN/'survey-runs.json',runs)
                assert r['status']=='PASS','Stop on failed measured configuration; preserve evidence'
        for definition in definitions:
            cell=[x for x in runs if x['dataset']==item['id'] and x['dimension']==d and x['method']==definition[0] and x['phase']=='measured']
            val=next(x for x in validations if x['dataset']==item['id'] and x['dimension']==d and x['method']==definition[0])
            summaries.append(dict(dataset=item['id'],dimension=d,method=definition[0],status=val['status'],scope='H0-H1' if d==2 else 'H0-H2',rounds=len(cell),median_s=statistics.median(x['elapsed_s'] for x in cell) if cell else '',min_s=min((x['elapsed_s'] for x in cell),default=''),max_s=max((x['elapsed_s'] for x in cell),default=''),median_rss_kib=statistics.median(x['rss_kib'] for x in cell) if cell else '',max_rss_kib=max((x['rss_kib'] for x in cell),default=''),validation_path=val['path']))
        table(RUN/'performance.csv',summaries);table(RUN/'survey-runs.csv',runs)
        progress('survey-'+item['id']+f'-d{d}','COMPLETE',f'{len(eligible)}/{len(definitions)} configurations passed; one warmup and five measured each')
    write(RUN/'survey-summary.json',summaries)

if __name__=='__main__':main()
