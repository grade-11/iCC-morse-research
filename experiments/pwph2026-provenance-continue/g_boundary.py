#!/usr/bin/env python3
"""Bounded old-grid instrumentation diagnostic near artificial RAM forecasts.

Four predeclared forecast limits, no model fitting and no grid0.1 G evaluation.
External process guard stays at the original ~6.05GiB.
"""
from matrix import *
from h0_ledger import ledger

def main():
    items={x['id']:x for x in setup()['cases']};fixed=json.loads((RUN/'correctness.json').read_text());existing=[]
    for r in json.loads((RUN/'g-replay.json').read_text()):
        original=next(x for x in fixed if x['version']=='D' and x['case']==r['case'] and x['dimension']==3 and x['threads']==r['threads'] and x['mode']=='ph')
        data=Path(items[r['case']]['file']);grid=[1.2,1.5,1.7,2.];text=exact_bounds((Path(r['G_path'])/'stdout.log').read_text(),grid)
        if r['verbose']:
            parsed=RUN/'G-existing-parser-inputs'/(Path(r['G_path']).name+'.log');parsed.parent.mkdir(exist_ok=True);parsed.write_text(text+'\niCC exit status: 0\n')
            meta=dict(dataset=str(data),dim=3,threads=r['threads'],breaks=[Decimal(format(x,'.17g')) for x in grid],alpha=Decimal('.85'))
            bars,_=reviewed.verbose_parser.parse_log(parsed,meta,'piecewise')
        else:bars,_,_=reviewed.compact(text,'piecewise',3)
        refs=reviewed.reference(RUN,data,3)
        good=bars==counters(original) and all(reviewed.projected(v)==refs[k] for k,v in bars.items())
        existing.append(dict(case=r['case'],threads=r['threads'],repeat=r['repeat'],verbose=r['verbose'],mask=r['mask'],status='PASS' if good else 'PH_MISMATCH',G_path=r['G_path'],baseline_path=original['path'],replay_equal=r['PH_equal'],diagrams={str(k):reviewed.serial(v) for k,v in bars.items()},scope='Independent PH validation of previously saved G runs; no new run/timing'))
    write(RUN/'G-existing-reference.json',existing);table(RUN/'G-existing-reference.csv',[{k:(''.join(map(str,r[k])) if k=='mask' else r[k]) for k in ('case','threads','repeat','verbose','mask','status','G_path','baseline_path','replay_equal','scope')} for r in existing])
    item=dict(next(x for x in setup()['cases'] if x['id']=='standard100'),breaks=[1.2,1.5,1.7,2.])
    model=RUN/'n200-reused/models/g-d3-ridge0.01.txt'
    plan=dict(at=now(),dataset=item['id'],model_sha256=sha(model),forecast_limits_mib=[12,16,32,48],repetitions=2,verbose=[False,True],threads=4,grid=item['breaks'],scope='Old-grid diagnostic stress ONLY. Limits bracket observed ~7-26MiB branch forecasts; no final-data/G0.1 evaluation. External timeout/RSS guard unchanged.')
    if (RUN/'G-boundary-plan.json').exists():
        previous=json.loads((RUN/'G-boundary-plan.json').read_text());assert previous['forecast_limits_mib']==plan['forecast_limits_mib'] and previous['model_sha256']==plan['model_sha256'];plan=previous
    else:write(RUN/'G-boundary-plan.json',plan)
    baseline=one('D',item,3,4,'ph',verbose=True,tag='-Gboundary-baseline');assert baseline['status']=='OK'
    baseline0=reviewed.msf((Path(baseline['path'])/'stdout.log').read_text())[0]['pairs'];refs=reviewed.reference(RUN,Path(item['file']),3);basebars=counters(baseline);rows=[]
    for cap in plan['forecast_limits_mib']:
      for rep in (1,2):
       for verbose in (False,True):
        r=one('D',item,3,4,extra=['--g-model',str(model),'--g-margin','.05','--g-memory-mib',str(cap)],verbose=verbose,tag=f'-Gboundary-cap{cap}-r{rep}-v{int(verbose)}')
        gtext=(Path(r['path'])/'stdout.log').read_text();states=reviewed.parse_states(gtext)
        if r['status']!='OK':
            assert 'G_MEMORY_LIMIT: neither forecast action fits headroom' in (Path(r['path'])/'stderr.log').read_text(),'Unexpected native error'
            rows.append(dict(case=item['id'],forecast_limit_mib=cap,repeat=rep,verbose=verbose,mask=[int(s['action']=='QUOTIENT') for s in states if s['action']!='FINAL'],trajectory_complete=False,PH_replay_equal=None,state_replay_equal=None,H0_ledger=dict(status='NOT_VERIFIED_INCOMPLETE'),PH_reference_status='NOT_VERIFIED_G_FORECAST_MEMORY_LIMIT',diagrams={},G_path=r['path'],replay_path='',state_records=states,G_decision_records=[],exit_status=r['exit_status']))
            write(RUN/'G-boundary.json',rows);progress(f'G-boundary-cap{cap}-r{rep}-v{int(verbose)}','G_FORECAST_MEMORY_LIMIT','Exit1 is an intentional controller forecast refusal, external RSS guard not triggered; partial PH not validated');continue
        mask=[int(s['action']=='QUOTIENT') for s in states[:-1]]
        decisions=[dict(field.split('=',1) for field in line.split()[1:]) for line in gtext.splitlines() if line.startswith('[g]')]
        q=one('D',item,3,4,extra=['--cost-telemetry','--quotient-mask',*map(str,mask)],verbose=True,tag=f'-Gboundary-replay-cap{cap}-r{rep}-v{int(verbose)}')
        assert q['status']=='OK';text=(Path(q['path'])/'stdout.log').read_text();rs=reviewed.parse_states(text)
        equal=counters(r)==counters(q);state_equal=[s['state_hash'] for s in states]==[s['state_hash'] for s in rs]
        if not (equal and state_equal):
            write(RUN/'G-boundary-failure.json',dict(G_path=r['path'],replay_path=q['path'],mask=mask,PH_equal=equal,state_equal=state_equal));progress('G-boundary-unexpected','PH_REPLAY_MISMATCH','Raw evidence retained; investigate before continuation')
            raise AssertionError('G/fixed-action replay differs')
        h0=ledger(reviewed.msf(text),q['windows'],baseline0,refs[0]);assert h0['status']=='PASS'
        ph=counters(r)==basebars and all(reviewed.projected(v)==refs[k] for k,v in counters(r).items())
        rows.append(dict(case=item['id'],forecast_limit_mib=cap,repeat=rep,verbose=verbose,mask=mask,trajectory_complete=True,PH_replay_equal=equal,state_replay_equal=state_equal,H0_ledger=h0,PH_reference_status='PASS' if ph else 'PH_MISMATCH',diagrams=r['diagrams'],G_path=r['path'],replay_path=q['path'],state_records=states,G_decision_records=decisions,exit_status=0))
        write(RUN/'G-boundary.json',rows)
    differences=[]
    for cap in plan['forecast_limits_mib']:
        selected=[x for x in rows if x['forecast_limit_mib']==cap]
        if len({tuple(x['mask']) for x in selected})>1:differences.append(dict(forecast_limit_mib=cap,observations=[{k:x[k] for k in ('repeat','verbose','mask','PH_reference_status')} for x in selected],PH_multiset_varies=len({json.dumps(x['diagrams'],sort_keys=True) for x in selected})>1))
    upstream=[]
    allq_failure=any(x['status']=='PH_MISMATCH' and all(x['mask']) and x['case']=='standard100' for x in existing) or any(x['trajectory_complete'] and x['PH_reference_status']=='PH_MISMATCH' and all(x['mask']) for x in rows)
    if allq_failure:
        for version in ('U','U17'):
            r=one(version,item,3,4,tag='-G-selected-allQ-grid-upstream')
            assert r['status']=='OK'
            oldbase=next(x for x in fixed if x['version']==version and x['case']=='standard100' and x['dimension']==3 and x['threads']==4 and x['mode']=='ph')
            bad=counters(r)!=counters(oldbase);assert bad
            if version=='U17':
                matching=[x for x in rows if x['trajectory_complete'] and all(x['mask'])]
                matching += [x for x in existing if x['case']=='standard100' and all(x['mask'])]
                for x in matching:assert x['diagrams']==r['diagrams']
            upstream.append(dict(version=version,status='PH_MISMATCH',result=r,baseline_path=oldbase['path'],scope='Upstream supports this actual allQ epsilon schedule without mask adaptation'))
    write(RUN/'G-selected-upstream.json',upstream)
    paired=[dict(forecast_limit_mib=cap,repeat=rep,mask_off=next(x['mask'] for x in rows if x['forecast_limit_mib']==cap and x['repeat']==rep and not x['verbose']),mask_on=next(x['mask'] for x in rows if x['forecast_limit_mib']==cap and x['repeat']==rep and x['verbose'])) for cap in plan['forecast_limits_mib'] for rep in (1,2)]
    paired=[x for x in paired if x['mask_off']!=x['mask_on']]
    table(RUN/'G-boundary.csv',[{k:(''.join(map(str,r[k])) if k=='mask' else r[k]) for k in ('case','forecast_limit_mib','repeat','verbose','mask','trajectory_complete','exit_status','PH_replay_equal','state_replay_equal','PH_reference_status','G_path','replay_path')} for r in rows])
    summary=dict(status='TRAJECTORY_DIFFERENCE_OBSERVED' if differences else 'NO_OBSERVED_TRAJECTORY_DIFFERENCE',runs=len(rows),complete_trajectories=sum(x['trajectory_complete'] for x in rows),forecast_refusals=sum(not x['trajectory_complete'] for x in rows),differences=differences,paired_verbosity_differences=paired,PH_reference_statuses=dict(Counter(x['PH_reference_status'] for x in rows)),existing52_PH_reference_statuses=dict(Counter(x['status'] for x in existing)),allQ_failure_upstream_reproduced=bool(upstream),mixed_mask_failure_verdict='NOT_EQUIVALENT_UPSTREAM_MASK; no claim of group causality or upstream reproduction for unsupported mixed masks',all_completed_replays_equal=True,H0_actual_ledger_pass=sum(x['trajectory_complete'] for x in rows),scope=plan['scope'],interpretation='Verbose also resets VmHWM per window and adds diagnostic work; a changed decision or refusal would not isolate pure wall-clock timing. Compare prefixes/refusals separately from complete trajectories. Fixed-mask replay isolates algorithmic output for completed trajectories; incomplete PH never called PASS.')
    write(RUN/'G-boundary-summary.json',summary);progress('G-boundary',summary['status'],str(summary));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
