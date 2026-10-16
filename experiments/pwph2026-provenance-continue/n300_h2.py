#!/usr/bin/env python3
"""Bounded n300 native PwPH H2 probes; disclose absent exact H2 baseline."""
from matrix import *
from h0_ledger import ledger
import statistics

def main():
    item=dict(id='standard300',file=str(RUN/'inputs/standard300.csv'))
    prior=json.loads((RUN/'survey-validation.json').read_text())
    base2=next(x for x in prior if x['dataset']=='standard300' and x['dimension']==2 and x['method']=='native-one')
    base3=next(x for x in prior if x['dataset']=='standard300' and x['dimension']==3 and x['method']=='native-one')
    assert base2['status']=='PASS'
    baseline1=counters(base2)[1];baseline0=reviewed.msf((Path(base2['path'])/'stdout.log').read_text())[0]['pairs']
    refs=reviewed.reference(RUN,Path(item['file']),3)
    baseline2=counters(base3).get(2) if base3['status']=='PASS' else None
    validations=[];runs=[];summaries=[];eligible=[]
    configs=[('native-author',[1.2,1.7,2.]),('native-prior-fixed',[1.2,2.])]
    for method,grid in configs:
        r=one('S',dict(item,breaks=grid),3,4,verbose=True,tag='-n300-H2-probe-'+method)
        if r['status']=='OK':
            bars=counters(r);rp={str(k):reviewed.projected(v)==refs[k] for k,v in bars.items()}
            h0=ledger(reviewed.msf((Path(r['path'])/'stdout.log').read_text()),r['windows'],baseline0,refs[0])
            exact1=bars[1]==baseline1;exact2=bars[2]==baseline2 if baseline2 is not None else 'NOT_VERIFIED_NO_COMPLETE_NATIVE_BASELINE'
            r.update(H0=h0['status'],H0_ledger=h0,H1_CPP_exact=exact1,H2_CPP_exact=exact2,Ripser_float32=rp)
            r['status']='PASS_H2_RIPSER_FLOAT32_ONLY' if all(rp.values()) and h0['status']=='PASS' and exact1 and baseline2 is None else 'PASS' if all(rp.values()) and h0['status']=='PASS' and exact1 and exact2 is True else 'PH_MISMATCH'
        r.update(dataset='standard300',method=method,dimension=3,baseline_d3_status=base3['status']);validations.append(r)
        write(RUN/'n300-h2-validation.json',validations)
        if r['status'] in ('PASS','PASS_H2_RIPSER_FLOAT32_ONLY'):eligible.append((method,grid,r))
    for rep in range(6):
        order=eligible[rep%len(eligible):]+eligible[:rep%len(eligible)] if eligible else []
        if rep%2:order=list(reversed(order))
        for method,grid,val in order:
            r=one('S',dict(item,breaks=grid),3,4,tag=f'-n300-H2-{method}-round{rep}')
            good=r['status']=='OK' and r['diagrams']==val['diagrams']
            runs.append(dict(dataset='standard300',method=method,dimension=3,round=rep,phase='warmup' if rep==0 else 'measured',status=val['status'] if good else r['status'] if r['status']!='OK' else 'PH_MISMATCH',elapsed_s=r.get('elapsed_s'),rss_kib=r.get('rss_kib'),path=r['path']))
            write(RUN/'n300-h2-runs.json',runs);assert good
    for val in validations:
        cell=[x for x in runs if x['method']==val['method'] and x['phase']=='measured']
        summaries.append(dict(dataset='standard300',dimension=3,method=val['method'],status=val['status'],scope='H0/H1 exact native and float32 Ripser; H2 '+('float32 Ripser only, exact CPP baseline incomplete' if baseline2 is None else 'exact native and float32 Ripser'),warmups=sum(x['method']==val['method'] and x['phase']=='warmup' for x in runs),rounds=len(cell),median_s=statistics.median(x['elapsed_s'] for x in cell) if cell else '',min_s=min((x['elapsed_s'] for x in cell),default=''),max_s=max((x['elapsed_s'] for x in cell),default=''),max_rss_kib=max((x['rss_kib'] for x in cell),default=''),probe_elapsed_s=val.get('elapsed_s'),probe_rss_kib=val.get('rss_kib'),validation_path=val['path']))
    table(RUN/'n300-h2-performance.csv',summaries);table(RUN/'n300-h2-runs.csv',runs)
    write(RUN/'n300-h2-summary.json',dict(baseline_d3_status=base3['status'],H2_exact_native_verified=baseline2 is not None,performance=summaries,scope='n300 intermediate extension, not a TableI S2 configuration or full-paper reproduction'))
    progress('n300-H2-bounded','COMPLETE','H2 native-baseline limitation explicit; only individually PH-validated quotient configurations measured')

if __name__=='__main__':main()
