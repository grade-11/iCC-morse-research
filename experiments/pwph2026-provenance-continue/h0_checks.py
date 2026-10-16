#!/usr/bin/env python3
from matrix import *
from h0_ledger import ledger

def check():
    rows=[];versions=json.loads((RUN/'versions.json').read_text());items={x['id']:x for x in setup()['cases']};baselines={}
    for id,item in items.items():
        if not item['small'] and id!='standard100':continue
        data=Path(item['file']);refs=reviewed.reference(RUN,data,3);binary=Path(versions['D']['binary'])
        bo,bc=execute('h0-direct-'+id,reviewed.command(binary,data,maxdim=1,verbose=True),[binary,data]);assert bc['status']=='OK';baselines[id]=reviewed.msf((bo/'stdout.log').read_text())[0]['pairs']
        for v in ('C','D_NO_CLOCK'):
            q=one(v,item,3,4,verbose=True,extra=['--cost-telemetry','--quotient-mask',*['1']*(len(item['breaks'])-1)] if v=='D_NO_CLOCK' else [],tag='-h0')
            text=(Path(q['path'])/'stdout.log').read_text();h=ledger(reviewed.msf(text),q['windows'],baselines[id],refs[0]);assert h['status']=='PASS';rows.append(dict(case=id,version=v,mode='fixed',dimension=3,threads=4,H0=h['status'],path=q['path'],ledger=h))
    for q in json.loads((RUN/'instrumentation.json').read_text()):
        if not q['verbose']:continue
        text=(Path(q['path'])/'stdout.log').read_text();refs=reviewed.reference(RUN,Path(items[q['case']]['file']),q['dimension']);h=ledger(reviewed.msf(text),q['windows'],baselines[q['case']],refs[0]);assert h['status']=='PASS';rows.append(dict(case=q['case'],version=q['version'],mode=q['setting'],dimension=q['dimension'],threads=q['threads'],H0=h['status'],path=q['path'],ledger=h))
    for q in json.loads((RUN/'g-replay.json').read_text()):
        item=dict(items[q['case']],breaks=[1.2,1.5,1.7,2.]);replay=one('D',item,3,q['threads'],verbose=True,extra=['--cost-telemetry','--quotient-mask',*map(str,q['mask'])],tag='-h0-Gmask'+''.join(map(str,q['mask'])))
        gs=reviewed.parse_states((Path(q['G_path'])/'stdout.log').read_text());text=(Path(replay['path'])/'stdout.log').read_text();rs=reviewed.parse_states(text);assert [s['state_hash'] for s in gs]==[s['state_hash'] for s in rs]
        h=ledger(reviewed.msf(text),replay['windows'],baselines[q['case']],reviewed.reference(RUN,Path(item['file']),3)[0]);assert h['status']=='PASS';rows.append(dict(case=q['case'],version='D',mode='G-actual-mask-replay',dimension=3,threads=q['threads'],H0=h['status'],path=replay['path'],G_path=q['G_path'],ledger=h))
    write(RUN/'h0-verification.json',rows);table(RUN/'h0-verification.csv',[{k:r[k] for k in ['case','version','mode','dimension','threads','H0','path']} for r in rows]);print('H0 actual MSF ledger PASS',len(rows))
    # Upstream + output-precision-only patch, one window => no quotient.
    item=dict(id='duplicates4',file=str(RUN/'inputs/duplicates4.csv'),breaks=[2.])
    zero=[]
    for dimension in (2,3):
        u=one('U17',item,dimension,1,'ph',tag='-zero');q=one('U17',item,dimension,1,tag='-zero')
        zero.append(dict(dimension=dimension,equal=counters(u)==counters(q),native_ph=u,native_piecewise_noquotient=q))
    write(RUN/'upstream-zero-exclusion.json',dict(upstream=versions['U']['commit'],precision_only=True,quotient_count=0,input_sha256=sha(item['file']),results=zero,conclusion='Known duplicate-point/zero-edge exclusion: record actual differences separately by dimension; no quotient occurs in the single-window piecewise run. No inference from the failed maxdim2 hypothesis.'))
    # Verify old native behavior was not changed by adding the opt-in engine.
    old=json.loads((RUN/'correctness.json').read_text());reg=[]
    for r in old:
        if r['version']!='D' or r['threads']!=4 or r['case'] not in ('old10','old11','new11','new14','standard100','pass32'):continue
        nowr=one('S',items[r['case']],r['dimension'],4,r['mode'],tag='-native-regression');equal=nowr['diagrams']==r['diagrams'];assert equal
        reg.append(dict(case=r['case'],dimension=r['dimension'],mode=r['mode'],native_unchanged=equal,old_status=r['status'],new_path=nowr['path'],old_path=r['path']))
    write(RUN/'stream-native-regression.json',reg);progress('H0-and-native-regression','PASS','Actual native H0 verified where exported; no native PH change from opt-in F2 engine; zero edge exclusion reproduced in U17 without quotient')

if __name__=='__main__':check()
