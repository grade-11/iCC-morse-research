#!/usr/bin/env python3
"""At most two new fixed-schedule failures, <=10 minutes each, H1 priority (then H2)."""
from study import *

def minimize(run,binary,source):
    folder=run/'new-shrunk'/source['dataset_id'];folder.mkdir(parents=True,exist_ok=True)
    if (folder/'result.json').exists():return json.loads((folder/'result.json').read_text())
    points=np.loadtxt(source['dataset'],delimiter=',');ids=list(range(len(points)));begin=time.monotonic();attempt=0;last=None;d=source['maxdim']
    def test(candidate):
        nonlocal attempt,last
        path=folder/f'candidate-{attempt:04d}.csv';attempt+=1;np.savetxt(path,points[candidate],delimiter=',',fmt='%.17g')
        bars=[]
        for q in (0,1):
            out,c=execute(run,f'new-shrink-{source["dataset_id"]}-{attempt:04d}-q{q}',command(binary,path,None if q==0 else source['breaks'],maxdim=d),[binary,path,Path(__file__)],timeout=max(.1,min(300,540-(time.monotonic()-begin))))
            if c['status']!='OK':return False
            bars.append(compact((out/'stdout.log').read_text(),'ph' if q==0 else 'piecewise',d)[0])
        mismatch=bars[0]!=bars[1]
        detail=dict(point_ids=candidate,mismatch=mismatch,missing={str(k):serial(bars[0][k]-bars[1][k]) for k in bars[0]},extra={str(k):serial(bars[1][k]-bars[0][k]) for k in bars[0]})
        write(path.with_suffix('.json'),detail)
        if mismatch:last=detail
        return mismatch
    reproduced=test(ids);granularity=2
    while reproduced and len(ids)>4 and time.monotonic()-begin<500:
        size=math.ceil(len(ids)/granularity);changed=False
        for k in range(0,len(ids),size):
            candidate=ids[:k]+ids[k+size:]
            if len(candidate)<4:continue
            if time.monotonic()-begin>=500:break
            if test(candidate):ids=candidate;granularity=max(2,granularity-1);changed=True;break
        if not changed:
            if granularity>=len(ids):break
            granularity=min(len(ids),granularity*2)
    retained=folder/'retained.csv';np.savetxt(retained,points[ids],delimiter=',',fmt='%.17g')
    check_result=check(run,'confirm-new-shrunk-'+source['dataset_id'],binary,retained,source['breaks'],maxdim=d)
    trace,tc=execute(run,'trace-new-shrunk-'+source['dataset_id'],command(binary,retained,source['breaks'],maxdim=d,verbose=True,extra=['--trace-pwph']),[binary,retained,Path(__file__)])
    result=dict(case=source['dataset_id'],source_config=source['case'],method=source['method'],grid=source['breaks'],dimension=d,status='PH_MISMATCH' if reproduced else 'NOT_REPRODUCED',exit_status=0,
                original_n=len(points),retained_n=len(ids),point_ids=ids,attempts=attempt,elapsed_s=time.monotonic()-begin,retained_sha256=sha(retained),difference=last,
                independent_ripser=check_result,trace_path=str(trace),claim='Reduced case, not a global minimum; investigation only, not performance')
    write(folder/'result.json',result);note(run,'new-shrunk-'+source['dataset_id'],status=result['status'],exit_status=0);return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);a=ap.parse_args();guard=runner_lock(a.run)
    reports=json.loads((a.run/'screen-dev.json').read_text())
    candidates=[x for x in reports if x['status']=='PH_MISMATCH' and x['method'] in ('original','sparse','early','dense')]
    candidates.sort(key=lambda x:(x['exact'].get('1',True),x['n'],x['method']!='original',x['dataset_id']))
    chosen=[]
    for x in candidates:
        if len(chosen)==2:break
        if not chosen or (x['geometry']!=chosen[0]['geometry'] and x['dataset_id']!=chosen[0]['dataset_id']):chosen.append(x)
    if len(chosen)==1:
        other=next((x for x in candidates if x['dataset_id']!=chosen[0]['dataset_id']),None)
        if other:chosen.append(other)
    write(a.run/'new-shrink-plan.json',[dict(dataset=x['dataset_id'],method=x['method'],dimension=x['maxdim']) for x in chosen])
    results=[minimize(a.run,a.binary.resolve(),x) for x in chosen];write(a.run/'new-shrink-summary.json',results)
