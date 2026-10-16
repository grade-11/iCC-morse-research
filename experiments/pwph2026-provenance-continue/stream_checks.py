#!/usr/bin/env python3
from matrix import *
from stream_reference import cells,StreamF2
import re,numpy as np

def stream(binary,item,grid,engine,threads=4,verbose=True,tag=''):
    extra=['--continue-engine',engine]+(['--quotient-mask',*['0']*(len(grid)-1)] if len(grid)>1 else [])
    out,c=execute('stream-'+item['id']+'-'+engine+'-t'+str(threads)+'-'+tag,reviewed.command(binary,Path(item['file']),grid,maxdim=2,threads=threads,verbose=verbose,extra=extra),[binary,Path(item['file']),BASE/'stream_checks.py'])
    if c['status']!='OK':return dict(status=c['status'],path=str(out),exit_status=c['exit_status'])
    text=(out/'stdout.log').read_text();start=text.index('window lower bound, window upper bound');filtered='\n'.join(l for l in text[start:].splitlines() if not l.startswith('[stream-'))
    bars,ignored,windows=reviewed.compact(exact_bounds(filtered,grid),'piecewise',2);h0={};states=[]
    for l in text.splitlines():
        if l.startswith('[stream-h0-pair]'):
            f=dict(t.split('=') for t in l.split()[1:]);h0.setdefault(f['eps'],Counter())[(Decimal(f['birth']),Decimal(f['death']))]+=1
        if l.startswith('[stream-state]'):states.append(dict(t.split('=') for t in l.split()[1:]))
    return dict(status='OK',path=str(out),bars=bars[1],h0=h0,states=states,elapsed_s=c['elapsed_monotonic_s'],rss_kib=int(c['resources']['max_rss_kib']))

def tests():
    versions=json.loads((RUN/'versions.json').read_text());binary=Path(versions['S']['binary']);items=[x for x in setup()['cases'] if x['small'] or x['id']=='standard100']
    for name in ['square','duplicates4','cutoff2','octahedron','disconnected']:items.append(dict(id=name,file=str(RUN/'inputs'/(name+'.csv')),breaks=[1.,1.5,2.]))
    decimal=RUN/'inputs/decimal-boundary.csv'
    np.savetxt(decimal,[[0,0],[.18,0],[.18,.24],[0,.24]],delimiter=',',fmt='%.17g')
    items.append(dict(id='decimal-boundary',file=str(decimal),breaks=[.2,.3,.4,2.]))
    rows=[]
    for item in items:
        baseline=one('S',item,2,4,'ph',tag='-stream-validation');want=counters(baseline)[1];refs=reviewed.reference(RUN,Path(item['file']),2)
        ho,hc=execute('stream-base-msf-'+item['id'],reviewed.command(binary,Path(item['file']),maxdim=1,verbose=True),[binary,Path(item['file'])]);base0=reviewed.msf((ho/'stdout.log').read_text())[0]['pairs']
        for t in (1,4):
            for name,grid in [('one',[2.]),('sparse',[.7,1.1,1.3,1.5,2.]),('grid01',[k/10 for k in range(1,21)])]:
                results={engine:stream(binary,item,grid,engine,t,tag=name) for engine in ['f2-stream','f2-rebuild']}
                for engine,r in results.items():
                    assert r['status']=='OK' and r['bars']==want,'Streaming C++ H1 differs from direct native'
                    assert r['h0']['2']==base0,'Streaming H0 differs from ACTUAL baseline MSF'
                    for eps,h0 in r['h0'].items():assert sum(h0.values())==int(np.loadtxt(item['file'],delimiter=',').shape[0]),'H0 conservation'
                    rp=reviewed.projected(r['bars'])==refs[1]
                    rows.append(dict(case=item['id'],threads=t,partition=name,engine=engine,status=('KNOWN_CUTOFF_EQUALITY' if not rp else 'PASS' if reviewed.projected(r['h0']['2'])==refs[0] else 'KNOWN_H0_ZERO_CONVENTION' if any(death==0 for birth,death in base0) else 'KNOWN_CUTOFF_EQUALITY'),H1_exact=True,H0_actual_MSF_exact=True,H0_Ripser_float32=reviewed.projected(r['h0']['2'])==refs[0],path=r['path'],columns_edges=int(r['states'][-1]['total_edges']),columns_triangles=int(r['states'][-1]['total_triangles'])))
                assert [s['state_hash'] for s in results['f2-stream']['states']]==[s['state_hash'] for s in results['f2-rebuild']['states']],'Retained state differs from fresh-prefix state'
                write(RUN/'stream-correctness.json',rows)
        # New mode MUST reject quotient and unsupported homology dimensions.
    errors=[];sample=items[0];data=Path(sample['file'])
    for name,extra,d in [('quotient',['--continue-engine','f2-stream','--quotient-mask','1'],2),('H2',['--continue-engine','f2-stream','--quotient-mask','0'],3),('G',['--continue-engine','f2-stream','--g-model',str(REPO/'experiments/pwph2026-stage2-pilot/v2/runs/v2-20261015T195640Z/models/global-d2.txt')],2)]:
        out,c=execute('stream-refuse-'+name,reviewed.command(binary,data,[1.2,2.],maxdim=d,extra=extra),[binary,data]);assert c['exit_status']!=0 and 'F2 stream/rebuild supports' in (out/'stderr.log').read_text();errors.append(dict(feature=name,status='EXPECTED_REJECTION',exit_status=c['exit_status'],path=str(out)))
    table(RUN/'stream-correctness.csv',rows);write(RUN/'stream-rejections.json',errors);progress('cpp-stream-tests','PASS_WITH_SCOPE','Explicit F2 engine H0/H1 no-quotient state vs rebuild vs native/ref verified; quotient/H2/G explicitly refused')
    print(json.dumps(dict(tests=len(rows),statuses=dict(Counter(x['status'] for x in rows)),rejections=errors),indent=2))

if __name__=='__main__':tests()
