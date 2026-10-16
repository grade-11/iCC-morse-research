#!/usr/bin/env python3
from matrix import *
from stream_checks import stream
from stream_reference import cells,StreamF2
import numpy as np,math

def run():
    planar=np.array([[0,0],[0,.4],[.4,0],[.4,.4],[1.4,0],[1.4,.4],[1.8,0],[1.8,.4]],float)
    s2=[]
    for x,y in planar:
        u=(x-.9)/1.5;v=(y-.2)/1.5;s2.append([u,v,math.sqrt(1-u*u-v*v)])
    p=RUN/'inputs/two-squares-s2-8.csv';np.savetxt(p,s2,delimiter=',',fmt='%.17g')
    write(RUN/'counterexample8-plan.json',dict(at=now(),data_sha256=sha(p),construction='Two .4 planar squares translated1.4 in x; center(.9,.2), divide1.5; lift(x,y,sqrt(1-x²-y²)) to unitS2. All8 points distinct; clustered deterministic, not uniform TableI data.',norms=[float(np.linalg.norm(x)) for x in s2],single_quotient=[.5,2.],dense_grid=[k/10 for k in range(1,21)],alpha=.85,max_simplex=2,field=2,metric='euclidean'))
    item=dict(id='two-squares-s2-8',file=str(p),breaks=[.5,2.]);refs=reviewed.reference(RUN,p,2);rows=[]
    for v in ['U','U17','B17','C','D','S']:
        for t in (1,4):
            base=one(v,item,2,t,'ph',tag='-S2-counterexample');b=counters(base)
            for name,grid in [('single',[.5,2.]),('grid01',[k/10 for k in range(1,21)]),('author',[1.2,1.7,2.])]:
                qq=dict(item,breaks=grid);q=one(v,qq,2,t,tag='-S2-'+name)
                same=counters(q)==b;ref='PRINTED_PRECISION_ONLY' if v=='U' else reviewed.projected(b[1])==refs[1]
                rows.append(dict(version=v,threads=t,schedule=name,baseline_matches_Ripser=ref,quotient_equal_baseline=same,baseline=base,quotient=q))
                write(RUN/'counterexample8.json',rows)
    binary=Path(json.loads((RUN/'versions.json').read_text())['S']['binary'])
    inc=stream(binary,item,[k/10 for k in range(1,21)],'f2-stream',4,tag='counterexample8');assert inc['status']=='OK' and reviewed.projected(inc['bars'])==refs[1]
    trace=one('S',item,2,4,verbose=True,extra=['--trace-pwph','--cost-telemetry'],tag='-trace-counterexample8')
    text=(Path(trace['path'])/'stdout.log').read_text();pvs={};traces={}
    for line in text.splitlines():
        if line.startswith('[pv-members]'):
            f=dict(x.split('=') for x in line.split()[1:]);pvs[int(f['label'])]=list(map(int,f['representatives'].strip(',').split(',')))
        elif line.startswith('[simplex]'):
            f=dict(x.split('=') for x in line.split()[1:]);traces.setdefault(float(f['eps_hi']),{})[tuple(sorted(map(int,f['labels'].strip(',').split(','))))]=float(f['weight'])
    from oracle import witness_complex,reduction
    pts=np.loadtxt(p,delimiter=',');windows=[]
    for eps,pv in [(.5,{}),(2.,pvs)]:
        geometric=witness_complex(pts,pv,eps,maxdim=2);actual=traces[eps]
        for simplex,w in geometric.items():
            if len(simplex)==1:actual[simplex]=w
        assert actual==geometric
        rb=reduction(actual);ob=reduction(witness_complex(pts,{},eps,maxdim=2));windows.append(dict(eps=eps,active_labels=sum(len(s)==1 for s in actual),CPP_simplex_weights_exact=True,quotient_H1=reviewed.serial(Counter({(Decimal(format(b,'.17g')),Decimal(format(d,'.17g'))):n for (b,d),n in rb[1].items()})),original_H1=reviewed.serial(Counter({(Decimal(format(b,'.17g')),Decimal(format(d,'.17g'))):n for (b,d),n in ob[1].items()}))))
    diameter=[max(np.linalg.norm(pts[a]-pts[b]) for a,b in __import__('itertools').combinations(vs,2)) for vs in pvs.values()]
    assert len(pvs)==2 and all(len(v)==4 for v in pvs.values()) and all(d<.85*.5 for d in diameter)
    summary=dict(status='UPSTREAM_H1_COUNTEREXAMPLE_CONFIRMED',input_sha256=sha(p),upstream='734be3dd396f71217b7693b6acb382b81d7db6a9',max_simplex=2,homology='H0/H1 configuration; H1 loses future interval',PV_representatives=pvs,PV_diameters=diameter,cap=.425,windows=windows,trace_path=trace['path'],stream_noquotient_path=inc['path'],stream_matches_Ripser=True,limits='Deterministic clustered S2 n8 counterexample; not a claim that all uniform TableI cases or all old failures share this cause.')
    write(RUN/'counterexample8-summary.json',summary);table(RUN/'counterexample8.csv',[{k:r[k] for k in ['version','threads','schedule','baseline_matches_Ripser','quotient_equal_baseline']} for r in rows])
    assert all(not r['quotient_equal_baseline'] for r in rows if r['schedule'] in ('single','grid01'))
    write(RUN/'G-GRID01-GATE.json',dict(status='NOT_RUN_FOUNDATION_PH_MISMATCH',at=now(),reason='New positive-distance unitS2 n8 fails H1 at max-simplex2 on untouched upstream and all group versions for single quotient and grid0.1. User requires related PH foundation checked before a new controller. No new G model fitted, no final dataset used. Continue with explicitly no-quotient checkpoint-cost survey and per-case validated n200/n300 resource study.',evidence='counterexample8-summary.json',reused_old_G='Only existing-grid instrumentation/replay diagnostic, not new-grid evaluation'))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':run()
