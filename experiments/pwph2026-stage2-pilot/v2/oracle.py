#!/usr/bin/env python3
"""Small-case F2 boundary reduction and brute-force joint-witness oracle.
Diagnostic only; not a substitute algorithm or an H0 implementation proof.
"""
from itertools import combinations,product
from collections import Counter
import re,math,json,argparse
from pathlib import Path
import numpy as np
from lab import *

def reduction(simplices,maxdim=3):
    cells=sorted(simplices,key=lambda s:(simplices[s],len(s),s));rank={s:i for i,s in enumerate(cells)}
    pivots={};births=set();paired=set();bars={d:Counter() for d in range(maxdim)}
    for j,s in enumerate(cells):
        col={rank[s[:k]+s[k+1:]] for k in range(len(s))} if len(s)>1 else set()
        while col and max(col) in pivots:col^=pivots[max(col)]
        if col:
            i=max(col);pivots[i]=col;paired.add(i);b=simplices[cells[i]];death=simplices[s];d=len(cells[i])-1
            if d<maxdim and b!=death:bars[d][b,death]+=1
        else:births.add(j)
    for i in births-paired:
        d=len(cells[i])-1
        if d<maxdim:bars[d][simplices[cells[i]],-1.]+=1
    return bars

def witness_complex(points,pvs,eps,maxdim=3):
    n=len(points);absorbed={v for vs in pvs.values() for v in vs}
    labels={i:[i] for i in range(n) if i not in absorbed};labels.update(pvs)
    dist={}
    for i,j in combinations(range(n),2):
        total=0.
        for a,b in zip(points[i],points[j]):delta=float(a)-float(b);total+=delta*delta
        dist[i,j]=dist[j,i]=math.sqrt(total)
    cells={(i,):0. for i in labels}
    for k in range(2,maxdim+2):
        for s in combinations(sorted(labels),k):
            weight=min(max(dist[a,b] for a,b in combinations(w,2)) for w in product(*(labels[i] for i in s)))
            if weight<eps:cells[s]=weight
    return cells

def investigate(run,binary):
    prior=BASE.parent/'runs/stage2-20261015T190630Z';results=[]
    for name,breaks in [('20261007-f15',[1.2,1.5,2]),('20261008-f4',[1.2,1.4,1.7,2])]:
        data=run/'data'/(name+'.csv');pts=np.loadtxt(data,delimiter=',')
        for t in (1,4):
            out,c=execute(run,f'trace-{name}-t{t}',command(binary,data,breaks,threads=t,verbose=True,extra=['--trace-pwph']),[binary,data,Path(__file__)])
            need(c['status']=='OK','Trace failed');text=(out/'stdout.log').read_text();states={};traces={};block=None;bars={}
            for line in text.splitlines():
                if line.startswith('[pv-members]'):
                    f=dict(x.split('=') for x in line.split()[1:]);states.setdefault(float(f['eps_hi']),{})[int(f['label'])]=list(map(int,f['representatives'].strip(',').split(',')))
                elif line.startswith('[simplex]'):
                    f=dict(x.split('=') for x in line.split()[1:]);eps=float(f['eps_hi']);s=tuple(sorted(map(int,f['labels'].strip(',').split(','))));traces.setdefault(eps,{})[s]=float(f['weight'])
                elif line.startswith('in eps range '):
                    m=re.search(r'in eps range (\S+)\s+(\S+)\s+dimension = (\d+)',line);block=float(m[2]),int(m[3])-1;bars.setdefault(block,Counter())
                elif block:
                    m=re.fullmatch(r'\s*\(([^,]+),\s*([^,]+)\)',line)
                    if m:bars[block][tuple(map(float,m.groups()))]+=1
            windows=[];pvs={};lo=0.
            for eps in breaks:
                geometric=witness_complex(pts,pvs,eps);actual=traces[float(eps)]
                # Add vertices, which are not emitted by C++ trace.
                for s,w in geometric.items():
                    if len(s)==1:actual[s]=w
                mismatch=[dict(simplex=s,cpp=actual.get(s),joint_witness=geometric.get(s)) for s in sorted(actual.keys()|geometric.keys()) if actual.get(s)!=geometric.get(s)]
                need(not mismatch,'Witness enumerator differs from brute force; inspect saved trace')
                qbars=reduction(actual);obars=reduction(witness_complex(pts,{},eps))
                row=dict(eps_hi=eps,eps_lo=lo,pvs=pvs,geometric_weight_mismatches=mismatch,dimensions={})
                for d in (1,2):
                    # Raw C++ frame emits only deaths in current window (plus live bars).
                    filt=lambda c:Counter({p:n for p,n in c.items() if p[1]==-1 or p[1]>=lo})
                    row['dimensions'][d]=dict(cpp_matches_boundary_reduction=bars[float(eps),d]==filt(qbars[d]),
                        quotient_vs_original=dict(missing=serial(filt(obars[d])-filt(qbars[d])),extra=serial(filt(qbars[d])-filt(obars[d]))))
                windows.append(row);pvs=states.get(float(eps),pvs);lo=eps
            result=dict(case=name,threads=t,binary_sha256=sha(binary),windows=windows,
                        trace_path=str(out),method='Brute joint-witness enumeration; explicit F2 boundary reduction on small traced complex')
            results.append(result);write(out/'oracle.json',result)
    write(run/'root-cause-oracle.json',results);print(json.dumps(results,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);a=ap.parse_args();investigate(a.run,a.binary.resolve())
