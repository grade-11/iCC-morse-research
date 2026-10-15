#!/usr/bin/env python3
"""Read-only brute-force analysis of at most two NEW reduced traces."""
from oracle import *

def analyze(run):
    results=[]
    for case in json.loads((run/'new-shrink-summary.json').read_text()):
        data=run/'new-shrunk'/case['case']/'retained.csv';pts=np.loadtxt(data,delimiter=',');text=(Path(case['trace_path'])/'stdout.log').read_text()
        states={};traces={};bars={};block=None
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
        rows=[];pvs={};lo=0.
        for eps in case['grid']:
            geometric=witness_complex(pts,pvs,eps,maxdim=case['dimension']);actual=traces[float(eps)]
            for s,w in geometric.items():
                if len(s)==1:actual[s]=w
            mismatches=[dict(simplex=s,cpp=actual.get(s),joint_witness=geometric.get(s)) for s in sorted(actual.keys()|geometric.keys()) if actual.get(s)!=geometric.get(s)]
            need(not mismatches,'New witness weights differ from oracle')
            qb=reduction(actual);ob=reduction(witness_complex(pts,{},eps,maxdim=case['dimension']));dims={}
            for d in range(1,case['dimension']):
                filt=lambda c:Counter({p:n for p,n in c.items() if p[1]==-1 or p[1]>=lo})
                match=bars[float(eps),d]==filt(qb[d]);need(match,'New C++ PH differs from boundary reduction')
                dims[str(d)]=dict(cpp_matches_boundary_reduction=match,missing=serial(filt(ob[d])-filt(qb[d])),extra=serial(filt(qb[d])-filt(ob[d])))
            rows.append(dict(eps_lo=lo,eps_hi=eps,active_labels=sum(len(s)==1 for s in actual),pvs=pvs,weight_mismatches=mismatches,dimensions=dims))
            pvs=states.get(float(eps),pvs);lo=eps
        results.append(dict(case=case['case'],retained_n=len(pts),windows=rows,scope='This reduced input only; not proof that full source dataset has the same cause'))
    write(run/'new-root-cause-oracle.json',results);print(json.dumps(results,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();analyze(a.run)
