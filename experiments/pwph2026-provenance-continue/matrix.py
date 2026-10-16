#!/usr/bin/env python3
import shutil,re
from workflow import *
sys.path.insert(0,str(REPO/'experiments/pwph2026-stage2-pilot/v2'))
import lab as reviewed
from collections import Counter
from decimal import Decimal

def setup():
    if (RUN/'PLAN.json').exists():return json.loads((RUN/'PLAN.json').read_text())
    old=REPO/'experiments/pwph2026-stage2-pilot/v2/runs/v2-20261015T195640Z';pilotold=REPO/'experiments/pwph2026-stage2-pilot/runs/stage2-20261015T190630Z'
    rows=[]
    def add(id,p,grid,small=False,known='PH_MISMATCH'):
        dest=RUN/'inputs'/(id+'.csv');dest.parent.mkdir(exist_ok=True);shutil.copy2(p,dest)
        rows.append(dict(id=id,file=str(dest),source=str(p),sha256=sha(dest),breaks=grid,small=small,prior_status=known))
    add('old11',old/'data/20261007-f15.csv',[1.2,1.5,2.],True)
    add('old10',old/'data/20261008-f4.csv',[1.2,1.4,1.7,2.],True)
    add('new11',old/'new-shrunk/uniform_s2-n64-s20261013/retained.csv',[.9,1.3,1.7,2.],True)
    add('new14',old/'new-shrunk/polar_s2-n100-s20261018/retained.csv',[1.2,1.7,2.],True)
    add('pass32',old/'data/uniform_s2-n32-s20261012.csv',[1.2,1.7,2.],True,'PASS')
    add('irregular5',REPO/'tests/data/irregular5.csv',[1.,1.5,2.],True,'Baseline CTest PASS; quotient tested here')
    add('standard100',pilotold/'data/S2_R3_n100_seed20261007.csv',[1.2,1.7,2.],False,'PASS')
    for seed in (20261007,20261008):
        for name,grid in [('f15',[1.2,1.5,2.]),('f4',[1.2,1.4,1.7,2.]),('f5',[1.2,1.4,1.6,1.8,2.])]:
            add(f'full-{seed}-{name}',pilotold/'data'/f'S2_R3_n100_seed{seed}.csv',grid)
    add('new-full64',old/'data/uniform_s2-n64-s20261013.csv',[.9,1.3,1.7,2.])
    add('new-full100',old/'data/polar_s2-n100-s20261018.csv',[1.2,1.7,2.])
    # Cache only references with the original complete config/version/checksum.
    shutil.copytree(old/'references',RUN/'references',dirs_exist_ok=True)
    plan=dict(created_at=now(),cases=rows,field=2,metric='euclidean',cutoff=2.,alpha=.85,
       dimensions='max-simplex2/3; output dimensions2/3 mean H1/H2',thread_rule='small1and4, full4',
       versions=['U','U17','B17','C','D'],reference_cache_source=str(old/'references'),
       U_precision='Unmodified output6sig: only internal printed precision comparison; supplement U17 output-only patch for exact Decimal/Ripser-float32 checks')
    write(RUN/'PLAN.json',plan);return plan

def exact_bounds(text,grid):
    if grid is None:return text
    bounds=[0.]+list(grid);formatted=[format(float(x),'.17g') for x in bounds];lines=text.splitlines();index=0
    for i,line in enumerate(lines):
        if line=='window lower bound, window upper bound':
            lo,hi=map(float,lines[i+1].split(','));expected=bounds[index:index+2]
            assert all(x in (y,float(format(y,'.6g'))) for x,y in zip((lo,hi),expected)),'Unexpected raw window metadata'
            lines[i+1]=formatted[index]+', '+formatted[index+1];index+=1
        elif line.startswith('  epsilon breaks: '):lines[i]='  epsilon breaks: '+' '.join(formatted[1:])
        elif line.startswith('in eps range '):
            match=re.match(r'in eps range (\S+)\s+(\S+)(\s+dimension = .*)',line)
            if match:
                found=next(j for j in range(len(grid)) if float(match[1]) in (bounds[j],float(format(bounds[j],'.6g'))) and float(match[2]) in (bounds[j+1],float(format(bounds[j+1],'.6g'))))
                lines[i]='in eps range '+formatted[found]+'  '+formatted[found+1]+match[3]
    return '\n'.join(lines)+'\n'

def parse_output(out,data,grid,d,threads,verbose=False,lossy=False):
    text=(out/'stdout.log').read_text()
    if not lossy:text=exact_bounds(text,grid)
    if not verbose:return reviewed.compact(text,'ph' if grid is None else 'piecewise',d,2)
    p=out/'parser-input.log';p.write_text(text+'\niCC exit status: 0\n')
    meta=dict(dataset=str(data),dim=d,threads=threads,breaks=[Decimal(str(x) if lossy else format(float(x),'.17g')) for x in (grid or [2.])],alpha=Decimal('.85') if grid else None)
    bars,ignored=reviewed.verbose_parser.parse_log(p,meta,'piecewise' if grid else 'ph')
    bound=lambda x: Decimal(str(x) if lossy else format(float(x),'.17g'))
    return bars,ignored,[(bound(a),bound(b)) for a,b in zip([0]+(grid or [2.])[:-1],grid or [2.])]

def one(version,item,d,threads,mode='quotient',extra=(),verbose=False,tag=''):
    versions=json.loads((RUN/'versions.json').read_text());binary=Path(versions[version]['binary']);data=Path(item['file']);grid=None if mode=='ph' else item['breaks'];key=f'{version}-{item["id"]}-d{d}-t{threads}-{mode}{tag}'
    argv=reviewed.command(binary,data,grid,threads=threads,maxdim=d,verbose=verbose,extra=extra)
    model_inputs=[Path(extra[extra.index('--g-model')+1])] if '--g-model' in extra else []
    out,c=execute(key,argv,[binary,data,BASE/'matrix.py',*model_inputs]);r=dict(version=version,case=item['id'],dimension=d,threads=threads,mode=mode,extra=list(extra),verbose=verbose,path=str(out),status=c['status'],exit_status=c['exit_status'],elapsed_s=c['elapsed_monotonic_s'],rss_kib=int(c.get('resources',{}).get('max_rss_kib',-1)),input_sha256=sha(data),binary_sha256=sha(binary))
    if c['status']=='OK':
        bars,ignored,windows=parse_output(out,data,grid,d,threads,verbose,lossy=version=='U')
        r.update(diagrams={str(k):reviewed.serial(v) for k,v in bars.items()},ignored_temporary_survivors=ignored,windows=[[str(a),str(b)] for a,b in windows]);r['H0']='NOT_EXPORTED' if version in ('U','U17','B17') else 'DIAGNOSTIC_NOT_REQUESTED'
    return r

def counters(r):return {int(k):Counter({(Decimal(x['birth']),Decimal(x['death'])):x['multiplicity'] for x in v}) for k,v in r['diagrams'].items()}

def matrix():
    plan=setup();rows=[]
    for version in plan['versions']:
        for item in plan['cases']:
            for d in ((2,3) if item['small'] or item['id']=='standard100' else (3,)):
                refs=reviewed.reference(RUN,Path(item['file']),d)
                for t in ((1,4) if item['small'] else (4,)):
                    base=one(version,item,d,t,'ph');rows.append(base)
                    q=one(version,item,d,t);rows.append(q)
                    if base['status']=='OK' and q['status']=='OK':
                        b=counters(base);x=counters(q);exact={str(k):b[k]==x[k] for k in b};different_counts={str(k):sum(b[k].values())!=sum(x[k].values()) for k in b}
                        base['status']='PRINTED_PRECISION_ONLY' if version=='U' else ('PASS' if all(reviewed.projected(b[k])==refs[k] for k in b) else 'PH_MISMATCH')
                        q.update(exact_CPP=exact,bar_count_diff=different_counts,baseline_path=base['path'])
                        if version=='U':q['status']='PH_MISMATCH' if any(different_counts.values()) else ('PRINTED_PRECISION_ONLY' if all(exact.values()) else 'OUTPUT_DIFFERENCE_NEEDS17')
                        else:q['status']='PASS' if all(exact.values()) and base['status']=='PASS' and all(reviewed.projected(x[k])==refs[k] for k in x) else 'PH_MISMATCH'
                        q['difference']={str(k):dict(missing=reviewed.serial(b[k]-x[k]),extra=reviewed.serial(x[k]-b[k])) for k in b}
                    for r in (base,q):write(Path(r['path'])/'ph-check.json',r)
                    write(RUN/'correctness.json',rows);progress(q['case']+'-'+version+f'-d{d}-t{t}',q['status'],'Exact CPP17, Ripserfloat32; U untouched output6sig explicitly qualified')
    table(RUN/'correctness.csv',[{k:r.get(k,'') for k in ('version','case','dimension','threads','mode','status','H0','elapsed_s','rss_kib','path')} for r in rows])
    print(json.dumps(dict(configurations=len(rows),statuses=dict(Counter(r['status'] for r in rows))),indent=2))

if __name__=='__main__':matrix()
