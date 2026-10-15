#!/usr/bin/env python3
"""Read-only S2/n100 archive comparison. Requires numpy and the reviewed benchmark_icc.py.
Exit 0: H1/H2 match after float32 projection; 1: mismatch; 2: invalid evidence.
No iCC/Ripser execution. H0 is not compared. JSON is written to stdout.
"""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile
import sys
import numpy as np

PARSER_SHA = '730a2aea8b946d2d04995e189647188115146bdcbf78f4024e8e6d219d275e72'

def need(ok, message):
    if not ok:
        raise ValueError(message)

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def unpack(path, dest):
    seen = set()
    total = 0
    with tarfile.open(path) as archive:
        for m in archive:
            p = PurePosixPath(m.name)
            need(not p.is_absolute() and '..' not in p.parts and len(p.parts)>0, 'Unsafe archive path')
            need(m.name not in seen, 'Duplicate archive member')
            seen.add(m.name)
            need(m.isdir() or m.isfile(), 'Links/special archive members unsupported')
            if m.isfile():
                total += m.size
                need(total <= 100_000_000, 'Archive exceeds 100 MB limit')
                target = dest.joinpath(*p.parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive.extractfile(m).read())
    roots = list(dest.iterdir())
    need(len(roots)==1 and roots[0].is_dir(), 'Expected one archive root')
    return roots[0]

def entries(c):
    return [dict(birth=str(b),death=str(d),multiplicity=n) for (b,d),n in sorted(c.items())]

def compare(args):
    need(sha(args.benchmark_script)==PARSER_SHA, 'Unreviewed benchmark parser; checksum differs')
    spec=importlib.util.spec_from_file_location('reviewed_benchmark',args.benchmark_script)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    with tempfile.TemporaryDirectory(prefix='icc-ripser-review-') as tmp:
        base=Path(tmp)
        batch=unpack(args.icc_archive,base/'icc')
        ref=unpack(args.ripser_archive,base/'ripser')
        manifest=json.loads((ref/'checksums.json').read_text())
        need(set(manifest)=={p.name for p in ref.iterdir() if p.is_file()}-{'checksums.json'},'Incomplete reference manifest')
        for n,h in manifest.items():
            need(Path(n).name==n and sha(ref/n)==h,'Reference checksum mismatch: '+n)
        ctx=json.loads((ref/'context.json').read_text())
        need(ctx['status']=='REFERENCE_CREATED' and ctx['exit_status']==0,'Reference unsuccessful')
        need(ctx['parameters']==dict(maxdim=2,thresh=2.0,coeff=2,metric='euclidean',distance_matrix=False,do_cocycles=False,n_perm=None),'Incompatible Ripser parameters')
        need(ctx['packages']['ripser']=='0.6.15','Unexpected Ripser version')
        need(sha(ref/'input.csv')==ctx['dataset_sha256']==ctx['dataset_sha256_after']==mod.DATA_SHA,'Dataset mismatch')
        need(sha(ref/'runner.py')==ctx['script_sha256'],'Runner provenance mismatch')
        benchmark=mod.analyze(batch)
        rows=list(csv.DictReader((ref/'intervals.csv').open()))
        need(all(r['homology_dimension'] in ('0','1','2') for r in rows),'Unexpected CSV dimension')
        reference={}
        with np.load(ref/'diagrams.npz',allow_pickle=False) as z:
            need(set(z.files)=={'H0','H1','H2'},'Unexpected NPZ dimensions')
            for dim in range(3):
                a=z[f'H{dim}']
                need(a.ndim==2 and a.shape[1]==2,'Invalid diagram shape')
                need(np.all(np.isfinite(a[:,0]) & (a[:,0]>=0) & (a[:,0]<=2)),'Invalid births')
                need(np.all((np.isfinite(a[:,1]) & (a[:,1]>=a[:,0]) & (a[:,1]<=2)) | np.isposinf(a[:,1])),'Invalid deaths')
                need(np.array_equal(a,a.astype(np.float32).astype(np.float64)),'Reference endpoints not float32 representable')
                reference[dim]=Counter(map(tuple,a.tolist()))
                c=Counter((float(r['birth']),float(r['death'])) for r in rows if int(r['homology_dimension'])==dim)
                need(c==reference[dim],'CSV/NPZ mismatch')
        results=[]
        for folder in sorted(batch.glob('round-*')):
            for mode in ('ph','piecewise'):
                diagrams,ignored,_=mod.parse_log(folder/mode/'program.log',mode)
                dimensions={}
                for dim in (1,2):
                    projected=Counter()
                    for (b,d),n in diagrams[dim].items():
                        projected[(float(np.float32(float(b))),float('inf') if d==-1 else float(np.float32(float(d))))]+=n
                    dimensions[f'H{dim}']=dict(count=sum(projected.values()),match=projected==reference[dim],missing=entries(reference[dim]-projected),extra=entries(projected-reference[dim]))
                results.append(dict(round=folder.name,mode=mode,ignored_intermediate_survivors=ignored,dimensions=dimensions))
        passed=benchmark['status']=='PASS' and all(d['match'] for r in results for d in r['dimensions'].values())
        return dict(status='PASS' if passed else 'FAIL',comparison='H1/H2 endpoint multisets after iCC float32 projection; no tolerance',H0='NOT_CHECKED',notes=['Dataset/binary provenance for iCC comes from recorded context; original binary is not included','This is not exact float64 equality, algorithmic proof, or full Stage 1 completion','Intermediate piecewise death=-1 excluded; final survivors mapped to inf at cutoff 2'],input_sha256={k:sha(v) for k,v in vars(args).items()},icc_head=benchmark['head'],ripser_head=ctx['head'],dataset_sha256=mod.DATA_SHA,icc_binary_sha256=benchmark['checksums']['./build/iCC'],icc_exact_decimal_status=benchmark['status'],icc_evidence_sha256=benchmark['input_sha256'],reference_counts={f'H{d}':sum(c.values()) for d,c in reference.items()},runs=results)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--icc-archive',type=Path,required=True)
    parser.add_argument('--ripser-archive',type=Path,required=True)
    parser.add_argument('--benchmark-script',type=Path,required=True)
    args=parser.parse_args()
    try:
        report=compare(args)
    except Exception as e:
        print(json.dumps(dict(status='ERROR',message=str(e)),indent=2))
        return 2
    print(json.dumps(report,indent=2))
    return 0 if report['status']=='PASS' else 1

if __name__=='__main__':
    sys.exit(main())
