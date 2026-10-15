#!/usr/bin/env python3
"""Validate one NEW binary sequentially; not a repeated performance benchmark.
Captures n100 verbose/nonverbose modes, known H0 fixtures and raw Ripser fixtures.
H0 diagnostics are the complement of iCC's actual MSF clearing mask; only the
first window (before quotient) is compared to the original n100 diagram.
"""
import argparse
from collections import Counter
from decimal import Decimal
import importlib.util
import json
from pathlib import Path
import re
import sys
import numpy as np
from ripser import ripser
from evidence_run import capture, sha, write
import compare_icc_ph as verbose
import benchmark_icc as compact

BASE = Path(__file__).resolve().parent
DATA = compact.DATA


def need(ok, why):
    if not ok:
        raise ValueError(why)


def meta(mode):
    return dict(dataset=DATA, dim=3, threads=4,
                breaks=[Decimal('2')] if mode=='ph' else [Decimal(x) for x in ('1.2','1.7','2')],
                alpha=Decimal('.85'))


def msf(text):
    frames=[]
    for line in text.splitlines():
        if line.startswith('[h0-msf] '):
            fields=dict(x.split('=') for x in line.split()[1:])
            frames.append(dict(metadata=fields, pairs=Counter()))
        elif line.startswith('[h0-msf-pair] '):
            fields=dict(x.split('=') for x in line.split()[1:])
            frames[-1]['pairs'][(Decimal(fields['birth']),Decimal(fields['death']))]+=1
    for f in frames:
        m, c = f['metadata'], f['pairs']
        finite=sum(n for (b,d),n in c.items() if d!=-1)
        alive=sum(n for (b,d),n in c.items() if d==-1)
        need(finite==int(m['merges']) and alive==int(m['components']), 'MSF count mismatch')
        need(finite+alive==int(m['vertices']), 'MSF vertex conservation failed')
    return frames


def projected(c):
    out=Counter()
    for (b,d),n in c.items():
        out[(float(np.float32(float(b))), float('inf') if d==-1 else float(np.float32(float(d))))]+=n
    return out


def serial(c):
    return [dict(birth=str(b),death=str(d),multiplicity=n) for (b,d),n in sorted(c.items())]


def run(binary, output):
    output.mkdir(parents=True, exist_ok=False)
    need(sha(DATA)==compact.DATA_SHA, 'Wrong n100 input')
    refdir=BASE/'runs/ripser-reference-noa_bg6w'
    manifest=json.loads((refdir/'checksums.json').read_text())
    need(all(sha(refdir/n)==h for n,h in manifest.items()), 'Reference manifest mismatch')
    with np.load(refdir/'diagrams.npz',allow_pickle=False) as z:
        refs={d:Counter(map(tuple,z[f'H{d}'].tolist())) for d in range(3)}
    reports={}
    diagrams={}
    texts={}
    for verbose_mode in (True,False):
        for mode in ('ph','piecewise'):
            name=mode+('-verbose' if verbose_mode else '-compact')
            cmd=compact.command(mode); cmd[0]=str(binary)
            if verbose_mode: cmd+=['--verbose']
            rc=capture(output/name,cmd,(binary,Path(DATA),Path(__file__),BASE/'evidence_run.py'),120)
            need(rc==0, f'{name} exit {rc}')
            # Derived parser input; raw stdout and actual exit code stay untouched.
            text=(output/name/'stdout.log').read_text()+f'\niCC exit status: {rc}\n'
            p=output/name/'parser-input.log'; p.write_text(text)
            if verbose_mode:
                diag,ignored=verbose.parse_log(p,meta(mode),mode)
            else:
                diag,ignored,_=compact.parse_log(p,mode)
            diagrams[name]=diag; texts[name]=text
            match={f'H{d}':projected(diag[d])==refs[d] for d in (1,2)}
            need(all(match.values()), f'{name} H1/H2 differs from Ripser')
            reports[name]=dict(exit_status=rc,counts={f'H{d}':sum(diag[d].values()) for d in (1,2)},
                               ripser_float32=match,ignored_intermediate_survivors=ignored)
    reference=diagrams['ph-verbose']
    need(all(d==reference for d in diagrams.values()), 'New modes differ in exact Decimal multisets')
    old=BASE/'runs/whole-run-monotonic-ma523yck/round-00/ph/program.log'
    old_diag,_,_=compact.parse_log(old,'ph')
    need(reference==old_diag,'H1/H2 changed relative to old binary')
    h0={m:msf(texts[m+'-verbose']) for m in ('ph','piecewise')}
    need(len(h0['ph'])==1 and len(h0['piecewise'])==3,'Missing H0 frames')
    first=h0['piecewise'][0]
    need(first['pairs']==h0['ph'][0]['pairs'],'Initial PwPH H0 differs from baseline')
    need(projected(first['pairs'])==refs[0], 'H0 differs from Ripser')
    need(all(f['metadata']['components']=='1' for m in h0.values() for f in m),
         'n100 connectivity criterion failed; do not stitch window-local H0')
    need(first['metadata']['vertices']=='100','Wrong original n100 vertex count')
    # n100 connectivity already resolves every finite H0 interval before the
    # first quotient. Later local bars are NOT concatenated into that diagram.
    filters=[]
    for line in texts['piecewise-verbose'].splitlines():
        if line.startswith('[pv-filters] '):
            row=dict(x.split('=') for x in line.split()[1:]); filters.append(row)
            v={k:int(x) for k,x in row.items() if k!='eps_hi'}
            need(v['raw']==v['rejected_pv']+v['rejected_protected']+v['label_clean'],'Label filter conservation')
            need(v['label_clean']==sum(v[k] for k in ('rejected_overlap','rejected_cardinality',
                 'rejected_diameter','rejected_separation','selected')), 'Sequential filter conservation')
            need(v['selected']==v['eligible_at_visit'], 'Eligibility mismatch')
    need(len(filters)==2,'Final window must not process supports')
    need([int(r['raw']) for r in filters]==[201,214] and
         [int(r['label_clean']) for r in filters]==[188,111] and
         [int(r['rejected_overlap']) for r in filters]==[178,106] and
         [int(r['selected']) for r in filters]==[10,5], 'PV counters changed')
    # Compare all deterministic old diagnostic lines: pairs/labels, cofacets,
    # support/PV statistics and FFI metrics. Exclude only clocks/resources and
    # newly added diagnostics, not invariant-bearing lines.
    oldtext=(BASE/'runs/repeat-verbose-TIVmZd/run-01/piecewise.log').read_text()
    def invariants(t):
        return [s.strip() for s in t.splitlines() if s.strip() and
                not s.startswith(('[h0-msf','[pv-filters]','[rebuild-timing]',
                    '[interface-timing]','[matching-timing]','[pv-timing]',
                    '[window-resources]','[witness-timing]','run time ='))]
    need(invariants(oldtext)==invariants(texts['piecewise-verbose']), 'Old/new verbose invariants changed')
    fixture_results=[]
    fixtures=[('multiplicity-disconnected', [[0,0],[1,0],[2,0],[10,0],[11,0]], 3.0,
               Counter({(Decimal(0),Decimal(1)):3,(Decimal(0),Decimal(-1)):2})),
              ('cutoff-boundary', [[0,0],[1,0],[3,0]], 1.0,
               Counter({(Decimal(0),Decimal(-1)):3})),
              ('coincident', [[0,0],[0,0],[1,0]], 2.0, None)]
    for name,pts,cutoff,expected in fixtures:
        data=output/(name+'.csv'); np.savetxt(data,pts,delimiter=',',fmt='%.17g')
        # Raw Ripser result: same coefficients, Euclidean metric, no permutation.
        rp=ripser(np.asarray(pts,dtype=float),maxdim=0,thresh=cutoff,coeff=2,
                  metric='euclidean',distance_matrix=False,do_cocycles=False,n_perm=None)['dgms'][0]
        np.savez_compressed(output/(name+'-ripser.npz'),H0=rp)
        ref=Counter(map(tuple,rp.tolist()))
        for mode in ('ph','piecewise'):
            cmd=[str(binary),'--tool',mode,'--file-name',str(data),'--max-dim','1',
                 '--threads','1','--verbose']
            cmd+=['--max-eps' if mode=='ph' else '--eps-breaks',str(cutoff)]
            if mode=='piecewise': cmd+=['--pv-cap-scale','0.85']
            rc=capture(output/(name+'-'+mode),cmd,(data,binary),30)
            need(rc==0,'Fixture command failed')
            frames=msf((output/(name+'-'+mode)/'stdout.log').read_text())
            need(len(frames)==1,'Fixture H0 missing')
            c=frames[0]['pairs']
            if expected is not None: need(c==expected,'Known H0 fixture mismatch')
            if name=='multiplicity-disconnected': need(projected(c)==ref,'Fixture Ripser mismatch')
            elif name=='cutoff-boundary':
                need(ref==Counter({(0.0,1.0):1,(0.0,float('inf')):2}), 'Unexpected Ripser cutoff convention')
                need(projected(c)!=ref,'Expected strict/inclusive cutoff difference missing')
            else:
                known = (Counter({(Decimal(0),Decimal(0)):1,(Decimal(0),Decimal(1)):1,
                                  (Decimal(0),Decimal(-1)):1}) if mode=='ph' else
                         Counter({(Decimal(0),Decimal(1)):2,(Decimal(0),Decimal(-1)):1}))
                need(c==known,'Coincident-point convention changed')
                need(ref==Counter({(0.0,1.0):1,(0.0,float('inf')):1}), 'Unexpected Ripser zero-length convention')
            fixture_results.append(dict(fixture=name,mode=mode,exit_status=rc,intervals=serial(c),
                ripser_intervals=[dict(birth=b,death=str(d),multiplicity=n) for (b,d),n in ref.items()],
                ripser_match=projected(c)==ref,
                interpretation='Expected cutoff/zero-edge convention difference' if name!='multiplicity-disconnected' else 'PASS'))
    report=dict(status='PASS',binary_sha256=sha(binary),dataset_sha256=sha(DATA),runs=reports,
        H1_H2_exact_new_vs_old='PASS',deterministic_piecewise_invariants='PASS',filters=filters,
        H0=dict(status='PASS_N100_INITIAL_MSF',scope='Original n100 H0 resolved before first quotient; window-local diagnostics only',
                baseline_piecewise_exact=True,ripser_float32=True,finite=99,surviving_at_cutoff=1,
                window_components={k:[int(f['metadata']['components']) for f in v] for k,v in h0.items()},
                original_intervals=serial(first['pairs']),
                limitations=['No general stitched H0 output through quotient',
                  'iCC uses distance < cutoff; Ripser includes equality',
                  'PwPH excludes zero-distance edges; Ripser omits zero-length H0 intervals']),fixtures=fixture_results)
    write(output/'report.json',report)
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    a=ap.parse_args()
    if a.output.exists():
        ap.error('Output already exists; choose a new evidence directory')
    try:
        report=run(a.binary.resolve(),a.output)
    except Exception as e:
        report=dict(status='ERROR',error=str(e))
        if a.output.is_dir(): write(a.output/'failure.json',report)
    print(json.dumps(report,indent=2))
    return 0 if report['status']=='PASS' else 2


if __name__=='__main__':
    sys.exit(main())
