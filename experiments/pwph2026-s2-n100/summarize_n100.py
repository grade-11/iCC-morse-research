#!/usr/bin/env python3
"""Read existing n100 evidence and write CSV tables. Never run PH/benchmarks."""
import argparse
import csv
import json
from pathlib import Path
import re
import statistics
import tarfile
from evidence_run import sha, write
import benchmark_icc as benchmark
import compare_icc_ph as verbose

BASE=Path(__file__).resolve().parent


def need(ok, why):
    if not ok: raise ValueError(why)


def csvout(path, rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(rows)


def verify_tree(archive, tree):
    evidence={}
    with tarfile.open(archive) as t:
        for m in t:
            if m.isfile():
                parts=Path(m.name).parts
                need('..' not in parts and not Path(m.name).is_absolute(),'Unsafe source path')
                rel=Path(*parts[1:]); p=tree/rel
                import hashlib
                h=hashlib.sha256(t.extractfile(m).read()).hexdigest()
                need(p.is_file() and sha(p)==h,'Folder differs from archive: '+str(p))
                evidence[str(p)]=h
    return evidence


def windows(text, run, source):
    rows=[]
    blocks=re.split(r'(?=^\[rebuild-timing\] stage=edges )',text,flags=re.M)[1:]
    need(len(blocks)==3,'Missing verbose windows')
    active,pvtotal=100,0
    timings=[]
    for block in blocks:
        header=re.search(r'eps_lo=(\S+) eps_hi=(\S+)',block)
        lo,hi=map(float,header.groups()); final=hi==2.0
        row=dict(source=source,run=run,warmup=run.endswith('00'),eps_lo=lo,eps_hi=hi,active_in=active)
        cofacets=re.findall(r'in eps range .*?dimension = (\d+)\s+cofacet num = (\d+)\s+facet num = (\d+)',block)
        need([int(d) for d,c,f in cofacets]==[2,3],'Missing cofacets')
        for dim,count,facets in cofacets:
            row['cofacets_dim'+dim]=int(count); row['facets_dim'+dim]=int(facets)
        if not final:
            raw=int(re.search(r'^raw pv support cofacets size = (\d+)',block,re.M)[1])
            clean=int(re.search(r'^pv support label sets with no pv/protected contents size = (\d+)',block,re.M)[1])
            overlap=int(re.search(r'^supports rejected by overlap with newly selected PVs = (\d+)',block,re.M)[1])
            change=re.search(r'new pv number = (\d+)\s+total pv number = (\d+)\s+active labels before quotient = (\d+)\s+active labels after quotient = (\d+)',block)
            new,pvtotal,before,after=map(int,change.groups())
            need(before==active,'Active labels discontinuity'); active=after
            row.update(active_before=before,active_after=after,raw_supports=raw,label_clean=clean,
                       rejected_overlap=overlap,new_pvs=new,total_pvs=pvtotal)
            diag=re.search(r'^\[pv-filters\] (.*)$',block,re.M)
            for k in ('rejected_pv','rejected_protected','rejected_cardinality','rejected_diameter','rejected_separation','eligible_at_visit'):
                row[k]=int(dict(x.split('=') for x in diag[1].split())[k]) if diag else 'NOT_MEASURED'
        else:
            row.update({k:'N/A' for k in ('active_before','active_after','raw_supports','label_clean',
                'rejected_overlap','rejected_pv','rejected_protected','rejected_cardinality',
                'rejected_diameter','rejected_separation','eligible_at_visit')})
            row.update(new_pvs=0,total_pvs=pvtotal)
        row.update(matching_ms=0.0,rebuild_ms=0.0,interface_ms=0.0,
                   witness_search_worker_ms=0.0,witness_facet_prepare_worker_ms=0.0,
                   witness_covt_prepare_worker_ms=0.0,witness_search_calls=0)
        for line in block.splitlines():
            tag=re.match(r'^\[([^]]+)\] (.+)',line)
            if not tag or tag[1] not in ('rebuild-timing','matching-timing','interface-timing',
                                       'witness-timing','pv-timing','window-resources'): continue
            group=tag[1]; fields=dict(x.split('=',1) for x in tag[2].split())
            for k,v in fields.items():
                if k.endswith(('_ms','_kib')) or k=='search_calls':
                    timings.append(dict(source=source,run=run,warmup=run.endswith('00'),eps_lo=lo,eps_hi=hi,
                        group=group,stage=fields.get('stage',''),dim=fields.get('dim',''),metric=k,value=float(v),
                        unit='KiB' if k.endswith('_kib') else 'calls' if k=='search_calls' else 'ms'))
            if group=='matching-timing': row['matching_ms']+=float(fields['matching_ms'])
            elif group=='rebuild-timing':
                row['rebuild_ms']+=sum(float(v) for k,v in fields.items() if k.endswith('_ms'))
            elif group=='interface-timing': row['interface_ms']+=float(fields['elapsed_ms'])
            elif group=='witness-timing':
                for k,target in [('search_worker_ms','witness_search_worker_ms'),
                   ('facet_prepare_worker_ms','witness_facet_prepare_worker_ms'),
                   ('covt_prepare_worker_ms','witness_covt_prepare_worker_ms'),('search_calls','witness_search_calls')]:
                    row[target]+=float(fields[k])
            elif group=='pv-timing':
                row['pv_selection_ms']=float(fields['selection_ms']); row['pv_commit_ms']=float(fields['commit_ms'])
            elif group=='window-resources':
                need(fields['memory_status']=='ok','Window RSS unavailable')
                row.update(elapsed_ms=float(fields['elapsed_ms']),rss_start_kib=int(fields['rss_start_kib']),
                    rss_end_kib=int(fields['rss_end_kib']),rss_peak_kernel_kib=int(fields['rss_peak_kernel_kib']))
        need('elapsed_ms' in row,'Missing window duration')
        if final: row.update(pv_selection_ms='N/A',pv_commit_ms='N/A')
        rows.append(row)
    return rows,timings


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dossier',required=True,type=Path)
    a=ap.parse_args(); out=a.dossier/'tables'; out.mkdir(exist_ok=False)
    batch=BASE/'runs/whole-run-monotonic-ma523yck'
    diagnostic=BASE/'runs/repeat-verbose-TIVmZd'
    verified={}
    for tree in (batch,diagnostic):
        verified.update(verify_tree(tree.with_suffix('.tar.gz'),tree))
    write(out/'source-files-sha256.json',verified)
    whole=benchmark.analyze(batch)
    need(whole['status']=='PASS','Old benchmark PH failed')
    write(out/'whole-run-analysis.json',whole)
    pairs=[]
    for row in whole['rows']:
        x=dict(round=row['round'],warmup=row['warmup'])
        for mode in ('ph','piecewise'):
            for k,v in row[mode].items():
                if isinstance(v,(int,float)): x[mode+'_'+k]=v
        x['paired_time_ratio_ph_over_piecewise']=row['ph']['elapsed_monotonic_s']/row['piecewise']['elapsed_monotonic_s']
        x['paired_peak_rss_ratio_ph_over_piecewise']=row['ph']['max_rss_kib']/row['piecewise']['max_rss_kib']
        pairs.append(x)
    csvout(out/'whole-run-paired.csv',pairs)
    summary=[]
    for mode,metrics in whole['summary'].items():
        for metric,st in metrics.items(): summary.append(dict(mode=mode,metric=metric,**st))
    csvout(out/'whole-run-summary.csv',summary)
    oldrows=[]; timings=[]
    baseline=BASE/'runs/comparison-baseline-rmfn6u'
    for run in sorted(diagnostic.glob('run-*')):
        ctx=(run/'context.txt').read_text()
        need('b1545e28453677d06e68280f311bcdf9ef0a40d9' in ctx and
             '04b34fef3cd88efe4415603411a569b9fc7ce918f3307c8cb0a47e72a20f0af9' in ctx,
             'Wrong diagnostic provenance')
        need(json.loads((run/'report.json').read_text())['status']=='PASS','Old report failed')
        meta=verbose.context(run/'context.txt','piecewise')
        ds,ignored=verbose.parse_log(run/'piecewise.log',meta,'piecewise')
        need([sum(ds[d].values()) for d in (1,2)]==[26,1],'Diagnostic pair counts differ')
        rows,long=windows((run/'piecewise.log').read_text(),run.name,'old-verbose-b1545e2')
        oldrows+=rows; timings+=long
    need(len(oldrows)==33,'Need 11 verbose rounds')
    csvout(out/'windows-old-verbose.csv',oldrows)
    newrows,newlong=windows((a.dossier/'validation-v2/piecewise-verbose/stdout.log').read_text(),'validation-once','new-validation-only')
    csvout(out/'windows-new-validation.csv',newrows)
    csvout(out/'timings-long.csv',timings+newlong)
    agg=[]
    metrics=('elapsed_ms','rss_peak_kernel_kib','matching_ms','rebuild_ms','interface_ms',
             'pv_selection_ms','pv_commit_ms','witness_search_worker_ms',
             'witness_facet_prepare_worker_ms','witness_covt_prepare_worker_ms')
    for hi in (1.2,1.7,2.0):
        subset=[r for r in oldrows if r['eps_hi']==hi and not r['warmup']]
        need(len(subset)==10,'Need 10 measured verbose rounds')
        for metric in metrics:
            vals=[r[metric] for r in subset if isinstance(r[metric],(int,float))]
            agg.append(dict(eps_hi=hi,metric=metric,n=len(vals),
                median=statistics.median(vals) if vals else 'N/A',
                minimum=min(vals) if vals else 'N/A',maximum=max(vals) if vals else 'N/A'))
    csvout(out/'windows-old-summary.csv',agg)
    report=dict(status='PASS',whole_run=whole['summary'],ratios=whole['ratios'],clock_warnings=whole['clock_warnings'],
        raw_archives={str(p):sha(p) for p in (batch.with_suffix('.tar.gz'),diagnostic.with_suffix('.tar.gz'))},
        old_head=whole['head'],old_binary_sha256=whole['checksums']['./build/iCC'],
        diagnostic_rows=len(oldrows),measured_verbose_rounds=10,new_validation_rows=len(newrows),
        notes=['All 11 raw rounds retained; round00 excluded from summaries',
        'Old verbose and nonverbose benchmark are separate workloads; new validation is not a benchmark',
        'Worker elapsed sums are neither CPU time nor wall time; do not add them to enclosing cofacet rebuild timers',
        'Rebuild sum covers separate edges/cofacets/active-facet stages; cofacets include witness work',
        'Interface sum includes FFI diagnostics; window elapsed includes stdout and other uninstrumented work',
        'Window RSS resets /proc/self/clear_refs=5 and measures VmHWM; not whole-run peak',
        'First-window witness worker fields 0 mean no PV witness calls; final supports/filters N/A',
        'Historical downstream rejection counters NOT_MEASURED; measured separately in new validation'])
    write(out/'report.json',report); print(json.dumps(report,indent=2))


if __name__=='__main__': main()
