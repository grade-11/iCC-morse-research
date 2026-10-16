#!/usr/bin/env python3
"""Audit persisted evidence without making new PH decisions or measurements."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def audit(run):
    plan=read(run/'PLAN.json');lock=read(run/'LOCK.json');checks=[]
    checks.append(dict(name='locked_model_and_code',PASS=sha(run/'model.json')==lock['model_sha256'] and all(sha(run/p)==h for p,h in lock['code'].items())))
    checks.append(dict(name='frozen_engine_source',PASS=all(sha(run/'engine-source'/p)==h for p,h in lock['engine']['source_files'].items())))
    checks.append(dict(name='measured_binary',PASS=sha(run/'build/iCC')==lock['engine']['sha256']))
    seeds=[x['seed'] for x in plan['datasets']];checks.append(dict(name='disjoint_seeds',PASS=len(seeds)==len(set(seeds)) and Counter(x['split'] for x in plan['datasets'])=={'dev':12,'selection':6,'final':6}))
    geometry=[]
    for item in plan['datasets']:
        p=run/'data'/(item['id']+'.csv');meta=read(p.with_suffix('.json'))
        pts=[list(map(float,l.split(','))) for l in p.read_text().splitlines()]
        distances=[math.sqrt(sum((a-b)**2 for a,b in zip(v,w))) for i,v in enumerate(pts) for w in pts[:i]]
        geometry.append(dict(dataset=item['id'],sha256=sha(p),matches_metadata=sha(p)==meta['sha256'],n=len(pts),
             minimum_distance=min(distances),maximum_distance=max(distances),exact_cutoff_equalities=sum(d==2 for d in distances),
             max_norm_squared_error=max(abs(sum(x*x for x in v)-1) for v in pts),
             generated_at=meta['created_at'],after_lock=meta['created_at']>lock['locked_at'] if item['split']=='final' else None))
    checks.append(dict(name='dataset_integrity_unit_sphere_positive_distances',PASS=all(x['matches_metadata'] and x['minimum_distance']>0 and x['max_norm_squared_error']<1e-12 and x['exact_cutoff_equalities']==0 for x in geometry)))
    checks.append(dict(name='final_generated_only_after_lock',PASS=lock['no_final_data_exists'] and all(x['after_lock'] for x in geometry if x['after_lock'] is not None)))
    contexts=[read(p) for p in (run/'cases').glob('*/context.json')]
    checks.append(dict(name='all_process_inputs_unchanged_during_process',PASS=all(c.get('inputs')==c.get('inputs_after') for c in contexts)))
    checks.append(dict(name='complete_exit_status_evidence',PASS=all('exit_status' in c and 'ended_at' in c for c in contexts)))
    rows=read(run/'final-timings.json');choices=read(run/'choices.json')
    checks.append(dict(name='final_five_methods_six_repeats',PASS=len(rows)==180 and all(sum(r['dataset']==x['dataset'] and r['method']==m for r in rows)==6 for x in choices for m in 'ABCDE')))
    checks.append(dict(name='final_choices_immutable_all_repeats',PASS=all(read(run/'choices'/(x['dataset']+f'-r{k}.json'))['schedule']==x['schedule'] for x in choices for k in range(6))))
    checks.append(dict(name='all_timed_methods_verbose_off',PASS=all('--verbose' not in c['argv'] for c in contexts if c['case'].startswith(('survey-','final-')))))
    checks.append(dict(name='exec_keeps_peak_RSS',PASS=read(run/'measurement-check.json')['PASS']))
    source_snapshots={sha(p):str(p.relative_to(run)) for p in run.rglob('*.py') if 'build' not in p.parts and 'third-party' not in p.parts and 'reproduction-check' not in p.parts}
    historical_driver_hashes=set(c['inputs'][str(Path(__file__).resolve().parent/'offline.py')] for c in contexts if str(Path(__file__).resolve().parent/'offline.py') in c['inputs'])
    checks.append(dict(name='historical_driver_sources_retained',PASS=all(h in source_snapshots for h in historical_driver_hashes),mapping={h:source_snapshots.get(h) for h in historical_driver_hashes}))
    final_checks=read(run/'final-correctness.json')
    result=dict(PASS=all(x['PASS'] for x in checks),checks=checks,geometry=geometry,raw_processes=len(contexts),
        final_PH_by_method=dict(Counter(r['status'] for r in rows)),unique_final_schedule_diagnostics=len(final_checks),
        final_zero_length_H1=sum(x['multiplicity'] for r in final_checks for x in r.get('bars',[]) if x['birth']==x['death']))
    (run/'AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(PASS=result['PASS'],checks=checks),indent=2));assert result['PASS']

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();audit(a.run.resolve())
