#!/usr/bin/env python3
"""Check evidence completeness without running PH/benchmarks again."""
import argparse
from lab import *

def audit(run):
    errors=[];contexts=[]
    for p in (run/'cases').glob('*/context.json'):
        c=json.loads(p.read_text());contexts.append(c)
        if 'exit_status' not in c:errors.append('Interrupted: '+str(p))
        if c.get('inputs_after')!=c.get('inputs'):errors.append('Input changed during process: '+str(p))
        if not all((p.parent/n).exists() for n in ('stdout.log','stderr.log','resources.txt')):errors.append('Missing raw evidence: '+str(p))
    final=json.loads((run/'final-validation.json').read_text());lock=json.loads((run/'LOCK-v2.json').read_text())
    for p,h in {**lock['source_hashes'],**lock['model_hashes']}.items():
        if sha(p)!=h:errors.append('Locked input changed: '+p)
    if sha(run/'build/iCC')!=lock['binary_sha256']:errors.append('Locked binary changed')
    final_ids={x['dataset_id'] for x in final};need(len(final_ids)>=8,'Fewer than eight final datasets')
    need(len(final)==len(final_ids)*2*4,'Missing final logical configurations')
    from datetime import datetime
    events=[json.loads(line) for line in (run/'index.jsonl').read_text().splitlines()]
    for dataset_id in final_ids:
        generated=[e for e in events if e.get('case')=='generate-'+dataset_id and e.get('status')=='DATA_CREATED']
        if len(generated)!=1 or datetime.fromisoformat(generated[0]['at'])<=datetime.fromisoformat(lock['locked_at']):
            errors.append('Final data creation not uniquely after lock: '+dataset_id)
    dataset_hashes=[]
    for p in (run/'data').glob('*.json'):
        m=json.loads(p.read_text());csv=p.with_suffix('.csv')
        # Metadata key schema is retained from the generator; validate any explicit checksum.
        expected=m.get('sha256',m.get('data_sha256'))
        if expected and sha(csv)!=expected:errors.append('Data checksum: '+str(csv))
        dataset_hashes.append(sha(csv))
    if len(dataset_hashes)!=len(set(dataset_hashes)):errors.append('Duplicate new datasets by byte hash')
    newpairs=json.loads((run/'calibration-global-dev.json').read_text())+json.loads((run/'calibration-global-validation.json').read_text())
    for x in newpairs:
        a=x['branches']['0']['state'];b=x['branches']['1']['state']
        if any(a[k]!=b[k] for k in ('state_hash','A','C','pv','raw','entries','eps','h')):errors.append('Global prefix differs: '+x['dataset'])
    measurements=[json.loads(p.read_text()) for p in (run/'cases').glob('measure-final-*/ph-check.json')]
    result=dict(at=now(),status='PASS' if not errors else 'ERROR',errors=errors,process_contexts=len(contexts),final_datasets=len(final_ids),
                final_logical_configurations=len(final),final_measurement_processes=len(measurements),
                non_OK_processes=[dict(case=c['case'],status=c['status'],exit_status=c.get('exit_status')) for c in contexts if c['status']!='OK'],
                note='PH failures remain failures; this audits provenance/completeness, not a general PH theorem.')
    write(run/'AUDIT.json',result);print(json.dumps(result,indent=2));need(not errors,'Evidence audit failed')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();audit(a.run)
