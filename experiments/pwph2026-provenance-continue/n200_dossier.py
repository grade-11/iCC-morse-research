#!/usr/bin/env python3
"""Reuse selected immutable n200 evidence, with original-manifest verification."""
from workflow import *
import shutil

def main():
    old=REPO/'experiments/pwph2026-stage2-pilot/v2/runs/v2-20261015T195640Z'
    target=RUN/'n200-reused';target.mkdir(exist_ok=True)
    dataset='uniform_s2-n200-s20261033';data=old/'data'/(dataset+'.csv');digest=sha(data)
    versions=json.loads((RUN/'versions.json').read_text());v=[x for x in json.loads((old/'final-validation.json').read_text()) if x['dataset_id']==dataset]
    assert len(v)==8 and all(x['status']=='PASS' and x['H0']=='PASS' and x['binary_sha256']==versions['D']['binary_sha256'] and x['dataset_sha256']==digest for x in v)
    manifest=json.loads((old/'MANIFEST.sha256.json').read_text());paths=set()
    def collect(p):
        if p.is_dir():paths.update(x for x in p.rglob('*') if x.is_file())
        elif p.is_file():paths.add(p)
    for p in old.glob('cases/*/context.json'):
        c=json.loads(p.read_text())
        if digest in c.get('inputs',{}).values() or str(data.relative_to(REPO)) in c.get('argv',[]) or str(data) in c.get('argv',[]):collect(p.parent)
    for p in old.glob('references/*.json'):
        m=json.loads(p.read_text())
        if m.get('config',{}).get('data_sha256')==digest:collect(p);collect(p.with_suffix('.npz'))
    for name in ('LOCK-v2.json','active-models.json','selected-models.json','global-models.json','initial-context.json','environment-context.json','budget.json','packaging-context.json'):
        collect(old/name)
    collect(old/'models');collect(data);collect(data.with_suffix('.json'))
    copied=[]
    for p in sorted(paths):
        rel=str(p.relative_to(old));assert manifest.get(rel)==sha(p),'Old manifest mismatch: '+rel
        dest=target/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
        copied.append(dict(original=str(p),copy=str(dest.relative_to(RUN)),sha256=sha(dest),verified_against_previous_manifest=True))
    measures=[x for x in csv.DictReader((old/'measurements.csv').open()) if x['dataset']==dataset]
    perf=[x for x in csv.DictReader((old/'performance.csv').open()) if x['dataset']==dataset]
    ratios=[x for x in csv.DictReader((old/'paired-ratios.csv').open()) if x.get('dataset')==dataset]
    write(target/'validation-selected.json',v);table(target/'measurements-selected.csv',measures);table(target/'performance-selected.csv',perf);table(target/'paired-ratios-selected.csv',ratios)
    write(target/'copy-map.json',copied)
    result=dict(status='REUSED_VERIFIED',dataset=dataset,input_sha256=digest,binary_version='D / 259d4bf (NOT new S)',binary_sha256=versions['D']['binary_sha256'],source_manifest=str(old/'MANIFEST.sha256.json'),source_manifest_sha256=sha(old/'MANIFEST.sha256.json'),configurations=len(v),validation_statuses=dict(__import__('collections').Counter(x['status'] for x in v)),warmups=sum(x['warmup']=='True' for x in measures),measured=sum(x['warmup']=='False' for x in measures),copied_files=len(copied),note='Only this dataset/configuration scope; no new-grid controller result, no full-paper reproduction. Original commands retain original paths; copy-map supplies relocations. Rebuild D using exported source and CLI11.')
    write(RUN/'n200-dossier.json',result);progress('n200-dossier','REUSED_VERIFIED',str(result));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
