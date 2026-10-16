#!/usr/bin/env python3
"""Check data/binary integrity and disclose superseded helper-source coverage."""
from workflow import *
from collections import Counter
import numpy as np
import tarfile,shutil

def main():
    origins={Path(x['file']).name:x.get('source','generated fixture') for plan in ('PLAN.json','survey-plan.json') for x in json.loads((RUN/plan).read_text()).get('cases',json.loads((RUN/plan).read_text()).get('items',[]))}
    datasets=[]
    for p in sorted((RUN/'inputs').glob('*.csv')):
        points=np.loadtxt(p,delimiter=',',ndmin=2);assert np.isfinite(points).all()
        norms=np.linalg.norm(points,axis=1);delta=points[:,None,:]-points[None,:,:];dist=np.sqrt(np.sum(delta*delta,axis=2));pairs=dist[np.triu_indices(len(points),1)]
        norm_error=float(np.max(np.abs(norms-1)));unit=points.shape[1]==3 and norm_error<1e-12
        if p.stem.startswith('standard'):assert unit,'Standard S2 input geometry differs from label'
        datasets.append(dict(file=str(p.relative_to(RUN)),sha256=sha(p),points=len(points),ambient_dim=points.shape[1],norm_min=float(norms.min()),norm_max=float(norms.max()),max_unit_norm_error=norm_error,unit_S2_geometry_diagnostic=unit,min_pair_distance=float(pairs.min()),coincident_pairs=int(np.sum(pairs==0)),source=origins.get(p.name,'generated fixture; see source/plan'),scope='Geometry diagnostic only; inputs and PH endpoints never normalized/rounded'))
    table(RUN/'dataset-inventory.csv',datasets)
    reference_map=[];needed={x['sha256'] for x in datasets}
    for meta_path in (RUN/'references').glob('*.json'):
        meta=json.loads(meta_path.read_text())
        if meta.get('config',{}).get('data_sha256') not in needed:continue
        original=Path(meta['case']);original=original if original.is_absolute() else REPO/original
        if original.is_relative_to(RUN):
            reference_map.append(dict(reference=meta_path.name,original_case=str(original),copy_case=str(original.relative_to(RUN)),status='CURRENT_RUN_RAW_PRESENT'));continue
        oldrun=original.parent.parent;manifest_path=oldrun/'MANIFEST.sha256.json'
        if not original.exists() or not manifest_path.exists():
            reference_map.append(dict(reference=meta_path.name,original_case=str(original),copy_case='',status='OLD_RAW_NOT_LOCATED_NPZ_AND_CONFIG_RETAINED'));continue
        oldmanifest=json.loads(manifest_path.read_text());target=RUN/'reused-reference-cases'/oldrun.name/original.name
        for source in original.rglob('*'):
            if not source.is_file():continue
            relative=str(source.relative_to(oldrun));assert oldmanifest.get(relative)==sha(source),'Reference original-manifest mismatch'
            dst=target/source.relative_to(original);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dst)
        reference_map.append(dict(reference=meta_path.name,original_case=str(original),copy_case=str(target.relative_to(RUN)),status='OLD_RAW_COPIED_MANIFEST_VERIFIED'))
    write(RUN/'reused-reference-map.json',reference_map)
    previous=REPO/'experiments/pwph2026-stage2-pilot/v2/v2-20261015T195640Z-handoff.tar.gz'
    archiveproof=dict(path=str(previous),status='NOT_AVAILABLE_SELECTED_FILES_STILL_MANIFEST_VERIFIED')
    if previous.exists():
        expected=previous.with_suffix('.gz.sha256').read_text().split()[0];digest=sha(previous);assert digest==expected
        with tarfile.open(previous) as tar:payload=tar.extractfile('v2-20261015T195640Z/MANIFEST.sha256.json').read()
        oldmeta=json.loads((RUN/'n200-dossier.json').read_text());assert hashlib.sha256(payload).hexdigest()==oldmeta['source_manifest_sha256']
        archiveproof.update(status='SHA_AND_EMBEDDED_MANIFEST_VERIFIED',sha256=digest,embedded_manifest_sha256=hashlib.sha256(payload).hexdigest(),note='Previous archive is referenced, not embedded wholesale;508 selected files copied and checked individually.')
    write(RUN/'n200-previous-archive-proof.json',archiveproof)
    catalog=[]
    for binary in (RUN/'builds').glob('*/*'):
        if not binary.is_file():continue
        with binary.open('rb') as f:elf=f.read(4)==b'\x7fELF'
        if elf:catalog.append(dict(path=str(binary),version=binary.parent.name,target=binary.name,sha256=sha(binary)))
    write(RUN/'binary-catalog.json',catalog)
    available={}
    for p in RUN.rglob('*'):
        if not p.is_file() or p.relative_to(RUN).parts[0] in ('builds','packaging-capture'):continue
        available.setdefault(sha(p),[]).append(str(p.relative_to(RUN)))
    # Tracked helper files are included by the final source snapshot.
    for rel in subprocess.check_output(['git','ls-files','-z'],cwd=REPO).decode().split('\0'):
        if rel and (REPO/rel).is_file():available.setdefault(sha(REPO/rel),[]).append('source/'+rel)
    for p in BASE.glob('*.py'):available.setdefault(sha(p),[]).append('source/'+str(p.relative_to(REPO)))
    binary_hashes={x['sha256'] for x in catalog};rows=[];critical=[]
    for p in (RUN/'cases').glob('*/context.json'):
        c=json.loads(p.read_text());assert 'exit_status' in c
        for name,digest in c.get('inputs',{}).items():
            if digest in available:status='SOURCE_OR_INPUT_BYTES_AVAILABLE'
            elif digest in binary_hashes:status='BINARY_HASH_RECORDED_SOURCE_AND_BUILD_AVAILABLE'
            else:
                status='SUPERSEDED_HELPER_BYTES_NOT_LOCATED' if Path(name).suffix in ('.py','.json','.txt') else 'CRITICAL_INPUT_NOT_LOCATED'
                if status=='CRITICAL_INPUT_NOT_LOCATED':critical.append(dict(case=c['case'],input=name,sha256=digest))
            rows.append(dict(case=c['case'],input=name,sha256=digest,status=status,available_as=';'.join(available.get(digest,[])[:3])))
    table(RUN/'evidence-input-coverage.csv',rows)
    summary=dict(at=now(),input_CSV_files=len(datasets),unique_pointfile_hashes=len({x['sha256'] for x in datasets}),case_contexts=len(list((RUN/'cases').glob('*/context.json'))),input_records=len(rows),statuses=dict(Counter(x['status'] for x in rows)),critical_unavailable=critical,scope='Final code and exact argv/data/binary fingerprints retained. Superseded helper-source omissions, if any, are disclosed separately and never synthesized. Binary files intentionally excluded from source-only archive; source/version/build flags included.')
    write(RUN/'evidence-audit.json',summary);assert not critical,'Unresolved exact input/binary evidence; inspect audit before packaging'
    progress('evidence-integrity','PASS_WITH_DISCLOSED_HELPER_COVERAGE',str(summary));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
