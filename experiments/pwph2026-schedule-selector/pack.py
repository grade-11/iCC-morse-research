#!/usr/bin/env python3
"""Package a frozen run, verify every tar member, and emit external SHA256 proof."""
import argparse
from datetime import datetime,timezone
import gzip
import hashlib
import json
from pathlib import Path
import tarfile
import time

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--archive',type=Path,required=True);a=p.parse_args()
    run=a.run.resolve();archive=a.archive.resolve();budget=json.loads((run/'budget.json').read_text());start=time.monotonic()
    assert not archive.exists(),'Never overwrite a delivered archive'
    def check_time():
        if time.monotonic()>=budget['deadline_monotonic']:raise TimeoutError('Original two-hour budget ended')
    def included(p):
        rel=p.relative_to(run)
        return (p.is_file() and not p.is_symlink() and not any(x in ('build','__pycache__','.git','tools-venv') for x in rel.parts)
                and rel.name not in ('manifest-sha256.json','runner.lock')
                and not str(rel).startswith(('reproduction-check/engine-source/','reproduction-check/third-party/')))
    paths=sorted(p for p in run.rglob('*') if included(p))
    files={}
    for f in paths:check_time();files[str(f.relative_to(run))]=sha(f)
    manifest=dict(created_at=datetime.now(timezone.utc).isoformat(),algorithm='SHA-256',files=files,
        exclusions=['manifest itself','build directories and binaries','venvs/cache/Git','duplicate reproduction engine/CLI11 source'],
        count=len(files),budget_start=budget['started_at'],deadline=budget['deadline'])
    mp=run/'manifest-sha256.json';mp.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    check_time();tmp=archive.with_suffix(archive.suffix+'.part')
    with tmp.open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',compresslevel=6,mtime=0) as gz:
            with tarfile.open(fileobj=gz,mode='w') as t:
                for f in [*paths,mp]:
                    check_time();t.add(f,arcname=run.name+'/'+str(f.relative_to(run)),recursive=False)
    verified=0
    with tarfile.open(tmp,'r:gz') as t:
        members=t.getmembers();assert len(members)==len(files)+1
        embedded=json.load(t.extractfile(run.name+'/manifest-sha256.json'));assert embedded==manifest
        for m in members:
            check_time();rel=str(Path(m.name).relative_to(run.name))
            assert m.isfile() and not Path(rel).is_absolute() and '..' not in Path(rel).parts
            if rel=='manifest-sha256.json':continue
            assert hashlib.sha256(t.extractfile(m).read()).hexdigest()==files[rel],rel
            verified+=1
    # Also ensure the workspace still matches the frozen manifest after packing.
    assert all(sha(run/name)==h for name,h in files.items())
    tmp.replace(archive);digest=sha(archive)
    archive.with_suffix(archive.suffix+'.sha256').write_text(digest+'  '+archive.name+'\n')
    result=dict(status='PASS',verified_at=datetime.now(timezone.utc).isoformat(),archive=str(archive),sha256=digest,
        manifest_sha256=sha(mp),verified_files=verified,archive_members=verified+1,archive_bytes=archive.stat().st_size,
        pack_elapsed_s=time.monotonic()-start,total_elapsed_from_original_start_s=time.monotonic()-budget['started_monotonic'],
        before_original_deadline=time.monotonic()<budget['deadline_monotonic'])
    archive.with_suffix(archive.suffix+'.verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
