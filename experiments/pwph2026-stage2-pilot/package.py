#!/usr/bin/env python3
"""Package a finished run; exclude all builds, binaries, venvs and caches."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
from pilot import BASE,REPO,sha,write,now

def package(run):
    snapshot=run/'source';snapshot.mkdir(exist_ok=False)
    tracked=subprocess.check_output(['git','ls-files','-z']).decode().split('\0')
    files=[Path(p) for p in tracked if p and not p.endswith(':Zone.Identifier')]
    files += [Path('deps/AdaptiveSchedule.hpp'),Path('Sphere gen/utilities.py')]
    files += [p.relative_to(REPO) for p in BASE.iterdir() if p.is_file() and p.suffix in ('.py','.md','.txt') and
              not p.name.startswith(('GIAO-','CURRENT-'))]
    for p in files:
        dest=snapshot/p;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    (run/'changes.patch').write_bytes(subprocess.check_output(['git','diff','c8cd0df34fe5e20d090994a53470312e069b1a24','--','main.cpp','deps','tests']))
    bc=run/'build-context';bc.mkdir(exist_ok=True)
    for folder in ('build','build-original'):
        for rel in ('CMakeCache.txt','CMakeFiles/iCC.dir/flags.make','CMakeFiles/icc_core.dir/flags.make','CMakeFiles/iCC.dir/link.txt'):
            p=run/folder/rel
            if p.exists():
                dest=bc/folder/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
    original=run/'n100-original-evidence';original.mkdir(exist_ok=True)
    n100=REPO/'experiments/pwph2026-s2-n100'
    for name in ('whole-run-monotonic-ma523yck.tar.gz','ripser-reference-noa_bg6w.tar.gz'):
        shutil.copy2(n100/'runs'/name,original/name)
    shutil.copy2(n100/'runs/repeat-verbose-TIVmZd/run-01/piecewise.log',original/'old-piecewise-verbose.log')
    (run/'packaging-context.json').write_text(json.dumps(dict(packaged_at=now(),
        head=subprocess.check_output(['git','rev-parse','HEAD']).decode().strip(),
        branch=subprocess.check_output(['git','branch','--show-current']).decode().strip(),
        note='Tracked source plus explicitly added task files; user utilities copied for generator provenance, not modified.',
        binary_sha256=sha(run/'build/iCC')),indent=2))
    def eligible(p):
        rel=p.relative_to(run)
        return (rel.parts[0] not in {'build','build-original','packaging-capture'} and
                not any(x in {'__pycache__','.cache','.venv'} for x in rel.parts) and
                not p.name.endswith(':Zone.Identifier'))
    paths=sorted(p for p in run.rglob('*') if p.is_file() and eligible(p) and p.name!='MANIFEST.sha256.json')
    manifest={str(p.relative_to(run)):sha(p) for p in paths}
    write(run/'MANIFEST.sha256.json',manifest)
    archive=BASE/(run.name+'-handoff.tar.gz')
    if archive.exists():raise RuntimeError('Archive exists; refuse overwrite')
    with tarfile.open(archive,'w:gz') as tar:
        for p in paths+[run/'MANIFEST.sha256.json']:tar.add(p,arcname=str(Path(run.name)/p.relative_to(run)),recursive=False)
    # Verify every archived member against the on-disk manifest, including no extras.
    import hashlib
    with tarfile.open(archive) as tar:
        expected={str(Path(run.name)/p) for p in manifest}|{str(Path(run.name)/'MANIFEST.sha256.json')}
        if set(tar.getnames())!=expected:raise RuntimeError('Archive membership mismatch')
        for name,h in manifest.items():
            if hashlib.sha256(tar.extractfile(str(Path(run.name)/name)).read()).hexdigest()!=h:
                raise RuntimeError('Archive checksum mismatch: '+name)
    (archive.with_suffix(archive.suffix+'.sha256')).write_text(sha(archive)+'  '+archive.name+'\n')
    print(json.dumps(dict(archive=str(archive),sha256=sha(archive),files=len(paths)+1,verified=True),indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',required=True,type=Path);a=ap.parse_args();package(a.run)
