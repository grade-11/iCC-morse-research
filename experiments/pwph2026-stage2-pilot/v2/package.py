#!/usr/bin/env python3
"""One self-contained evidence archive; no binaries/builds/venvs/caches."""
import argparse,hashlib,shutil,subprocess,tarfile
from lab import *

def package(run,prepare_only=False):
    snapshot=run/'source';snapshot.mkdir(exist_ok=True)
    tracked=subprocess.check_output(['git','ls-files','-z'],cwd=REPO).decode().split('\0')
    paths={Path(p) for p in tracked if p and not p.endswith(':Zone.Identifier')}
    paths.add(Path('Sphere gen/utilities.py'))
    paths.update(p.relative_to(REPO) for p in BASE.iterdir() if p.is_file() and p.suffix in ('.py','.md','.txt'))
    for rel in sorted(paths):
        src=REPO/rel
        if not src.is_file():continue
        dest=snapshot/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    cli=REPO/'build/_deps/cli11-src'
    cli_paths=subprocess.check_output(['git','ls-files','-z'],cwd=cli).decode().split('\0')
    for name in cli_paths:
        if not name or not (cli/name).is_file():continue
        dest=run/'third-party/CLI11'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(cli/name,dest)
    write(run/'CLI11-provenance.json',dict(head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=cli).decode().strip(),
          note='Tracked dependency source only, copied from existing FetchContent checkout; no .git/build/binaries.'))
    (run/'changes.patch').write_bytes(subprocess.check_output(['git','diff','dd92867ab1cb3a999fbaea05c4003ed2de53bc7b','--','main.cpp','deps','tests','experiments/pwph2026-stage2-pilot/v2'],cwd=REPO))
    for folder in ('build','build-trace'):
        for rel in ('CMakeCache.txt','CMakeFiles/iCC.dir/flags.make','CMakeFiles/icc_core.dir/flags.make','CMakeFiles/iCC.dir/link.txt'):
            src=run/folder/rel
            if src.exists():
                dest=run/'build-context'/folder/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    provenance=run/'previous-pilot';provenance.mkdir(exist_ok=True)
    old_archive=BASE.parent/'stage2-20261015T190630Z-handoff.tar.gz'
    need(sha(old_archive)=='5b4cc2e313c5d75d88175e023e842d1ef810e1505dfc8b740ea762b674b18ae3','Previous archive changed')
    shutil.copy2(old_archive,provenance/old_archive.name)
    write(run/'packaging-context.json',dict(at=now(),head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO).decode().strip(),
          binary_sha256=sha(run/'build/iCC'),note='Previous pilot archive embedded for old inputs/dependencies; source snapshot includes task files and unchanged generator utilities.'))
    if prepare_only:
        print('Source/dependency snapshot prepared; no archive created.');return
    def eligible(p):
        rel=p.relative_to(run)
        return rel.parts[0] not in {'build','build-trace','packaging-capture'} and not any(s in {'__pycache__','.cache','.venv'} for s in rel.parts) and not p.name.endswith(':Zone.Identifier')
    files=sorted(p for p in run.rglob('*') if p.is_file() and eligible(p) and p.name!='MANIFEST.sha256.json')
    manifest={str(p.relative_to(run)):sha(p) for p in files};write(run/'MANIFEST.sha256.json',manifest)
    archive=BASE/(run.name+'-handoff.tar.gz');need(not archive.exists(),'Refuse overwrite archive')
    with tarfile.open(archive,'w:gz') as tar:
        for p in files+[run/'MANIFEST.sha256.json']:tar.add(p,arcname=str(Path(run.name)/p.relative_to(run)),recursive=False)
    with tarfile.open(archive) as tar:
        expected={str(Path(run.name)/p) for p in manifest}|{str(Path(run.name)/'MANIFEST.sha256.json')}
        need(set(tar.getnames())==expected,'Archive membership differs')
        for rel,h in manifest.items():
            payload=tar.extractfile(str(Path(run.name)/rel)).read()
            need(not payload.startswith(b'\x7fELF'),'ELF binary in archive: '+rel)
            need(hashlib.sha256(payload).hexdigest()==h,'Archive checksum: '+rel)
    archive.with_suffix('.gz.sha256').write_text(sha(archive)+'  '+archive.name+'\n')
    print(json.dumps(dict(archive=str(archive),sha256=sha(archive),files=len(files)+1,verified=True),indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--prepare-only',action='store_true');a=ap.parse_args();package(a.run,a.prepare_only)
