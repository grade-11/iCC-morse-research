#!/usr/bin/env python3
"""Source-only, offline-buildable handoff with every member verified."""
from workflow import *
import argparse,shutil,tarfile

def snapshot():
    destroot=RUN/'source';destroot.mkdir(exist_ok=True)
    names=subprocess.check_output(['git','ls-files','-z'],cwd=REPO).decode().split('\0')
    paths={Path(n) for n in names if n and not n.endswith(':Zone.Identifier')}
    paths.add(Path('Sphere gen/utilities.py'))
    paths.update(p.relative_to(REPO) for p in BASE.iterdir() if p.is_file() and (p.suffix in ('.py','.md','.txt') or p.name=='.gitignore'))
    for rel in sorted(paths):
        src=REPO/rel
        if src.is_file():
            dst=destroot/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    versions=json.loads((RUN/'versions.json').read_text())
    for rel,digest in versions['S']['source_hashes'].items():assert sha(destroot/rel)==digest,'Compiled S source drift: '+rel
    cli=REPO/'build/_deps/cli11-src'
    cli_names=subprocess.check_output(['git','ls-files','-z'],cwd=cli).decode().split('\0')
    for name in cli_names:
        if name and (cli/name).is_file():
            dst=RUN/'third-party/CLI11'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(cli/name,dst)
    write(RUN/'CLI11-provenance.json',dict(head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=cli).decode().strip(),source=str(cli),note='Tracked dependency source only, from existing checkout; no network needed for reproduction.'))
    for build in (RUN/'builds').iterdir():
        for pattern in ('CMakeCache.txt','CMakeFiles/*/CMakeCXXCompiler.cmake','CMakeFiles/*.dir/flags.make','CMakeFiles/*.dir/link.txt','Testing/Temporary/*.log'):
            for src in build.glob(pattern):
                dst=RUN/'build-context'/build.name/src.relative_to(build);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    result=subprocess.run(['git','diff','259d4bf','HEAD','--','main.cpp','deps','tests','CMakeLists.txt','README.md','experiments/pwph2026-provenance-continue'],cwd=REPO,capture_output=True)
    assert result.returncode==0;(RUN/'changes.patch').write_bytes(result.stdout)
    # Full small historical tree for inspection; no claim that it is executable-equivalent.
    commit=subprocess.check_output(['git','rev-parse','044f443'],cwd=REPO).decode().strip()
    for name in subprocess.check_output(['git','ls-tree','-r','--name-only',commit],cwd=REPO).decode().splitlines():
        content=subprocess.check_output(['git','show',commit+':'+name],cwd=REPO)
        dst=RUN/'versions/historical-2025'/name;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(content)
    write(RUN/'packaging-context.json',dict(at=now(),head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO).decode().strip(),binary_S_sha256=versions['S']['binary_sha256'],snapshot_cpp_hashes_match_compiled_S=True,excluded='builds, ELF executables, virtualenvs/cache, mutable packaging-capture/locks',manifest='All eligible raw results, inputs, source, patches, selected prior n200 evidence and CLI11; SHA checked for every archive member.'))

def pack():
    archive=BASE/(RUN.name+'-handoff.tar.gz');assert not archive.exists(),'Refuse overwriting archive'
    def eligible(p):
        rel=p.relative_to(RUN)
        return rel.parts[0] not in {'builds','packaging-capture'} and not any(x in {'__pycache__','.cache','.venv','.git'} for x in rel.parts) and not p.name.endswith((':Zone.Identifier','.lock')) and p.name!='MANIFEST.sha256.json'
    files=sorted(p for p in RUN.rglob('*') if p.is_file() and eligible(p))
    manifest={str(p.relative_to(RUN)):sha(p) for p in files};write(RUN/'MANIFEST.sha256.json',manifest)
    with tarfile.open(archive,'w:gz',dereference=True) as tar:
        for p in files+[RUN/'MANIFEST.sha256.json']:tar.add(p,arcname=str(Path(RUN.name)/p.relative_to(RUN)),recursive=False)
    with tarfile.open(archive) as tar:
        expected={str(Path(RUN.name)/p) for p in manifest}|{str(Path(RUN.name)/'MANIFEST.sha256.json')}
        assert set(tar.getnames())==expected,'Archive membership mismatch'
        for rel,digest in manifest.items():
            payload=tar.extractfile(str(Path(RUN.name)/rel)).read()
            assert not payload.startswith(b'\x7fELF'),'ELF executable included: '+rel
            assert hashlib.sha256(payload).hexdigest()==digest,'Member checksum mismatch: '+rel
    digest=sha(archive);archive.with_suffix('.gz.sha256').write_text(digest+'  '+archive.name+'\n')
    print(json.dumps(dict(archive=str(archive),sha256=digest,files=len(files)+1,bytes=archive.stat().st_size,all_members_verified=True),indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prepare-only',action='store_true');ap.add_argument('--pack-only',action='store_true');args=ap.parse_args()
    if not args.pack_only:snapshot()
    if not args.prepare_only:pack()
