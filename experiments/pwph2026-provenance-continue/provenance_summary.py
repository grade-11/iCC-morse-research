#!/usr/bin/env python3
from matrix import *
import subprocess,difflib

def main():
    versions=json.loads((RUN/'versions.json').read_text());rows=[]
    descriptions={
      'U':'Untouched upstream remote HEAD/main/parallel verified directly; exact paper artifact NOT pinned',
      'U17':'U + output precision ONLY; algorithm unchanged; never call untouched',
      'B17':'Group pre instrumentation/controller; precision + steady-clock changes',
      'C':'Group instrumentation before G',
      'D':'Frozen G v2 current at task start; fixed actions and old-grid G/replay',
      'D_NO_CLOCK':'D diagnostic fixed-action copy; clocks return zero; NOT G/performance eligible',
      'D_NO_INTERNAL_RAM':'D diagnostic fixed-action copy; internal RAM sampling disabled, external monitor active',
      'S':'Opt-in explicit F2 H0/H1 no-quotient stream/rebuild; native algorithm unchanged'}
    for key,v in versions.items():
        rows.append(dict(version=key,commit=v['commit'],description=descriptions[key],binary_sha256=v.get('binary_sha256',''),paper_artifact='UNPINNED',compatible_2026_multiwindow=True,source=v['source']))
    for label,short,description in [('group-start','f410424','First group commit adds sphere generator; C++ kernel identical to U, not separately rebuilt'),('group-data-baseline','3c59f45','Group dataset baseline before C++ edits; C++ kernel identical to U, not separately rebuilt')]:
        full=subprocess.check_output(['git','rev-parse',short],cwd=REPO).decode().strip()
        assert not subprocess.check_output(['git','diff','734be3d',short,'--','main.cpp','deps','CMakeLists.txt'],cwd=REPO)
        rows.append(dict(version=label,commit=full,description=description,binary_sha256='',paper_artifact='UNPINNED',compatible_2026_multiwindow=True,source='Git commit '+full+'; kernel source exported as U'))
    # Export historical source as evidence, do not adapt its entrypoint into upstream.
    commit=subprocess.check_output(['git','rev-parse','044f443'],cwd=REPO).decode().strip()
    for rel in ('README.md','main.cpp','CMakeLists.txt','deps/QuotientAndExpand.cpp'):
        result=subprocess.run(['git','show',commit+':'+rel],cwd=REPO,capture_output=True)
        dest=RUN/'versions/historical-2025'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        if result.returncode==0:dest.write_bytes(result.stdout)
        else:write(dest.with_suffix(dest.suffix+'.unavailable.json'),dict(argv=['git','show',commit+':'+rel],exit_status=result.returncode,stderr=result.stderr.decode()))
    rows.append(dict(version='historical-2025',commit=commit,description='Historical snapshot, not proven paper commit; hardcoded entrypoint and different Q&E/merged supports; NOT_EQUIVALENT',binary_sha256='',paper_artifact='UNPINNED',compatible_2026_multiwindow=False,source=str(RUN/'versions/historical-2025')))
    table(RUN/'provenance.csv',rows)
    u=Path(versions['U']['source']);p=Path(versions['U17']['source'])
    (RUN/'versions/U17-output-only.patch').write_text(''.join(difflib.unified_diff((u/'deps/Pipelines.cpp').read_text().splitlines(True),(p/'deps/Pipelines.cpp').read_text().splitlines(True),fromfile='a/deps/Pipelines.cpp',tofile='b/deps/Pipelines.cpp')))
    history={}
    for label,args in {
      'remotes':['remote','-v'],
      'upstream_ancestry':['log','--reverse','--format=%H %aI %s','734be3d..259d4bf'],
      'base_cpp_diff':['diff','734be3d','3c59f45','--','main.cpp','deps','CMakeLists.txt'],
      'B17_cpp_diff':['diff','734be3d','36b9db8','--','main.cpp','deps','CMakeLists.txt'],
      'C_cpp_diff':['diff','36b9db8','c8cd0df','--','main.cpp','deps','CMakeLists.txt'],
      'D_cpp_diff':['diff','c8cd0df','259d4bf','--','main.cpp','deps','CMakeLists.txt'],
      'S_cpp_diff':['diff','259d4bf','50636a5','--','main.cpp','deps','tests','CMakeLists.txt']}.items():
        result=subprocess.run(['git',*args],cwd=REPO,capture_output=True)
        dest=RUN/'provenance-details'/(label+'.txt');dest.parent.mkdir(exist_ok=True);dest.write_bytes(result.stdout)
        history[label]=dict(argv=['git',*args],exit_status=result.returncode,stderr=result.stderr.decode(),path=str(dest),sha256=sha(dest))
    write(RUN/'provenance-git-commands.json',history)
    original=json.loads((RUN/'correctness.json').read_text());lookup={(r['version'],r['case'],r['dimension'],r['threads'],r['mode']):r for r in original};comparisons=[]
    for r in original:
        if r['version']!='U17':continue
        key=(r['case'],r['dimension'],r['threads'],r['mode']);b=counters(r)
        for v in ('B17','C','D'):
            other=lookup[(v,*key)];equal=counters(other)==b;assert equal
            comparisons.append(dict(case=r['case'],dimension=r['dimension'],threads=r['threads'],mode=r['mode'],comparison='U17 vs '+v,status='EXACT_MULTISET_EQUAL',left_path=r['path'],right_path=other['path']))
        u6=lookup[('U',*key)]
        projected={d:Counter() for d in b}
        for d,pairs in b.items():
            for (a,z),n in pairs.items():projected[d][(Decimal(format(float(a),'.6g')),Decimal(format(float(z),'.6g')))]+=n
        equal=projected==counters(u6);assert equal,'Unmodified U differs beyond six significant output digits'
        comparisons.append(dict(case=r['case'],dimension=r['dimension'],threads=r['threads'],mode=r['mode'],comparison='U untouched vs U17 projected6sig',status='PRINTED_PRECISION_MULTISET_EQUAL',left_path=u6['path'],right_path=r['path']))
    for r in original:
        if r['threads']!=1:continue
        other=lookup[(r['version'],r['case'],r['dimension'],4,r['mode'])];assert counters(other)==counters(r)
        comparisons.append(dict(case=r['case'],dimension=r['dimension'],threads='1 vs4',mode=r['mode'],comparison=r['version']+' threads',status='EXACT_AT_AVAILABLE_OUTPUT_PRECISION',left_path=r['path'],right_path=other['path']))
    table(RUN/'cross-version.csv',comparisons)
    verdicts=[]
    for r in original:
        if r['version']!='U17' or r['mode']!='quotient':continue
        verdicts.append(dict(case=r['case'],dimension=r['dimension'],threads=r['threads'],status=r['status'],verdict='Reproduced on upstream734 (U untouched + U17 precision-only corroboration)' if r['status']=='PH_MISMATCH' else 'PASS on tested case only',H0='UPSTREAM_NOT_EXPORTED',upstream_U_path=lookup[('U',r['case'],r['dimension'],r['threads'],r['mode'])]['path'],upstream_U17_path=r['path']))
    for r in json.loads((RUN/'counterexample8.json').read_text()):
        if r['version']=='U17':verdicts.append(dict(case='S2-n8-'+r['schedule'],dimension=2,threads=r['threads'],status='PASS' if r['quotient_equal_baseline'] else 'PH_MISMATCH',verdict='Reproduced on upstream734 (H1 future interval loss)' if not r['quotient_equal_baseline'] else 'PASS on this late schedule only',H0='UPSTREAM_NOT_EXPORTED',upstream_U_path=next(x['quotient']['path'] for x in json.loads((RUN/'counterexample8.json').read_text()) if x['version']=='U' and x['threads']==r['threads'] and x['schedule']==r['schedule']),upstream_U17_path=r['quotient']['path']))
    table(RUN/'case-verdicts.csv',verdicts)
    write(RUN/'cross-version-summary.json',dict(comparisons=len(comparisons),all_equal=True,group_first_bad_change='NONE_OBSERVED_ON_TESTED_FIXED_CONFIGURATIONS',bisect='NOT_REQUIRED: no group PASS-to-FAIL pair found',upstream_H0='NOT_EXPORTED; external UF cannot prove upstream execution path',paper_artifacts='UNPINNED for both2025 and2026; do not identify repository currentHEAD with exact publication artifact',limits='Counterexamples demonstrate these configurations, not a universal conclusion about every failure or paper statement.'))
    progress('provenance-summary','COMPLETE','Direct U/U17/B17/C/D matrix and thread comparison saved; paper commits unpinned, historical2025 not equivalent')

if __name__=='__main__':main()
