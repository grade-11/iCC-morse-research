#!/usr/bin/env python3
"""Build/test the source-only snapshot with the bundled CLI11 dependency."""
from stream_checks import *

def main():
    source=RUN/'source';build=RUN/'builds/OFFLINE_SOURCE';steps=[]
    commands=[('configure',['cmake','-S',str(source),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DBUILD_TESTING=ON','-DCMAKE_DISABLE_FIND_PACKAGE_CLI11=ON','-DFETCHCONTENT_FULLY_DISCONNECTED=ON','-DFETCHCONTENT_SOURCE_DIR_CLI11='+str(RUN/'third-party/CLI11')]),('build',['cmake','--build',str(build),'-j4']),('ctest',['ctest','--test-dir',str(build),'--output-on-failure'])]
    for name,argv in commands:
        out,c=execute('offline-source-'+name,argv,[source/'CMakeLists.txt',source/'deps/StreamingF2.cpp',RUN/'third-party/CLI11/CMakeLists.txt'])
        steps.append(dict(step=name,status=c['status'],exit_status=c['exit_status'],path=str(out)));write(RUN/'offline-source-build.json',steps);assert c['status']=='OK'
    data=RUN/'inputs/two-squares-s2-8.csv';item=dict(id='two-squares-s2-8',file=str(data));binary=build/'iCC'
    original=next(r for r in json.loads((RUN/'counterexample8.json').read_text()) if r['version']=='S' and r['threads']==4 and r['schedule']=='single')
    checks=[]
    for mode,grid,want in [('ph',None,counters(original['baseline'])),('quotient',[.5,2.],counters(original['quotient']))]:
        out,c=execute('offline-source-n8-'+mode,reviewed.command(binary,data,grid,maxdim=2,threads=4),[binary,data]);assert c['status']=='OK'
        bars,_,_=parse_output(out,data,grid,2,4);assert bars==want
        checks.append(dict(mode=mode,status='EXPECTED_OUTPUT_MATCH',path=str(out),scope='Native quotient retains its known mismatch; not silently repaired'))
    inc=stream(binary,item,[.5,2.],'f2-stream',4,True,'offline-source-smoke');assert inc['status']=='OK' and inc['bars']==counters(original['baseline'])[1]
    checks.append(dict(mode='f2-stream',status='EXACT_NATIVE_BASELINE_MATCH',path=inc['path']))
    write(RUN/'offline-smoke.json',dict(status='PASS',source=str(source),dependency='Bundled CLI11 tracked source; fully disconnected CMake',binary=str(binary),binary_sha256=sha(binary),steps=steps,checks=checks,performance='Not a benchmark; no timing reused in performance tables'))
    progress('offline-source-smoke','PASS','Source-only offline build +4 CTests +n8 native/quotient/F2 behavior verified')

if __name__=='__main__':main()
