#!/usr/bin/env python3
import tarfile,shutil
from workflow import *

def prepare():
    versions={}
    for label,commit in [('U','734be3d'),('B17','36b9db8'),('C','c8cd0df'),('D','259d4bf')]:
        src=RUN/'versions'/label;src.mkdir(parents=True,exist_ok=True);archive=RUN/'versions'/(label+'.tar')
        o,c=execute('export-'+label,['git','archive','--format=tar','--output='+str(archive),commit],[REPO/'CMakeLists.txt'])
        assert c['status']=='OK'
        with tarfile.open(archive) as t:t.extractall(src,filter='data')
        archive.unlink() # transient Git transport, source tree preserved below
        versions[label]=dict(commit=subprocess.check_output(['git','rev-parse',commit],text=True).strip(),source=str(src),changes='NONE')
    src=RUN/'versions/U17';shutil.copytree(RUN/'versions/U',src,dirs_exist_ok=True)
    p=src/'deps/Pipelines.cpp';s=p.read_text().replace('#include <iostream>','#include <iostream>\n#include <iomanip>');old='std::cout << pair.facet_weight << ", " << pair.cofacet_weight;';assert old in s;s=s.replace(old,'const auto precision = std::cout.precision();\n    std::cout << std::setprecision(17) << pair.facet_weight << ", " << pair.cofacet_weight;\n    std::cout.precision(precision);');p.write_text(s)
    versions['U17']=dict(commit=versions['U']['commit'],source=str(src),changes='Output precision only, 6 to17 significant digits; no algorithm/timing changes')
    for label,v in versions.items():
        src=Path(v['source']);v['source_hashes']={str(p.relative_to(src)):sha(p) for p in src.rglob('*') if p.is_file() and ('deps' in p.parts or p.name in ('main.cpp','CMakeLists.txt'))}
        build=RUN/'builds'/label;v['binary']=str(build/'iCC');write(RUN/'versions.json',versions)
        for kind,argv in [('configure',['cmake','-S',str(src),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DBUILD_TESTING=ON','-DFETCHCONTENT_SOURCE_DIR_CLI11='+str(REPO/'build/_deps/cli11-src'),'-DFETCHCONTENT_FULLY_DISCONNECTED=ON']),('build',['cmake','--build',str(build),'-j4'])]:
            out,c=execute(kind+'-'+label,argv,[src/'CMakeLists.txt',src/'deps/Pipelines.cpp'],300);print(label,kind,c['status'],flush=True)
            if c['status']!='OK':raise SystemExit(1)
        v['binary_sha256']=sha(build/'iCC');write(RUN/'versions.json',versions)
    progress('prepare-versions','PASS','Five isolated builds complete; U untouched, U17 precision-only, B17 pre-instrumentation, C pre-controller instrumentation, D current')

if __name__=='__main__':prepare()
