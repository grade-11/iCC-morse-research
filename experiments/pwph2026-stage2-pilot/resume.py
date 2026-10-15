#!/usr/bin/env python3
"""Resume this pilot within its original six-hour budget, preserving evidence."""
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import sys
from pilot import BASE,sha,now
sys.path.insert(0,str(BASE.parent/'pwph2026-s2-n100'))
from evidence_run import capture

run=Path((BASE/'CURRENT-RUN.txt').read_text().strip())
context=json.loads((run/'initial-context.json').read_text())
remaining=6*3600-(datetime.now(timezone.utc)-datetime.fromisoformat(context['started_at'])).total_seconds()-600
if remaining<=0:raise SystemExit('Original budget ended; inspect completed evidence. A new experiment budget requires a new authorized run.')
lock=json.loads((run/'LOCK.json').read_text());binary=(run/'build/iCC').resolve()
if sha(binary)!=lock['binary_sha256'] or sha(BASE/'pilot.py')!=lock['pilot_sha256']:
    raise SystemExit('Binary/harness changed; do not mix batches. Inspect and create a new versioned run.')
os.environ.update(MPLCONFIGDIR='/tmp/icc-stage2-mpl',OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',OMP_THREAD_LIMIT='4')
code='''import sys;from pathlib import Path;sys.path.insert(0,sys.argv[1]);from pilot import batch,DEV,HELD
run=Path(sys.argv[2]);binary=Path(sys.argv[3])
for n in (100,200):
 for seeds in (DEV,HELD):
  methods=['ph','f2','f15','n100','f4','f5','adaptive'] if seeds==DEV else ['ph','f2','n100','adaptive']
  batch(run,binary,n,seeds,methods,5)
'''
dest=run/('resume-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
rc=capture(dest,[sys.executable,'-c',code,str(BASE),str(run),str(binary)],
           [binary,BASE/'pilot.py',run/'LOCK.json'],remaining)
print(json.dumps(dict(evidence=str(dest),exit_status=rc)))
raise SystemExit(rc)
