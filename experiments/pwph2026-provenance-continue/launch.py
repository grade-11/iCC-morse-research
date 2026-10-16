#!/usr/bin/env python3
from workflow import *
from evidence_run import capture
name=sys.argv[1];script=BASE/(name+'.py')
assert script.is_file() and script.parent==BASE
lock=(RUN/'runner.lock').open('a+')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
for previous in RUN.glob('stage-*/context.json'):
    if 'exit_status' not in json.loads(previous.read_text()):
        raise SystemExit('REFUSED: unfinished stage metadata; inspect before resume: '+str(previous))
if (RUN/'completion.json').exists() and not (name=='package' and '--pack-only' in sys.argv):
    raise SystemExit('REFUSED: completed run; no new work/budget without a new explicit task')
left=remaining(pack=name in ('package','finalize'));assert left>0,'Original budget expired'
stamp=str(time.monotonic_ns())
out=RUN/'packaging-capture'/stamp if name=='package' and '--pack-only' in sys.argv else RUN/('finalization-capture-'+stamp) if name=='finalize' else RUN/('stage-'+name+'-'+stamp)
lock.seek(0);lock.truncate();lock.write(json.dumps(dict(pid=os.getpid(),stage=name,at=now())));lock.flush()
rc=capture(out,[sys.executable,str(script),*sys.argv[2:]],[Path(__file__),script,BASE/'workflow.py',BASE/'matrix.py'],max(1,left+10))
print(json.dumps(dict(stage=name,exit_status=rc,path=str(out))));raise SystemExit(rc)
