#!/usr/bin/env python3
import sys,os,json,time,subprocess,hashlib,fcntl,csv
from pathlib import Path
BASE=Path(__file__).resolve().parent;REPO=BASE.parents[1]
sys.path.insert(0,str(REPO/'experiments/pwph2026-stage2-pilot'))
import pilot
from evidence_run import sha,now
RUN=REPO/(BASE/'CURRENT-RUN.txt').read_text().strip()
def write(p,obj):pilot.write(p,obj)
def budget():return json.loads((RUN/'budget.json').read_text())
def remaining(pack=False):
 b=budget();assert Path('/proc/sys/kernel/random/boot_id').read_text().strip()==b['boot_id'];return b['deadline_monotonic' if pack else 'experiments_stop_monotonic']-time.monotonic()
def execute(name,argv,inputs=(),timeout=300):
 assert remaining()>0,'Experiment budget expired'
 return pilot.execute(RUN,name,argv,[*inputs,Path(__file__)],min(timeout,remaining()))
def progress(stage,status,detail=''):
 pilot.index(RUN,dict(case=stage,status=status,exit_status=0,at=now(),detail=detail))
def table(path,rows):
 if not rows:return
 with Path(path).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
