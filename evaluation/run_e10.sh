#!/bin/bash
# E10 after the lock: generation, preparation, timing, correctness check, locked analysis.
# Progress in run_e10.status; logs in generate.log, prepare.log, timing.log, check.log, analyze.log.
P=/home/alex_1709/.venvs/icc-ripser-0.6.15-VnIso6/bin/python
cd /home/alex_1709/projects/iCC-team/experiments/pwph-e10-even-ladder-confirm-20261105
online() { timeout 60 powershell.exe -NoProfile -Command "(Get-CimInstance -Namespace root/wmi -ClassName BatteryStatus).PowerOnline" 2>/dev/null | tr -d '\r' | grep -q True; }
for step in generate prepare time check analyze; do
  if [ $step = time ]; then
    until online; do sleep 60; done
    echo "AC power online $(date -Is)" >> run_e10.status
    ( while [ ! -f TIMING-END.json ]; do echo "$(date -Is) $(online && echo AC || echo BATTERY)" >> power.log; sleep 60; done ) &
  fi
  echo "$step started $(date -Is)" >> run_e10.status
  $P e10_$step.py > $step.log 2>&1 || { echo "$step FAILED $(date -Is)" >> run_e10.status; exit 1; }
  echo "$step done $(date -Is)" >> run_e10.status
done
