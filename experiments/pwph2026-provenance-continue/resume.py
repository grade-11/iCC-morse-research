#!/usr/bin/env python3
"""Resume within the original budget; never reset clock or overwrite raw cases."""
from workflow import *
import argparse

MARKERS={'matrix':'correctness.csv','instrumentation':'g-replay.csv',
 'h0_checks':'stream-native-regression.json','stream_checks':'stream-rejections.json',
 'counterexample8':'G-GRID01-GATE.json','survey':'survey-summary.json',
 'n200_dossier':'n200-dossier.json','provenance_summary':'cross-version-summary.json',
 'grid_foundation':'grid01-expanded-upstream.json','n300_h2':'n300-h2-summary.json',
 'g_boundary':'G-boundary-summary.json','bounds_h0_audit':'metadata-H0-boundary-audit.json',
 'evidence_audit':'evidence-audit.json'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=MARKERS);args=ap.parse_args()
    b=budget();sameboot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()==b['boot_id']
    status=dict(run=str(RUN),at=now(),same_boot=sameboot,remaining_experiment_s=remaining() if sameboot else None,completed=(RUN/'completion.json').exists(),stages={s:(RUN/m).exists() for s,m in MARKERS.items()})
    print(json.dumps(status,indent=2))
    if not args.stage:return
    if status['stages'][args.stage]:print('SKIP_COMPLETED_STAGE: '+args.stage);return
    if status['completed'] or not sameboot or status['remaining_experiment_s']<=0:raise SystemExit('REFUSED: completed run, expired original budget, or different boot; no automatic new budget')
    result=subprocess.run([sys.executable,str(BASE/'launch.py'),args.stage]);raise SystemExit(result.returncode)

if __name__=='__main__':main()
