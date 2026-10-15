#!/usr/bin/env python3
"""Supplement historical reports, preserving the original H0 verdicts/raw logs."""
from h0_ledger import *

def enhance(run):
    audit_saved(run)
    results={x['path']:x for x in json.loads((run/'h0-ledger-audit.json').read_text())}
    for split in ('dev','validation'):
        target=run/f'screen-h0-enhanced-{split}.json'
        if target.exists():continue
        rows=json.loads((run/f'screen-{split}.json').read_text())
        for row in rows:
            if row['path'] not in results:continue
            h=results[row['path']];row.update(original_H0=row['H0'],h0_ledger=h,h0_verifier_sha256=sha(BASE/'h0_ledger.py'))
            if h['status']=='PASS':
                row['H0']='PASS'
                if row['status']=='H0_NOT_VERIFIED':row['status']='PASS'
        write(target,rows)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();guard=runner_lock(a.run);enhance(a.run)
