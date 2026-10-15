#!/usr/bin/env python3
"""Final retained inputs vs fresh Ripser; no performance conclusions."""
import argparse
import json
from pathlib import Path
from pilot import execute,reference,parse_compact,command,projected,serial,write,index,sha

def verify(run,binary):
    for folder in sorted((run/'shrunk').iterdir()):
        result=json.loads((folder/'result.json').read_text());data=folder/'retained.csv'
        # Unique reference filename for each retained point set.
        named=folder/(folder.name+'-retained.csv')
        if not named.exists():named.write_bytes(data.read_bytes())
        refs=reference(run,named);diags={}
        for method in ('ph',result['method']):
            out,c=execute(run,'confirm-shrunk-'+folder.name+'-'+method,command(binary,named,method),[binary,named,Path(__file__)])
            if c['status']!='OK':raise RuntimeError('Retained case process failed: '+str(out))
            diags[method]=parse_compact((out/'stdout.log').read_text(),'ph' if method=='ph' else 'piecewise')[0]
        baseline=diags['ph'];piece=diags[result['method']]
        report=dict(case='confirm-shrunk-'+folder.name,status='PH_MISMATCH' if baseline!=piece else 'NOT_REPRODUCED',exit_status=0,
                    purpose='Confirm final retained data and classify against independent fresh Ripser; not a retry for performance',
                    data_sha256=sha(named),binary_sha256=sha(binary),
                    baseline_ripser_float32={str(d):projected(baseline[d])==refs[d] for d in (1,2)},
                    piecewise_ripser_float32={str(d):projected(piece[d])==refs[d] for d in (1,2)},
                    difference={str(d):dict(missing=serial(baseline[d]-piece[d]),extra=serial(piece[d]-baseline[d])) for d in (1,2)},
                    H0='NOT_USED_FOR_PERFORMANCE; deletion investigation of H1/H2 only')
        write(folder/'independent-confirmation.json',report);index(run,report)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',required=True,type=Path);ap.add_argument('--binary',required=True,type=Path)
    a=ap.parse_args();verify(a.run,a.binary.resolve())
