#!/usr/bin/env python3
"""Exact argv thresholds in the H0 ledger as well as the PH parser."""
from matrix import *
from h0_ledger import ledger
import numpy as np

def main():
    data=RUN/'inputs/decimal-H0-boundary.csv';np.savetxt(data,[[0.,0.],[.3,0.],[0.,1.]],delimiter=',',fmt='%.17g')
    item=dict(id='decimal-H0-boundary',file=str(data),breaks=[.2,.3,.4,2.])
    base=one('S',item,2,4,'ph',verbose=True,tag='-H0-metadata-audit')
    q=one('S',item,2,4,verbose=True,extra=['--quotient-mask','0','0','0'],tag='-H0-metadata-audit')
    baseline=reviewed.msf((Path(base['path'])/'stdout.log').read_text())[0]['pairs'];refs=reviewed.reference(RUN,data,2)
    frames=reviewed.msf((Path(q['path'])/'stdout.log').read_text())
    oldwindows=[[str(a),str(b)] for a,b in zip([0]+item['breaks'][:-1],item['breaks'])]
    old=ledger(frames,oldwindows,baseline,refs[0]);correct=ledger(frames,q['windows'],baseline,refs[0])
    assert old['status']=='H0_NOT_VERIFIED' and correct['status']=='PASS'
    assert counters(base)==counters(q)
    checked=0
    for row in json.loads((RUN/'h0-verification.json').read_text()):
        c=json.loads((Path(row['path'])/'context.json').read_text());argv=c['argv'];start=argv.index('--eps-breaks')+1;grid=[]
        for v in argv[start:]:
            if v.startswith('--'):break
            grid.append(float(v))
        windows=[[format(a,'.17g'),format(b,'.17g')] for a,b in zip([0.]+grid[:-1],grid)]
        known=Counter({(Decimal(v['birth']),Decimal(v['death'])):v['multiplicity'] for v in row['ledger']['intervals']})
        reference=reviewed.projected(known)
        result=ledger(reviewed.msf((Path(row['path'])/'stdout.log').read_text()),windows,known,reference)
        assert result['status']=='PASS';checked+=1
    write(RUN/'metadata-H0-boundary-audit.json',dict(status='PASS',fixture_sha256=sha(data),raw_logs_unchanged=True,pair_endpoints_unchanged=True,old_shortest_decimal_bounds=old,exact_argv_binary64_bounds=correct,previous_H0_checks_reaudited_unchanged=checked,baseline_path=base['path'],piecewise_path=q['path'],scope='Parser/ledger metadata only, no native algorithm or input modification. Old results retained; source before fix saved.'))
    progress('H0-boundary-metadata','PASS','Exact argv binary64 bounds apply to verbose-return windows too; decimal edge death at0.3 regression; all210 previous H0 ledger checks unchanged')

if __name__=='__main__':main()
