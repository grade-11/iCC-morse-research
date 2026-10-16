#!/usr/bin/env python3
"""Corroborate newly observed fixed-grid n200 failure, not a G evaluation."""
from matrix import *

def main():
    rows=[]
    for n in (200,300):
      item=dict(id=f'standard{n}',file=str(RUN/'inputs'/f'standard{n}.csv'),breaks=[k/10 for k in range(1,21)])
      refs=reviewed.reference(RUN,Path(item['file']),2)
      for version in ('U','U17','D'):
        base=one(version,item,2,4,'ph',tag='-grid-foundation');q=one(version,item,2,4,tag='-grid-foundation')
        assert base['status']==q['status']=='OK'
        b=counters(base);x=counters(q)
        rows.append(dict(dataset=item['id'],version=version,baseline=base,quotient=q,exact_equal=b==x,status='PH_MISMATCH' if b!=x else 'PASS',reference='PRINTED_PRECISION_ONLY' if version=='U' else all(reviewed.projected(v)==refs[k] for k,v in b.items()),difference={str(k):dict(missing=reviewed.serial(b[k]-x[k]),extra=reviewed.serial(x[k]-b[k])) for k in b}))
    assert all(x['status']=='PH_MISMATCH' for x in rows)
    write(RUN/'grid01-expanded-upstream.json',rows);progress('grid01-expanded-upstream','PH_MISMATCH','Fixed all-quotient grid0.1 loses H1 on standard n200/n300 seed20261007, reproduced U/U17/D; no new G fit/evaluation')
    print('Upstream fixed-grid0.1 n200/n300 PH_MISMATCH corroborated')

if __name__=='__main__':main()
