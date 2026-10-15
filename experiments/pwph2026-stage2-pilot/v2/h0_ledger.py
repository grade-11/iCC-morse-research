#!/usr/bin/env python3
"""Conditional H0 endpoint ledger from the ACTUAL per-window MSF masks.

With positive distances and 0<alpha<1, every accepted PV is a connected clique
at formation. Collapsing connected regions preserves components, and min-weight
quotient edges preserve all later component adjacencies. Thus later MSF deaths
>= previous strict cutoff are exactly the new component mergers, not all local
MSF bars. Each window checks component conservation, then the global ledger is
compared exactly to baseline MSF and float32 Ripser. No independent union-find
is used as proof of the PwPH execution path.
"""
from lab import *

def ledger(frames,windows,baseline,reference,alpha=.85):
    if not 0<alpha<1:return dict(status='H0_NOT_VERIFIED',reason='Proof scope requires strict alpha<1')
    if any(d==0 for b,d in baseline):return dict(status='KNOWN_DIFFERENCE',reason='Zero-distance/zero-length H0 convention; excluded from positive-distance proof')
    need(len(frames)==len(windows),'H0 frames/window mismatch');out=Counter();previous=None;n=int(frames[0]['metadata']['vertices']);steps=[]
    for i,(f,bounds) in enumerate(zip(frames,windows)):
        lo,hi=map(Decimal,bounds);components=int(f['metadata']['components'])
        new=Counter({(b,d):count for (b,d),count in f['pairs'].items() if d!=-1 and (i==0 or d>=lo)})
        count=sum(new.values());expected=n-components if i==0 else previous-components
        if count!=expected:return dict(status='H0_NOT_VERIFIED',reason='Actual MSF component conservation failed',window=i,new=count,expected=expected)
        out.update(new);steps.append(dict(window=i,lo=str(lo),hi=str(hi),components=components,new_mergers=count));previous=components
    out[(Decimal(0),Decimal(-1))]+=previous
    need(sum(out.values())==n,'H0 global count failed')
    exact=out==baseline;rp=projected(out)==reference
    return dict(status='PASS' if exact and rp else 'H0_NOT_VERIFIED',baseline_exact=exact,ripser_float32=rp,steps=steps,intervals=serial(out),
                scope='Positive distances, alpha in (0,1), existing disjoint/frozen PV rules; endpoint ledger only, not a general H1/H2 repair')

def audit_saved(run):
    reports=[]
    for split in ('dev','validation'):
        reports+=json.loads((run/f'screen-{split}.json').read_text())
    results=[]
    for r in reports:
        if 'diagnostic' not in r:continue
        frames=msf((Path(r['diagnostic'])/'stdout.log').read_text());data=Path(r['dataset'])
        matches=[]
        for p in (run/'cases').glob('base-msf-'+r['dataset_sha256'][:12]+'-*/context.json'):
            c=json.loads(p.read_text())
            if r['binary_sha256'] in c['inputs'].values() and c['status']=='OK':matches.append(p.parent)
        need(matches,'No actual baseline MSF');baseline=msf((matches[0]/'stdout.log').read_text())[0]['pairs']
        refs=reference(run,data,r['maxdim'],r['cutoff'])
        result=ledger(frames,r['windows'],baseline,refs[0]);result.update(case=r['case'],path=r['path'],previous_H0=r['H0'])
        results.append(result)
    write(run/'h0-ledger-audit.json',results)
    print(json.dumps(dict(configurations=len(results),statuses=dict(Counter(x['status'] for x in results)),newly_verified=sum(x['status']=='PASS' and x['previous_H0']!='PASS' for x in results)),indent=2))

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();audit_saved(a.run)
