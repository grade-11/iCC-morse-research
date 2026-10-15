#!/usr/bin/env python3
"""Run two actual RSS caps only if a real pre-final state predicts changed action."""
from study import *

def probe(run,binary):
    reports=json.loads((run/'screen-dev.json').read_text())+json.loads((run/'screen-global-dev.json').read_text());candidates=[]
    for report in reports:
        if report['method'] not in GRID or report['status'] not in ('PASS','H0_NOT_VERIFIED','PH_MISMATCH'):continue
        text=(Path(report['path'])/'stdout.log').read_text();gs=[dict(x.split('=') for x in l.split()[1:]) for l in text.splitlines() if l.startswith('[g]')]
        for idx,g in enumerate(gs):
            negative=g['action']=='DEFER' and g['reason']=='G-below-margin' and float(g['G_s'])<-.0005
            positive=idx==0 and g['action']=='QUOTIENT' and float(g['G_s'])>float(g['threshold_s'])+.0005
            if not (negative or positive):continue
            if any(p['reason'] in ('memory-constraint','memory-reject-quotient') for p in gs[:idx]):continue
            previous=max([float(p['mem_quotient_kib'] if p['action']=='QUOTIENT' else p['mem_continue_kib']) for p in gs[:idx]]+[0])
            floor=max(previous,float(g['mem_quotient_kib'] if negative else g['mem_continue_kib']));gap=float(g['mem_continue_kib'] if negative else g['mem_quotient_kib'])-floor
            if gap>512:candidates.append((gap,report,g,idx,floor))
    if not candidates:
        write(run/'ram-study.json',dict(status='NOT_RUN_NO_NONTRIVIAL_CANDIDATE',reason='No negative-G checkpoint with >512KiB memory gap and an unchanged feasible prefix; do not run budgets that only change labels',unit_memory_branch_test='CTest verified action switch at two budgets'));return
    gap,report,g,idx,floor=max(candidates,key=lambda p:(p[3]==0,p[0]))
    low=(floor+.25*gap)/.7/1024
    budgets={'low':low,'high':max(2*low,512.,1.2*float(g['mem_continue_kib'])/.7/1024)};root=run/'ram-study';root.mkdir(exist_ok=True);results=[]
    write(run/'ram-study-plan.json',dict(at=now(),dataset=report['dataset_id'],method=report['method'],dimension=report['maxdim'],budgets_mib=budgets,
        checkpoint_index=idx,expected='Two predicted actions under different RAM constraints (positive-G quotient may be rejected at low RAM); stop after next checkpoint, no full-PH/performance claim'))
    breaks,extra=method_options(run,report['method'],report['maxdim']);data=Path(report['dataset'])
    for label,mib in budgets.items():
        sub=root/label;sub.mkdir(exist_ok=True);b=budget(run);b['memory_limit_bytes']=min(b['memory_limit_bytes'],int(mib*1024**2));write(sub/'budget.json',b);write(sub/'initial-context.json',b)
        out,c=execute(sub,'budget-probe',command(binary,data,breaks,maxdim=report['maxdim'],extra=extra+['--g-memory-mib',str(mib),'--cost-telemetry','--stop-after-window',str(idx+2)]),[binary,data,Path(__file__)])
        states=parse_states((out/'stdout.log').read_text());results.append(dict(budget=label,budget_mib=mib,status=c['status'],exit_status=c['exit_status'],
            action=states[idx]['action'] if len(states)>idx else None,state_hash=states[idx]['state_hash'] if len(states)>idx else None,path=str(out),peak_tree_bytes=c['tree_rss_sampled_bytes'],PH_scope='PARTIAL_PREFIX_ONLY'))
    write(run/'ram-study.json',dict(status='COMPLETED',dataset=report['dataset_id'],dimension=report['maxdim'],results=results,
          action_changed=len({x['action'] for x in results})==2 and all(x['action'] for x in results),
          equivalent_initial_state=results[0]['state_hash']==results[1]['state_hash'] and results[0]['state_hash'] is not None))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);a=ap.parse_args();guard=runner_lock(a.run);probe(a.run,a.binary.resolve())
