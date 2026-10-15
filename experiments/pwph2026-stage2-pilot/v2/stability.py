#!/usr/bin/env python3
"""Two targeted validation pairs, two repeat labels, to quantify timing noise."""
from study import *

def stability(run,binary):
    metrics=json.loads((run/'estimator-validation.json').read_text());pairs=json.loads((run/'calibration-validation.json').read_text());selected=json.loads((run/'active-models.json').read_text());chosen=[];results=[]
    for d in (2,3):
        lam=json.loads((run/'models'/(selected[str(d)]+'.json')).read_text())['lambda_ridge']
        rows=next(m['rows'] for m in metrics if m['dimension']==d and m['ridge']==lam)
        worst=max(rows,key=lambda r:r['absolute_error_s']);original=next(p for p in pairs if p['dimension']==d and p['dataset']==worst['dataset'] and p['pair']==worst['pair'])
        chosen.append(dict(dimension=d,dataset=original['dataset'],pair=original['pair'],reason='Largest validation absolute G error for selected ridge; repeat for noise diagnosis, no refit'))
        item=next(x for x in plan(run)['validation'] if x['id']==original['dataset']);data=dataset(run,item);cp=plan(run)['calibration'][original['pair']];gains=[original['G_actual_s']]
        whole_gains=[original['branches']['0']['whole_run_s']-original['branches']['1']['whole_run_s']]
        for repeat in (1,2):
            costs={};whole={}
            for q in ([0,1] if repeat==1 else [1,0]):
                mask=cp['prefix']+[q]+[0]*(len(cp['grid'])-len(cp['prefix'])-2)
                out,c=execute(run,f'stability-{item["id"]}-d{d}-p{original["pair"]}-r{repeat}-q{q}',command(binary,data,cp['grid'],maxdim=d,
                    extra=['--cost-telemetry','--quotient-mask',*map(str,mask),'--stop-after-window',str(cp['index']+2)]),[binary,data,Path(__file__)])
                need(c['status']=='OK','Stability process failed');states=parse_states((out/'stdout.log').read_text());s=states[cp['index']]
                need(s['state_hash']==original['state']['state_hash'],'Stability prefix changed')
                costs[q]=float(s['q_s'])+float(states[cp['index']+1]['work_s'])
                whole[q]=c['elapsed_monotonic_s']
            gains.append(costs[0]-costs[1])
            whole_gains.append(whole[0]-whole[1])
        results.append(dict(**chosen[-1],G_actual_repeats_s=gains,G_pred_s=worst['G_pred_s'],mean_G_actual_s=statistics.mean(gains),min_G_actual_s=min(gains),max_G_actual_s=max(gains),stdev_G_actual_s=statistics.stdev(gains),
                            whole_process_replay_delta_s=whole_gains,whole_process_note='Includes startup and prefix timing noise; not an exact conditional cost from S'))
    write(run/'stability-plan.json',chosen);write(run/'stability.json',results);print(json.dumps(results,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);a=ap.parse_args();guard=runner_lock(a.run);stability(a.run,a.binary.resolve())
