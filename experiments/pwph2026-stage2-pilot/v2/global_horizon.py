#!/usr/bin/env python3
"""One additional pre-final hypothesis: G to the terminal horizon, two windows.
Collect matching counterfactual labels instead of extrapolating from short h.
"""
from study import *

def collect(run,binary,split):
    rows=[]
    for item in plan(run)[split]:
        data=dataset(run,item)
        for d in (2,3):
            branches={}
            for q in ([0,1] if item['seed']%2 else [1,0]):
                out,c=execute(run,f'global-cal-{split}-{item["id"]}-d{d}-q{q}',command(binary,data,[1.2,2.],maxdim=d,
                    extra=['--cost-telemetry','--quotient-mask',str(q)]),[binary,data,Path(__file__)])
                need(c['status']=='OK','Global calibration failed');s=parse_states((out/'stdout.log').read_text())
                branches[q]=dict(state=s[0],next=s[1],path=str(out),whole_run_s=c['elapsed_monotonic_s'])
            need(branches[0]['state']['state_hash']==branches[1]['state']['state_hash'],'Global replay prefix differs')
            costs={q:float(branches[q]['state']['q_s'])+float(branches[q]['next']['work_s']) for q in (0,1)}
            rows.append(dict(dataset=item['id'],split=split,dimension=d,pair='global',state=branches[0]['state'],branches=branches,
                             measured_continue_s=costs[0],measured_quotient_s=costs[1],G_actual_s=costs[0]-costs[1],prefix_equivalent=True,
                             PH_status='COST_PAIR_NOT_A_PH_PASS'))
    write(run/f'calibration-global-{split}.json',rows)

def fit_global(run):
    train=json.loads((run/'calibration-dev.json').read_text())+json.loads((run/'calibration-global-dev.json').read_text())
    validation=json.loads((run/'calibration-validation.json').read_text())+json.loads((run/'calibration-global-validation.json').read_text())
    selected={};metrics=[]
    for d in (2,3):
        X=[];R=[];F=[];M=[];QX=[];QY=[]
        for p in train:
            if p['dimension']!=d:continue
            for q in (0,1):
                s=p['state'];b=p['branches'][str(q)];X.append(features(s,d,q));R.append(float(b['next']['rebuild_s']))
                F.append(float(b['next']['residual_s'])+(float(b['state']['q_s']) if q==0 else 0));M.append(float(b['next']['peak_kib']))
                if q:QX.append([1,float(s['raw'])/100,float(s['entries'])/1000]);QY.append(float(b['state']['q_s']))
        model=dict(R=ridge(X,R,.01).tolist(),F=ridge(X,F,.01).tolist(),memory=ridge(X,M,.01).tolist(),Q=ridge(QX,QY,.01).tolist(),
                   dimension=d,lambda_ridge=.01,max_A=max(float(p['state']['A']) for p in train if p['dimension']==d),
                   max_C=max(float(p['state']['C']) for p in train if p['dimension']==d),purpose='Terminal-horizon variant, calibrated including terminal branches; parameters fixed from previous validation')
        name=f'global-d{d}';write(run/'models'/(name+'.json'),model)
        (run/'models'/(name+'.txt')).write_text(f'PWPH_G_V1 {d} {model["max_A"]:.17g} {model["max_C"]:.17g}\n'+'\n'.join(' '.join(format(x,'.17g') for x in model[k]) for k in ('R','F','memory','Q'))+'\n');selected[str(d)]=name
        rows=[]
        for p in validation:
            if p['dimension']!=d:continue
            c=predict(model,p['state'],d,0);q=predict(model,p['state'],d,1);g=c['total']-q['total'];actual=p['G_actual_s']
            rows.append(dict(dataset=p['dataset'],pair=str(p['pair']),dimension=d,ridge=.01,G_pred_s=g,G_actual_s=actual,absolute_error_s=abs(g-actual),correct_sign=(g>0)==(actual>0),strong=abs(actual)>.0005,continue_pred_s=c['total'],continue_actual_s=p['measured_continue_s'],quotient_pred_s=q['total'],quotient_actual_s=p['measured_quotient_s']))
        globalrows=[p for p in rows if p['pair']=='global']
        metrics.append(dict(dimension=d,model=name,G_MAE_s=statistics.mean(x['absolute_error_s'] for x in rows),sign_accuracy=statistics.mean(x['correct_sign'] for x in rows),
                            terminal_G_MAE_s=statistics.mean(x['absolute_error_s'] for x in globalrows),terminal_sign_accuracy=statistics.mean(x['correct_sign'] for x in globalrows),samples=len(rows),terminal_samples=len(globalrows),rows=rows))
    write(run/'global-models.json',selected);write(run/'estimator-global-validation.json',metrics);print(json.dumps([{k:v for k,v in x.items() if k!='rows'} for x in metrics],indent=2))

def screen_global(run,binary,split):
    models=json.loads((run/'global-models.json').read_text());results=[]
    for item in plan(run)[split]:
        data=dataset(run,item)
        for d in (2,3):
            model=run/'models'/(models[str(d)]+'.txt')
            r=check(run,f'screen-global-{split}-{item["id"]}-d{d}',binary,data,[1.2,2.],maxdim=d,extra=['--g-model',str(model),'--g-margin','.05'])
            r.update(method='gs5',split=split,n=item['n'],geometry=item['geometry'],dataset_id=item['id']);results.append(r);write(Path(r['path'])/'ph-check.json',r)
    write(run/f'screen-global-{split}.json',results)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['collect','fit','screen']);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--split',default='dev',choices=['dev','validation']);a=ap.parse_args();guard=runner_lock(a.run)
    if a.stage=='collect':collect(a.run,a.binary.resolve(),a.split)
    elif a.stage=='fit':fit_global(a.run)
    else:screen_global(a.run,a.binary.resolve(),a.split)
