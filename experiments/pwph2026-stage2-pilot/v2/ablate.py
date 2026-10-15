#!/usr/bin/env python3
"""One predeclared feature-group ablation; validation selection, never final."""
from study import *

def ablate(run):
    train=json.loads((run/'calibration-dev.json').read_text());val=json.loads((run/'calibration-validation.json').read_text())
    original=json.loads((run/'selected-models.json').read_text());metrics=json.loads((run/'estimator-validation.json').read_text());selected={};outputs=[]
    for d in (2,3):
        name=original[str(d)];base=json.loads((run/'models'/(name+'.json')).read_text());lam=base['lambda_ridge']
        X=[];R=[];F=[];M=[]
        for pair in train:
            if pair['dimension']!=d:continue
            for q in (0,1):
                x=features(pair['state'],d,q);x[4]=x[6]=0;X.append(x);b=pair['branches'][str(q)]
                R.append(float(b['next']['rebuild_s']));F.append(float(b['next']['residual_s'])+(float(b['state']['q_s']) if q==0 else 0));M.append(float(b['next']['peak_kib']))
        model=dict(base,R=ridge(X,R,lam).tolist(),F=ridge(X,F,lam).tolist(),memory=ridge(X,M,lam).tolist(),feature_ablation='Drop PV-fraction*work and action*raw-support-count; keep action*work and action intercept')
        pred=[]
        for pair in val:
            if pair['dimension']!=d:continue
            c=predict(model,pair['state'],d,0);q=predict(model,pair['state'],d,1);g=c['total']-q['total'];actual=pair['G_actual_s']
            pred.append(dict(dataset=pair['dataset'],pair=pair['pair'],dimension=d,ridge=lam,G_pred_s=g,G_actual_s=actual,absolute_error_s=abs(g-actual),correct_sign=(g>0)==(actual>0),strong=abs(actual)>.0005,continue_pred_s=c['total'],continue_actual_s=pair['measured_continue_s'],quotient_pred_s=q['total'],quotient_actual_s=pair['measured_quotient_s']))
        mae=statistics.mean(x['absolute_error_s'] for x in pred)
        denom=statistics.mean(x['measured_continue_s']+x['measured_quotient_s'] for x in val if x['dimension']==d)
        score=mae/denom;full=next(m for m in metrics if m['dimension']==d and m['ridge']==lam)
        newname=name+'-lean';write(run/'models'/(newname+'.json'),model)
        (run/'models'/(newname+'.txt')).write_text(f'PWPH_G_V1 {d} {model["max_A"]:.17g} {model["max_C"]:.17g}\n'+'\n'.join(' '.join(format(x,'.17g') for x in model[k]) for k in ('R','F','memory','Q'))+'\n')
        selected[str(d)]=newname if score<full['normalized_G_MAE'] else name
        outputs.append(dict(dimension=d,ridge=lam,G_MAE_s=mae,normalized_G_MAE=score,sign_accuracy=statistics.mean(x['correct_sign'] for x in pred),
                            strong_sign_accuracy=statistics.mean(x['correct_sign'] for x in pred if x['strong']) if any(x['strong'] for x in pred) else None,
                            selected=selected[str(d)],baseline_full_model=name,full_normalized_G_MAE=full['normalized_G_MAE'],rows=pred))
    write(run/'estimator-ablation.json',outputs);write(run/'active-models.json',selected)
    print(json.dumps(dict(active=selected,metrics=[{k:v for k,v in m.items() if k!='rows'} for m in outputs]),indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();ablate(a.run)
