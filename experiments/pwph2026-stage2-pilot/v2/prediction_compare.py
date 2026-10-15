#!/usr/bin/env python3
"""Compare already fitted estimators on the same validation labels; no refit."""
from study import *

def compare(run):
    pairs=json.loads((run/'calibration-validation.json').read_text())+json.loads((run/'calibration-global-validation.json').read_text());out=[]
    for d in (2,3):
        for name in (json.loads((run/'active-models.json').read_text())[str(d)],'global-d'+str(d)):
            model=json.loads((run/'models'/(name+'.json')).read_text());rows=[]
            for p in pairs:
                if p['dimension']!=d:continue
                g=predict(model,p['state'],d,0)['total']-predict(model,p['state'],d,1)['total'];actual=p['G_actual_s']
                rows.append(dict(dataset=p['dataset'],pair=str(p['pair']),G_actual_s=actual,G_pred_s=g,absolute_error_s=abs(g-actual),correct_sign=(g>0)==(actual>0)))
            for scope,subset in [('all',rows),('short',[r for r in rows if r['pair']!='global']),('terminal',[r for r in rows if r['pair']=='global'])]:
                out.append(dict(dimension=d,model=name,scope=scope,samples=len(subset),G_MAE_s=statistics.mean(r['absolute_error_s'] for r in subset),sign_accuracy=statistics.mean(r['correct_sign'] for r in subset),always_positive_accuracy=statistics.mean(r['G_actual_s']>0 for r in subset),rows=subset))
    write(run/'estimator-common-validation.json',out)
    print(json.dumps([{k:v for k,v in r.items() if k!='rows'} for r in out],indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();compare(a.run)
