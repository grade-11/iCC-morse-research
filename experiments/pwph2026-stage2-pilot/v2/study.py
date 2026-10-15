#!/usr/bin/env python3
"""Predeclared follow-up study; no final data are generated before LOCK-v2.json."""
import argparse,csv,math,statistics,os
from itertools import product
from lab import *

FIXED={'original':[1.2,1.7,2.], 'sparse':[1.2,2.], 'early':[.9,1.3,1.7,2.],
       'dense':[1.2,1.4,1.6,1.8,2.], 'single':[2.]}
GRID={'gs5':[1.2,2.], 'g0':[1.2,1.5,1.7,2.], 'g5':[1.2,1.5,1.7,2.],
      'gd5':[1.2,1.4,1.6,1.8,2.]}
V1=[1.2,1.4,1.5,1.6,1.7,1.8,2.]

def plan(run):
    p=run/'PLAN.json'
    if p.exists():return json.loads(p.read_text())
    def entry(split,n,geo,seed):return dict(split=split,n=n,geometry=geo,seed=seed,id=f'{geo}-n{n}-s{seed}')
    dev=[entry('dev',n,g,20261012+i+4*k) for i,n in enumerate([32,64,100,200]) for g in ['uniform_s2','polar_s2','multiscale'] for k in (0,1)]
    val=[entry('validation',n,g,s) for n,g,s in [(32,'uniform_s2',20261020),(64,'uniform_s2',20261021),(100,'uniform_s2',20261022),(200,'uniform_s2',20261023),(64,'polar_s2',20261020),(100,'polar_s2',20261021),(64,'multiscale',20261022),(100,'multiscale',20261023)]]
    final=[entry('final',n,g,s) for n,g,s in [(32,'uniform_s2',20261030),(64,'uniform_s2',20261031),(100,'uniform_s2',20261032),(200,'uniform_s2',20261033),(64,'polar_s2',20261034),(100,'polar_s2',20261035),(64,'multiscale',20261036),(100,'multiscale',20261037)]]
    obj=dict(created_at=now(),dev=dev,validation=val,final=final,dimensions=[2,3],
             fixed=FIXED,G_candidates=GRID,ridge_candidates=[.01,1.],
             selection='First maximize exact H1(+H2) coverage, then H0 verified coverage, then geometric mean runtime; per dimension and size small<=64 / large>=100.',
             calibration=[dict(grid=[.9,1.3,1.7,2.],index=0,prefix=[]),dict(grid=[.9,1.3,1.7,2.],index=1,prefix=[1]),dict(grid=[1.2,1.7,2.],index=0,prefix=[]),dict(grid=[1.2,1.5,1.7,2.],index=2,prefix=[1,1])])
    write(p,obj);return obj

def dataset(run,item):
    p=run/'data'/(item['id']+'.csv');meta=p.with_suffix('.json')
    if p.exists():need(json.loads(meta.read_text())['sha256']==sha(p),'Dataset changed');return p
    if item['split']=='final':need((run/'LOCK-v2.json').exists(),'Final data before design lock forbidden')
    sys.path.insert(0,str(REPO/'Sphere gen'));from dSphere_gen import dSphere_uniformRandom
    np.random.seed(item['seed']);n=item['n'];g=item['geometry'];description=''
    if g=='uniform_s2':
        pts=dSphere_uniformRandom(dimension=3,numPoints=n,r1=1.,r2=1.,trials=0,selection='min')
        description='Existing Muller generator, dimension=3, radii=1, trials=0, selection=min'
    elif g=='polar_s2':
        k=int(.7*n);polar=np.random.normal(size=(k,3))*.18+np.array([0,0,1.])
        polar/=np.linalg.norm(polar,axis=1)[:,None]
        uniform=dSphere_uniformRandom(dimension=3,numPoints=n-k,r1=1.,r2=1.,trials=0,selection='min')
        pts=np.concatenate([polar,uniform]);description='70% N((0,0,1),.18^2 I) normalized per point + 30% existing uniform S2; no global rescaling'
    else:
        centers=np.array([[-1.,0,0],[.5,.2,0],[0,.9,.1]]);group=np.random.choice(3,size=n,p=[.5,.3,.2])
        pts=centers[group]+np.random.normal(size=(n,3))*np.array([.035,.15,.35])[group,None]
        diameter=max(np.linalg.norm(pts[i]-pts[j]) for i in range(n) for j in range(i))
        pts*=1.8/diameter;description='3 Gaussian clusters centers=(-1,0,0),(.5,.2,0),(0,.9,.1); weights=.5/.3/.2; sigma=.035/.15/.35; global max pair distance normalized to1.8'
    p.parent.mkdir(exist_ok=True)
    with p.open('w',newline='') as f:csv.writer(f).writerows(pts)
    write(meta,dict(**item,sha256=sha(p),distribution=description,numpy=np.__version__,source_sha256=sha(__file__),
                   generator_sha256=sha(REPO/'Sphere gen/dSphere_gen.py'),utilities_sha256=sha(REPO/'Sphere gen/utilities.py'),author_data=False))
    note(run,'generate-'+item['id'],status='DATA_CREATED',exit_status=0,data_sha256=sha(p));return p

def calibration(run,binary,split):
    design=plan(run);rows=[]
    for item in design[split]:
        data=dataset(run,item)
        for d in design['dimensions']:
            for j,cp in enumerate(design['calibration']):
                branch={}
                # Alternate branch order; no performance comparison of Python wrappers.
                for q in ([0,1] if (j+item['seed'])%2==0 else [1,0]):
                    mask=cp['prefix']+[q]+[0]*(len(cp['grid'])-len(cp['prefix'])-2)
                    out,c=execute(run,f'cal-{split}-{item["id"]}-d{d}-p{j}-q{q}',command(binary,data,cp['grid'],maxdim=d,
                        extra=['--cost-telemetry','--quotient-mask',*map(str,mask),'--stop-after-window',str(cp['index']+2)]),[binary,data,Path(__file__)])
                    need(c['status']=='OK','Calibration process failed')
                    states=parse_states((out/'stdout.log').read_text());s=states[cp['index']];nxt=states[cp['index']+1]
                    branch[q]=dict(state=s,next=nxt,path=str(out),whole_run_s=c['elapsed_monotonic_s'])
                s0,s1=branch[0]['state'],branch[1]['state'];need(s0['state_hash']==s1['state_hash'],'Non-equivalent replay prefix')
                need(all(s0[k]==s1[k] for k in ('A','C','pv','raw','entries','eps','h')),'Prefix observables differ')
                measured={q:float(branch[q]['state']['q_s'])+float(branch[q]['next']['work_s']) for q in (0,1)}
                row=dict(dataset=item['id'],split=split,dimension=d,pair=j,state=s0,branches=branch,
                         measured_continue_s=measured[0],measured_quotient_s=measured[1],G_actual_s=measured[0]-measured[1],
                         prefix_equivalent=True,PH_status='COST_ONLY_PARTIAL_PREFIX_NOT_A_PH_PASS')
                rows.append(row)
                dest=run/'calibration';dest.mkdir(exist_ok=True);write(dest/f'{split}-{item["id"]}-d{d}-p{j}.json',row)
                note(run,'cal-pair-'+item['id']+f'-d{d}-p{j}',status='COST_PAIR_EQUIVALENT',exit_status=0)
    write(run/f'calibration-{split}.json',rows)

def features(s,d,q):
    a,c,pv,raw,eps,h=[float(s[k]) for k in ('A','C','pv','raw','eps','h')];growth=math.log(h/eps);w=c/1e5
    return [1,w,w*growth,(a/100)**(d+1)*growth,pv/max(1,a)*w,q*w,q*raw/100,q]
def ridge(X,y,lam):
    X=np.array(X,float);y=np.array(y,float);mean=X.mean(0);sd=X.std(0);mean[0]=0;sd[0]=1;sd[sd<1e-12]=1
    Z=(X-mean)/sd;reg=np.eye(X.shape[1])*lam;reg[0,0]=0
    beta=np.linalg.solve(Z.T@Z+reg,Z.T@y);coef=beta/sd;coef[0]-=np.dot(coef[1:],mean[1:]);return coef
def predict(model,s,d,q):
    x=np.array(features(s,d,q));r=max(0,float(x@model['R']));f=max(0,float(x@model['F']));mem=max(float(s['peak_kib']),float(x@model['memory']))
    select=max(0,float(np.array([1,float(s['raw'])/100,float(s['entries'])/1000])@model['Q'])) if q else 0
    return dict(R=r,F=f,Q=select,memory=mem,total=r+f+select)
def fit(run):
    dev=json.loads((run/'calibration-dev.json').read_text());val=json.loads((run/'calibration-validation.json').read_text());allmetrics=[];selected={}
    dest=run/'models';dest.mkdir(exist_ok=True)
    for d in (2,3):
        train=[p for p in dev if p['dimension']==d];validation=[p for p in val if p['dimension']==d];options=[]
        for lam in (.01,1.):
            X=[];R=[];F=[];M=[];QX=[];QY=[]
            for p in train:
                for q in (0,1):
                    s=p['state'];b=p['branches'][str(q)];X.append(features(s,d,q));R.append(float(b['next']['rebuild_s']))
                    F.append(float(b['next']['residual_s'])+(float(b['state']['q_s']) if q==0 else 0));M.append(float(b['next']['peak_kib']))
                    if q:QX.append([1,float(s['raw'])/100,float(s['entries'])/1000]);QY.append(float(b['state']['q_s']))
            model=dict(R=ridge(X,R,lam).tolist(),F=ridge(X,F,lam).tolist(),memory=ridge(X,M,lam).tolist(),Q=ridge(QX,QY,lam).tolist(),
                       dimension=d,lambda_ridge=lam,max_A=max(float(p['state']['A']) for p in train),max_C=max(float(p['state']['C']) for p in train))
            predrows=[]
            for p in validation:
                c=predict(model,p['state'],d,0);q=predict(model,p['state'],d,1);g=c['total']-q['total'];actual=p['G_actual_s']
                predrows.append(dict(dataset=p['dataset'],pair=p['pair'],dimension=d,ridge=lam,G_pred_s=g,G_actual_s=actual,
                                     absolute_error_s=abs(g-actual),correct_sign=(g>0)==(actual>0),strong=abs(actual)>.0005,
                                     continue_pred_s=c['total'],continue_actual_s=p['measured_continue_s'],quotient_pred_s=q['total'],quotient_actual_s=p['measured_quotient_s']))
            mae=statistics.mean(p['absolute_error_s'] for p in predrows)
            score=mae/statistics.mean(p['measured_continue_s']+p['measured_quotient_s'] for p in validation)
            metrics=dict(dimension=d,ridge=lam,G_MAE_s=mae,normalized_G_MAE=score,sign_accuracy=statistics.mean(p['correct_sign'] for p in predrows),
                         strong_sign_accuracy=statistics.mean(p['correct_sign'] for p in predrows if p['strong']) if any(p['strong'] for p in predrows) else None,
                         samples=len(predrows),rows=predrows)
            modelname=f'g-d{d}-ridge{lam}';write(dest/(modelname+'.json'),model)
            text=f'PWPH_G_V1 {d} {model["max_A"]:.17g} {model["max_C"]:.17g}\n'+'\n'.join(' '.join(format(x,'.17g') for x in model[k]) for k in ('R','F','memory','Q'))+'\n'
            (dest/(modelname+'.txt')).write_text(text);options.append((score,modelname));allmetrics.append(metrics)
        selected[str(d)]=min(options)[1]
    write(run/'estimator-validation.json',allmetrics);write(run/'selected-models.json',selected)
    print(json.dumps(dict(selected=selected,metrics=[{k:v for k,v in m.items() if k!='rows'} for m in allmetrics]),indent=2))

def method_options(run,method,d,lock=None):
    if method in FIXED:return FIXED[method],[]
    if method=='v1':return V1,['--adaptive']
    pointer=run/('global-models.json' if method=='gs5' else 'active-models.json')
    models=json.loads((pointer if pointer.exists() else run/'selected-models.json').read_text());path=run/'models'/(models[str(d)]+'.txt')
    return GRID[method],['--g-model',str(path),'--g-margin','0' if method.endswith('0') else '.05']

def screen(run,binary,split):
    results=[]
    for item in plan(run)[split]:
        data=dataset(run,item)
        for d in (2,3):
            for method in [*FIXED,'v1',*GRID]:
                breaks,extra=method_options(run,method,d)
                r=check(run,f'screen-{split}-{item["id"]}-d{d}-{method}',binary,data,breaks,maxdim=d,extra=extra)
                r.update(method=method,split=split,n=item['n'],geometry=item['geometry'],dataset_id=item['id'])
                results.append(r);write(Path(r['path'])/'ph-check.json',r)
    write(run/f'screen-{split}.json',results)

def lock_design(run,binary):
    need(not (run/'LOCK-v2.json').exists(),'Design already locked')
    reports=sum([json.loads((run/f'screen-h0-enhanced-{s}.json').read_text())+json.loads((run/f'screen-global-{s}.json').read_text()) for s in ('dev','validation')],[])
    best={};scores={}
    for d in (2,3):
        for size in ('small','large'):
            group=[r for r in reports if r['maxdim']==d and ('small' if r['n']<=64 else 'large')==size]
            def choose(methods):
                ranked=[]
                for m in methods:
                    rs=[r for r in group if r['method']==m];valid=[r for r in rs if r['status'] in ('PASS','H0_NOT_VERIFIED')]
                    times=[json.loads((Path(r['path'])/'context.json').read_text())['elapsed_monotonic_s'] for r in valid]
                    score=(-len(valid),-sum(r['H0']=='PASS' for r in valid),math.exp(statistics.mean(math.log(t) for t in times)) if times else float('inf'))
                    ranked.append((score,m));scores[f'{d}-{size}-{m}']=dict(PH_positive_dimensions=len(valid),H0_verified=sum(r['H0']=='PASS' for r in valid),attempted=len(rs),geom_time_s=score[2])
                return min(ranked)[1]
            best[f'{d}-{size}']=dict(fixed=choose(FIXED),G=choose(GRID))
    modelpaths=list((run/'models').glob('*.txt'))
    obj=dict(locked_at=now(),lock_monotonic=time.monotonic(),binary_sha256=sha(binary),source_hashes={str(p):sha(p) for p in [Path(__file__),BASE/'lab.py',BASE/'h0_ledger.py',REPO/'deps/GCostController.hpp',run/'active-models.json',run/'global-models.json',run/'PLAN.json',run/'PLAN-AMENDMENT-2.json']},
             selected=best,global_estimators=json.loads((run/'global-models.json').read_text()),selected_estimators=json.loads((run/'active-models.json').read_text()),model_hashes={str(p):sha(p) for p in modelpaths},scores=scores,grids=GRID,
             rules=plan(run)['selection'],online_update='Common time scale EMA .2; clamp [.25,4]; no online counterfactual execution',final_plan=plan(run)['final'])
    write(run/'LOCK-v2.json',obj);print(json.dumps(obj,indent=2))

def final(run,binary):
    lock=json.loads((run/'LOCK-v2.json').read_text());need(sha(binary)==lock['binary_sha256'],'Locked binary changed')
    for path,digest in {**lock['source_hashes'],**lock['model_hashes']}.items():need(sha(path)==digest,'Locked source/model changed: '+path)
    results=[]
    for item in plan(run)['final']:
        data=dataset(run,item)
        for d in (2,3):
            group=lock['selected'][f'{d}-'+('small' if item['n']<=64 else 'large')]
            methods={'original':'original','best_fixed':group['fixed'],'v1':'v1','G':group['G']}
            checks={};seen={}
            for label,method in methods.items():
                if method in seen:
                    r=dict(checks[seen[method]],method=label,alias_of=seen[method]);checks[label]=r;results.append(r);continue
                breaks,extra=method_options(run,method,d)
                r=check(run,f'final-{item["id"]}-d{d}-{label}',binary,data,breaks,maxdim=d,extra=extra)
                r.update(method=label,implementation=method,split='final',n=item['n'],geometry=item['geometry'],dataset_id=item['id']);checks[label]=r;results.append(r);write(Path(r['path'])/'ph-check.json',r);seen[method]=label
            eligible=[m for m in methods if 'alias_of' not in checks[m] and checks[m]['status'] in ('PASS','H0_NOT_VERIFIED')];blocked=set()
            for round in range(6):
                order=eligible[round%len(eligible):]+eligible[:round%len(eligible)] if eligible else []
                if round%2:order.reverse()
                for label in order:
                    if label in blocked:continue
                    method=methods[label];breaks,extra=method_options(run,method,d)
                    out,c=execute(run,f'measure-final-{item["id"]}-d{d}-{label}-r{round}',command(binary,data,breaks,maxdim=d,extra=extra),[binary,data,Path(__file__),run/'LOCK-v2.json'])
                    ph=dict(status=c['status'],exit_status=c['exit_status'],round=round,warmup=round==0,method=label,implementation=method,dataset_id=item['id'],n=item['n'],dimension=d,
                            H0=checks[label]['H0'],scope='H0-H'+str(d-1) if checks[label]['H0']=='PASS' else 'H1'+('-H2' if d==3 else '')+'; H0_NOT_VERIFIED')
                    if c['status']=='OK':
                        diag,ign,windows=compact((out/'stdout.log').read_text(),'piecewise',d)
                        expected={int(k):Counter({(Decimal(row['birth']),Decimal(row['death'])):row['multiplicity'] for row in v}) for k,v in checks[label]['diagrams'].items()}
                        ph.update(status='PASS' if diag==expected else 'PH_MISMATCH',windows=[[str(a),str(b)] for a,b in windows])
                    if ph['status']=='PASS' and '--g-model' in extra:
                        observed_states=parse_states((out/'stdout.log').read_text())
                        validated_states=parse_states((Path(checks[label]['path'])/'stdout.log').read_text())
                        identity=lambda ss:[(s['state_hash'],s['action']) for s in ss]
                        if identity(observed_states)==identity(validated_states):
                            ph['H0_verification']='Same actual logical states/actions as validated MSF trajectory'
                        else:
                            mask=[int(s['action']=='QUOTIENT') for s in observed_states[:-1]]
                            replay_name=f'measure-trajectory-{item["id"]}-d{d}-mask'+''.join(map(str,mask))
                            replay=check(run,replay_name,binary,data,[float(b) for a,b in windows],maxdim=d,extra=['--quotient-mask',*map(str,mask)])
                            replay_states=parse_states((Path(replay['path'])/'stdout.log').read_text())
                            need(identity(observed_states)==identity(replay_states),'Measured G trajectory replay differs')
                            ph.update(H0=replay['H0'],H0_verification='Actual mask replay with state identity and MSF ledger',trajectory_check=replay['path'])
                            if replay['status']!='PASS':ph['status']=replay['status']
                        ph['scope']='H0-H'+str(d-1) if ph['H0']=='PASS' else 'H1'+('-H2' if d==3 else '')+'; H0_NOT_VERIFIED'
                    write(out/'ph-check.json',ph);note(run,ph['dataset_id']+f'-{label}-r{round}',**{k:v for k,v in ph.items() if k!='dataset_id'})
                    if ph['status']!='PASS':blocked.add(label)
    write(run/'final-validation.json',results)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['plan','calibrate','fit','screen','lock','final']);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--split',choices=['dev','validation'],default='dev')
    a=ap.parse_args();guard=runner_lock(a.run);plan(a.run)
    if a.stage=='calibrate':calibration(a.run,a.binary.resolve(),a.split)
    elif a.stage=='fit':fit(a.run)
    elif a.stage=='screen':screen(a.run,a.binary.resolve(),a.split)
    elif a.stage=='lock':lock_design(a.run,a.binary.resolve())
    elif a.stage=='final':final(a.run,a.binary.resolve())
