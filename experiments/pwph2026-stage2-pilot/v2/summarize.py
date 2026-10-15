#!/usr/bin/env python3
"""Read saved study evidence, including failures and unverified H0."""
import argparse,csv,re,statistics as st
from collections import defaultdict,Counter
from lab import *

def table(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def summarize(run):
    checks=[]
    for split in ('dev','validation'):
        enhanced=run/f'screen-h0-enhanced-{split}.json'
        p=enhanced if enhanced.exists() else run/f'screen-{split}.json'
        if p.exists():checks+=json.loads(p.read_text())
        p=run/f'screen-global-{split}.json'
        if p.exists():checks+=json.loads(p.read_text())
    if (run/'final-validation.json').exists():checks+=json.loads((run/'final-validation.json').read_text())
    flat=[dict(dataset=x['dataset_id'],split=x['split'],n=x['n'],geometry=x['geometry'],dimension=x['maxdim'],method=x['method'],
               status=x['status'],H0=x.get('H0','H0_NOT_VERIFIED'),alias_of=x.get('alias_of',''),path=x['path']) for x in checks]
    table(run/'correctness.csv',flat);write(run/'correctness.json',checks)
    measurements=[];groups=defaultdict(list)
    for p in sorted((run/'cases').glob('measure-final-*/ph-check.json')):
        ph=json.loads(p.read_text());c=json.loads((p.parent/'context.json').read_text())
        row=dict(dataset=ph['dataset_id'],dimension=ph['dimension'],n=ph['n'],method=ph['method'],round=ph['round'],warmup=ph['warmup'],status=ph['status'],
                 scope=ph['scope'],exit_status=c['exit_status'],elapsed_s=c['elapsed_monotonic_s'],rss_kib=int(c.get('resources',{}).get('max_rss_kib',-1)),path=str(p.parent))
        measurements.append(row)
        if ph['round']>0:groups[row['dataset'],row['dimension'],row['method']].append(row)
    summary=[]
    for key,rows in sorted(groups.items()):
        ok=len(rows)==5 and {r['round'] for r in rows}==set(range(1,6)) and all(r['status']=='PASS' for r in rows)
        t=[r['elapsed_s'] for r in rows];rss=[r['rss_kib'] for r in rows]
        summary.append(dict(dataset=key[0],dimension=key[1],method=key[2],status='PASS' if ok else 'INCOMPLETE_OR_FAILED',scope=rows[0]['scope'],
                            median_s=st.median(t),min_s=min(t),max_s=max(t),median_rss_kib=st.median(rss),min_rss_kib=min(rss),max_rss_kib=max(rss),rounds=len(rows)))
    # Alias identical original/best fixed configs without pretending new processes ran.
    for x in checks:
        if x['split']=='final' and x.get('alias_of'):
            key=x['dataset_id'],x['maxdim'],x['alias_of']
            if key in groups:
                groups[x['dataset_id'],x['maxdim'],x['method']]=groups[key]
                src=next(s for s in summary if (s['dataset'],s['dimension'],s['method'])==key)
                summary.append(dict(src,method=x['method']))
    pairs=[]
    for x in checks:
        if x['split']!='final' or x['method']!='G':continue
        for other in ('original','best_fixed','v1'):
            k=x['dataset_id'],x['maxdim'];a=groups.get((*k,'G'),[]);b=groups.get((*k,other),[])
            if len(a)!=5 or len(b)!=5 or any(r['status']!='PASS' for r in a+b):continue
            ar={r['round']:r for r in a};t=[r['elapsed_s']/ar[r['round']]['elapsed_s'] for r in b];m=[r['rss_kib']/ar[r['round']]['rss_kib'] for r in b]
            pairs.append(dict(dataset=k[0],dimension=k[1],reference=other,ratio='reference/G; >1 favors G',median_time_ratio=st.median(t),min_time_ratio=min(t),max_time_ratio=max(t),median_rss_ratio=st.median(m),min_rss_ratio=min(m),max_rss_ratio=max(m)))
    table(run/'measurements.csv',measurements);table(run/'performance.csv',summary);table(run/'paired-ratios.csv',pairs)
    aggregate=[]
    for dimension in (2,3):
        for ref in ('original','best_fixed','v1'):
            rs=[x for x in pairs if x['dimension']==dimension and x['reference']==ref]
            if rs:aggregate.append(dict(dimension=dimension,reference=ref,datasets=len(rs),
                    geometric_mean_time_ratio=float(np.exp(st.mean(np.log(x['median_time_ratio']) for x in rs))),
                    geometric_mean_rss_ratio=float(np.exp(st.mean(np.log(x['median_rss_ratio']) for x in rs))),
                    aggregation='GM across datasets of median paired reference/G ratios; >1 favors G'))
    table(run/'aggregate-ratios.csv',aggregate)
    metas=[json.loads(p.read_text()) for p in (run/'data').glob('*.json')]
    refs=[p for p in (run/'references').glob('*.json')]
    process_states=Counter();elapsed=0
    offline=Counter();offline_seconds=Counter()
    for p in (run/'cases').glob('*/context.json'):
        c=json.loads(p.read_text());process_states[c.get('status','INTERRUPTED')]+=1;elapsed+=c.get('elapsed_monotonic_s',0)
        if c.get('case','').startswith(('cal-','global-cal-','stability-')):
            kind='stability' if c['case'].startswith('stability-') else 'offline_branch_calibration'
            offline[kind]+=1;offline_seconds[kind]+=c.get('elapsed_monotonic_s',0)
    wrappers=[]
    for p in run.rglob('context.json'):
        parts=p.relative_to(run).parts
        if 'cases' in parts or parts[0] in {'build','build-trace','source','source-controller-initial','third-party'}:continue
        c=json.loads(p.read_text())
        if 'input_sha256' in c and (len(c.get('argv',[]))<2 or Path(c['argv'][1]).name not in {'summarize.py','audit.py','plots.py','report.py','package.py','resume.py'}):
            wrappers.append(dict(path=str(p),exit_status=c.get('exit_status'),elapsed_s=c.get('elapsed_monotonic_s')))
    supplemental={}
    for name in ('reproduction.json','regression-v2.json','defer-regression.json','new-shrink-summary.json'):
        if (run/name).exists():
            rows=json.loads((run/name).read_text());supplemental[name]=dict(configurations=len(rows),statuses=dict(Counter(x['status'] for x in rows)))
    ram_contexts=[json.loads(p.read_text()) for p in (run/'ram-study').glob('*/cases/*/context.json')]
    inventory=dict(unique_new_datasets=len(metas),datasets_by_split=dict(Counter(m['split'] for m in metas)),datasets_by_geometry=dict(Counter(m['geometry'] for m in metas)),
                   validation_configurations=len(flat),validation_statuses=dict(Counter(r['status'] for r in flat)),
                   unique_executed_validation_configurations=sum(not r['alias_of'] for r in flat),
                   final_measurement_processes=len(measurements),warmups=sum(r['warmup'] for r in measurements),measured=sum(not r['warmup'] for r in measurements),
                   measurement_statuses=dict(Counter(r['status'] for r in measurements)),process_states=dict(process_states),process_elapsed_sum_s=elapsed,
                   budget_elapsed_s=time.monotonic()-budget(run)['start_monotonic'],references=len(refs),
                   offline_processes=dict(offline),offline_elapsed_sum_s=dict(offline_seconds),
                   supplemental_checks=supplemental,RAM_partial_process_statuses=dict(Counter(x['status'] for x in ram_contexts)),
                   capture_wrappers_count=len(wrappers),capture_wrapper_exit_statuses=dict(Counter(str(x['exit_status']) for x in wrappers)),
                   capture_wrapper_failures=[x for x in wrappers if x['exit_status']!=0],
                   known_old_datasets='2 reduced inputs + six configurations on two old n100 CSVs + standard n100 regression; never counted as new datasets',
                   note='COST_PAIR_EQUIVALENT labels are not PH PASS; H0_NOT_VERIFIED not combined with exact H0-H1 claims')
    write(run/'inventory.json',inventory);print(json.dumps(inventory,indent=2))
    inventory_rows=[]
    statuses=('PASS','PH_MISMATCH','H0_NOT_VERIFIED','TIMEOUT','MEMORY_LIMIT','ERROR')
    for split in ('dev','validation','final'):
        subset=[r for r in flat if r['split']==split];counts=Counter(r['status'] for r in subset)
        inventory_rows.append(dict(split=split,unique_new_datasets=sum(m['split']==split for m in metas),logical_configurations=len(subset),
              distinct_executions=sum(not r['alias_of'] for r in subset),**{s:counts[s] for s in statuses}))
    table(run/'inventory.csv',inventory_rows)
    header=['split','unique_new_datasets','logical_configurations','distinct_executions',*statuses]
    md=['# Kiểm kê — các batch research mới','', '| '+' | '.join(header)+' |','|'+'|'.join('---' for _ in header)+'|']
    md += ['| '+' | '.join(str(r[k]) for k in header)+' |' for r in inventory_rows]
    md += ['',f"Final: {inventory['warmups']} warmup +{inventory['measured']} measured processes; {inventory['measurement_statuses']}.",
           f"Process states (main cases): {dict(process_states)}. RAM partial processes: {inventory['RAM_partial_process_statuses']}.",
           f"Offline calibration/stability: {dict(offline)}; elapsed sums {dict(offline_seconds)} seconds. These labels are NOT PH PASS.",
           f"Capture/build/orchestration wrappers: {inventory['capture_wrapper_exit_statuses']}; failures in inventory.json and TOOLING-ERRORS.json.",
           f"Wall-clock used through this inventory: {inventory['budget_elapsed_s']/60:.2f}/240 min; do not sum nested wrapper and child elapsed times.",
           'Old seed07–11/reduced inputs/standard n100 and synthetic convention fixtures are regression only, excluded from the new-dataset count.',
           'Point-deletion candidates are shrink attempts, excluded from research-dataset count. New reduced inputs are confirmations, not new held-out datasets.',
           f"Supplemental checks (separate from main table): {supplemental}.",
           'Initial 18 H0_NOT_VERIFIED reports remain unchanged; conditional-ledger supplements verify them. No H1/H2 failure is removed by that update.']
    (run/'INVENTORY.md').write_text('\n'.join(md)+'\n')
    pred=[];prediction_summary=[]
    for filename in ('estimator-validation.json','estimator-ablation.json','estimator-global-validation.json'):
        if not (run/filename).exists():continue
        for m in json.loads((run/filename).read_text()):
            rows=m['rows'];positive=sum(r['G_actual_s']>0 for r in rows)
            confusion=Counter(('actual_positive' if r['G_actual_s']>0 else 'actual_nonpositive')+'__'+('pred_positive' if r['G_pred_s']>0 else 'pred_nonpositive') for r in rows)
            prediction_summary.append(dict(source=filename,dimension=m['dimension'],ridge=m.get('ridge',.01),samples=len(rows),G_MAE_s=st.mean(r['absolute_error_s'] for r in rows),
                    sign_accuracy=st.mean(r['correct_sign'] for r in rows),always_positive_accuracy=positive/len(rows),majority_sign_accuracy=max(positive,len(rows)-positive)/len(rows),confusion=dict(confusion)))
            for row in rows:
                pred.append(dict(source=filename,dataset=row['dataset'],pair=row['pair'],dimension=row['dimension'],ridge=m.get('ridge',.01),
                                 G_pred_s=row['G_pred_s'],G_actual_s=row['G_actual_s'],absolute_error_s=row['absolute_error_s'],correct_sign=row['correct_sign'],strong=row['strong']))
    table(run/'G-predictions.csv',pred)
    write(run/'G-prediction-summary.json',prediction_summary)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();summarize(a.run)
