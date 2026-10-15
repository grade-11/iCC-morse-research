#!/usr/bin/env python3
"""Standalone scientific figures from saved tables, no additional experiments."""
from summarize import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,FormatStrFormatter,FuncFormatter

def plots(run):
    rows=[x for x in csv.DictReader((run/'paired-ratios.csv').open()) if int(x['dimension'])==2 and x['reference']=='best_fixed']
    fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True)
    for ax,kind in zip(axes,('time','rss')):
        mid=np.array([float(x[f'median_{kind}_ratio']) for x in rows]);low=np.array([float(x[f'min_{kind}_ratio']) for x in rows]);high=np.array([float(x[f'max_{kind}_ratio']) for x in rows]);y=np.arange(len(rows))
        ax.errorbar(mid,y,xerr=[mid-low,high-mid],fmt='o',capsize=3);ax.axvline(1,color='.5',linestyle='--');ax.xaxis.set_major_locator(MaxNLocator(5));ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'));ax.grid(axis='x',alpha=.2)
        ax.set_xlabel('Fixed / G '+('elapsed time' if kind=='time' else 'peak RSS')+' ratio\n>1 favors G; median and min-max across five pairs')
    axes[0].set_yticks(range(len(rows)),[x['dataset'].replace('202610','') for x in rows]);axes[0].invert_yaxis()
    fig.suptitle('Final H0-H1: eligible complete pairs only; best fixed selected before final');fig.tight_layout();fig.savefig(run/'final-paired-ratios.png',dpi=180);fig.savefig(run/'final-paired-ratios.pdf');plt.close(fig)
    metrics=json.loads((run/'estimator-common-validation.json').read_text());fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,d in zip(axes,(2,3)):
        bound=.001
        for m in [x for x in metrics if x['dimension']==d and x['scope']=='all']:
            a=np.array([r['G_actual_s'] for r in m['rows']]);p=np.array([r['G_pred_s'] for r in m['rows']]);bound=max(bound,max(abs(a)),max(abs(p)));ax.scatter(a,p,label=m['model'],alpha=.65,s=20)
        ax.plot([-bound,bound],[-bound,bound],color='.5',linestyle='--');ax.axhline(0,color='.8');ax.axvline(0,color='.8');ax.set_xscale('symlog',linthresh=.0005);ax.set_yscale('symlog',linthresh=.0005)
        ticks=[-.1,-.001,0,.001,.1] if d==2 else [-10,-.1,-.001,0,.001,.1,10]
        for axis in (ax.xaxis,ax.yaxis):
            axis.set_ticks(ticks);axis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
        ax.tick_params(labelsize=8)
        ax.set_title(f'Simplex max-dim {d}, same 40 validation pairs');ax.set_xlabel('Measured G phase difference (s)');ax.set_ylabel('Predicted G (s)');ax.legend(fontsize=8)
    fig.suptitle('Offline paired phase costs; not exact whole-process counterfactual costs');fig.tight_layout();fig.savefig(run/'G-validation.png',dpi=180);fig.savefig(run/'G-validation.pdf');plt.close(fig)
    print('Saved PNG/PDF performance and prediction figures.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args();plots(a.run)
