#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
rows=list(csv.DictReader((a.run/'performance.csv').open()));ids=list(dict.fromkeys(r['dataset'] for r in rows))
fig,axes=plt.subplots(1,2,figsize=(13,5));colors=['#657b83','#d69e2e','#2f855a','#111111','#805ad5']
labels=['A: không quotient','B: lịch tác giả','C: cố định tốt','D: G đủ chi phí','E: chạy lại lịch G']
for ax,metric,lo,hi,ylabel in [(axes[0],'median_s','min_s','max_s','Thời gian toàn tiến trình (giây)'),(axes[1],'median_rss_kib','min_rss_kib','max_rss_kib','Bộ nhớ thường trú cao nhất (MiB)')]:
    for j,m in enumerate('ABCDE'):
        selected=[next(r for r in rows if r['dataset']==id and r['method']==m) for id in ids]
        scale=1024 if metric.endswith('kib') else 1
        y=np.array([float(r[metric])/scale for r in selected]);low=np.array([float(r[lo])/scale for r in selected]);high=np.array([float(r[hi])/scale for r in selected])
        ax.bar(np.arange(len(ids))+(j-2)*.15,y,.14,label=labels[j],color=colors[j],yerr=np.stack([y-low,high-y]),capsize=2,error_kw={'lw':.65})
    ax.set_xticks(np.arange(len(ids)),['64 đều','64 cụm','100 đều','100 cụm','200 đều','200 cụm'],rotation=25,ha='right');ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
axes[0].legend(fontsize=8);fig.suptitle('S² đơn vị, H0–H1, F2, cutoff 2; trung vị và khoảng nhỏ nhất–lớn nhất, 5 lượt')
fig.tight_layout();fig.savefig(a.run/'performance.png',dpi=180);fig.savefig(a.run/'performance.pdf');plt.close(fig)
