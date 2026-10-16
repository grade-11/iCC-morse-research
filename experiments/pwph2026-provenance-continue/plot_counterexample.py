#!/usr/bin/env python3
"""Scientific figure from the verified S2 n8 evidence; no new PH claim."""
from workflow import *
import numpy as np,itertools
os.environ.setdefault('MPLCONFIGDIR','/tmp/icc-provenance-plot')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def main():
    summary=json.loads((RUN/'counterexample8-summary.json').read_text());data=RUN/'inputs/two-squares-s2-8.csv';assert sha(data)==summary['input_sha256']
    points=np.loadtxt(data,delimiter=',');eps=.68;groups=[summary['PV_representatives'][str(i)] for i in (8,9)];colors=['#2563a6','#ba4a34'];bridges=[]
    intervals=summary['windows'][-1]['original_H1'];live=[x for x in intervals if float(x['birth'])<eps<float(x['death'])];assert len(live)==1 and live[0]['multiplicity']==1
    fig,axes=plt.subplots(1,3,figsize=(11.4,3.5),gridspec_kw=dict(width_ratios=[1.4,1,1.2]));ax=axes[0]
    for vertices,color in zip(groups,colors):
        center=points[vertices,:2].mean(axis=0);ordered=sorted(vertices,key=lambda i:np.arctan2(points[i,1]-center[1],points[i,0]-center[0]))
        ax.add_patch(Polygon(points[ordered,:2],facecolor=color,alpha=.13,edgecolor='none'))
        for i,j in itertools.combinations(vertices,2):
            assert np.linalg.norm(points[i]-points[j])<eps;ax.plot(points[[i,j],0],points[[i,j],1],color=color,lw=1,alpha=.45)
        ax.scatter(points[vertices,0],points[vertices,1],s=40,color=color,zorder=3)
        for i in vertices:ax.annotate(str(i),points[i,:2],xytext=(0,7),textcoords='offset points',ha='center',fontsize=9)
    for i in groups[0]:
        for j in groups[1]:
            if np.linalg.norm(points[i]-points[j])<eps:
                bridges.append([i,j]);ax.plot(points[[i,j],0],points[[i,j],1],color='#24805a',lw=2.4)
    assert len(bridges)==2
    ax.set_title('Original: two bridges, one H1 class',fontsize=10);ax.set_aspect('equal');ax.set_xlim(-.73,.73);ax.set_ylim(-.28,.30);ax.set_xlabel('XY projection; edges use 3D distance',fontsize=8);ax.set_yticks([]);ax.set_xticks([])
    ax=axes[1];ax.plot([0,1],[0,0],color='#52606d',lw=2.4)
    for x,color,label in zip([0,1],colors,['PV8','PV9']):ax.scatter([x],[0],s=150,color=color,zorder=3);ax.annotate(label,(x,0),xytext=(0,15),textcoords='offset points',ha='center',fontsize=10)
    ax.text(.5,-.24,'Two labels, one simple edge\nH1 = 0',ha='center',va='center',fontsize=10);ax.set_xlim(-.3,1.3);ax.set_ylim(-.5,.5);ax.axis('off');ax.set_title('Quotient viewed at epsilon = 0.68\n(PVs formed at 0.5)',fontsize=10)
    ax=axes[2];birth,death=float(live[0]['birth']),float(live[0]['death']);ax.plot([birth,death],[1,1],lw=5,color='#24805a',solid_capstyle='butt');ax.scatter([(birth+death)/2],[0],marker='x',s=50,color='#ba4a34');ax.text((birth+death)/2,.13,'missing',ha='center',fontsize=9,color='#ba4a34');ax.axvline(.5,color='#697586',ls='--',lw=1);ax.text(.505,1.32,'quotient',fontsize=8);ax.set_yticks([0,1],['PwPH','Original']);ax.set_xlim(.47,.76);ax.set_ylim(-.4,1.6);ax.set_xlabel('epsilon');ax.set_title('Future interval lost',fontsize=10);ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Verified unit-S2 counterexample: n=8, F2, alpha=0.85, max-simplex=2',fontsize=12)
    fig.text(.5,.015,f'At epsilon=0.68: original H1 rank 1, quotient H1 rank 0. Lost interval ({birth:.9f}, {death:.9f}).',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.08,1,.91]);outputs=[]
    for ext in ('png','pdf'):
        dst=RUN/('counterexample8.'+ext);fig.savefig(dst,dpi=180,bbox_inches='tight');outputs.append(dict(file=dst.name,sha256=sha(dst)))
    plt.close(fig);write(RUN/'counterexample8-figure.json',dict(input_sha256=sha(data),source_sha256=sha(__file__),epsilon=eps,bridges_original=bridges,original_H1_rank=1,quotient_H1_rank=0,verified_evidence='counterexample8-summary.json',outputs=outputs,scope='Illustration of already verified case only; projected geometry, original 3D distances used.'))
    progress('counterexample-figure','COMPLETE','PNG/PDF generated from verified3D input; no new PH/model experiment')

if __name__=='__main__':main()
