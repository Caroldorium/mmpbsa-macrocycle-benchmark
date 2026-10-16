#!/usr/bin/env python3
"""SO3LR-SF vs MM/PBSA vs FEP+ on macrocyclic peptides, under the manuscript's estimators."""
import csv, itertools, json, random
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

import pathlib
HERE = pathlib.Path(__file__).resolve().parent
R = str(HERE.parent / 'results' / 'zha25_cohort_scores.csv')
SO3 = str(HERE / 'so3lr_sf_scores.csv')
so3={x['slug']:float(x['binding_kcal_mol']) for x in csv.DictReader(open(SO3))
     if x.get('status')=='ok' and x.get('binding_kcal_mol')}
rows=[]
for x in csv.DictReader(open(R)):
    if x.get('status')!='ok': continue
    try:
        if float(x['pb_2_VDWAALS'])>0: continue
        fep=x['published_fepplus_dG_kcal_mol']
        rows.append(dict(t=x['target'],dg=float(x['experimental_dG_kcal_mol']),
            nheavy=float(x['heavy_atoms']),h=-float(x['heavy_atoms']),pb=float(x['pb_2']),
            gb=float(x['gb2_1']),so3=so3.get(x['slug']),fep=float(fep) if fep else None,
            cens=x['censored'].strip().lower() in ('true','1')))
    except (ValueError,KeyError): continue
rows=[r for r in rows if r['so3'] is not None]
T=sorted({r['t'] for r in rows}); MARGIN=0.5
KEYS=[('pb','MM/PBSA','#B51F58'),('gb','MM/GBSA','#E69F00'),('so3','SO3LR-SF','#2E8B57'),
      ('fep','FEP+','#4C78A8'),('h','heavy atoms','#888888')]

def acc(sub,key,band):
    ok=tot=0
    for a,b in itertools.combinations(sub,2):
        if a['t']!=b['t'] or abs(a['dg']-b['dg'])<MARGIN or abs(a['nheavy']-b['nheavy'])>band: continue
        if a[key] is None or b[key] is None or a[key]==b[key]: continue
        tot+=1
        if ((-a[key])>(-b[key]))==(a['dg']<b['dg']): ok+=1
    return ok,tot
def macro(sub,key,band):
    v=[o/n for o,n in (acc([r for r in sub if r['t']==t],key,band) for t in T) if n>=5]
    return np.mean(v) if v else None
rng=random.Random(11)
def ci(pw,key,band,n=1500):
    m=macro(pw,key,band)
    if m is None: return None
    byt={t:[r for r in pw if r['t']==t] for t in T}; bs=[]
    for _ in range(n):
        rs=[]
        for t in T:
            s=byt[t]; rs+=[s[rng.randrange(len(s))] for _ in s]
        v=macro(rs,key,band)
        if v is not None: bs.append(v)
    return m,np.percentile(bs,2.5),np.percentile(bs,97.5)

pw=[r for r in rows if not r['cens']]
fig,axes=plt.subplots(1,3,figsize=(15.5,5.2))

# panel A: macro Spearman per target
ax=axes[0]; w=0.16; xs=np.arange(len(T))
for i,(k,lab,c) in enumerate(KEYS):
    v=[spearmanr([r[k] for r in rows if r['t']==t],[r['dg'] for r in rows if r['t']==t]).correlation for t in T]
    ax.bar(xs+(i-2)*w,v,w,color=c,label=lab,edgecolor='white',lw=.6)
ax.axhline(0,color='k',lw=.8); ax.set_xticks(xs)
ax.set_xticklabels([t.replace('_','/') for t in T],rotation=20,ha='right',fontsize=9)
ax.set_ylabel('Spearman $\\rho$ vs experiment\n(higher = better ranking)')
ax.set_title('A  Within-receptor ranking',loc='left',fontweight='bold')
ax.legend(fontsize=8,ncol=2,frameon=False); ax.set_ylim(-0.45,1.05)

# panel B: macro mean
ax=axes[1]
ms=[np.mean([spearmanr([r[k] for r in rows if r['t']==t],
    [r['dg'] for r in rows if r['t']==t]).correlation for t in T]) for k,_,_ in KEYS]
ax.barh(range(len(KEYS)),ms,color=[c for _,_,c in KEYS],edgecolor='white')
for i,m in enumerate(ms): ax.text(m+.012,i,f'{m:.3f}',va='center',fontsize=9)
ax.set_yticks(range(len(KEYS))); ax.set_yticklabels([l for _,l,_ in KEYS]); ax.invert_yaxis()
ax.set_xlabel('macro Spearman $\\rho$'); ax.set_xlim(0,1.0)
ax.axvline(ms[4],color='#888',ls=':',lw=1.2)
ax.text(ms[4]+0.015,4.42,'heavy-atom null',fontsize=8,color='#555',ha='left',style='italic')
ax.set_title('B  Averaged over the five targets',loc='left',fontweight='bold')

# panel C: size-matched pairwise accuracy
ax=axes[2]; w=0.36
for j,band in enumerate((0,3)):
    for i,(k,lab,c) in enumerate(KEYS):
        r=ci(pw,k,band)
        if r is None: continue
        m,lo,hi=r
        ax.bar(i+(j-0.5)*w,m,w,color=c,alpha=1.0 if band==3 else 0.52,
               edgecolor='white',lw=.6,yerr=[[m-lo],[hi-m]],capsize=2.5,
               error_kw=dict(lw=1,ecolor='#333'))
ax.axhline(0.5,color='k',ls='--',lw=1)
ax.text(len(KEYS)-0.4,0.507,'chance',fontsize=8,ha='right')
ax.set_xticks(range(len(KEYS))); ax.set_xticklabels([l for _,l,_ in KEYS],rotation=20,ha='right',fontsize=9)
ax.set_ylabel('macro pairwise accuracy'); ax.set_ylim(0.4,1.0)
ax.set_title('C  Size-matched pairs\n(pale $|\\Delta|$=0, solid $|\\Delta|\\leq$3)',loc='left',fontweight='bold')
for a in axes: a.spines[['top','right']].set_visible(False)
fig.suptitle('SO3LR-SF on macrocyclic peptides: the small-molecule advantage does not transfer',
             fontsize=13,fontweight='bold',x=0.01,ha='left')
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig('so3lr_comparison.png',dpi=200)
print('wrote so3lr_comparison.png')
