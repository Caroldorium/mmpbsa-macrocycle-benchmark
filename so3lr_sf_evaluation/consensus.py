#!/usr/bin/env python3
"""Do MM/PBSA and SO3LR-SF errors combine?

s3.6.1 observes the two methods fail on different targets and calls it worth following up.
This tests it: if the errors are orthogonal, a consensus should beat both components.  Ranks are
z-scored WITHIN each receptor before combining, since the two scores are on different scales.
Reported under the manuscript's own estimators, with the same compound-level bootstrap.
"""
import csv, itertools, random
import numpy as np
from scipy.stats import spearmanr, rankdata

R='/tmp/mmpbsa-macrocycle-benchmark/results/zha25_cohort_scores.csv'
so3={x['slug']:float(x['binding_kcal_mol']) for x in csv.DictReader(open('so3lr_full.csv'))
     if x.get('status')=='ok' and x.get('binding_kcal_mol')}
rows=[]
for x in csv.DictReader(open(R)):
    if x.get('status')!='ok': continue
    try:
        if float(x['pb_2_VDWAALS'])>0: continue
        fep=x['published_fepplus_dG_kcal_mol']
        rows.append(dict(t=x['target'],dg=float(x['experimental_dG_kcal_mol']),
            nheavy=float(x['heavy_atoms']),h=-float(x['heavy_atoms']),pb=float(x['pb_2']),
            gb=float(x['gb2_1']),so3=so3[x['slug']],fep=float(fep) if fep else None,
            cens=x['censored'].strip().lower() in ('true','1')))
    except (ValueError,KeyError): continue
T=sorted({r['t'] for r in rows}); MARGIN=0.5

# within-receptor rank-z, then consensus
for t in T:
    sub=[r for r in rows if r['t']==t]; n=len(sub)
    for k in ('pb','so3','gb'):
        z=(rankdata([r[k] for r in sub])-(n+1)/2)/max(np.std(rankdata([r[k] for r in sub])),1e-9)
        for r,v in zip(sub,z): r['z_'+k]=v
    for r in sub:
        r['cons']=(r['z_pb']+r['z_so3'])/2          # equal-weight MM/PBSA + SO3LR-SF
        r['cons_gb']=(r['z_pb']+r['z_gb'])/2        # control: two MM methods

print("=== error orthogonality ===")
for a,b in (('pb','so3'),('pb','gb'),('so3','gb')):
    rr=[np.mean([spearmanr([x[k] for x in rows if x['t']==t],
        [x['dg'] for x in rows if x['t']==t]).correlation for t in T]) for k in (a,b)]
    # correlation of residual errors, within target
    es=[]
    for t in T:
        sub=[r for r in rows if r['t']==t]
        ra=rankdata([r[a] for r in sub]); rb=rankdata([r[b] for r in sub]); rd=rankdata([r['dg'] for r in sub])
        es.append(spearmanr(ra-rd, rb-rd).correlation)
    print(f"  {a:4s} vs {b:4s}: rank-error correlation {np.mean(es):+.3f}   (component rho {rr[0]:.3f}, {rr[1]:.3f})")

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
pw=[r for r in rows if not r['cens']]
def ci(key,band,n=1500):
    m=macro(pw,key,band)
    byt={t:[r for r in pw if r['t']==t] for t in T}; bs=[]
    for _ in range(n):
        rs=[]
        for t in T:
            s=byt[t]; rs+=[s[rng.randrange(len(s))] for _ in s]
        v=macro(rs,key,band)
        if v is not None: bs.append(v)
    bs=np.array(bs); return m,np.percentile(bs,2.5),np.percentile(bs,97.5)

KEYS=[('pb','MM/PBSA'),('so3','SO3LR-SF'),('cons','consensus PBSA+SO3LR'),
      ('cons_gb','control: PBSA+GBSA'),('fep','FEP+')]
print("\n=== macro Spearman ===")
for k,lab in KEYS:
    v=[spearmanr([r[k] for r in rows if r['t']==t],[r['dg'] for r in rows if r['t']==t]).correlation for t in T]
    print(f"  {lab:<24}{np.mean(v):>7.3f}   " + " ".join(f"{t[:4]}={x:+.2f}" for t,x in zip(T,v)))
for band in (0,3):
    print(f"\n=== pairwise accuracy, band {band} ===")
    for k,lab in KEYS:
        m,lo,hi=ci(k,band)
        print(f"  {lab:<24}{m:>7.3f}  [{lo:.3f}, {hi:.3f}]")
