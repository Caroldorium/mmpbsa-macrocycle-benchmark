#!/usr/bin/env python3
"""Macro size-matched pairwise accuracy -- the manuscript's estimator.

Pairs are formed WITHIN a target only (a KRAS-vs-PCSK9 comparison is meaningless), per-target
accuracy is computed, and the headline is the unweighted mean across targets.  Intervals
resample COMPOUNDS within each target.
"""
import csv, itertools, json, random, sys
import numpy as np
from scipy.stats import spearmanr

COHORT, SO3LR, MARGIN = sys.argv[1], sys.argv[2], 0.5
so3 = {x['slug']: float(x['binding_kcal_mol']) for x in csv.DictReader(open(SO3LR))
       if x.get('status')=='ok' and x.get('binding_kcal_mol')}
rows=[]
for x in csv.DictReader(open(COHORT)):
    if x.get('status')!='ok': continue
    try:
        if float(x['pb_2_VDWAALS'])>0: continue
        fep=x['published_fepplus_dG_kcal_mol']
        rows.append(dict(slug=x['slug'], t=x['target'], nheavy=float(x['heavy_atoms']),
                         h=-float(x['heavy_atoms']), dg=float(x['experimental_dG_kcal_mol']),
                         pb=float(x['pb_2']), so3=so3.get(x['slug']),
                         fep=float(fep) if fep else None,
                         cens=x['censored'].strip().lower() in ('true','1')))
    except (ValueError,KeyError): continue

KEYS=[('pb','MM/PBSA'),('so3','SO3LR-SF'),('fep','FEP+'),('h','-heavy atoms')]
pw=[r for r in rows if not r['cens'] and r['so3'] is not None]
TGT=sorted({r['t'] for r in pw})
print(f"compounds in pairwise set: {len(pw)}   targets: {len(TGT)}")

def acc(sub,key,band):
    ok=tot=0
    for a,b in itertools.combinations(sub,2):
        if a['t']!=b['t']: continue                       # WITHIN TARGET ONLY
        if abs(a['dg']-b['dg'])<MARGIN or abs(a['nheavy']-b['nheavy'])>band: continue
        if a[key] is None or b[key] is None or a[key]==b[key]: continue
        tot+=1
        if ((-a[key])>(-b[key]))==(a['dg']<b['dg']): ok+=1
    return ok,tot

rng=random.Random(11)
def macro(sub,key,band):
    per={}
    for t in TGT:
        o,n=acc([r for r in sub if r['t']==t],key,band)
        if n>=5: per[t]=o/n
    return (np.mean(list(per.values())) if per else None), per, sum(
        acc([r for r in sub if r['t']==t],key,band)[1] for t in TGT)

def ci(key,band,n=2000):
    m,per,tot=macro(pw,key,band)
    if m is None: return None
    byt={t:[r for r in pw if r['t']==t] for t in TGT}
    bs=[]
    for _ in range(n):
        rs=[]
        for t in TGT:                                     # resample compounds within each target
            s=byt[t]; rs+= [s[rng.randrange(len(s))] for _ in s]
        v,_,_=macro(rs,key,band)
        if v is not None: bs.append(v)
    bs=np.array(bs)
    return m,tot,per,np.percentile(bs,2.5),np.percentile(bs,97.5),float((bs>0.5).mean())

out={}
for band in (0,3):
    print(f"\n=== macro pairwise accuracy, |dHeavy| <= {band}, margin {MARGIN} ===")
    print(f"{'method':<16}{'macro':>8}{'pairs':>7}{'95% CI':>18}{'P(>0.5)':>9}   per-target")
    for k,lab in KEYS:
        r=ci(k,band)
        if r is None: print(f"{lab:<16}{'(too few pairs)':>33}"); continue
        m,tot,per,lo,hi,mass=r
        out[f'band{band}_{k}']=dict(macro=m,pairs=tot,lo=lo,hi=hi,mass=mass,per_target=per)
        pt=" ".join(f"{t.split('_')[0][:4]}={v:.2f}" for t,v in sorted(per.items()))
        print(f"{lab:<16}{m:>8.3f}{tot:>7}{f'[{lo:.3f}, {hi:.3f}]':>18}{mass:>9.3f}   {pt}")
json.dump(out,open('so3lr_macro.json','w'),indent=1,default=float)
print("\nwrote so3lr_macro.json")
