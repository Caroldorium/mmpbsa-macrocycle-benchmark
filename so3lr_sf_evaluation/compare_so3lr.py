#!/usr/bin/env python3
"""SO3LR-SF vs MM/PBSA vs FEP+ on the macrocycle cohort, under the manuscript's own estimators.

Two estimators, both from the paper:
  1. within-receptor (macro) Spearman -- pooling across targets is confounded by size
  2. size-matched pairwise accuracy at |dHeavy| = 0 and <= 3, bootstrap over COMPOUNDS
"""
import csv, itertools, json, os, random, sys
import numpy as np
from scipy.stats import spearmanr

COHORT = sys.argv[1]
SO3LR  = sys.argv[2]
MARGIN = 0.5

so3 = {}
for x in csv.DictReader(open(SO3LR)):
    if x.get('status') == 'ok' and x.get('binding_kcal_mol'):
        try: so3[x['slug']] = float(x['binding_kcal_mol'])
        except ValueError: pass

rows = []
for x in csv.DictReader(open(COHORT)):
    if x.get('status') != 'ok': continue
    try:
        pb   = float(x['pb_2'])
        vdw  = float(x['pb_2_VDWAALS'])
        dg   = float(x['experimental_dG_kcal_mol'])
        hv   = float(x['heavy_atoms'])
    except (ValueError, KeyError): continue
    if vdw > 0: continue                       # MCL_1__cpd_55, unrelieved clash
    fep = x['published_fepplus_dG_kcal_mol']
    # every predictor is expressed dG-like: MORE NEGATIVE = PREDICTED TIGHTER.
    # heavy-atom count is therefore negated -- more atoms predicts tighter binding.
    rows.append(dict(slug=x['slug'], t=x['target'], h=-hv, nheavy=hv, dg=dg, pb=pb,
                     so3=so3.get(x['slug']),
                     fep=float(fep) if fep else None,
                     cens=x['censored'].strip().lower() in ('true','1')))

have = [r for r in rows if r['so3'] is not None]
print(f"cohort rows (status ok, VDW<=0): {len(rows)}")
print(f"  with an SO3LR-SF score:        {len(have)}")
print(f"  missing:                       {len(rows)-len(have)}")

# ---------- 1. within-receptor Spearman ----------
KEYS = [('pb','MM/PBSA (pb_2)'), ('so3','SO3LR-SF'), ('fep','FEP+ (published)'), ('h','-heavy atoms')]
print("\n=== within-receptor Spearman (rho vs experimental dG; POSITIVE = correct ranking) ===")
print(f"{'target':<12}{'n':>5}" + "".join(f"{lab:>20}" for _,lab in KEYS))
macro = {k: [] for k,_ in KEYS}
for t in sorted({r['t'] for r in have}):
    sub = [r for r in have if r['t']==t]
    line = f"{t:<12}{len(sub):>5}"
    for k,_ in KEYS:
        v = [(r[k], r['dg']) for r in sub if r[k] is not None]
        if len(v) < 5: line += f"{'--':>20}"; continue
        rho = spearmanr([a for a,_ in v], [b for _,b in v]).correlation
        macro[k].append(rho); line += f"{rho:>20.3f}"
    print(line)
print(f"{'MACRO MEAN':<12}{'':>5}" + "".join(
    f"{np.mean(macro[k]):>20.3f}" if macro[k] else f"{'--':>20}" for k,_ in KEYS))

# ---------- 2. size-matched pairwise accuracy ----------
pw = [r for r in have if not r['cens']]
print(f"\n=== size-matched pairwise accuracy (margin {MARGIN} kcal/mol, n={len(pw)}) ===")

def acc(sub, key, band):
    ok = tot = 0
    for a,b in itertools.combinations(sub, 2):
        if abs(a['dg']-b['dg']) < MARGIN or abs(a['nheavy']-b['nheavy']) > band: continue
        if a[key] is None or b[key] is None or a[key]==b[key]: continue
        tot += 1
        if ((-a[key]) > (-b[key])) == (a['dg'] < b['dg']): ok += 1
    return ok, tot

rng = random.Random(11)
def ci(sub, key, band, n=2000):
    o,t = acc(sub, key, band)
    if t < 15: return None
    bs = []
    for _ in range(n):
        rs = [sub[rng.randrange(len(sub))] for _ in sub]
        o2,t2 = acc(rs, key, band)
        if t2 >= 5: bs.append(o2/t2)
    bs = np.array(bs)
    return o/t, t, np.percentile(bs,2.5), np.percentile(bs,97.5), float((bs>0.5).mean())

out = {}
for band in (0,3):
    print(f"\n-- band |dHeavy| <= {band} --")
    print(f"{'method':<22}{'acc':>8}{'pairs':>8}{'95% CI':>20}{'P(>chance)':>13}")
    for k,lab in KEYS:
        r = ci(pw, k, band)
        if r is None: print(f"{lab:<22}{'(too few pairs)':>49}"); continue
        a,t,lo,hi,m = r
        out[f'band{band}_{k}'] = dict(acc=a, pairs=t, lo=lo, hi=hi, mass=m)
        print(f"{lab:<22}{a:>8.3f}{t:>8}{f'[{lo:.3f}, {hi:.3f}]':>20}{m:>13.3f}")

print("\n-- per target, band 3 --")
print(f"{'target':<12}{'n':>5}" + "".join(f"{lab:>22}" for _,lab in KEYS))
for t in sorted({r['t'] for r in pw}):
    sub = [r for r in pw if r['t']==t]
    line = f"{t:<12}{len(sub):>5}"
    for k,_ in KEYS:
        r = ci(sub, k, 3)
        line += f"{'--':>22}" if r is None else f"{f'{r[0]:.3f} [{r[2]:.2f},{r[3]:.2f}]':>22}"
    print(line)

json.dump(out, open('so3lr_comparison.json','w'), indent=1)
print("\nwrote so3lr_comparison.json")
