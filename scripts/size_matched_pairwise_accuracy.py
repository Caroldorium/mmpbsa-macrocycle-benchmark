#!/usr/bin/env python3
"""Pairwise accuracy among size-matched peptides.

The comparison that matters for lead optimisation: given two analogues of the same size,
how often does a method pick the tighter binder?  A size descriptor cannot answer at all
when the sizes are equal, which is the point -- it is a proxy for buried area, not a rival
scoring function.
"""
import csv, itertools, os, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = os.path.expanduser('~/Desktop/Lab_computation_results.nosync/projects/'
                       'transferable-peptide-affinity-ranker-20260824/goals/'
                       'zha25-mmpbsa-replication-20260924/grand_results')
OUT = os.path.expanduser('~/Downloads/lpta_methods_20260927/size_matched_pairwise_accuracy.png')

rows = []
for x in csv.DictReader(open(f'{G}/zha25_cohort_scores.csv')):
    if x.get('status') != 'ok':
        continue
    try:
        if float(x['pb_2_VDWAALS']) > 0:          # unrelieved clash, pre-registered screen
            continue
        if x['censored'].strip().lower() in ('true', '1'):
            continue
        rows.append(dict(t=x['target'], h=float(x['heavy_atoms']),
                         dg=float(x['experimental_dG_kcal_mol']), pb=float(x['pb_2']),
                         fep=float(x['published_fepplus_dG_kcal_mol'])
                         if x['published_fepplus_dG_kcal_mol'] else None))
    except (ValueError, KeyError):
        continue

def flags(sub, band, margin=0.5):
    out = []
    for a, b in itertools.combinations(sub, 2):
        if abs(a['dg'] - b['dg']) < margin:
            continue
        if band is not None and abs(a['h'] - b['h']) > band:
            continue
        rec = {}
        for key in ('pb', 'fep'):
            if a[key] is None or b[key] is None:
                rec[key] = None; continue
            pa, pv = -a[key], -b[key]          # more negative = tighter
            rec[key] = None if pa == pv else int((pa > pv) == (a['dg'] < b['dg']))
        out.append(rec)
    return out

rng = random.Random(7)
def ci(v, n=2000):
    bs = [np.mean([v[rng.randrange(len(v))] for _ in v]) for _ in range(n)]
    return np.mean(v), np.percentile(bs, 2.5), np.percentile(bs, 97.5)

targets = sorted(set(r['t'] for r in rows))
fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5), sharey=True)
for ax, band, title in zip(axes, (0, 3),
                           ("Identical heavy-atom count", "Within 3 heavy atoms")):
    xs = np.arange(len(targets) + 1)
    for off, key, col, lab in ((-0.19, 'pb', '#B51F58', 'MM/PBSA'),
                               (0.19, 'fep', '#4C78A8', 'FEP+ (published)')):
        m, lo, hi, per = [], [], [], {}
        for t in targets:
            v = [x[key] for x in flags([r for r in rows if r['t'] == t], band) if x[key] is not None]
            a, l, h = ci(v) if len(v) >= 15 else (np.nan, np.nan, np.nan)
            per[t] = a; m.append(a); lo.append(a - l); hi.append(h - a)
        mv = [v for v in per.values() if not np.isnan(v)]
        bs = [np.mean([mv[rng.randrange(len(mv))] for _ in mv]) for _ in range(5000)]
        mac = np.mean(mv)
        m.append(mac); lo.append(mac - np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5) - mac)
        ax.bar(xs + off, m, 0.36, yerr=[lo, hi], capsize=3, color=col, label=lab,
               edgecolor='white', linewidth=0.6)
    ax.axhline(0.5, color='#444', lw=1.2, ls='--')
    ax.text(len(targets) + 0.52, 0.505, 'chance', fontsize=8, color='#444', ha='right')
    ax.set_xticks(xs); ax.set_xticklabels(targets + ['MACRO'], rotation=20, ha='right', fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.set_ylim(0.25, 1.0); ax.spines[['top', 'right']].set_visible(False)
axes[0].set_ylabel('pairwise accuracy\n(picks the tighter binder)', fontsize=10)
axes[0].legend(frameon=False, fontsize=9, loc='upper left')
fig.suptitle('Among peptides of the same size, a heavy-atom count has no information — '
             'MM/PBSA still ranks above chance', fontsize=11.5, y=0.99)
fig.text(0.5, 0.005, 'Awoonor-Williams et al. cohort, 220 compounds after excluding one '
         'unrelieved clash and all censored labels · decisive pairs only (|ΔG difference| '
         '≥ 0.5 kcal/mol) · bars are 95% bootstrap intervals',
         ha='center', fontsize=7.5, color='#5E6873')
fig.tight_layout(rect=[0, 0.035, 1, 0.94])
fig.savefig(OUT, dpi=300)
print("wrote", OUT)
