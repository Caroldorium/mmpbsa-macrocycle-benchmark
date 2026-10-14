#!/usr/bin/env python3
"""Pairwise accuracy among size-matched peptides.

The comparison that matters for lead optimisation: given two analogues of the same size, how
often does a method pick the tighter binder?  A size descriptor cannot answer at all when the
sizes are equal, which is the point -- it is a proxy for buried interface area, not a rival
scoring function.

Intervals resample COMPOUNDS, not pairs.  Pairs drawn from n compounds share compounds and are
not independent, so a pair-level bootstrap badly understates the uncertainty: on the exact-match
band it gives [0.530, 0.693] where the compound-level interval is [0.472, 0.720].
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
MARGIN = 0.5

rows = []
for x in csv.DictReader(open(f'{G}/zha25_cohort_scores.csv')):
    if x.get('status') != 'ok':
        continue
    try:
        if float(x['pb_2_VDWAALS']) > 0 or x['censored'].strip().lower() in ('true', '1'):
            continue
        rows.append(dict(t=x['target'], h=float(x['heavy_atoms']),
                         dg=float(x['experimental_dG_kcal_mol']), pb=float(x['pb_2']),
                         fep=float(x['published_fepplus_dG_kcal_mol'])
                         if x['published_fepplus_dG_kcal_mol'] else None))
    except (ValueError, KeyError):
        continue

def acc(sub, key, band):
    ok = tot = 0
    for a, b in itertools.combinations(sub, 2):
        if abs(a['dg'] - b['dg']) < MARGIN or abs(a['h'] - b['h']) > band:
            continue
        if a[key] is None or b[key] is None or a[key] == b[key]:
            continue
        tot += 1
        if ((-a[key]) > (-b[key])) == (a['dg'] < b['dg']):
            ok += 1
    return ok, tot

rng = random.Random(11)

def ci(sub, key, band, n=2000):
    o, t = acc(sub, key, band)
    if t < 15:
        return np.nan, np.nan, np.nan, []
    bs = []
    for _ in range(n):
        rs = [sub[rng.randrange(len(sub))] for _ in sub]
        o2, t2 = acc(rs, key, band)
        if t2 >= 5:
            bs.append(o2 / t2)
    return o / t, np.percentile(bs, 2.5), np.percentile(bs, 97.5), bs

targets = sorted({r['t'] for r in rows})
fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0), sharey=True)
for ax, band, title in zip(axes, (0, 3),
                           ("Identical heavy-atom count\n(strictest; fewest pairs)",
                            "Within 3 heavy atoms\n(better powered)")):
    xs = np.arange(len(targets) + 1)
    for off, key, col, lab in ((-0.19, 'pb', '#B51F58', 'MM/PBSA'),
                               (0.19, 'fep', '#4C78A8', 'FEP+ (published)')):
        m, lo, hi, per, draws = [], [], [], {}, {}
        for t in targets:
            sub = [r for r in rows if r['t'] == t]
            a, l, h, bs = ci(sub, key, band)
            if np.isnan(a):
                m.append(np.nan); lo.append(0); hi.append(0); continue
            per[t] = a; draws[t] = bs
            m.append(a); lo.append(a - l); hi.append(h - a)
        ts = list(draws)
        mac = float(np.mean([per[t] for t in ts]))
        mbs = [np.mean([draws[t][rng.randrange(len(draws[t]))] for t in ts]) for _ in range(5000)]
        m.append(mac); lo.append(mac - np.percentile(mbs, 2.5)); hi.append(np.percentile(mbs, 97.5) - mac)
        ax.bar(xs + off, m, 0.36, yerr=[lo, hi], capsize=3, color=col, label=lab,
               edgecolor='white', linewidth=0.6)
        if key == 'pb':
            ax.text(len(targets) + off, min(mac + hi[-1] + 0.03, 0.97),
                    f"P>chance\n{np.mean(np.array(mbs) > 0.5):.2f}",
                    ha='center', va='bottom', fontsize=7.8, color='#B51F58', linespacing=0.95)
    ax.axhline(0.5, color='#444', lw=1.2, ls='--')
    ax.text(len(targets) + 0.52, 0.513, 'chance', fontsize=8, color='#444', ha='right')
    ax.set_xticks(xs); ax.set_xticklabels(targets + ['MACRO'], rotation=20, ha='right', fontsize=9)
    ax.set_title(title, fontsize=10.5)
    ax.set_ylim(0.20, 1.06); ax.spines[['top', 'right']].set_visible(False)
axes[0].set_ylabel('pairwise accuracy\n(picks the tighter binder)', fontsize=10)
axes[0].legend(frameon=False, fontsize=9, loc='upper left')
fig.suptitle('A heavy-atom count has no information among equally sized peptides.\n'
             'MM/PBSA ranks above chance — decisively once the window admits enough pairs.',
             fontsize=11.5, y=0.995, linespacing=1.35)
fig.text(0.5, 0.005, 'Awoonor-Williams et al. cohort, 220 compounds after excluding one unrelieved '
         'clash and all censored labels · decisive pairs only (|ΔG difference| ≥ 0.5 '
         'kcal/mol) · bars are 95% bootstrap intervals over COMPOUNDS, not pairs',
         ha='center', fontsize=7.3, color='#5E6873')
fig.tight_layout(rect=[0, 0.045, 1, 0.885])
fig.savefig(OUT, dpi=300)
print("wrote", OUT)
