#!/usr/bin/env python3
"""Why size matching is necessary: affinity and the score both track heavy-atom count."""
import csv, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

G = os.path.expanduser('~/Desktop/Lab_computation_results.nosync/projects/'
                       'transferable-peptide-affinity-ranker-20260824/goals/'
                       'zha25-mmpbsa-replication-20260924/grand_results')
OUT = os.path.expanduser('~/Downloads/lpta_methods_20260927/fig2_size_confound.png')

rows = []
for x in csv.DictReader(open(f'{G}/zha25_cohort_scores.csv')):
    if x.get('status') != 'ok':
        continue
    try:
        if float(x['pb_2_VDWAALS']) > 0 or x['censored'].strip().lower() in ('true', '1'):
            continue
        rows.append((x['target'], float(x['heavy_atoms']),
                     float(x['experimental_dG_kcal_mol']), float(x['pb_2'])))
    except (ValueError, KeyError):
        continue
targets = sorted({r[0] for r in rows})
cols = dict(zip(targets, ['#4C78A8', '#F58518', '#54A24B', '#B51F58', '#9D7BD8']))

fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.1))

# A: affinity genuinely tracks size
for t in targets:
    s = [r for r in rows if r[0] == t]
    ax[0].scatter([r[1] for r in s], [r[2] for r in s], s=16, alpha=.75,
                  color=cols[t], label=f'{t} (n={len(s)})', edgecolor='none')
rho = np.mean([spearmanr([r[1] for r in rows if r[0] == t],
                         [r[2] for r in rows if r[0] == t]).statistic for t in targets])
ax[0].set_xlabel('ligand heavy atoms'); ax[0].set_ylabel('experimental $\\Delta G$ (kcal/mol)')
ax[0].set_title(f'A  Affinity tracks size\nmean within-target $\\rho$ = {rho:+.2f}', fontsize=10.5)
ax[0].legend(frameon=False, fontsize=7.2, loc='lower left')

# B: the score tracks size too
for t in targets:
    s = [r for r in rows if r[0] == t]
    ax[1].scatter([r[1] for r in s], [r[3] for r in s], s=16, alpha=.75,
                  color=cols[t], edgecolor='none')
rho2 = np.mean([spearmanr([r[1] for r in rows if r[0] == t],
                          [r[3] for r in rows if r[0] == t]).statistic for t in targets])
ax[1].set_xlabel('ligand heavy atoms'); ax[1].set_ylabel('MM/PBSA score (kcal/mol)')
ax[1].set_title(f'B  So does the score\nmean within-target $\\rho$ = {rho2:+.2f}', fontsize=10.5)

# C: size distributions overlap enough to match on
for i, t in enumerate(targets):
    h = [r[1] for r in rows if r[0] == t]
    ax[2].scatter(h, np.full(len(h), i) + np.random.default_rng(i).normal(0, .07, len(h)),
                  s=14, alpha=.65, color=cols[t], edgecolor='none')
ax[2].set_yticks(range(len(targets)))
ax[2].set_yticklabels(targets, fontsize=8.5)
ax[2].set_xlabel('ligand heavy atoms')
ax[2].set_title('C  Enough compounds share a size\nto compare them directly', fontsize=10.5)
for a in ax:
    a.spines[['top', 'right']].set_visible(False)
fig.suptitle('Size is a real driver of peptide binding, not a nuisance — which is exactly why '
             'it must be held fixed to test a scoring function', fontsize=11.5, y=1.0)
fig.text(0.5, 0.003, 'Awoonor-Williams et al. cohort, 220 compounds after excluding one '
         'unrelieved clash and all censored labels', ha='center', fontsize=7.5, color='#5E6873')
fig.tight_layout(rect=[0, 0.03, 1, 0.93])
fig.savefig(OUT, dpi=300)
print('wrote', OUT, f'| panel A rho={rho:+.3f}  panel B rho={rho2:+.3f}')
