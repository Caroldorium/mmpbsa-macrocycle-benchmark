#!/usr/bin/env python3
"""Consensus scoring was tested and does not work: the errors are not orthogonal.

The argument for combining two scoring functions is that their mistakes are independent, so
averaging cancels them.  That is testable, and here it fails.  Within each receptor, rank the
compounds by each method and by experiment; the two methods' rank ERRORS correlate at +0.714.
They are wrong about the same compounds, in the same direction, so a rank-z consensus is worse
than the better component alone -- 0.445 against MM/PBSA's 0.545.

The MM/PBSA + MM/GBSA pair is the positive control for the measurement: two variants of the same
energy function should be maximally redundant, and they are (+0.920).  SO3LR-SF sits closer to
that control than to independence.

Numbers come from so3lr_sf_evaluation/consensus.py, which prints them all.
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = 'figures/consensus_fails.png'

PAIRS = [('MM/PBSA\nvs SO3LR-SF', 0.714, '#B51F58'),
         ('MM/PBSA\nvs MM/GBSA', 0.920, '#8C8C8C'),
         ('SO3LR-SF\nvs MM/GBSA', 0.780, '#8C8C8C')]
METHODS = [('MM/PBSA', 0.545, '#B51F58'),
           ('consensus\nPBSA + SO3LR', 0.445, '#D98CA8'),
           ('control\nPBSA + GBSA', 0.513, '#8C8C8C'),
           ('MM/GBSA', 0.478, '#8C8C8C'),
           ('SO3LR-SF', 0.292, '#8C8C8C')]
NULL = 0.487

fig, ax = plt.subplots(1, 2, figsize=(13.0, 5.2),
                       gridspec_kw=dict(width_ratios=[1.0, 1.35]))

x = np.arange(len(PAIRS))
ax[0].bar(x, [v for _, v, _ in PAIRS], 0.54,
          color=[c for _, _, c in PAIRS], edgecolor='white')
for xi, (_, v, _) in enumerate(PAIRS):
    ax[0].text(xi, v + 0.022, f'{v:+.3f}', ha='center', fontsize=10.5, fontweight='bold')
ax[0].axhline(0, color='#333', lw=1.0)
ax[0].text(1.0, -0.075, 'errors that cancelled under averaging would sit near zero',
           ha='center', va='center', fontsize=8.8, color='#2F6B3C', style='italic')
ax[0].set_xticks(x); ax[0].set_xticklabels([n for n, _, _ in PAIRS], fontsize=9.5)
ax[0].set_ylabel('within-receptor rank-error correlation')
ax[0].set_ylim(-0.14, 1.04)
ax[0].spines['bottom'].set_position(('data', 0.0))
ax[0].tick_params(axis='x', length=0, pad=34)
ax[0].set_title('A   The two methods are wrong about\nthe same compounds', fontsize=11)
ax[0].spines[['top', 'right']].set_visible(False)

x2 = np.arange(len(METHODS))
ax[1].bar(x2, [v for _, v, _ in METHODS], 0.56,
          color=[c for _, _, c in METHODS], edgecolor='white')
for xi, (_, v, _) in enumerate(METHODS):
    ax[1].text(xi, v + 0.014, f'{v:.3f}', ha='center', fontsize=10, zorder=5,
               bbox=dict(fc='white', ec='none', pad=1.2))
ax[1].axhline(NULL, color='#2F6B3C', ls='--', lw=1.4, zorder=1)
ax[1].text(len(METHODS) - 0.42, NULL + 0.008, f'heavy-atom null  {NULL:.3f}',
           ha='right', fontsize=8.6, color='#2F6B3C')
ax[1].set_xticks(x2); ax[1].set_xticklabels([n for n, _, _ in METHODS], fontsize=9.3)
ax[1].set_ylabel('macro Spearman ρ against experiment')
ax[1].set_ylim(0, 0.63)
ax[1].set_title('B   So the consensus loses to MM/PBSA alone\n'
                '-- and so does the two-MM control', fontsize=11)
ax[1].spines[['top', 'right']].set_visible(False)

fig.suptitle('Combining the two scoring functions was tested, not assumed -- and it does '
             'not help.', fontsize=12.5, y=0.985)
fig.text(0.5, 0.012, '234 compounds, 5 targets · ranks z-scored within receptor before '
         'combining · same estimator and bootstrap as Table 5 · '
         'so3lr_sf_evaluation/consensus.py', ha='center', fontsize=7.4, color='#5E6873')
fig.tight_layout(rect=[0, 0.04, 1, 0.925])
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
