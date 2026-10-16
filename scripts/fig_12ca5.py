#!/usr/bin/env python3
"""The sealed external test: on a sixth target, the score read pose quality instead of chemistry.

17 peptidomimetics against the anti-haemagglutinin antibody 12CA5, affinities sealed before any
prediction.  These compounds were threaded onto a template they do not match, and it shows: the
pose builder left 0-44 clashes per complex and backbone RMSD up to 4.72 A.

Panel A is the diagnostic.  Three of the built poses carry positive van der Waals terms and are
removed by the screen of Section 2.2 (drawn as crosses, not used in the correlation).  Across the
13 that survive, the MM/PBSA score tracks the builder's CLASH COUNT at rho = -0.891 -- the score
is reading how badly the pose was built, not how well the compound binds.

Panel B is Table 6: every method, including the size baselines, lands below chance on this
target, and the size-affinity relationship reverses sign relative to the benchmark.  Those values
come from the full 17-compound panel; the deposited score file holds the subset that could be
scored, so they are quoted from the manuscript rather than recomputed here.
"""
import csv, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

SRC = sys.argv[1] if len(sys.argv) > 1 else 'results/zha25_12ca5_built_pose_scores.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'figures/external_12ca5.png'

rows = list(csv.DictReader(open(SRC)))
keep = [r for r in rows if float(r['vdw']) <= 0]
drop = [r for r in rows if float(r['vdw']) > 0]
pb = np.array([float(r['pb_2']) for r in keep])
cl = np.array([float(r['clash']) for r in keep])
rho = spearmanr(pb, cl).correlation

fig, ax = plt.subplots(1, 2, figsize=(13.2, 5.2),
                       gridspec_kw=dict(width_ratios=[1.12, 1.0]))

col = ['#B51F58' if r['cls'] == 'binder' else '#4C78A8' for r in keep]
ax[0].scatter(cl, pb, c=col, s=62, zorder=3, edgecolor='white', linewidth=0.8)
z = np.polyfit(cl, pb, 1)
xs = np.linspace(cl.min(), cl.max(), 50)
ax[0].plot(xs, np.polyval(z, xs), '--', color='#555', lw=1.3, zorder=2)
tr = ax[0].get_xaxis_transform()
for r in drop:
    ax[0].scatter([float(r['clash'])], [0.965], marker='x', s=70, c='#999',
                  transform=tr, clip_on=False, zorder=4)
if drop:
    ax[0].text(0.02, 0.03, f'✕  {len(drop)} poses carried positive VDW and were excluded by the\n'
               f'    Section 2.2 screen — marked along the top, no usable score',
               transform=ax[0].transAxes, ha='left', va='bottom', fontsize=8.2, color='#777')
ax[0].set_xlabel('clashes below 1.8 Å left by the pose builder')
ax[0].set_ylabel('MM/PBSA pb_2 (kcal/mol)')
ax[0].set_title(f'A   The score tracks pose quality, not binding\n'
                f'$\\rho$ = {rho:+.3f} against clash count, n = {len(keep)}', fontsize=11)
ax[0].scatter([], [], c='#B51F58', s=62, label='binder')
ax[0].scatter([], [], c='#4C78A8', s=62, label='non-binder')
ax[0].legend(frameon=False, fontsize=9, loc='upper right')
ax[0].spines[['top', 'right']].set_visible(False)

# Table 6, full 17-compound panel
T6 = [('heavy-atom\ncount', -0.750, 0.257), ('molecular\nweight', -0.750, 0.271),
      ('chemical model\n(pAffinity)', -0.107, 0.214), ('chemical model\n(rank score)', -0.306, 0.200)]
x = np.arange(len(T6))
ax[1].bar(x - 0.19, [t[1] for t in T6], 0.36, color='#B51F58', edgecolor='white',
          label='$\\rho$ vs p$K_D$ (7 binders)')
ax[1].bar(x + 0.19, [t[2] for t in T6], 0.36, color='#8C8C8C', edgecolor='white',
          label='binder / non-binder AUC')
ax[1].axhline(0.5, color='#2F6B3C', ls='--', lw=1.4)
ax[1].text(len(T6) - 0.45, 0.53, 'chance AUC', ha='right', fontsize=8.6, color='#2F6B3C')
ax[1].axhline(0, color='#333', lw=1.0)
for xi, t in enumerate(T6):
    ax[1].text(xi - 0.19, t[1] - 0.075, f'{t[1]:+.3f}', ha='center', fontsize=8.8, color='#B51F58')
    ax[1].text(xi + 0.19, t[2] + 0.03, f'{t[2]:.3f}', ha='center', fontsize=8.8, color='#555')
ax[1].set_xticks(x); ax[1].set_xticklabels([t[0] for t in T6], fontsize=8.8)
ax[1].set_ylim(-1.0, 0.78)
ax[1].set_ylabel('correlation  /  AUC')
ax[1].set_title('B   Every method lands below chance,\nand size flips sign', fontsize=11)
ax[1].legend(frameon=False, fontsize=8.6, loc='upper left')
ax[1].spines[['top', 'right']].set_visible(False)

fig.suptitle('A sixth target, held entirely outside development, with affinities sealed before '
             'prediction — everything failed.', fontsize=12.5, y=0.985)
fig.text(0.5, 0.012, '17 peptidomimetics vs anti-haemagglutinin antibody 12CA5, biolayer '
         'interferometry, 7 binders (0.21–77 nM) and 10 non-binders · '
         'heavy-atom count correlates −0.750 here against +0.49 to +0.61 on the benchmark',
         ha='center', fontsize=7.4, color='#5E6873')
fig.tight_layout(rect=[0, 0.04, 1, 0.925])
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
print(f'  kept {len(keep)} of {len(rows)} after the positive-VDW screen; rho vs clash {rho:+.3f}')
