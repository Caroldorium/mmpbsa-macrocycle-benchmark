#!/usr/bin/env python3
"""Longer sampling does not help: replica noise exceeds the between-compound signal.

40 KRAS compounds, four independent 10 ns explicit-solvent replicas each, scored over the
inclusive 2-10 ns window at 81 frames per replica -- 324 frames per compound.  The question is
not whether the ensemble score is noisy in absolute terms but whether it is noisier than the
spread it has to resolve.  It is: the median replicate-to-replicate standard deviation is larger
than the standard deviation BETWEEN compounds, so which replicas you happened to run moves a
compound further than real chemistry does.

Panel A shows that directly -- each compound's four replica values against its own median, with
the between-compound spread drawn as the band they would have to beat.  Panel B gives the two
standard deviations side by side for the whole panel and for each half of it, so the conclusion
can be seen not to depend on the subset.
"""
import csv, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC = sys.argv[1] if len(sys.argv) > 1 else (
    '/Users/daddycl/Desktop/Lab_computation_results.nosync/projects/'
    'transferable-peptide-affinity-ranker-20260824/goals/kras-md-metric-correlation-20260917/'
    'grand_results/kras_all40_md_score_and_stability_metrics.csv')
OUT = sys.argv[2] if len(sys.argv) > 2 else 'figures/sampling_noise.png'

rows = [r for r in csv.DictReader(open(SRC)) if r['status'] == 'PASS']
med = np.array([float(r['median_score_kcal_mol']) for r in rows])
sd = np.array([float(r['replica_score_stdev_kcal_mol']) for r in rows])
reps = [[float(v) for v in r['replica_values_kcal_mol'].split(';')] for r in rows]
role = np.array([r['panel_role'] for r in rows])

def stats(mask):
    return float(np.median(sd[mask])), float(np.std(med[mask]))

groups = [('all 40', np.ones(len(rows), bool)),
          ('confirmatory\n(n = 20)', role == 'CONFIRMATORY'),
          ('primary\n(n = 20)', role == 'PRIMARY')]

fig, ax = plt.subplots(1, 2, figsize=(13.4, 5.2),
                       gridspec_kw=dict(width_ratios=[1.75, 1.0]))

# --- A: every replica value, compounds ordered by their own median ---
order = np.argsort(med)
x = np.arange(len(rows))
between = np.std(med)
ax[0].axhspan(-between, between, color='#4C78A8', alpha=0.16, zorder=0,
              label=f'between-compound SD  ±{between:.2f}')
for xi, i in enumerate(order):
    d = np.array(reps[i]) - med[i]
    ax[0].plot([xi] * len(d), d, 'o', ms=4.2, color='#B51F58', alpha=0.72,
               zorder=3, mew=0)
ax[0].axhline(0, color='#333', lw=1.0, zorder=2)
ax[0].set_xlabel('KRAS compounds, ordered by their own median score')
ax[0].set_ylabel('replica score minus that compound’s median\n(kcal/mol)')
ax[0].set_title('A   Each compound’s four replicas scatter further\n'
                'than the compounds scatter from each other', fontsize=11)
ax[0].legend(frameon=False, fontsize=9, loc='upper left')
ax[0].spines[['top', 'right']].set_visible(False)

# --- B: the two standard deviations, per subset ---
xs = np.arange(len(groups))
noise = [stats(m)[0] for _, m in groups]
signal = [stats(m)[1] for _, m in groups]
ax[1].bar(xs - 0.19, noise, 0.36, color='#B51F58', edgecolor='white',
          label='replica noise (median SD)')
ax[1].bar(xs + 0.19, signal, 0.36, color='#4C78A8', edgecolor='white',
          label='between-compound signal (SD)')
for xi, (n, s) in enumerate(zip(noise, signal)):
    ax[1].text(xi - 0.19, n + 0.12, f'{n:.2f}', ha='center', fontsize=9, color='#B51F58')
    ax[1].text(xi + 0.19, s + 0.12, f'{s:.2f}', ha='center', fontsize=9, color='#4C78A8')
    ax[1].text(xi, max(n, s) + 0.75, f'ratio {n / s:.2f}', ha='center', fontsize=9.5,
               color='#1F2933', fontweight='bold')
ax[1].set_xticks(xs); ax[1].set_xticklabels([g for g, _ in groups], fontsize=9.5)
ax[1].set_ylabel('kcal/mol')
ax[1].set_ylim(0, max(noise) * 1.52)
ax[1].set_title('B   Noise exceeds signal in the whole panel\nand in each half of it', fontsize=11)
ax[1].legend(frameon=False, fontsize=9, loc='upper right')
ax[1].spines[['top', 'right']].set_visible(False)

fig.suptitle('Four 10 ns replicas per compound -- 324 frames -- do not make the ensemble '
             'score reproducible at the resolution a ranking needs.',
             fontsize=12.5, y=0.985)
fig.text(0.5, 0.012, f'{len(rows)} KRAS compounds, 4 independent 10 ns replicas, inclusive '
         '2-10 ns window, 81 frames per replica · '
         'kras_all40_md_score_and_stability_metrics.csv',
         ha='center', fontsize=7.4, color='#5E6873')
fig.tight_layout(rect=[0, 0.04, 1, 0.925])
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
for g, m in groups:
    n, s = stats(m)
    print(f'  {g.splitlines()[0]:14s} noise {n:.2f}  signal {s:.2f}  ratio {n/s:.2f}')
