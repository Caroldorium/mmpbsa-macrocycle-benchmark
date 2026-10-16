#!/usr/bin/env python3
"""Why the long-MD stability result was retired: it did not survive longer sampling.

The July 2026 screen ran 3 replicas x 10 ns at 300 K and concluded that three-residue ncAA swaps
stabilise the backbone relative to their canonical parent.  The August 2026 rerun -- 4 replicas
x 250 ns at 350 K, the same systems -- reverses that ranking.  S06, the most stable system in
the short run, becomes the least stable; S01, the canonical parent, becomes the best; and no
three-swap beats its parent any more.

This is the figure behind the decision to stop ranking designs by long explicit-solvent
simulation.  It is kept as history, not as a current result: the source repository marks the
10 ns files "historical context [that] must not replace the current conclusion"
(ncflow-3-swap/results/00_READ_ME_FIRST.md).

Panel B gives the current binding energies for all 13 systems on one consistent run, which is
the number that actually motivates the LptA campaign: thanatin still beats everything designed.
"""
import csv, json, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/Users/daddycl/Desktop/Lab_computation_results.nosync/projects'
LEG = f'{BASE}/untitled folder/ncflow-3-swap/results/10_legacy_10ns_all_systems.tsv'
CUR = f'{BASE}/molecular-dynamics/results/01_all_13_structural_summary.tsv'
BND = f'{BASE}/molecular-dynamics/results/02_all_13_binding_energy_summary.tsv'
OUT = sys.argv[1] if len(sys.argv) > 1 else 'figures/md_legacy_vs_current.png'

leg = {r['figure_id']: r for r in csv.DictReader(open(LEG), delimiter='\t')}
cur = {r['system_id'].split('_')[0]: r for r in csv.DictReader(open(CUR), delimiter='\t')}
shared = [k for k in sorted(leg) if k in cur]
M1 = [k for k in shared if k <= 'S06']          # the M1/M10 family the deck's claim was about

fig, ax = plt.subplots(1, 2, figsize=(13.6, 5.4),
                       gridspec_kw=dict(width_ratios=[1.0, 1.32]))

# --- A: paired slope chart, backbone RMSD ---
for k in M1:
    a = float(leg[k]['median_replica_backbone_rmsd_A'])
    b = float(cur[k]['median_peptide_backbone_rmsd_A'])
    kind = cur[k]['kind']
    hot = k == 'S06'
    col = '#B51F58' if hot else ('#2F6B3C' if 'canonical' in kind.lower() else '#4C78A8')
    ax[0].plot([0, 1], [a, b], '-o', color=col, lw=2.2 if hot else 1.4,
               ms=7 if hot else 5.5, alpha=1.0 if hot else 0.8, zorder=3 if hot else 2)
    ax[0].text(-0.045, a, k, ha='right', va='center', fontsize=9, color=col)
    ax[0].text(1.045, b, k, ha='left', va='center', fontsize=9, color=col,
               fontweight='bold' if hot else 'normal')
ax[0].set_xlim(-0.33, 1.33); ax[0].set_xticks([0, 1])
ax[0].set_xticklabels(['3 × 10 ns, 300 K\n(Jul 2026 -- the deck’s numbers)',
                       '4 × 250 ns, 350 K\n(Aug 2026 -- current)'], fontsize=9.2)
ax[0].set_ylabel('median peptide backbone RMSD (Å)\nlower is more stable')
ax[0].set_title('A   Longer sampling inverts the ranking', fontsize=11)
ax[0].spines[['top', 'right']].set_visible(False)
ax[0].annotate('S06 goes from most stable\nto least stable',
               xy=(1.0, float(cur['S06']['median_peptide_backbone_rmsd_A'])),
               xytext=(0.30, 4.6), fontsize=8.8, color='#B51F58',
               arrowprops=dict(arrowstyle='->', color='#B51F58', lw=1.2))

# --- B: current binding energies, all 13 ---
bnd = [r for r in csv.DictReader(open(BND), delimiter='\t') if 'mature' in r['window']]
bnd.sort(key=lambda r: float(r['mean_across_replica_means_kcal_mol']))
names = [r['system_id'].split('_', 1)[1].replace('_', ' ') for r in bnd]
vals = np.array([float(r['mean_across_replica_means_kcal_mol']) for r in bnd])
sds = np.array([float(r['replica_SD_kcal_mol']) for r in bnd])
cols = ['#B51F58' if 'THANATIN' in r['system_id'] else '#8C8C8C' for r in bnd]
y = np.arange(len(bnd))
ax[1].barh(y, vals, xerr=sds, color=cols, edgecolor='white',
           error_kw=dict(ecolor='#555', lw=1.0, capsize=2.5))
ax[1].axvline(0, color='#333', lw=1.0)
ax[1].set_yticks(y); ax[1].set_yticklabels(names, fontsize=8.4)
ax[1].invert_yaxis()
ax[1].set_xlabel('MM/PBSA, mature 20-250 ns window (kcal/mol) · more negative is tighter')
ax[1].set_title('B   On the current run, thanatin still beats\nevery design we made', fontsize=11)
ax[1].spines[['top', 'right']].set_visible(False)

fig.suptitle('The 10 ns stability result did not survive longer sampling -- which is why '
             'designs are no longer ranked this way.', fontsize=12.5, y=0.985)
fig.text(0.5, 0.012, 'Left: 12 shared systems, M1/M10 family shown · '
         'Right: 13 systems, 4 replicas × 250 ns at 350 K, error bars are replica SD · '
         'thanatin derivative-7, scored separately on the same protocol, reaches −36.04 '
         '± 6.75 kcal/mol',
         ha='center', fontsize=7.4, color='#5E6873')
fig.tight_layout(rect=[0, 0.045, 1, 0.925])
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
for k in M1:
    print(f"  {k}: {float(leg[k]['median_replica_backbone_rmsd_A']):.3f} -> "
          f"{float(cur[k]['median_peptide_backbone_rmsd_A']):.3f} A")
