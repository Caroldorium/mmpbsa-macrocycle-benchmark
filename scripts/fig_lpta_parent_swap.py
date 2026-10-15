#!/usr/bin/env python3
"""Figure 6 -- the size confound on the LptA campaign itself.

Each of the 96 non-canonical swaps is paired with the canonical parent it was built from, so the
comparison is matched by construction: one residue differs.  Panel A asks whether the swap scores
better than its parent, panel B whether MM/PBSA and Rosetta agree on WHICH swap helps, and panel C
whether that agreement survives removing the change in buried interface area.

The partial correlation in panel C is computed on ranks (Spearman partial): invert the 3x3 rank
correlation matrix of dMM/PBSA, dRosetta and d(buried dSASA) and read off the dMM/PBSA-dRosetta
element.  No experimental label exists for any of these designs; nothing here is a skill estimate.
"""
import csv, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC = sys.argv[1] if len(sys.argv) > 1 else 'results/zha25_lpta_parent_swap_delta.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'figures/zha25_lpta_parent_swap.png'

rows = list(csv.DictReader(open(SRC)))
col = lambda k: np.array([float(r[k]) for r in rows])
d03, d14, d19 = col('d_ff03'), col('d_ff14SB'), col('d_ff19SB')
dros, dsasa = col('d_rosetta'), col('d_dSASA')
nparent = len({r['parent'] for r in rows})

def rank(a):
    order = np.argsort(a, kind='mergesort')
    r = np.empty(len(a), float); r[order] = np.arange(len(a), dtype=float)
    for v in np.unique(a):                      # average ties
        m = a == v
        if m.sum() > 1: r[m] = r[m].mean()
    return r

def rho(a, b):
    return float(np.corrcoef(rank(a), rank(b))[0, 1])

def partial(a, b, c):
    """Spearman partial correlation of a and b given c."""
    P = np.linalg.inv(np.corrcoef(np.vstack([rank(a), rank(b), rank(c)])))
    return float(-P[0, 1] / np.sqrt(P[0, 0] * P[1, 1]))

fig, ax = plt.subplots(1, 3, figsize=(19.0, 5.8))

better, med = int((d03 < 0).sum()), float(np.median(d03))
ax[0].hist(d03, bins=24, color='#4C78A8', edgecolor='white', linewidth=0.5)
ax[0].axvline(med, color='#C0392B', ls='--', lw=2, label=f'median {med:.2f} kcal/mol')
ax[0].axvline(0, color='#333', lw=1.2)
ax[0].set_xlabel('Δ MM/PBSA, swap minus its canonical parent (kcal/mol)')
ax[0].set_ylabel('ncAA swaps')
ax[0].set_title(f'A   {better} of {len(rows)} swaps score better than their parent\n'
                f'MM/PBSA ff03, {len(rows)} matched pairs across {nparent} parents', fontsize=11)
ax[0].legend(frameon=False, fontsize=9, loc='upper left')

sc = ax[1].scatter(dros, d03, c=dsasa, cmap='viridis', s=34, edgecolor='none')
ax[1].axhline(0, color='#888', lw=0.9); ax[1].axvline(0, color='#888', lw=0.9)
ax[1].set_xlabel('Δ Rosetta dG_separated (REU)')
ax[1].set_ylabel('Δ MM/PBSA (kcal/mol)')
ax[1].set_title(f'B   but they barely agree on WHICH swap helps\n'
                f'$\\rho$ = {rho(d03, dros):+.3f}; colour shows both track buried area', fontsize=11)
fig.colorbar(sc, ax=ax[1], label='Δ buried dSASA ($\\AA^2$)')

labels = ['ff03', 'ff14SB', 'ff19SB']
raw = [rho(d, dros) for d in (d03, d14, d19)]
ctl = [partial(d, dros, dsasa) for d in (d03, d14, d19)]
x = np.arange(3)
ax[2].bar(x - 0.19, raw, 0.36, color='#9EC4E3', label='raw', edgecolor='white')
ax[2].bar(x + 0.19, ctl, 0.36, color='#2B6CB0', label='buried area controlled', edgecolor='white')
ax[2].axhline(0, color='#333', lw=1.0)
ax[2].set_xticks(x); ax[2].set_xticklabels(labels)
ax[2].set_ylim(0, 0.42)
ax[2].set_ylabel('$\\rho$ between ΔMM/PBSA and ΔRosetta')
ax[2].set_title('C   remove the size change and the agreement\nis gone in every force field', fontsize=11)
ax[2].legend(frameon=False, fontsize=9)

fig.suptitle('Zha25 MM/PBSA on the 20 canonical LptA binders and their 96 ncAA swaps '
             '— matched parent-to-swap pairs.  Swaps bury '
             f'{np.mean(dsasa):+.0f} $\\AA^2$ more on average,\n'
             'and both scoring functions reward that; neither has an experimental label to '
             'check against.', fontsize=13, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.86])
fig.savefig(OUT, dpi=200)
print('wrote', OUT)
for lab, r, c in zip(labels, raw, ctl):
    print(f'  {lab}: raw rho {r:+.3f}  controlled {c:+.3f}')
print(f'  {better}/{len(rows)} better than parent, median {med:.2f}, mean dSASA {np.mean(dsasa):+.1f}')
