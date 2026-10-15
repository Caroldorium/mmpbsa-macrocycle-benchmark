#!/usr/bin/env python3
"""Figure 3 -- per-target parity: raw MM/PBSA, per-family calibrated MM/PBSA, published FEP+.

Three rows over the five FEP+ families.  Row 1 plots the raw score, whose scale is ~7x the
experimental range, so y=x is meaningless there.  Row 2 applies a per-family linear calibration
(dG ~ a + b*score) and plots against y=x, which is what makes the RMSE interpretable -- and the
fitted slopes near 0.1-0.2 are the point: the raw score moves 5-10 kcal/mol for every 1 kcal/mol
of real affinity.  Row 3 is published FEP+, which needs no fitting.

Censored measurements (assay ceiling) are drawn as open circles and kept in the correlations,
matching the cohort definition used throughout.
"""
import csv, os, sys
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

SRC = sys.argv[1] if len(sys.argv) > 1 else 'results/zha25_cohort_scores.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'figures/zha25_fepplus_parity_5families.png'
SCORE, VDW = 'pb_2', 'pb_2_VDWAALS'

rows, dropped = [], []
for x in csv.DictReader(open(SRC)):
    if x.get('status') != 'ok': continue
    try:
        v = float(x[VDW]); s = float(x[SCORE]); dg = float(x['experimental_dG_kcal_mol'])
    except (ValueError, KeyError): continue
    if v > 0:                                   # unrelieved clash: excluded, reported in title
        dropped.append((x['slug'], v)); continue
    fep = x['published_fepplus_dG_kcal_mol']
    rows.append(dict(t=x['target'], dg=dg, s=s, fep=float(fep) if fep else None,
                     cens=x['censored'].strip().lower() in ('true', '1')))

T = sorted({r['t'] for r in rows})
fig, axes = plt.subplots(3, len(T), figsize=(4.0*len(T), 8.6))
macro = {'s': [], 'fep': []}; rmse = {'cal': [], 'fep': []}

for j, t in enumerate(T):
    sub = [r for r in rows if r['t'] == t]
    dg = np.array([r['dg'] for r in sub]); sc = np.array([r['s'] for r in sub])
    cen = np.array([r['cens'] for r in sub])
    fo = [r for r in sub if r['fep'] is not None]
    fdg = np.array([r['dg'] for r in fo]); fv = np.array([r['fep'] for r in fo])
    fcen = np.array([r['cens'] for r in fo])

    rho = spearmanr(sc, dg).correlation; R = pearsonr(sc, dg)[0]
    macro['s'].append(rho)
    b, a = np.polyfit(sc, dg, 1)                        # dG ~ a + b*score
    cal = a + b*sc
    e_cal = float(np.sqrt(np.mean((cal-dg)**2))); rmse['cal'].append(e_cal)
    frho = spearmanr(fv, fdg).correlation; macro['fep'].append(frho)
    e_fep = float(np.sqrt(np.mean((fv-fdg)**2))); rmse['fep'].append(e_fep)

    def draw(ax, x, y, m, col):
        ax.scatter(x[~m], y[~m], s=22, c=col, alpha=.85, lw=0)
        ax.scatter(x[m], y[m], s=30, facecolors='none', edgecolors=col, lw=1.1)
    def box(ax, txt):
        ax.text(.03,.97,txt, transform=ax.transAxes, va='top', ha='left', fontsize=7.4,
                family='monospace', bbox=dict(fc='white', ec='0.6', lw=.6, pad=3))
    def diag(ax):
        lo=min(ax.get_xlim()[0],ax.get_ylim()[0]); hi=max(ax.get_xlim()[1],ax.get_ylim()[1])
        ax.plot([lo,hi],[lo,hi],ls='--',c='0.45',lw=1); ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)

    ax = axes[0, j]; draw(ax, dg, sc, cen, '#C0504D')
    xs = np.linspace(dg.min(), dg.max(), 2); ax.plot(xs, np.polyval(np.polyfit(dg, sc, 1), xs), c='#C0504D', lw=1.3)
    box(ax, f"rho = {rho:+.3f}\nR   = {R:+.3f}\noffset {np.mean(sc)-np.mean(dg):.1f}")
    ax.set_title(f"{t.replace('_','/')}  (n={len(sub)})", fontweight='bold')

    ax = axes[1, j]; draw(ax, dg, cal, cen, '#4C78A8')
    box(ax, f"rho = {rho:+.3f}\nRMSE {e_cal:.2f}\nfitted slope\n  {b:.3f}  (=1/{1/b:.1f})"); diag(ax)

    ax = axes[2, j]; draw(ax, fdg, fv, fcen, '#4F8F4F')
    box(ax, f"rho = {frho:+.3f}\nRMSE {e_fep:.2f}\n(no calibration)"); diag(ax)
    ax.set_xlabel('experimental $\\Delta G$ (kcal/mol)')

axes[0,0].set_ylabel('RAW MM/PBSA $\\Delta G_{calc}$ (kcal/mol)\ny=x is off-scale: no entropy term')
axes[1,0].set_ylabel('MM/PBSA after a per-family\nlinear calibration (kcal/mol)')
axes[2,0].set_ylabel('published FEP+ $\\Delta G$\n(kcal/mol)')
from matplotlib.lines import Line2D
axes[2,-1].legend(handles=[Line2D([],[],marker='o',ls='',mfc='0.3',mec='0.3',ms=5,label='measured'),
                           Line2D([],[],marker='o',ls='',mfc='none',mec='0.3',ms=5.5,label='censored (assay ceiling)')],
                  loc='lower right', fontsize=7, frameon=False, handletextpad=.3)
d = f"; {dropped[0][0]} excluded: unrelieved clash (VDW +{dropped[0][1]:,.0f} kcal/mol)" if dropped else ""
fig.suptitle(f"Minimised-structure MM/PBSA (PB, $\\epsilon_{{in}}$=2, ff03, no entropy) vs published FEP+ "
             f"— 5 FEP+ families, n={len(rows)}\n"
             f"macro-mean Spearman   MM/PBSA {np.mean(macro['s']):+.3f}    FEP+ {np.mean(macro['fep']):+.3f}"
             f"      accuracy   MM/PBSA {np.mean(rmse['cal']):.2f} kcal/mol only after a per-family fit"
             f"    FEP+ {np.mean(rmse['fep']):.2f} kcal/mol with no fitting{d}", fontsize=9.5)
fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig(OUT, dpi=170)
print(f"wrote {OUT}")
print(f"  macro Spearman  MM/PBSA {np.mean(macro['s']):+.3f}   FEP+ {np.mean(macro['fep']):+.3f}")
print(f"  macro RMSE      calibrated {np.mean(rmse['cal']):.2f}   FEP+ {np.mean(rmse['fep']):.2f}")
for t,a,b_,c,d_ in zip(T,macro['s'],rmse['cal'],macro['fep'],rmse['fep']):
    print(f"    {t:<12} rho {a:+.3f}  calRMSE {b_:.2f}   FEP rho {c:+.3f}  FEP RMSE {d_:.2f}")
