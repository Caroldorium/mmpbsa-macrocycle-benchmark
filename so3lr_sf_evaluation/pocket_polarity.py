#!/usr/bin/env python3
"""Pocket-based polar SASA ratio -- the descriptor SO3LR-SF's authors report (r=0.81 vs skill).

The manuscript's s3.4 tested an INTERFACE polar fraction (polar dSASA / total dSASA on burial).
That is a different quantity from theirs: polar SASA / total SASA of the BINDING POCKET,
computed on the apo-like receptor surface rather than on what gets buried.  This recomputes
their version so the manuscript's "polarity does not explain where the method works" claim is
tested against the descriptor it is actually about.

Shrake-Rupley, implemented here so the script needs only numpy.
"""
import glob, json, os, sys
import numpy as np

PROBE, N_SPHERE, POCKET_CUT = 1.4, 256, 5.0
# Bondi radii
RAD = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'P':1.80,'F':1.47,
       'CL':1.75,'BR':1.85,'I':1.98,'H':1.20}
POLAR = {'N','O','S'}          # polar heavy atoms; S included as in common SASA conventions

def sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2*i/n)
    th  = np.pi * (1 + 5**0.5) * i
    return np.c_[np.cos(th)*np.sin(phi), np.sin(th)*np.sin(phi), np.cos(phi)]

SPH = sphere(N_SPHERE)

def read_pdb(path):
    xyz, el = [], []
    for ln in open(path):
        if not ln.startswith(('ATOM','HETATM')): continue
        e = ln[76:78].strip().upper()
        if not e:
            e = ''.join(c for c in ln[12:16].strip() if c.isalpha())[:1].upper()
        if e == 'H': continue                       # heavy atoms only
        xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
        el.append(e)
    return np.array(xyz), np.array(el)

def sasa(xyz, el, subset=None):
    """Per-atom SASA of `xyz` (all atoms occlude). Returns array over `subset` indices."""
    r = np.array([RAD.get(e, 1.70) for e in el]) + PROBE
    idx = np.arange(len(xyz)) if subset is None else np.asarray(subset)
    out = np.zeros(len(idx))
    for k, i in enumerate(idx):
        pts = xyz[i] + r[i]*SPH
        d2 = ((pts[:,None,:] - xyz[None,:,:])**2).sum(-1)
        d2[:, i] = np.inf
        out[k] = 4*np.pi*r[i]**2 * (d2 >= (r**2)[None,:]).all(1).mean()
    return out

rows = []
for d in sorted(glob.glob(sys.argv[1] + '/*')):
    rp, lp = f'{d}/receptor.pdb', f'{d}/ligand.pdb'
    if not (os.path.exists(rp) and os.path.exists(lp)): continue
    try:
        rx, re_ = read_pdb(rp)
        lx, _   = read_pdb(lp)
        if not len(rx) or not len(lx): continue
        # pocket = receptor heavy atoms within POCKET_CUT of any ligand heavy atom
        dmin = np.sqrt(((rx[:,None,:]-lx[None,:,:])**2).sum(-1)).min(1)
        pocket = np.where(dmin <= POCKET_CUT)[0]
        if len(pocket) < 10: continue
        a = sasa(rx, re_, pocket)                   # receptor alone, pocket atoms
        tot = a.sum()
        pol = a[np.isin(re_[pocket], list(POLAR))].sum()
        rows.append(dict(slug=os.path.basename(d), n_pocket=int(len(pocket)),
                         pocket_sasa=round(float(tot),1),
                         pocket_polar_sasa=round(float(pol),1),
                         polar_ratio=round(float(pol/tot),4) if tot else None))
        print(rows[-1], flush=True)
    except Exception as e:
        print(f'{os.path.basename(d)}: FAILED {e}', flush=True)

import csv
with open(sys.argv[2],'w',newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f'\nwrote {sys.argv[2]}  n={len(rows)}')
