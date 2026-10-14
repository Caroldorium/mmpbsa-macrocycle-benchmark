#!/usr/bin/env python3
"""Does the GPU minimiser substitute for the published pmemd.cuda step?

Compares MM/PBSA scores from sander (the published SD+CG protocol) and from
OpenMM (same force field and implicit solvent on the GPU) on every Dataset I
complex where both arms ran.  What matters is not the absolute agreement but
whether the two arms give the same correlation with pKd, since that is the only
quantity any conclusion rests on.
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr


def load(d, key="pdb"):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "*", "result.json"))):
        try:
            r = json.load(open(f))
            if r.get("status") == "ok":
                rows.append(r)
        except Exception:
            pass
    return pd.DataFrame(rows)


def main():
    a = load(sys.argv[1] if len(sys.argv) > 1 else "results_dataset1_ff03")
    b = load(sys.argv[2] if len(sys.argv) > 2 else "results_dataset1_ff03_omm")
    d1 = pd.read_csv(sys.argv[3] if len(sys.argv) > 3 else "../MM_PBGBSA-CP/dataset1.csv")
    d1 = d1.rename(columns={"PDB code": "pdb"})[["pdb", "pKd"]]
    if a.empty or b.empty:
        print("one arm has no results yet"); return
    cols = [c for c in a.columns if c.split("_")[0] in ("pb", "pb2", "gb1", "gb2", "gb5")
            and c.rsplit("_", 1)[-1] in ("1", "2", "4", "6", "8") and c in b.columns]
    m = a[["pdb"] + cols].merge(b[["pdb"] + cols], on="pdb", suffixes=("_sander", "_omm"))
    m = m.merge(d1, on="pdb")
    print(f"complexes with both minimisers: {len(m)}")
    if len(m) < 4:
        print("too few for a comparison"); return
    out = []
    for c in cols:
        x, y = m[f"{c}_sander"].astype(float), m[f"{c}_omm"].astype(float)
        k = x.notna() & y.notna()
        if k.sum() < 4:
            continue
        r_s = pearsonr(x[k], m.pKd[k])[0]
        r_o = pearsonr(y[k], m.pKd[k])[0]
        out.append(dict(arm=c, n=int(k.sum()),
                        mean_diff=round(float((y[k] - x[k]).mean()), 2),
                        sd_diff=round(float((y[k] - x[k]).std(ddof=1)), 2),
                        r_between=round(pearsonr(x[k], y[k])[0], 3),
                        rho_between=round(spearmanr(x[k], y[k])[0], 3),
                        Rp_sander=round(r_s, 3), Rp_openmm=round(r_o, 3),
                        Rp_shift=round(r_o - r_s, 3)))
    T = pd.DataFrame(out)
    print(T.to_string(index=False))
    T.to_csv("zha25_minimiser_comparison.csv", index=False)
    print(f"\nacross arms: median |Rp shift| = {T.Rp_shift.abs().median():.3f}, "
          f"max = {T.Rp_shift.abs().max():.3f}; "
          f"median rank agreement between minimisers rho = {T.rho_between.median():.3f}")


if __name__ == "__main__":
    main()
