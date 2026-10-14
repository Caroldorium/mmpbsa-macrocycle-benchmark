#!/usr/bin/env python3
"""What Zha25's estimator says about the IL23R benchmark.

The IL23R patent digest has no structures, so the MM/PBSA protocol cannot be run on
it.  What can be done, and is the part of the question the data supports, is to score
the models this project already fitted with the estimator the paper reports -- a
Pearson correlation pooled over every system in the evaluation -- and put it next to
the within-series estimator the IL23R benchmark reports.
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr

SRC = sys.argv[1]
OUT = sys.argv[2]
RNG = np.random.default_rng(20260924)
NBOOT = 4000


def pooled(s):
    return pearsonr(s.prediction, s.experimental_pic50)[0] if s.prediction.nunique() > 1 else np.nan


def fold_pooled(s):
    """The benchmark's own estimator: one Spearman over each fold's whole test set."""
    v = [spearmanr(g.prediction, g.experimental_pic50)[0] for _, g in s.groupby("fold")
         if g.prediction.nunique() > 1 and g.experimental_pic50.nunique() > 1]
    return float(np.mean(v)) if v else np.nan


def within(s, min_n=3):
    """Series weighted equally inside a fold, then folds equally -- the estimator the
    de novo benchmark in this project uses, and the one a lead-optimisation campaign
    actually experiences."""
    per_fold, nser = [], 0
    for _, f in s.groupby("fold"):
        rs = [spearmanr(g.prediction, g.experimental_pic50)[0]
              for _, g in f.groupby("family_id")
              if len(g) >= min_n and g.prediction.nunique() > 1
              and g.experimental_pic50.nunique() > 1]
        if rs:
            per_fold.append(float(np.mean(rs))); nser += len(rs)
    return (float(np.mean(per_fold)) if per_fold else np.nan), nser


def boot_pooled(s, n=NBOOT):
    """Cluster the resample on series, the way the pooled statistic's dependence runs."""
    fams = s.family_id.unique()
    idx = {f: s.index[s.family_id == f].to_numpy() for f in fams}
    v = []
    for _ in range(n):
        sel = np.concatenate([idx[f] for f in RNG.choice(fams, len(fams), replace=True)])
        g = s.loc[sel]
        if g.prediction.nunique() > 2:
            v.append(pearsonr(g.prediction, g.experimental_pic50)[0])
    v = np.array([x for x in v if np.isfinite(x)])
    return np.percentile(v, [2.5, 97.5]) if len(v) > 100 else (np.nan, np.nan)


def main():
    d = pd.read_csv(os.path.join(SRC, "patent_v8_out_of_fold_predictions.csv"))
    d = d[d.label_type == "exact"].dropna(subset=["experimental_pic50", "prediction"])
    desc = pd.read_csv(os.path.join(SRC, "patent_v8_monomer_descriptor_sums.csv"))
    rows = []
    for (proto, arm), s in d.groupby(["protocol", "arm"]):
        if arm == "Training mean":
            continue
        w, k = within(s)
        lo, hi = boot_pooled(s)
        rows.append(dict(protocol=proto, arm=arm, n=len(s), series=k,
                         pooled_Rp=round(pooled(s), 3),
                         pooled_lo=round(lo, 3), pooled_hi=round(hi, 3),
                         benchmark_pooled_rho=round(fold_pooled(s), 3),
                         within_series_rho=round(w, 3),
                         inflation=round(abs(fold_pooled(s)) - abs(w), 3)))
    # the same two estimators for a molecular-size control on the same entities
    base = d[["entity_id", "family_id", "protocol", "experimental_pic50"]].drop_duplicates()
    b = base.merge(desc, on="entity_id", how="left").dropna(
        subset=["sum_free_monomer_heavy_atom_count"])
    for proto, s in b.groupby("protocol"):
        s = s.rename(columns={"sum_free_monomer_heavy_atom_count": "prediction"})
        s = s.assign(fold=0)
        w, k = within(s)
        rows.append(dict(protocol=proto, arm="heavy atoms (size control)", n=len(s), series=k,
                         pooled_Rp=round(pooled(s), 3), pooled_lo=None, pooled_hi=None,
                         benchmark_pooled_rho=round(fold_pooled(s), 3),
                         within_series_rho=round(w, 3),
                         inflation=round(abs(fold_pooled(s)) - abs(w), 3)))
    T = pd.DataFrame(rows)
    for p in T.protocol.unique():
        print(f"\n=== {p} ===")
        print(T[T.protocol == p].drop(columns="protocol")
              .sort_values("pooled_Rp", key=abs, ascending=False).to_string(index=False))
    T.to_csv(os.path.join(OUT, "zha25_il23r_estimator_comparison.csv"), index=False)
    print(f"\nwrote {OUT}/zha25_il23r_estimator_comparison.csv")


if __name__ == "__main__":
    main()
