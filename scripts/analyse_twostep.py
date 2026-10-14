#!/usr/bin/env python3
"""The two-step workflow on the cohort, reported the way Zha25 report it.

Their Figure 6 contrasts three things: the docking program's own top pose (−0.316),
the crystal complex as an upper limit (−0.813), and the two-step workflow (−0.732).
The equivalents here are a candidate pose picked at random, the deposited pose scored
directly, and the pose their step 1 selects then rescored by their step 2.

Both estimators are reported for every arm: their pooled Pearson over all systems, and
the within-family macro Spearman this project deploys against.
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr

TS = sys.argv[1] if len(sys.argv) > 1 else "results_cohort_twostep"
SINGLE = sys.argv[2] if len(sys.argv) > 2 else "results_cohort_ff03_omm"
MAN = sys.argv[3] if len(sys.argv) > 3 else "cohort_manifest.csv"
OUT = sys.argv[4] if len(sys.argv) > 4 else "."
RNG = np.random.default_rng(20260924)
NBOOT = int(os.environ.get("NBOOT", 4000))


def load(d, keep=None):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "*", "result.json"))):
        try:
            r = json.load(open(f))
        except Exception:
            continue
        if r.get("status") != "ok":
            continue
        rows.append({k: v for k, v in r.items() if keep is None or k in keep})
    return pd.DataFrame(rows)


def macro(df, col, min_n=3):
    per = [(t, len(s), spearmanr(s[col], s.experimental_dG_kcal_mol)[0])
           for t, s in df.groupby("target")
           if len(s) >= min_n and s[col].nunique() > 1]
    return (float(np.mean([p[2] for p in per])) if per else np.nan), per


def fam_boot(df, col, n=NBOOT):
    fams = df.target.unique()
    idx = {t: df.index[df.target == t].to_numpy() for t in fams}
    v = []
    for _ in range(n):
        rs = []
        for t in RNG.choice(fams, len(fams), replace=True):
            s = df.loc[RNG.choice(idx[t], len(idx[t]), replace=True)]
            if s[col].nunique() > 1:
                rs.append(spearmanr(s[col], s.experimental_dG_kcal_mol)[0])
        if rs:
            v.append(np.mean(rs))
    v = np.array([x for x in v if np.isfinite(x)])
    return np.percentile(v, [2.5, 97.5]) if len(v) > 100 else (np.nan, np.nan)


def main():
    man = pd.read_csv(MAN)
    ts = load(TS)
    sp = load(SINGLE)
    if ts.empty:
        print("no two-step results yet"); return
    sp = sp[["slug", "pb_2"]].rename(columns={"pb_2": "deposited_pose"})
    d = ts.merge(sp, on="slug", how="left").merge(man, on="slug", how="left")
    d = d[d.experimental_dG_kcal_mol.notna()].copy()
    print(f"compounds with a two-step result and a label: {len(d)} / {d.target.nunique()} families")

    print("\n=== what step 1 actually selected ===")
    print(f"  candidates per compound: median {d.n_candidates.median():.0f}")
    print(f"  decoy set ligand RMSD after minimisation: "
          f"{d.decoy_rmsd_min.median():.2f} to {d.decoy_rmsd_max.median():.2f} A (medians)")
    dep = (d.selected_pose == 0).mean()
    print(f"  deposited pose selected as top-1: {100*dep:.1f}% of compounds")
    print(f"  deposited pose's rank among candidates: median {d.native_rank.median():.0f}")
    print(f"  selected pose RMSD from deposited: median {d.selected_rmsd.median():.2f} A, "
          f"90th pct {d.selected_rmsd.quantile(.9):.2f} A")
    near = (d.selected_rmsd <= 2.5).mean()
    print(f"  selected pose within 2.5 A (their near-native criterion): {100*near:.1f}%")

    arms = [("deposited_pose", "deposited pose scored directly (their crystal upper limit)"),
            ("two_step_score", "two-step: step-1 selection then step-2 rescoring"),
            ("random_pose_score", "a candidate pose picked at random (control)"),
            ("heavy_atoms", "molecular size control"),
            ("published_fepplus_dG_kcal_mol", "published FEP+")]
    rows = []
    for c, label in arms:
        if c not in d.columns:
            continue
        s = d[d[c].notna()]
        if len(s) < 20:
            continue
        m, per = macro(s, c)
        lo, hi = fam_boot(s, c)
        rows.append(dict(arm=label, n=len(s),
                         pooled_Rp=round(pearsonr(s[c], s.experimental_dG_kcal_mol)[0], 3),
                         pooled_rho=round(spearmanr(s[c], s.experimental_dG_kcal_mol)[0], 3),
                         macro_rho=round(m, 3), lo=round(lo, 3), hi=round(hi, 3),
                         families=len(per)))
    T = pd.DataFrame(rows)
    print("\n=== both estimators, all arms ===")
    print(T.to_string(index=False))
    T.to_csv(os.path.join(OUT, "zha25_cohort_twostep_metrics.csv"), index=False)

    print("\n=== per family, two-step against the deposited pose ===")
    _, pt = macro(d[d.two_step_score.notna()], "two_step_score")
    _, pd_ = macro(d[d.deposited_pose.notna()], "deposited_pose")
    dd = {p[0]: p[2] for p in pd_}
    for t, n, r in pt:
        print(f"  {t:<12} n={n:3d}  two-step {r:+.3f}   deposited {dd.get(t, np.nan):+.3f}")
    d.drop(columns=[c for c in ("candidates",) if c in d.columns]).to_csv(
        os.path.join(OUT, "zha25_cohort_twostep_scores.csv"), index=False)


if __name__ == "__main__":
    main()
