#!/usr/bin/env python3
"""Zha25 MM/PBSA protocol on the five-family FEP+ cohort.

Reported under both estimators: the paper's pooled cross-system Pearson, and the
within-family macro Spearman this project deploys against.  Every arm is shown
next to the molecular-size control and the published FEP+ comparator.
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr

RES = sys.argv[1] if len(sys.argv) > 1 else "results_cohort"
MAN = sys.argv[2] if len(sys.argv) > 2 else "cohort_manifest.csv"
OUT = sys.argv[3] if len(sys.argv) > 3 else "."
RNG = np.random.default_rng(20260924)
NBOOT = int(os.environ.get("NBOOT", 20000))


def load():
    rows = []
    for f in sorted(glob.glob(os.path.join(RES, "*", "result.json"))):
        try:
            rows.append(json.load(open(f)))
        except Exception:
            pass
    r = pd.DataFrame(rows)
    man = pd.read_csv(MAN)
    m = r.merge(man, on="slug", how="right")
    return m[m.experimental_dG_kcal_mol.notna()].copy()


def macro_rho(df, col, min_n=3):
    per = []
    for t, s in df.groupby("target"):
        s = s[s[col].notna()]
        if len(s) < min_n or s[col].nunique() < 2:
            continue
        per.append((t, len(s), spearmanr(s[col], s.experimental_dG_kcal_mol)[0]))
    if not per:
        return np.nan, per
    return float(np.mean([p[2] for p in per])), per


def family_boot(df, col, other=None, n=NBOOT):
    """Resample families; within a resampled family, resample compounds."""
    fams = df.target.unique()
    idx = {t: df.index[df.target == t].to_numpy() for t in fams}
    out = []
    for _ in range(n):
        pick = RNG.choice(fams, size=len(fams), replace=True)
        rs = []
        rs2 = []
        for t in pick:
            sel = RNG.choice(idx[t], size=len(idx[t]), replace=True)
            s = df.loc[sel]
            if s[col].nunique() < 2:
                continue
            rs.append(spearmanr(s[col], s.experimental_dG_kcal_mol)[0])
            if other is not None and s[other].nunique() > 1:
                rs2.append(spearmanr(s[other], s.experimental_dG_kcal_mol)[0])
        if not rs:
            continue
        if other is None:
            out.append(np.mean(rs))
        elif rs2:
            out.append(abs(np.mean(rs)) - abs(np.mean(rs2)))
    v = np.array([x for x in out if np.isfinite(x)])
    return np.percentile(v, [2.5, 97.5]) if len(v) > 100 else (np.nan, np.nan)


def top_k_hit(df, col, k=5):
    """Fraction of families whose true best compound is in the model's top k.

    The ranking direction is taken from the sign of the within-family
    correlation, so a descriptor that predicts affinity inversely -- molecular
    size does -- is still given its best shot rather than its worst."""
    sign = np.sign(macro_rho(df, col)[0]) or 1.0
    hits, n = 0, 0
    for t, s in df.groupby("target"):
        s = s[s[col].notna()]
        if len(s) <= k:
            continue
        n += 1
        order = s.assign(_k=sign * s[col]).sort_values("_k", ascending=True)
        # eleven of the 45 KRAS labels sit on a censored assay ceiling and are
        # numerically identical, so "the best compound" is a set, not a row
        best = set(s.index[s.experimental_dG_kcal_mol == s.experimental_dG_kcal_mol.min()])
        if best & set(order.head(k).index):
            hits += 1
    return hits / n if n else np.nan, n


def main():
    m = load()
    built = m[m.get("status", pd.Series(dtype=object)) == "ok"] if "status" in m else m
    print(f"cohort rows with a label: {len(m)}   scored: {len(built)}   families: {m.target.nunique()}")
    if "status" in m:
        print("build status:", m.status.value_counts(dropna=False).to_dict())
    if "pose_drift_A" in m:
        d = pd.to_numeric(m.pose_drift_A, errors="coerce")
        print(f"pose drift after parameterisation: max {d.max():.4f} A (must be ~0)")
    ok = built.copy()
    if not len(ok):
        print("nothing scored yet"); return
    score_cols = sorted([c for c in ok.columns
                         if c.split("_")[0] in ("pb", "pb2", "gb1", "gb2", "gb5")
                         and c.rsplit("_", 1)[-1] in ("1", "2", "4", "6", "8")])
    ok.to_csv(os.path.join(OUT, "zha25_cohort_scores.csv"), index=False)

    arms = score_cols + ["heavy_atoms", "published_fepplus_dG_kcal_mol"]
    rows = []
    for c in arms:
        s = ok[ok[c].notna()]
        if len(s) < 20:
            continue
        pooled_r = pearsonr(s[c], s.experimental_dG_kcal_mol)[0]
        pooled_rho = spearmanr(s[c], s.experimental_dG_kcal_mol)[0]
        mac, per = macro_rho(s, c)
        lo, hi = family_boot(s, c)
        t5, n5 = top_k_hit(s, c, 5)
        rows.append(dict(arm=c, n=len(s), pooled_Rp=round(pooled_r, 3),
                         pooled_rho=round(pooled_rho, 3), macro_rho=round(mac, 3),
                         macro_lo=round(lo, 3), macro_hi=round(hi, 3),
                         top5_best_recovered=round(t5, 3) if np.isfinite(t5) else None,
                         families=len(per)))
    if not rows:
        print("\nno arm has enough scored compounds yet"); return
    T = pd.DataFrame(rows).sort_values("macro_rho", key=abs, ascending=False)
    print("\n=== all arms, both estimators ===")
    print(T.to_string(index=False))
    T.to_csv(os.path.join(OUT, "zha25_cohort_metrics.csv"), index=False)

    if "censored" in ok.columns:
        nc = ok[ok.censored != True]
        print(f"\n=== repeated with the {int((ok.censored == True).sum())} censored labels dropped "
              f"({len(nc)} compounds) ===")
        for c in ["pb_2", "heavy_atoms", "published_fepplus_dG_kcal_mol"]:
            if c not in nc.columns:
                continue
            sub = nc[nc[c].notna()]
            if len(sub) < 10:
                continue
            mm, per = macro_rho(sub, c)
            print(f"  {c:<32} macro rho={mm:+.3f} over {len(per)} families")

    phys = [r for r in rows if r["arm"].startswith(("pb", "gb"))]
    if phys:
        # pb_2 is the paper's own optimum -- MM/PBSA, LCPO nonpolar, interior
        # dielectric 2 -- chosen by them on their data, so it is the one arm here
        # that is not selected on this cohort.  Everything else is exploratory.
        designated = next((r for r in rows if r["arm"] == "pb_2"), None)
        if designated:
            print(f"\ndesignated arm (the paper's optimum, MM/PBSA epsin=2): "
                  f"macro rho={designated['macro_rho']:+.3f} "
                  f"[{designated['macro_lo']:+.3f}, {designated['macro_hi']:+.3f}]  "
                  f"pooled Rp={designated['pooled_Rp']:+.3f}")
        best = max(phys, key=lambda r: abs(r["macro_rho"]))
        print(f"best arm chosen on this cohort ({len(phys)} tried): {best['arm']}  "
              f"macro rho={best['macro_rho']:+.3f}  -- selected here, so optimistic")
        s = ok[ok[best["arm"]].notna()]
        lo, hi = family_boot(s, best["arm"], other="heavy_atoms")
        print(f"|macro rho| advantage over heavy atoms: "
              f"{abs(best['macro_rho'])-abs(macro_rho(s,'heavy_atoms')[0]):+.3f} "
              f"(95% CI {lo:+.3f} to {hi:+.3f})")

        # ---- which term carries the signal
        print("\n=== energy decomposition of the best arm ===")
        base = best["arm"]
        terms = [c for c in ok.columns if c.startswith(base + "_")]
        for t in ["heavy_atoms"] + terms:
            if t not in ok.columns:
                continue
            ss = ok[ok[t].notna()]
            if len(ss) < 20 or ss[t].nunique() < 3:
                continue
            mm, _ = macro_rho(ss, t)
            print(f"  {t:<28} macro rho={mm:+.3f}   pooled Rp={pearsonr(ss[t],ss.experimental_dG_kcal_mol)[0]:+.3f}")

        # ---- size partialled out of the physics score, within each family
        print("\n=== physics score with heavy atoms partialled out, within family ===")
        import numpy.linalg as la
        per = []
        for t, ssub in ok.groupby("target"):
            ssub = ssub[ssub[base].notna()]
            if len(ssub) < 5:
                continue
            X = np.column_stack([np.ones(len(ssub)), ssub.heavy_atoms.to_numpy()])
            res_s = ssub[base].to_numpy() - X @ la.lstsq(X, ssub[base].to_numpy(), rcond=None)[0]
            res_y = (ssub.experimental_dG_kcal_mol.to_numpy()
                     - X @ la.lstsq(X, ssub.experimental_dG_kcal_mol.to_numpy(), rcond=None)[0])
            per.append((t, len(ssub), spearmanr(res_s, res_y)[0]))
        if per:
            print(f"  partial macro rho = {np.mean([p[2] for p in per]):+.3f}")
            for t, n, r in per:
                print(f"      {t:<12} n={n:3d}  {r:+.3f}")
        print("\nper-family Spearman:")
        _, per = macro_rho(s, best["arm"])
        _, perh = macro_rho(s, "heavy_atoms")
        hh = {p[0]: p[2] for p in perh}
        for t, n, r in per:
            print(f"  {t:<12} n={n:3d}  MM/PBSA {r:+.3f}   heavy atoms {hh.get(t,float('nan')):+.3f}")


if __name__ == "__main__":
    main()
