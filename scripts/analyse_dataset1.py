#!/usr/bin/env python3
"""Analysis of the Zha25 Dataset I replication.

Reproduces the published Table 1 implicit-solvent column, then applies the three
controls the paper does not: an interface-size baseline on the same complexes,
receptor-group clustering (50 complexes are 21 receptors), and a within-receptor
ranking estimator.
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr

RES = sys.argv[1] if len(sys.argv) > 1 else "results_dataset1"
REPO = sys.argv[2] if len(sys.argv) > 2 else "../MM_PBGBSA-CP"
GEOM = sys.argv[3] if len(sys.argv) > 3 else "dataset1_geom.csv"
OUT = sys.argv[4] if len(sys.argv) > 4 else "."
RNG = np.random.default_rng(20260924)
NBOOT = int(os.environ.get("NBOOT", 20000))

SIZE_FEATURES = ["pep_heavy", "pep_res", "dSASA", "dSASA_apolar", "contacts_4A",
                 "contacts_5A", "iface_pep_atoms", "iface_prot_atoms"]


def load():
    rows = []
    for f in sorted(glob.glob(os.path.join(RES, "*", "result.json"))):
        try:
            rows.append(json.load(open(f)))
        except Exception:
            pass
    r = pd.DataFrame(rows)
    d1 = pd.read_csv(os.path.join(REPO, "dataset1.csv")).rename(
        columns={"PDB code": "pdb", "Peptide length": "pep_len",
                 "Polar interface (%)": "polar_pct", "Chain type": "chain_type"})
    g = pd.read_csv(GEOM)
    m = r.merge(d1[["pdb", "pep_len", "polar_pct", "pKd", "chain_type"]], on="pdb")
    m = m.merge(g, on="pdb", suffixes=("", "_g"))
    d2 = set(pd.read_csv(os.path.join(REPO, "dataset2.csv"))["PDB code"])
    m["in_decoy"] = m.pdb.isin(d2)
    grp = json.load(open("dataset1_receptor_groups.json"))
    lut = {p: int(k) for k, v in grp.items() for p in v}
    m["receptor_group"] = m.pdb.map(lut)
    return m


def cluster_boot_ci(df, col, stat, n=NBOOT):
    """Bootstrap over receptor groups, not complexes: 50 complexes are 21 receptors."""
    groups = df.receptor_group.unique()
    idx = {g: df.index[df.receptor_group == g].to_numpy() for g in groups}
    vals = []
    for _ in range(n):
        pick = RNG.choice(groups, size=len(groups), replace=True)
        sel = np.concatenate([idx[g] for g in pick])
        sub = df.loc[sel]
        if sub[col].nunique() < 3:
            continue
        try:
            vals.append(stat(sub[col].to_numpy(), sub.pKd.to_numpy()))
        except Exception:
            pass
    v = np.array([x for x in vals if np.isfinite(x)])
    if len(v) < 100:
        return (np.nan, np.nan)
    return tuple(np.percentile(v, [2.5, 97.5]))


def rp(x, y):
    return pearsonr(x, y)[0]


def within_receptor_macro(df, col, min_n=3):
    """One Spearman per receptor group, groups weighted equally."""
    out = []
    for g, sub in df.groupby("receptor_group"):
        if len(sub) < min_n or sub[col].nunique() < 2 or sub.pKd.nunique() < 2:
            continue
        out.append((g, len(sub), spearmanr(sub[col], sub.pKd)[0]))
    if not out:
        return np.nan, 0, []
    return float(np.mean([o[2] for o in out])), len(out), out


def main():
    m = load()
    ok = m[m.status == "ok"].copy()
    print(f"complexes with a score: {len(ok)}/{len(m)}   receptor groups: {ok.receptor_group.nunique()}")
    score_cols = sorted([c for c in ok.columns
                         if c.split("_")[0] in ("pb", "pb2", "gb1", "gb2", "gb5")
                         and c.split("_")[-1] in ("1", "2", "4", "6", "8")])
    ok.to_csv(os.path.join(OUT, "zha25_dataset1_scores.csv"), index=False)

    # ---- 1. reproduction: the published Table 1 implicit-solvent column
    print("\n=== Table 1 replication (implicit solvent, ff03): Rp against pKd ===")
    grid = []
    for c in score_cols:
        s = ok[c].astype(float)
        keep = s.notna()
        if keep.sum() < 10:
            continue
        meth, indi = c.rsplit("_", 1)
        r = rp(s[keep], ok.pKd[keep])
        rho = spearmanr(s[keep], ok.pKd[keep])[0]
        lo, hi = cluster_boot_ci(ok[keep], c, rp)
        grid.append(dict(method=meth, epsin=int(indi), n=int(keep.sum()),
                         Rp=round(r, 3), rho=round(rho, 3),
                         Rp_lo=round(lo, 3), Rp_hi=round(hi, 3)))
    G = pd.DataFrame(grid)
    if not G.empty:
        print(G.pivot(index="epsin", columns="method", values="Rp").to_string())
        G.to_csv(os.path.join(OUT, "zha25_dataset1_grid.csv"), index=False)

    # ---- 2. size control on the same complexes
    print("\n=== interface-size control (same 50 complexes) ===")
    ctrl = []
    for f in SIZE_FEATURES + ["pep_len", "polar_pct"]:
        r = rp(ok[f], ok.pKd)
        lo, hi = cluster_boot_ci(ok, f, rp)
        ctrl.append(dict(descriptor=f, Rp=round(r, 3), Rp_lo=round(lo, 3), Rp_hi=round(hi, 3),
                         rho=round(spearmanr(ok[f], ok.pKd)[0], 3)))
    C = pd.DataFrame(ctrl).sort_values("Rp", key=abs, ascending=False)
    print(C.to_string(index=False))
    C.to_csv(os.path.join(OUT, "zha25_dataset1_size_control.csv"), index=False)

    # ---- 3. paired contrast: best physics arm minus best size descriptor
    if not G.empty:
        # pb_2 is the paper's own optimum, chosen by them on their data; everything
        # else in the grid is being read after the fact and is optimistic
        des = G[(G.method == "pb") & (G.epsin == 2)]
        if len(des):
            print(f"\ndesignated arm (the paper's optimum, MM/PBSA epsin=2): "
                  f"Rp = {des.iloc[0].Rp:+.3f}  [{des.iloc[0].Rp_lo:+.3f}, {des.iloc[0].Rp_hi:+.3f}]"
                  f"   published: -0.545")
        best = G.loc[G.Rp.abs().idxmax()]
        bc = f"{best.method}_{best.epsin}"
        bs = C.iloc[0].descriptor
        print(f"best arm read off this grid: {bc}  Rp={best.Rp}  -- selected here, so optimistic")
        print(f"best size descriptor: {bs}  Rp={C.iloc[0].Rp}")
        groups = ok.receptor_group.unique()
        idx = {g: ok.index[ok.receptor_group == g].to_numpy() for g in groups}
        diffs = []
        for _ in range(NBOOT):
            pick = RNG.choice(groups, size=len(groups), replace=True)
            sub = ok.loc[np.concatenate([idx[g] for g in pick])]
            if sub[bc].nunique() < 3:
                continue
            try:
                diffs.append(abs(rp(sub[bc], sub.pKd)) - abs(rp(sub[bs], sub.pKd)))
            except Exception:
                pass
        d = np.array(diffs)
        print(f"|Rp| advantage of physics over size: {abs(best.Rp)-abs(C.iloc[0].Rp):+.3f} "
              f"(95% CI {np.percentile(d,2.5):+.3f} to {np.percentile(d,97.5):+.3f}, "
              f"P(>0)={np.mean(d>0):.3f})")

        # ---- 4. within-receptor ranking, the deployment estimator
        print("\n=== within-receptor Spearman, receptor groups weighted equally ===")
        for col, label in (("pb_2", "MM/PBSA epsin=2"), (bc, f"best arm {bc}"), (bs, f"size ({bs})")):
            mac, k, per = within_receptor_macro(ok, col)
            print(f"  {label:<24} macro rho = {mac:+.3f} over {k} receptor groups")
            for g, n, r_ in sorted(per):
                print(f"      group {g:2d} (n={n}): {r_:+.3f}")

        # ---- 5. the decoy subset the two-step workflow was reported on
        sub = ok[ok.in_decoy]
        print(f"\n=== Dataset II overlap (n={len(sub)}) ===")
        print(f"  pb_2: Rp={rp(sub['pb_2'],sub.pKd):+.3f}   {bc}: Rp={rp(sub[bc],sub.pKd):+.3f}   "
              f"{bs}: Rp={rp(sub[bs],sub.pKd):+.3f}")

        # ---- 6. high / low polarity split
        print("\n=== polarity split at 0.6 (their Fig. 1) ===")
        for label, s in (("high polarity", ok[ok.polar_pct > 60]),
                         ("low polarity", ok[ok.polar_pct <= 60])):
            if len(s) > 4:
                print(f"  {label:<15} n={len(s):2d}  pb_2: Rp={rp(s['pb_2'],s.pKd):+.3f}   "
                      f"{bc}: Rp={rp(s[bc],s.pKd):+.3f}   {bs}: Rp={rp(s[bs],s.pKd):+.3f}")


if __name__ == "__main__":
    main()
