#!/usr/bin/env python3
"""Build the executed notebook for the Zha25 replication goal."""
import json, sys, nbformat as nbf

def md(t): return nbf.v4.new_markdown_cell(t)
def co(t): return nbf.v4.new_code_cell(t)

cells = []

cells.append(md("""# Replicating Zhao et al. 2025 fine-tuned MM/PBSA on our benchmark

Zhao H, Huang J, Weng G, Jiang D, Hu R, Kang Y, Hou T.
*Improving the predictive performance of binding affinities and poses for
protein–cyclic peptide complexes through fine-tuned MM/PBSA(GBSA)-based methods.*
Briefings in Bioinformatics 2025;26(6):bbaf632. Data and script:
https://github.com/huifengzhao/MM_PBGBSA-CP

The paper reports Pearson **Rp = −0.545** between MM/PBSA and experimental pKd over
50 protein–cyclic-peptide crystal complexes, and **−0.732** for a two-step workflow on
a 24-complex decoy subset. This project measured MM/PBSA at Spearman −0.101 on its
20-compound KRAS panel. This notebook rebuilds their protocol, runs the controls the
paper does not, and applies the same protocol to the five-family FEP+ cohort.

Everything here is retrospective scoring of existing structures. No model is promoted."""))

cells.append(co('''import glob, json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

RES  = os.environ.get("ZHA_RESULTS", ".")
REPO = os.environ.get("ZHA_REPO", ".")
RNG = np.random.default_rng(20260924)\nNBOOT = int(os.environ.get("ZHA_NBOOT", 20000))

def load(arm):
    rows = []
    for f in sorted(glob.glob(os.path.join(RES, f"results_{arm}", "*", "result.json"))):
        try:
            r = json.load(open(f))
            if r.get("status") == "ok": rows.append(r)
        except Exception: pass
    return pd.DataFrame(rows)

def score_cols(df):
    return sorted([c for c in df.columns
                   if c.split("_")[0] in ("pb","pb2","gb1","gb2","gb5")
                   and c.rsplit("_",1)[-1] in ("1","2","4","6","8")])

d1_omm    = load("dataset1_ff03_omm")
d1_sander = load("dataset1_ff03")
cohort    = load("cohort_ff03_omm")
print(f"Dataset I, GPU minimiser : {len(d1_omm)} complexes")
print(f"Dataset I, sander        : {len(d1_sander)} complexes (validation sample)")
print(f"five-family cohort       : {len(cohort)} compounds")'''))

cells.append(md("""## 1. Does the published protocol reproduce?

Their optimum is MM/PBSA with interior dielectric 2 on structures energy-minimised in
implicit solvent under ff03, entropy excluded. The grid below is the same
force field over five interior dielectrics and five solvation treatments, scored on the
same 50 deposited crystal complexes, next to the published Table 1 values."""))

cells.append(co('''meta = pd.read_csv(os.path.join(REPO, "dataset1.csv")).rename(
    columns={"PDB code":"pdb","Peptide length":"pep_len","Polar interface (%)":"polar_pct",
             "Chain type":"chain_type"})
geom = pd.read_csv(os.path.join(RES, "dataset1_geom.csv"))
groups = json.load(open(os.path.join(RES, "dataset1_receptor_groups.json")))
lut = {p:int(k) for k,v in groups.items() for p in v}

# result.json also carries pep_res, so keep the geometry file's copy and drop the other
D = d1_omm.drop(columns=[c for c in ("pep_res",) if c in d1_omm.columns])
D = D.merge(meta[["pdb","pep_len","polar_pct","pKd","chain_type"]], on="pdb").merge(geom, on="pdb")
D["receptor_group"] = D.pdb.map(lut)
D["in_decoy"] = D.pdb.isin(set(pd.read_csv(os.path.join(REPO,"dataset2.csv"))["PDB code"]))

rows = []
for c in score_cols(D):
    m, e = c.rsplit("_", 1)
    rows.append(dict(method=m, epsin=int(e), Rp=round(pearsonr(D[c], D.pKd)[0], 3)))
grid = pd.DataFrame(rows).pivot(index="epsin", columns="method", values="Rp")
pub = pd.read_csv(os.path.join(RES, "zha25_published_table1_ff03.csv"))
pub = pub[pub.solvent=="implicit"].pivot(index="epsin", columns="method", values="published_Rp")
print("this replication (n=%d):" % len(D)); print(grid.to_string())
print("\\npublished Table 1, ff03, implicit solvent (n=50):"); print(pub.to_string())
print("\\ndifference (replication minus published):")
print((grid[["gb1","gb2","gb5","pb"]] - pub[["gb1","gb2","gb5","pb"]]).round(3).to_string())'''))

cells.append(md("""## 2. Is the GPU minimiser a valid substitute for the published pmemd.cuda step?

The released protocol minimises with `pmemd.cuda`; the conda AmberTools build ships no
pmemd, and `sander` needs hours per complex at the protocol's infinite cutoff. The
substitute runs the same force field and implicit solvent model on the GPU, capped at
500 L-BFGS iterations — a cap calibrated so the structural drift matches what sander's
5000 SD+CG cycles actually produce."""))

cells.append(co('''cols = [c for c in score_cols(d1_omm) if c in d1_sander.columns]
M = d1_omm[["pdb"]+cols].merge(d1_sander[["pdb"]+cols], on="pdb", suffixes=("_gpu","_sander"))
print(f"complexes run both ways: {len(M)}")
out = []
for c in cols:
    x, y = M[f"{c}_sander"].astype(float), M[f"{c}_gpu"].astype(float)
    k = x.notna() & y.notna()
    if k.sum() < 4: continue
    out.append(dict(arm=c, n=int(k.sum()), mean_diff=round(float((y[k]-x[k]).mean()),2),
                    sd_diff=round(float((y[k]-x[k]).std(ddof=1)),2),
                    rank_agreement=round(spearmanr(x[k], y[k])[0],3)))
V = pd.DataFrame(out)
print(V.to_string(index=False))
print(f"\\nmedian rank agreement across arms: {V.rank_agreement.median():+.3f}")'''))

cells.append(md("""## 3. The controls the paper does not run

An entropy-free MM/PBSA score sums van der Waals, electrostatic and nonpolar surface
terms that are all extensive in contact area, so the first question is whether it beats
a descriptor of interface size. The second is that the 50 complexes are not 50
independent systems: at 90% sequence identity they are 21 receptors, and six receptors
supply 32 of them. The third is that a cross-system Pearson is not the question a lead
optimisation campaign asks, which is how well candidates rank against one receptor."""))

cells.append(co('''SIZE = ["pep_heavy","pep_res","dSASA","dSASA_apolar","contacts_4A","contacts_5A",
        "iface_pep_atoms","iface_prot_atoms","pep_len","polar_pct"]
ctrl = pd.DataFrame([dict(descriptor=f, Rp=round(pearsonr(D[f],D.pKd)[0],3),
                          rho=round(spearmanr(D[f],D.pKd)[0],3)) for f in SIZE])
ctrl = ctrl.sort_values("Rp", key=abs, ascending=False)
print("interface-size control on the same complexes:"); print(ctrl.to_string(index=False))

DESIG = "pb_2"   # the paper's own optimum, not selected here
best_size = ctrl.iloc[0].descriptor
print(f"\\ndesignated physics arm {DESIG}: Rp = {pearsonr(D[DESIG],D.pKd)[0]:+.3f}")
print(f"best size descriptor {best_size}: Rp = {ctrl.iloc[0].Rp:+.3f}")

gs = D.receptor_group.unique()
idx = {g: D.index[D.receptor_group==g].to_numpy() for g in gs}
diffs = []
for _ in range(NBOOT):
    sub = D.loc[np.concatenate([idx[g] for g in RNG.choice(gs, len(gs), replace=True)])]
    if sub[DESIG].nunique() < 3: continue
    diffs.append(abs(pearsonr(sub[DESIG],sub.pKd)[0]) - abs(pearsonr(sub[best_size],sub.pKd)[0]))
diffs = np.array([d for d in diffs if np.isfinite(d)])
print(f"|Rp| advantage of physics over size, bootstrapped over the {len(gs)} receptor groups: "
      f"{np.mean(diffs):+.3f} (95% CI {np.percentile(diffs,2.5):+.3f} to {np.percentile(diffs,97.5):+.3f})")

print("\\nwithin-receptor Spearman, receptor groups weighted equally:")
for col,label in ((DESIG,"MM/PBSA"),(best_size,f"size ({best_size})")):
    per = [(g,len(s),spearmanr(s[col],s.pKd)[0]) for g,s in D.groupby("receptor_group")
           if len(s)>=3 and s[col].nunique()>1 and s.pKd.nunique()>1]
    if per:
        print(f"  {label:<26} macro rho = {np.mean([p[2] for p in per]):+.3f} over {len(per)} groups")
        for g,n,r in per: print(f"      group {g:2d} (n={n}): {r:+.3f}")

sub = D[D.in_decoy]
print(f"\\nDataset II overlap, where the two-step workflow reports -0.732 (n={len(sub)}):")
print(f"  {DESIG}: {pearsonr(sub[DESIG],sub.pKd)[0]:+.3f}   peptide heavy atoms: "
      f"{pearsonr(sub.pep_heavy,sub.pKd)[0]:+.3f}")'''))

cells.append(md("""## 4. The same protocol on the five-family FEP+ cohort

241 macrocyclic peptides against Cyclin A/B, JAK2, KRAS, MCL-1 and PCSK9, 235 with an
experimental ΔG. Reported under both estimators: the paper's pooled cross-system
Pearson, and the within-family macro Spearman this project deploys against. Ligands are
parameterised with GAFF2 and AM1-BCC because they carry non-canonical residues; the AM1
step is a single point so the bound pose is preserved."""))

cells.append(co('''man = pd.read_csv(os.path.join(RES, "cohort_manifest.csv"))
C = cohort.merge(man, on="slug").query("experimental_dG_kcal_mol == experimental_dG_kcal_mol").copy()
print(f"scored with a label: {len(C)}   families: {C.target.nunique()}")
print(f"max pose drift after parameterisation: {pd.to_numeric(C.pose_drift_A,errors='coerce').max():.4f} A")

def macro(df, col, min_n=3):
    per = [(t,len(s),spearmanr(s[col],s.experimental_dG_kcal_mol)[0])
           for t,s in df.groupby("target") if len(s)>=min_n and s[col].nunique()>1]
    return (float(np.mean([p[2] for p in per])) if per else np.nan), per

def fam_boot(df, col, other=None, n=NBOOT):
    fams = df.target.unique(); idx = {t: df.index[df.target==t].to_numpy() for t in fams}
    out = []
    for _ in range(n):
        a, b = [], []
        for t in RNG.choice(fams, len(fams), replace=True):
            s = df.loc[RNG.choice(idx[t], len(idx[t]), replace=True)]
            if s[col].nunique() < 2: continue
            a.append(spearmanr(s[col], s.experimental_dG_kcal_mol)[0])
            if other is not None and s[other].nunique() > 1:
                b.append(spearmanr(s[other], s.experimental_dG_kcal_mol)[0])
        if not a: continue
        out.append(np.mean(a) if other is None else abs(np.mean(a)) - abs(np.mean(b or [0])))
    v = np.array([x for x in out if np.isfinite(x)])
    return np.percentile(v, [2.5, 97.5]) if len(v) > 100 else (np.nan, np.nan)

arms = score_cols(C) + ["heavy_atoms", "published_fepplus_dG_kcal_mol"]
crows = []
for c in arms:
    s = C[C[c].notna()]
    if len(s) < 20: continue
    mac, per = macro(s, c); lo, hi = fam_boot(s, c)
    crows.append(dict(arm=c, n=len(s), pooled_Rp=round(pearsonr(s[c],s.experimental_dG_kcal_mol)[0],3),
                     macro_rho=round(mac,3), lo=round(lo,3), hi=round(hi,3), families=len(per)))
T = pd.DataFrame(crows).sort_values("macro_rho", key=abs, ascending=False)
print(); print(T.to_string(index=False))'''))

cells.append(co('''DES = "pb_2"
s = C[C[DES].notna()]
mac, per = macro(s, DES); mh, perh = macro(s, "heavy_atoms")
lo, hi = fam_boot(s, DES, other="heavy_atoms")
print(f"designated arm {DES}: macro rho {mac:+.3f}; heavy atoms {mh:+.3f}")
print(f"|advantage| over size: {abs(mac)-abs(mh):+.3f} (95% CI {lo:+.3f} to {hi:+.3f})")
print("\\nper family:")
hh = {p[0]: p[2] for p in perh}
for t,n,r in per: print(f"  {t:<12} n={n:3d}   MM/PBSA {r:+.3f}   heavy atoms {hh.get(t,np.nan):+.3f}")

if "censored" in C.columns:
    nc = C[C.censored != True]
    print(f"\\nwith the {int((C.censored==True).sum())} censored labels dropped ({len(nc)} compounds):")
    for c in [DES, "heavy_atoms", "published_fepplus_dG_kcal_mol"]:
        sub = nc[nc[c].notna()]
        print(f"  {c:<32} macro rho {macro(sub,c)[0]:+.3f}")

print("\\nenergy decomposition of the designated arm:")
for t in [c for c in C.columns if c.startswith(DES+"_")]:
    sub = C[C[t].notna()]
    if len(sub) > 20 and sub[t].nunique() > 3:
        print(f"  {t:<24} macro rho {macro(sub,t)[0]:+.3f}")'''))

cells.append(co('''fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

g = grid[["gb1","gb2","gb5","pb"]]
for m, mk in zip(["pb","gb1","gb2","gb5"], ["o","s","^","v"]):
    ax[0].plot(g.index, g[m], mk+"-", label=f"{m} (this work)")
    ax[0].plot(pub.index, pub[m], mk+"--", alpha=.45, label=f"{m} (published)")
ax[0].axhline(0, color="k", lw=.6)
ax[0].set_xlabel("interior dielectric"); ax[0].set_ylabel("Rp vs pKd")
ax[0].set_title(f"Dataset I, ff03, implicit (n={len(D)})"); ax[0].legend(fontsize=6, ncol=2)

names = [r["arm"] for r in crows]; vals = [r["macro_rho"] for r in crows]
o = np.argsort([-abs(v) for v in vals])
cols_ = ["tab:green" if names[i]=="published_fepplus_dG_kcal_mol" else
         "tab:grey" if names[i]=="heavy_atoms" else
         "tab:red" if names[i]==DES else "tab:blue" for i in o]
ax[1].barh(range(len(o)), [abs(vals[i]) for i in o], color=cols_)
ax[1].set_yticks(range(len(o))); ax[1].set_yticklabels([names[i] for i in o], fontsize=7)
ax[1].invert_yaxis(); ax[1].set_xlabel("|within-family macro Spearman|")
ax[1].set_title(f"five-family cohort (n={len(C)})")

fams = [p[0] for p in per]
x = np.arange(len(fams)); w = .38
ax[2].bar(x-w/2, [abs(p[2]) for p in per], w, label="MM/PBSA", color="tab:red")
ax[2].bar(x+w/2, [abs(hh.get(f, np.nan)) for f in fams], w, label="heavy atoms", color="tab:grey")
ax[2].set_xticks(x); ax[2].set_xticklabels(fams, rotation=30, ha="right", fontsize=7)
ax[2].set_ylabel("|Spearman|"); ax[2].set_title("per receptor family"); ax[2].legend(fontsize=7)

plt.tight_layout()
plt.savefig(os.path.join(RES, "zha25_replication_summary.png"), dpi=150)
print("wrote zha25_replication_summary.png")'''))

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"name":"python3","display_name":"Python 3","language":"python"}
nbf.write(nb, sys.argv[1] if len(sys.argv)>1 else "zha25_replication.ipynb")
print("wrote notebook")
