#!/usr/bin/env python3
"""Place each 12ca5 peptidomimetic in the antibody groove, motif-blind.

Six of the seven binders in this panel carry the D-Y-A epitope anchor and only one of
the ten non-binders does, so threading onto the HA epitope would hand the binders a
good pose by construction and the comparison would be circular.  Every ligand is
instead embedded onto the *peptide backbone* of peptidomimetic 16 -- the one compound
in the panel with an experimentally modelled complex.  The backbone is common to all
seventeen, binder and non-binder alike, and carries no information about the epitope.

The template backbone is read from PDB atom names, which are unambiguous; the target
backbone is found by SMARTS and ordered by walking the peptide bonds.
"""
import argparse, csv, os, sys
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from scipy.spatial import cKDTree
RDLogger.DisableLog("rdApp.*")

# any amide nitrogen, including the N-substituted residues (Pro, Tic, Pip, Thp, Hyp)
BB = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])")


def split_pdb(pdb, lig_chain="C"):
    rec, lig = [], []
    for l in open(pdb, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            (lig if l[21] == lig_chain else rec).append(l)
    return rec, lig


def template_backbone(lig_lines):
    """N, CA, C of each ligand residue, in residue order, from PDB atom names."""
    res = []
    for l in lig_lines:
        key = (l[21], l[22:27])
        nm = l[12:16].strip()
        if nm in ("N", "CA", "C"):
            if not res or res[-1][0] != key:
                res.append((key, {}))
            res[-1][1][nm] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    out = []
    for key, d in res:
        if all(k in d for k in ("N", "CA", "C")):
            out.append((d["N"], d["CA"], d["C"]))
    return out


def target_backbone(mol):
    """Backbone (N, CA, C) triples ordered along the chain."""
    hits = mol.GetSubstructMatches(BB)
    if not hits:
        return []
    units = {h[1]: (h[0], h[1], h[2]) for h in hits}          # keyed on CA
    cto = {}                                                   # carbonyl C -> next N
    for ca, (n, _, c) in units.items():
        for nb in mol.GetAtomWithIdx(c).GetNeighbors():
            if nb.GetSymbol() == "N" and nb.GetIdx() != n:
                cto[c] = nb.GetIdx()
    nxt = {}
    ncmap = {v[0]: k for k, v in units.items()}                # N -> CA
    for ca, (_, _, c) in units.items():
        if c in cto and cto[c] in ncmap:
            nxt[ca] = ncmap[cto[c]]
    # take the longest chain reachable, not the first: a missed peptide bond at an
    # N-substituted residue otherwise truncates the backbone and leaves different
    # compounds tethered over different numbers of residues
    best = []
    for start in units:
        chain, seen, ca = [], set(), start
        while ca is not None and ca not in seen:
            seen.add(ca); chain.append(units[ca]); ca = nxt.get(ca)
        if len(chain) > len(best):
            best = chain
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--graphs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nconf", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rec_lines, lig_lines = split_pdb(a.template)
    open(os.path.join(a.out, "receptor.pdb"), "w").writelines(rec_lines + ["END\n"])
    tb = template_backbone(lig_lines)
    rec_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in rec_lines])
    print(f"template: {len(tb)} backbone residues, receptor {len(rec_lines)} atoms", flush=True)

    report = []
    for r in csv.DictReader(open(a.graphs), delimiter="\t"):
        sid, num = r["system_id"], r["peptide_number"]
        m = Chem.MolFromSmiles(r["ligand_smiles"])
        if m is None:
            report.append(dict(system_id=sid, status="smiles_failed")); print(sid, "SMILES failed", flush=True); continue
        m = Chem.AddHs(m)
        qb = target_backbone(m)
        k = min(len(qb), len(tb))
        if k < 4:
            report.append(dict(system_id=sid, status="backbone_too_short", n=k))
            print(sid, "backbone too short", k, flush=True); continue
        # tether the ligand's first k backbone residues to the template's first k
        cmap = {}
        for i in range(k):
            for j in range(3):
                cmap[qb[i][j]] = Chem.rdGeometry.Point3D(*tb[i][j])
        tgt = np.array([[cmap[i].x, cmap[i].y, cmap[i].z] for i in cmap])
        order = list(cmap)
        best = None
        for seed in range(a.nconf):
            mm = Chem.Mol(m)
            # embed freely -- constrained distance geometry does not converge on a
            # 250-atom peptidomimetic -- then move it onto the template backbone as a
            # rigid body and only then relax against the tethers
            if AllChem.EmbedMolecule(mm, randomSeed=0xf00d + seed, useRandomCoords=True,
                                     maxAttempts=20) < 0:
                continue
            AllChem.UFFOptimizeMolecule(mm, maxIters=200)
            conf = mm.GetConformer()
            P = np.array([list(conf.GetAtomPosition(i)) for i in order])
            pc, tc = P.mean(0), tgt.mean(0)
            U, S, Vt = np.linalg.svd((P - pc).T @ (tgt - tc))
            d = np.sign(np.linalg.det(Vt.T @ U.T))
            R = Vt.T @ np.diag([1, 1, d]) @ U.T
            allp = np.array([list(conf.GetAtomPosition(i)) for i in range(mm.GetNumAtoms())])
            allp = (allp - pc) @ R.T + tc
            for i, p in enumerate(allp):
                conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(*p))
            ff = AllChem.UFFGetMoleculeForceField(mm)
            for idx, pt in cmap.items():
                p = ff.AddExtraPoint(pt.x, pt.y, pt.z, fixed=True) - 1
                ff.AddDistanceConstraint(p, idx, 0, 0.5, 100.0)
            ff.Initialize()          # extra points resize the position vector
            ff.Minimize(maxIts=600)
            conf = mm.GetConformer()
            xyz = np.array([list(conf.GetAtomPosition(i)) for i in range(mm.GetNumAtoms())])
            heavy = np.array([at.GetAtomicNum() > 1 for at in mm.GetAtoms()])
            nclash = int((cKDTree(rec_xyz).query_ball_point(xyz[heavy], 1.8, return_length=True) > 0).sum())
            bbrms = float(np.sqrt(np.mean([np.sum((xyz[i] - np.array([cmap[i].x, cmap[i].y, cmap[i].z])) ** 2)
                                           for i in cmap])))
            if best is None or nclash < best[0]:
                best = (nclash, bbrms, mm)
        if best is None:
            report.append(dict(system_id=sid, status="embed_failed")); print(sid, "embed failed", flush=True); continue
        nclash, bbrms, mm = best
        path = os.path.join(a.out, f"{sid.replace(':', '_')}.sdf")
        w = Chem.SDWriter(path); w.write(mm); w.close()
        hv = sum(1 for at in mm.GetAtoms() if at.GetAtomicNum() > 1)
        report.append(dict(system_id=sid, peptide_number=num, status="ok", backbone_residues_tethered=k,
                           backbone_rmsd_A=round(bbrms, 2), heavy_atoms=hv,
                           formal_charge=Chem.GetFormalCharge(mm), clashes_1p8A=nclash, sdf=path))
        print(f"{sid} pep{num}: {k} backbone residues tethered, bb RMSD {bbrms:.2f} A, "
              f"{hv} heavy, q {Chem.GetFormalCharge(mm):+d}, clashes {nclash}", flush=True)
    with open(os.path.join(a.out, "pose_build_report.csv"), "w", newline="") as fh:
        keys = sorted({k for d in report for k in d})
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(report)
    print(f"\nbuilt {sum(1 for d in report if d.get('status') == 'ok')}/{len(report)}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
