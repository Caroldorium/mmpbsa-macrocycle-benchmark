#!/usr/bin/env python3
"""Build one Amber dual-topology TI edge (ligand A -> ligand B) for a FEP+ cohort family.

Reuses the GAFF2/AM1-BCC ligand parameters the MM/PBSA cohort run already produced
(lig.mol2 + lig.frcmod per compound, ligand in its bound pose).  Receptor coordinates are
identical for every compound in a receptor group (verified: CA RMSD 0.000 A), so the two
ligands are already in one frame and no alignment is needed.

Amber's softcore TI requires the common-core atoms to be the FIRST atoms of each TI region
and to correspond in order between the two regions.  Independently built mol2 files do not
satisfy that, so both ligands are renumbered here from the MCS mapping.
"""
import argparse, json, os, shutil, sys
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS
RDLogger.DisableLog("rdApp.*")


def read_mol2(path):
    sec = None
    head, atoms, bonds, tail = [], [], [], []
    for l in open(path, errors="ignore"):
        if l.startswith("@<TRIPOS>"):
            sec = l.strip()
            continue
        if sec == "@<TRIPOS>ATOM":
            atoms.append(l.rstrip("\n"))
        elif sec == "@<TRIPOS>BOND":
            bonds.append(l.rstrip("\n"))
        elif atoms:
            tail.append(l)
        else:
            head.append(l)
    return head, [a for a in atoms if a.strip()], [b for b in bonds if b.strip()], tail


def write_mol2(path, name, atoms, bonds, resname):
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{name}\n{len(atoms):5d}{len(bonds):6d}     1     0     0\n"
                "SMALL\nbcc\n\n@<TRIPOS>ATOM\n")
        for i, a in enumerate(atoms):
            t = a.split()
            f.write(f"{i+1:7d} {t[1]:<8s}{float(t[2]):10.4f}{float(t[3]):10.4f}{float(t[4]):10.4f} "
                    f"{t[5]:<8s}{1:3d} {resname:<8s}{float(t[8]):10.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for j, b in enumerate(bonds):
            t = b.split()
            f.write(f"{j+1:6d}{int(t[1]):6d}{int(t[2]):6d} {t[3]}\n")
        # tleap's mol2 reader expects the substructure record and errors at EOF without it
        f.write("@<TRIPOS>SUBSTRUCTURE\n")
        f.write(f"     1 {resname:<8s}    1 TEMP              0 ****  ****    0 ROOT\n")


def reorder(atoms, bonds, order):
    """order: old 1-based indices, listed in their new position."""
    newpos = {old: i + 1 for i, old in enumerate(order)}
    a = [atoms[o - 1] for o in order]
    b = []
    for l in bonds:
        t = l.split()
        i, j = int(t[1]), int(t[2])
        if i in newpos and j in newpos:
            b.append(f"  1 {newpos[i]} {newpos[j]} {t[3]}")
    return a, b


def is_h(line):
    return line.split()[5].lower().startswith("h")


def attached_h(atoms, bonds, idx):
    out = []
    for l in bonds:
        t = l.split()
        i, j = int(t[1]), int(t[2])
        for u, v in ((i, j), (j, i)):
            if u == idx and is_h(atoms[v - 1]):
                out.append(v)
    return sorted(set(out))


def build_order(atoms, bonds, hv, match):
    """Common block first: each core heavy atom immediately followed by its hydrogens."""
    seen, order = set(), []
    for k in match:
        c = hv[k]
        if c in seen:
            continue
        order.append(c)
        seen.add(c)
        for h in attached_h(atoms, bonds, c):
            if h not in seen:
                order.append(h)
                seen.add(h)
    rest = [i for i in range(1, len(atoms) + 1) if i not in seen]
    return order, len(order), rest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slug_a")
    ap.add_argument("slug_b")
    ap.add_argument("outdir")
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--prep", required=True)
    ap.add_argument("--mcs-timeout", type=int, default=150)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    def sdf(slug):
        f = f"{a.prep}/{slug}/ligand_v2000.sdf"
        return Chem.MolFromMolFile(f, sanitize=True) or Chem.MolFromMolFile(f, sanitize=False)

    A, B = sdf(a.slug_a), sdf(a.slug_b)
    if A is None or B is None:
        print(json.dumps(dict(status="fail", reason="sdf unreadable")))
        sys.exit(2)
    res = rdFMCS.FindMCS([A, B], timeout=a.mcs_timeout, ringMatchesRingOnly=True,
                         completeRingsOnly=False, matchValences=False,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareAny)
    patt = Chem.MolFromSmarts(res.smartsString) if res.smartsString else None
    mA = A.GetSubstructMatch(patt) if patt else ()
    mB = B.GetSubstructMatch(patt) if patt else ()
    if not mA or not mB or len(mA) != len(mB):
        print(json.dumps(dict(status="fail", reason="no usable MCS match")))
        sys.exit(3)

    hA, atA, bdA, _ = read_mol2(f"{a.cohort}/{a.slug_a}/lig.mol2")
    hB, atB, bdB, _ = read_mol2(f"{a.cohort}/{a.slug_b}/lig.mol2")
    hvA = [i + 1 for i, l in enumerate(atA) if not is_h(l)]
    hvB = [i + 1 for i, l in enumerate(atB) if not is_h(l)]
    if len(hvA) != A.GetNumHeavyAtoms() or len(hvB) != B.GetNumHeavyAtoms():
        print(json.dumps(dict(status="fail", reason="heavy-atom count mismatch between sdf and mol2",
                              sdf_a=A.GetNumHeavyAtoms(), mol2_a=len(hvA),
                              sdf_b=B.GetNumHeavyAtoms(), mol2_b=len(hvB))))
        sys.exit(4)

    ordA, ncA, restA = build_order(atA, bdA, hvA, mA)
    ordB, ncB, restB = build_order(atB, bdB, hvB, mB)
    if ncA != ncB:
        # drop core atoms whose hydrogen counts differ; otherwise the leading blocks of the
        # two TI regions do not correspond atom for atom, which Amber requires
        keep = [k for k in range(len(mA))
                if len(attached_h(atA, bdA, hvA[mA[k]])) == len(attached_h(atB, bdB, hvB[mB[k]]))]
        mA = [mA[k] for k in keep]
        mB = [mB[k] for k in keep]
        ordA, ncA, restA = build_order(atA, bdA, hvA, mA)
        ordB, ncB, restB = build_order(atB, bdB, hvB, mB)
    if ncA != ncB or ncA == 0:
        print(json.dumps(dict(status="fail", reason="common blocks differ", nc_a=ncA, nc_b=ncB)))
        sys.exit(5)

    nA, nbA = reorder(atA, bdA, ordA + restA)
    nB, nbB = reorder(atB, bdB, ordB + restB)
    write_mol2(f"{a.outdir}/lga.mol2", "LGA", nA, nbA, "LGA")
    write_mol2(f"{a.outdir}/lgb.mol2", "LGB", nB, nbB, "LGB")
    shutil.copy(f"{a.cohort}/{a.slug_a}/lig.frcmod", f"{a.outdir}/lga.frcmod")
    shutil.copy(f"{a.cohort}/{a.slug_b}/lig.frcmod", f"{a.outdir}/lgb.frcmod")
    shutil.copy(f"{a.cohort}/{a.slug_a}/rec_clean.pdb", f"{a.outdir}/rec.pdb")

    scA = [nA[i].split()[1] for i in range(ncA, len(nA))]
    scB = [nB[i].split()[1] for i in range(ncB, len(nB))]
    info = dict(status="ok", slug_a=a.slug_a, slug_b=a.slug_b,
                n_common=ncA, n_atoms_a=len(nA), n_atoms_b=len(nB),
                n_soft_a=len(scA), n_soft_b=len(scB),
                timask1=":LGA", timask2=":LGB",
                scmask1=(":LGA@" + ",".join(scA)) if scA else "",
                scmask2=(":LGB@" + ",".join(scB)) if scB else "",
                mcs_converged=bool(not res.canceled), mcs_heavy_atoms=res.numAtoms)
    json.dump(info, open(f"{a.outdir}/edge.json", "w"), indent=1)
    print(json.dumps(info))


if __name__ == "__main__":
    main()
