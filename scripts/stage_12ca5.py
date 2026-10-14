#!/usr/bin/env python3
"""Stage 12ca5 complexes for the Zha25 pipeline.

Two arms.  `modelled` is the one compound in the panel with an experimentally
modelled complex, peptidomimetic 16, plus the native HA peptide, taken as deposited.
`built` is all seventeen placed by the backbone-tethered builder.  The second arm is
exploratory and its own diagnostic is reported alongside it: how far the builder puts
peptidomimetic 16 from the complex it was templated on.
"""
import argparse, csv, glob, os, sys
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")


def write_complex(rec_lines, mol, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    conf = mol.GetConformer()
    lines = list(rec_lines)
    for i, at in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        el = at.GetSymbol()
        lines.append("HETATM%5d %-4s LIG X   1    %8.3f%8.3f%8.3f  1.00 20.00          %2s\n"
                     % (i + 1, (el + str(i + 1))[:4], p.x, p.y, p.z, el.rjust(2)))
    open(os.path.join(out_dir, "native_complex.pdb"), "w").writelines(lines + ["END\n"])
    w = Chem.SDWriter(os.path.join(out_dir, "ligand_v2000.sdf")); w.write(mol); w.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poses", required=True)
    ap.add_argument("--graphs", required=True)
    ap.add_argument("--modelled", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rec_lines = [l for l in open(os.path.join(a.poses, "receptor.pdb"), errors="ignore")
                 if l.startswith(("ATOM", "HETATM"))]
    graphs = {r["system_id"]: r for r in csv.DictReader(open(a.graphs), delimiter="\t")}

    # arm 1: the modelled peptidomimetic 16 complex, bond orders from its own SMILES
    lig_pdb = [l for l in open(a.modelled, errors="ignore")
               if l.startswith(("ATOM", "HETATM")) and l[21] == "C"]
    ref = Chem.MolFromPDBBlock("".join(lig_pdb) + "END\n", removeHs=False,
                               sanitize=False, proximityBonding=True)
    smi = Chem.MolFromSmiles(graphs["12ca5:16"]["ligand_smiles"])
    ok16 = False
    if ref is not None and smi is not None:
        try:
            ref = AllChem.AssignBondOrdersFromTemplate(smi, ref)
            ok16 = True
        except Exception as e:
            print("bond-order assignment on the modelled ligand failed:", e)
    if ok16:
        write_complex(rec_lines, ref, os.path.join(a.out, "modelled", "12ca5_16_modelled"))
        print(f"modelled arm: peptidomimetic 16, {ref.GetNumAtoms()} ligand atoms")

    # arm 2: every built pose, plus the same diagnostic for 16
    rows = []
    for f in sorted(glob.glob(os.path.join(a.poses, "*.sdf"))):
        sid = os.path.basename(f)[:-4]
        m = Chem.SDMolSupplier(f, removeHs=False)[0]
        if m is None:
            print("unreadable", sid); continue
        write_complex(rec_lines, m, os.path.join(a.out, "built", sid))
        rows.append(sid)
    print(f"built arm: {len(rows)} complexes staged")

    # how far does the builder put 16 from the complex it was templated on?
    if ok16:
        b = Chem.SDMolSupplier(os.path.join(a.poses, "12ca5_16.sdf"), removeHs=False)[0]
        if b is not None:
            mt = b.GetSubstructMatch(Chem.RemoveHs(ref)) or None
            rc = ref.GetConformer(); bc = b.GetConformer()
            heavy_ref = [i for i, at in enumerate(ref.GetAtoms()) if at.GetAtomicNum() > 1]
            if mt and len(mt) == len(heavy_ref):
                d = [np.linalg.norm(np.array(list(bc.GetAtomPosition(mt[k])))
                                    - np.array(list(rc.GetAtomPosition(i))))
                     for k, i in enumerate(heavy_ref)]
                print(f"DIAGNOSTIC: the builder places peptidomimetic 16 at "
                      f"{np.sqrt(np.mean(np.square(d))):.2f} A heavy-atom RMSD from its own "
                      f"modelled complex")
            else:
                print("DIAGNOSTIC: could not map the built pose onto the modelled ligand "
                      "one-to-one; backbone RMSD from the build report stands")


if __name__ == "__main__":
    sys.exit(main())
