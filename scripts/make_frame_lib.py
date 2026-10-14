#!/usr/bin/env python3
"""Manufacture curated-library-shaped inputs for the thanatin derivative-7 frame residues.

The D7 peptide frame is itself non-canonical: EU0, HYP, LE1, DAB and 4FO sit alongside the
swapped ncAA, and none of them are in the NCFLOW curated library.  EU0, LE1 and 4FO have a
Rosetta .params from the campaign; HYP and DAB have nothing, so their connectivity is derived
from the geometry -- both are fully saturated (4-hydroxyproline and 2,4-diaminobutyric acid),
so the only multiple bond is the backbone carbonyl, which is assigned by giving order 2 to a
C-O bond whose oxygen has a single heavy neighbour.

Writes <out>/<CODE>/{<CODE>.params,<CODE>_rotamer.pdb} so ncaa_param.py runs on them unchanged.
"""
import argparse, math, os, sys

COV = {"H":0.31,"C":0.76,"N":0.71,"O":0.66,"S":1.05,"P":1.07,"F":0.57,
       "CL":1.02,"BR":1.20,"I":1.39,"SE":1.20}

def residue_atoms(pdb, resname):
    out=[]; key=None
    for l in open(pdb, errors="ignore"):
        if not l.startswith(("ATOM","HETATM")): continue
        if l[17:20].strip()!=resname: continue
        k=(l[21], l[22:27])
        if key is None: key=k
        if k!=key: break                      # first instance only
        el=(l[76:78].strip() or l[12:16].strip()[0]).strip().upper()
        out.append((l[12:16].strip(), el,
                    (float(l[30:38]), float(l[38:46]), float(l[46:54])), l))
    return out

def infer_bonds(atoms):
    n=len(atoms); bonds=[]; nb={i:[] for i in range(n)}
    for i in range(n):
        for j in range(i+1,n):
            ri=COV.get(atoms[i][1],0.75); rj=COV.get(atoms[j][1],0.75)
            d=math.dist(atoms[i][2],atoms[j][2])
            if d < ri+rj+0.45:
                if atoms[i][1]=="H" and atoms[j][1]=="H": continue
                bonds.append([i,j,"1"]); nb[i].append(j); nb[j].append(i)
    heavy_nb={i:[j for j in nb[i] if atoms[j][1]!="H"] for i in range(n)}
    for b in bonds:
        i,j=b[0],b[1]
        for a,o in ((i,j),(j,i)):
            # a carbonyl oxygen has exactly one neighbour of any kind; a hydroxyl
            # oxygen also has one heavy neighbour but carries a hydrogen as well
            if atoms[a][1]=="O" and len(nb[a])==1 and atoms[o][1]=="C":
                b[2]="2"
    return bonds

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("resname"); ap.add_argument("pdb"); ap.add_argument("outdir")
    ap.add_argument("--params", default=None, help="Rosetta .params to take names and bonds from")
    a=ap.parse_args()
    atoms=residue_atoms(a.pdb, a.resname)
    if not atoms: sys.exit(f"{a.resname} not found in {a.pdb}")
    d=os.path.join(a.outdir, a.resname); os.makedirs(d, exist_ok=True)
    def nn(s): return s[1:]+s[0] if s and s[0].isdigit() else s
    if a.params:
        praw=open(a.params, errors="ignore").read().splitlines()
        pnames=[l.split()[1] for l in praw if l.startswith("ATOM ")]
        present={nn(x[0]) for x in atoms}
        missing=[x for x in pnames if nn(x) not in present]
        if missing:
            # LE1's thiol hydrogen HG is absent from the structure because SG is in the
            # disulfide -- the CYX case.  Drop it and every bond that referenced it.
            print(f"  {a.resname}: dropping {missing} (absent from the structure)")
        keep={nn(x) for x in pnames} - {nn(x) for x in missing}
        # the .params and the structure list the same atoms in different orders (4FO puts the
        # HZ hydrogens before HA); coordinates are matched positionally downstream, so emit the
        # ATOM records in the structure's order
        byname={nn(l.split()[1]): l for l in praw if l.startswith("ATOM ")}
        order=[nn(x[0]) for x in atoms]
        if set(order)!=set(byname)-{nn(x) for x in missing}:
            sys.exit(f"{a.resname}: atom sets differ: only in params "
                     f"{sorted(set(byname)-set(order)-{nn(x) for x in missing})}, "
                     f"only in structure {sorted(set(order)-set(byname))}")
        out=[l for l in praw if not l.startswith(("ATOM ","BOND"))]
        head=[l for l in out if not l.startswith(("ICOOR","NBR","CHI","PROPERTIES"))]
        body=[byname[n] for n in order]
        bonds=[l for l in praw if l.startswith("BOND")
               and nn(l.split()[1]) in keep and nn(l.split()[2]) in keep]
        src="\n".join(head+body+bonds)+"\n"
        if len(body)!=len(atoms):
            sys.exit(f"{a.resname}: {len(body)} params atoms vs {len(atoms)} structure atoms")
    else:
        bonds=infer_bonds(atoms)
        lines=[f"NAME {a.resname}", f"IO_STRING {a.resname} X", "TYPE POLYMER", "AA UNK"]
        for nm,el,_,_ in atoms: lines.append(f"ATOM  {nm:<4s} X    X    0.000")
        for i,j,o in bonds: lines.append(f"BOND_TYPE  {atoms[i][0]:<4s} {atoms[j][0]:<4s} {o}")
        src="\n".join(lines)+"\n"
        print(f"  {a.resname}: {len(atoms)} atoms, {len(bonds)} bonds inferred from geometry "
              f"({sum(1 for b in bonds if b[2]=='2')} double)")
    open(os.path.join(d, f"{a.resname}.params"),"w").write(src)
    with open(os.path.join(d, f"{a.resname}_rotamer.pdb"),"w") as f:
        for _,_,_,raw in atoms: f.write(raw)
        f.write("END\n")
    print(f"  wrote {d}")

if __name__=="__main__": main()
