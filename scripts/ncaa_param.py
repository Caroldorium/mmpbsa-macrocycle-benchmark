#!/usr/bin/env python3
"""Amber parameterisation of non-canonical alpha amino acids for the Zha25 MM/PBSA pipeline.

Input   <code>_rotamer.pdb from the curated NCFLOW library: a Rosetta polymer residue with
        proper backbone atom names (N, CA, C, O, H, HA) and one or more stacked conformers,
        plus <code>.params for connectivity.
Output  <code>.prep   Amber residue, HEAD_NAME N / TAIL_NAME C, ready for `loadAmberPrep`
        <code>.frcmod missing GAFF2 parameters

Route   first conformer -> ACE/NME capping -> antechamber AM1-BCC on the neutral capped
        dipeptide -> prepgen (omits the caps, redistributes their charge) -> parmchk2.

The caps are the point.  AM1-BCC on a free amino acid (NH2 / COOH) fits backbone charges in
an environment the residue is never actually in; capping fits them in a peptide-like one.
Cap atoms are named CY/OY/CAY/HY* and NZ/HZ/CAZ/HZ* so they cannot collide with the
residue's own N, C and O.
"""
import argparse, json, math, os, re, subprocess, sys

GAFF2_ELEMENTS = {"H","C","N","O","S","P","F","CL","BR","I"}

CAP_ATOMS = ["CY","OY","CAY","HY1","HY2","HY3","NZ","HZ","CAZ","HZ1","HZ2","HZ3"]

def v(a,b): return [b[i]-a[i] for i in range(3)]
def add(a,b): return [a[i]+b[i] for i in range(3)]
def scale(a,s): return [x*s for x in a]
def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def norm(a):
    n=math.sqrt(dot(a,a)); return [x/n for x in a] if n>1e-9 else [0.,0.,0.]

def place(a, b, c, dist, angle_deg, dih_deg):
    """Position an atom at `dist` from c, angle(new,c,b), dihedral(new,c,b,a) -- NeRF."""
    ang=math.radians(angle_deg); dih=math.radians(dih_deg)
    bc=norm(v(b,c)); ab=v(a,b)
    n=norm(cross(ab,bc)); m=cross(n,bc)
    d=[-dist*math.cos(ang), dist*math.sin(ang)*math.cos(dih), dist*math.sin(ang)*math.sin(dih)]
    return [c[i] + d[0]*bc[i] + d[1]*m[i] + d[2]*n[i] for i in range(3)]

def norm_name(nm):
    """Rosetta writes PDB v2 hydrogen names with the digit first (1HG2); prepgen cannot parse
    those ("Unrecognized case-sensitive atomic symbol"). Move the digit to the end, v3 style.
    Only the template's own hydrogen names are affected -- the pipeline's pdb4amber pass
    strips hydrogens and tleap rebuilds them from this template, so nothing downstream
    depends on the old spelling."""
    return nm[1:]+nm[0] if nm and nm[0].isdigit() else nm

def first_conformer(pdb, params):
    """Atom names come from the .params ATOM records, coordinates from the first conformer of
    the rotamer PDB, matched positionally (verified to correspond for all 165 residues).

    The two files disagree on names for 9 halogenated residues: the PDB truncates CLM1 to CLM
    and, in Q2S, gives BOTH chlorines the same name CLM.  Taking names from the PDB leaves
    those atoms unbonded -- the .params BOND records name them CLM1/CLM2 -- and prepgen then
    silently discards them, emitting an atom with a nan internal coordinate.  The .params
    names are unique and are what the BOND records use, so they are the authority."""
    names=[norm_name(l[5:9].strip()) for l in open(params, errors="ignore") if l.startswith("ATOM ")]
    coords=[]
    for l in open(pdb, errors="ignore"):
        if not l.startswith(("ATOM","HETATM")): continue
        if len(coords)>=len(names): break
        el=(l[76:78].strip() or l[12:16].strip()[0]).strip().capitalize()
        coords.append((el,[float(l[30:38]),float(l[38:46]),float(l[46:54])]))
    if len(coords)!=len(names):
        raise RuntimeError(f"{len(names)} atoms in params but {len(coords)} coordinates in the rotamer PDB")
    if len(set(names))!=len(names):
        raise RuntimeError("duplicate atom names in the params file")
    return [{"name":n,"el":e,"xyz":x} for n,(e,x) in zip(names,coords)]


ORDER={"1":"1","2":"2","3":"3","ARO":"ar"}

def params_bonds(pf):
    """(a, b, order) with the order the library recorded -- antechamber's own bond-order
    perception fails on hypervalent sulfur (sulfonamides) and on iodoaryls, so we supply it."""
    bonds=[]
    for l in open(pf, errors="ignore"):
        if l.startswith("BOND"):
            p=l.split()
            if len(p)>=3:
                o=ORDER.get(p[3].upper(),"1") if len(p)>=4 else "1"
                bonds.append((norm_name(p[1]),norm_name(p[2]),o))
    return bonds

CAP_SPEC=[("CY","C"),("OY","O"),("CAY","C"),("HY1","H"),("HY2","H"),("HY3","H"),
          ("NZ","N"),("HZ","H"),("CAZ","C"),("HZ1","H"),("HZ2","H"),("HZ3","H")]

def unique_cap_names(existing):
    """ACE/NME cap names clash with side-chain atoms all through this library -- a Lys-like
    residue already owns NZ/HZ1, an aromatic one owns HZ.  Pick a free name per residue."""
    used=set(existing); out={}
    for pref,el in CAP_SPEC:
        nm=pref; k=0
        while nm in used:
            k+=1; nm=f"{el}{89+k}"
        used.add(nm); out[pref]=nm
    return out

def bisector_site(centre, neighbours, dist):
    """Place a substituent on an sp2 centre along the direction that bisects away from the
    bonds it already has.  Using a fixed dihedral instead puts the NME nitrogen on top of the
    residue's own carbonyl oxygen whenever psi happens to sit near 180 degrees."""
    d=[0.,0.,0.]
    for nb in neighbours:
        u=norm(v(centre,nb)); d=[d[i]-u[i] for i in range(3)]
    d=norm(d)
    if dot(d,d)<0.5:                      # neighbours were collinear: pick any perpendicular
        ref=norm(v(centre,neighbours[0]))
        alt=[1.,0.,0.] if abs(ref[0])<0.9 else [0.,1.,0.]
        d=norm(cross(ref,alt))
    return [centre[i]+dist*d[i] for i in range(3)]

def build_capped(atoms, bonds=None, cap_h=None):
    pos={a["name"]:a["xyz"] for a in atoms}
    for need in ("N","CA","C","O"):
        if need not in pos: raise RuntimeError(f"missing backbone atom {need}")
    N,CA,C,O = pos["N"],pos["CA"],pos["C"],pos["O"]
    cn=unique_cap_names(pos.keys())
    new=[]
    # ACE on the N terminus.  N already carries CA and (usually) H; an N-alkylated residue
    # has no amide H, so fall back to CA alone.
    # N's real neighbours, from the bond records: a proline-like N has CD instead of an amide
    # hydrogen, and taking only CA would place the cap collinear with CA-N, which makes every
    # downstream internal coordinate degenerate (two methyl hydrogens land on the same point).
    nbr={}
    for x,y,_ in (bonds or []):
        nbr.setdefault(x,[]).append(y); nbr.setdefault(y,[]).append(x)
    n_nb=[pos[b] for b in nbr.get("N",[]) if b in pos] or [CA]
    CY  = bisector_site(N, n_nb, 1.335)
    OY  = place(CA, N, CY, 1.229, 122.9,   0.0)
    CAY = place(CA, N, CY, 1.508, 116.6, 180.0)
    new += [(cn["CY"],"C",CY),(cn["OY"],"O",OY),(cn["CAY"],"C",CAY)]
    for i,d in enumerate((60.,180.,300.)):
        new.append((cn[f"HY{i+1}"],"H", place(N, CY, CAY, 1.090, 109.5, d)))
    # NME on the C terminus: C carries CA and O
    c_nb=[pos[b] for b in nbr.get("C",[]) if b in pos] or [CA,O]
    NZ  = bisector_site(C, c_nb, 1.335)
    HZ  = place(CA, C, NZ, 1.010, 119.8,   0.0)
    CAZ = place(CA, C, NZ, 1.448, 121.9, 180.0)
    new += [(cn["NZ"],"N",NZ),(cn["HZ"],"H",HZ),(cn["CAZ"],"C",CAZ)]
    for i,d in enumerate((60.,180.,300.)):
        new.append((cn[f"HZ{i+1}"],"H", place(C, NZ, CAZ, 1.090, 109.5, d)))
    extra=[]
    if cap_h:
        # LE1's SG carries no hydrogen because it sits in a disulfide -- correct for the
        # residue, but an open valence that AM1 cannot treat.  Add the hydrogen for the charge
        # fit only and omit it in prepgen, exactly as the ACE/NME caps are handled.
        if cap_h not in pos: raise RuntimeError(f"--cap-h atom {cap_h} not in the residue")
        nbrs=[pos[b] for b in {x for bpair in [] for x in bpair}] or []
        host=pos[cap_h]
        near=sorted(((math.dist(host,a["xyz"]), a) for a in atoms if a["name"]!=cap_h),
                    key=lambda t:t[0])[:1]
        hname=f"H{cap_h}"[:4]
        k=0
        while hname in pos:
            k+=1; hname=f"HX{k}"
        extra=[{"name":hname,"el":"H","xyz":bisector_site(host,[near[0][1]["xyz"]],1.34)}]
        cn[hname]=hname
    out = atoms + [{"name":n,"el":e,"xyz":x} for n,e,x in new] + extra
    worst=(None,None,9e9)
    for i in range(len(out)):
        for j in range(i+1,len(out)):
            dd=math.dist(out[i]["xyz"], out[j]["xyz"])
            if dd<worst[2]: worst=(out[i]["name"],out[j]["name"],dd)
    if worst[2]<0.9:
        raise RuntimeError(f"cap clash {worst[0]}-{worst[1]} at {worst[2]:.2f} A")
    return out, worst, cn, [e["name"] for e in extra]

SYBYL={"C":"C.3","N":"N.3","O":"O.3","S":"S.3","P":"P.3","H":"H",
       "F":"F","Cl":"Cl","Br":"Br","I":"I","Se":"Se"}

def write_mol2(atoms, bonds, capbonds, path, resname):
    """Atoms are written under element-safe placeholder names and the real names are mapped
    back into the .ac afterwards.  Residue atom names collide with element symbols all over
    this library -- every residue's alpha carbon is CA (calcium) and 61 aromatic carbons are
    named CL* (chlorine) -- and antechamber infers the element from the name, so writing the
    real names makes it silently build a different molecule.  Returns the placeholder order."""
    idx={a["name"]:i+1 for i,a in enumerate(atoms)}
    allb=list(bonds)+capbonds
    allb=[(x,y,o) for x,y,o in allb if x in idx and y in idx]
    with open(path,"w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{resname}\n{len(atoms):5d}{len(allb):6d}     1     0     0\n"
                "SMALL\nbcc\n\n@<TRIPOS>ATOM\n")
        counts={}; placeholders=[]
        for i,a in enumerate(atoms):
            el=a["el"]; ty=SYBYL.get(el, el)
            counts[el]=counts.get(el,0)+1
            # uppercase: prepgen parses atomic symbols case-sensitively and silently DROPS a
            # two-letter element it cannot match ("Cl1"), leaving a malformed atom with a nan
            # internal coordinate that only surfaces later as a tleap fatal error.
            ph=f"{el.upper()}{counts[el]}"
            placeholders.append(ph)
            f.write(f"{i+1:7d} {ph:<8s}{a['xyz'][0]:10.4f}{a['xyz'][1]:10.4f}"
                    f"{a['xyz'][2]:10.4f} {ty:<8s}1 {resname:<8s}0.000000\n")
        f.write("@<TRIPOS>BOND\n")
        for j,(x,y,o) in enumerate(allb):
            f.write(f"{j+1:6d}{idx[x]:6d}{idx[y]:6d} {o}\n")
    return placeholders

# ff14SB protein types for the shared backbone.  antechamber types the whole molecule with
# GAFF2, which leaves the backbone as ns/c/c3; tleap then has no bond, angle or dihedral
# parameters across the peptide bond to a neighbouring standard residue ("Could not find bond
# parameter for atom types: C - ns").  The side chain stays GAFF2, the backbone does not.
# CA is typed CT, not ff14SB's protein-specific CX: these residues are scored under ff03,
# whose alpha carbon is CT and which has no CX at all ("could not find vdW parameters for
# type (CX)").  CT exists in ff03, ff99, ff02 and ff14SB, so one prep serves every arm.
BACKBONE_TYPE={"N":"N","H":"H","CA":"CT","HA":"H1","C":"C","O":"O"}

# ff03/ff99 call the alpha carbon CT; ff14SB and ff19SB introduced CX for it.  A prep must
# be typed for the force field it will be scored under, or the residue gets generic sp3
# torsions instead of the protein ones.
PROTEIN_FF={"ff03": ("CT","oldff/leaprc.ff03"),
            "ff14SB": ("CX","leaprc.protein.ff14SB"),
            "ff19SB": ("CX","leaprc.protein.ff19SB")}

def retype_backbone(ac_path, real_names, placeholders, cap=None, ca_type="CT"):
    """Rewrite the backbone atoms' GAFF2 types to their protein equivalents, in the .ac.

    The ACE/NME cap atoms are retyped too.  They stand in for the neighbouring residues, so
    leaving them GAFF2-typed makes parmchk2 generate the cross terms against the wrong types:
    a proline-like residue whose side chain is bonded to N needs an angle C-N-c5, and without
    this it only ever gets c-N-c5, which tleap cannot use."""
    ph={r:p for r,p in zip(real_names, placeholders)}
    bb=dict(BACKBONE_TYPE); bb["CA"]=ca_type
    want={ph[r]:t for r,t in bb.items() if r in ph}
    for nm,t in (cap or {}).items():
        if nm in ph: want[ph[nm]]=t
    done={}; out=[]
    for l in open(ac_path, errors="ignore"):
        if l.startswith("ATOM"):
            nm=l[12:16].strip()
            if nm in want:
                old=l[72:].strip()
                l=l[:72]+f"{want[nm]:>2s}\n"
                done[nm]=(old,want[nm])
        out.append(l)
    open(ac_path,"w").write("".join(out))
    return {r:done[ph[r]] for r in BACKBONE_TYPE if r in ph and ph[r] in done}

def restore_prep_names(prep_path, mapping):
    """Rewrite placeholder atom names back to the library's real names in the .prep.

    prepgen runs on placeholders on purpose: it rejects real names whose leading characters
    are not an element symbol (an iodine named IZ2), and the ACE/NME cap names collide with
    side-chain atoms.  The rename has to reach the LOOP and IMPROPER sections as well as the
    atom block -- a ring-closure bond left pointing at a placeholder makes tleap drop the
    residue ("LOOP atom C6 not found - bond not formed").  Atom *type* columns are never
    touched, because a placeholder like H1 is also a real ff14SB type name."""
    # antechamber rewrites a two-letter placeholder to proper element case (CL1 -> Cl1), so
    # the mapping has to be matched case-insensitively or the halogens keep their placeholder
    # names and tleap reports "does not have a type" for them.
    mapping={k.upper(): v for k, v in mapping.items()}
    out=[]; section=None
    for raw in open(prep_path, errors="ignore"):
        l=raw.rstrip("\n"); up=l.strip().upper()
        if up in ("LOOP","IMPROPER","DONE","CHARGE"):
            section=up; out.append(l+"\n"); continue
        if not l.strip():
            section=None; out.append(l+"\n"); continue
        if section in ("LOOP","IMPROPER"):
            toks=l.split()
            if toks:
                new_l=l
                for t in toks:
                    if t.upper() in mapping:
                        new_l=re.sub(rf"(?<![A-Za-z0-9]){re.escape(t)}(?![A-Za-z0-9])",
                                     mapping[t.upper()], new_l, count=1)
                l=new_l
            out.append(l+"\n"); continue
        m=re.match(r"^(\s*\d+\s+)(\S+)(.*)$", l)
        if m and m.group(2).upper() in mapping:
            old=m.group(2); newn=mapping[old.upper()]; rest=m.group(3)
            if len(newn)<=len(old):
                l=m.group(1)+newn+" "*(len(old)-len(newn))+rest
            else:
                extra=len(newn)-len(old)
                rest=rest[extra:] if rest[:extra].strip()=="" else rest
                l=m.group(1)+newn+rest
        out.append(l+"\n")
    open(prep_path,"w").write("".join(out))


def clean_frcmod(path):
    """Drop anything parmchk2 invented for the ff14SB backbone types.

    parmchk2 is told the molecule is GAFF2, so it does not recognise N/CX/C/O/H/H1 and
    writes MASS and NONBON entries for them ("same as c3").  Loading that frcmod would
    silently redefine ff14SB's own backbone.  Terms that MIX a backbone type with a GAFF2
    side-chain type are genuinely missing from both force fields and are kept."""
    if not os.path.exists(path): return 0
    protein={"N","CT","CX","C","O","H","H1"}
    sec=None; out=[]; dropped=0
    for raw in open(path, errors="ignore"):
        l=raw.rstrip("\n"); up=l.strip().upper()
        if up in ("MASS","BOND","ANGLE","DIHE","IMPROPER","NONBON"):
            sec=up; out.append(l+"\n"); continue
        if not l.strip():
            out.append(l+"\n"); continue
        if sec in ("MASS","NONBON"):
            t=l.split()[0] if l.split() else ""
            if t in protein: dropped+=1; continue
        elif sec in ("BOND","ANGLE","DIHE","IMPROPER"):
            types=[x for x in re.split(r"[-\s]+", l[:11 if sec=="BOND" else 14]) if x]
            if types and all(t in protein for t in types): dropped+=1; continue
        out.append(l+"\n")
    open(path,"w").write("".join(out))
    return dropped


def write_mc(path, ph, caps):
    """ph maps real name -> placeholder; the mainchain file must speak placeholders too."""
    with open(path,"w") as f:
        f.write(f"HEAD_NAME {ph['N']}\nTAIL_NAME {ph['C']}\nMAIN_CHAIN {ph['CA']}\n")
        for a in caps: f.write(f"OMIT_NAME {ph[a]}\n")
        f.write("PRE_HEAD_TYPE C\nPOST_TAIL_TYPE N\nCHARGE 0.0\n")

def sh(cmd, cwd, log):
    with open(log,"a") as lg:
        lg.write(f"\n$ {cmd}\n"); lg.flush()
        return subprocess.run(cmd, shell=True, cwd=cwd, stdout=lg, stderr=subprocess.STDOUT).returncode

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("code"); ap.add_argument("libdir"); ap.add_argument("outdir")
    ap.add_argument("--atype", default="gaff2")
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--protein-ff", default="ff03", choices=sorted(PROTEIN_FF),
                    help="force field this prep will be scored under")
    ap.add_argument("--net-charge", type=int, default=0,
                    help="formal charge; a side-chain NH3+ makes the residue +1")
    ap.add_argument("--cap-h", default=None,
                    help="add a hydrogen on this atom for the charge fit, then omit it "
                         "(for a sulfur that is cross-linked in the real structure)")
    a=ap.parse_args()
    global CA_TYPE, VAL_LEAPRC
    CA_TYPE, VAL_LEAPRC = PROTEIN_FF[a.protein_ff]
    os.makedirs(a.outdir, exist_ok=True)
    log=os.path.join(a.outdir,"build.log"); open(log,"w").close()
    res={"code":a.code,"atype":a.atype,"net_charge":a.net_charge,
     "protein_ff":a.protein_ff,"ca_type":CA_TYPE}
    try:
        atoms=first_conformer(f"{a.libdir}/{a.code}/{a.code}_rotamer.pdb",
                              f"{a.libdir}/{a.code}/{a.code}.params")
        res["n_residue_atoms"]=len(atoms)
        res["elements"]=sorted({x["el"].upper() for x in atoms})
        bonds=params_bonds(f"{a.libdir}/{a.code}/{a.code}.params")
        # an atom with no bond is silently dropped by prepgen; refuse to build one
        bonded={x for b in bonds for x in b[:2]}
        orphan=[x["name"] for x in atoms if x["name"] not in bonded]
        if orphan: raise RuntimeError(f"atoms with no bond record: {orphan}")
        capped,worst,cn,extra_omit=build_capped(atoms, bonds, a.cap_h)
        res["n_capped_atoms"]=len(capped)
        res["closest_contact_A"]=round(worst[2],3)
        capbonds=[(x, a.cap_h, "1") for x in extra_omit]+[(cn["CY"],"N","1"),(cn["CY"],cn["OY"],"2"),(cn["CY"],cn["CAY"],"1"),
                  (cn["CAY"],cn["HY1"],"1"),(cn["CAY"],cn["HY2"],"1"),(cn["CAY"],cn["HY3"],"1"),
                  ("C",cn["NZ"],"1"),(cn["NZ"],cn["HZ"],"1"),(cn["NZ"],cn["CAZ"],"1"),
                  (cn["CAZ"],cn["HZ1"],"1"),(cn["CAZ"],cn["HZ2"],"1"),(cn["CAZ"],cn["HZ3"],"1")]
        placeholders=write_mol2(capped, bonds, capbonds, f"{a.outdir}/capped.mol2", a.code)
        real=[x["name"] for x in capped]
        ph=dict(zip(real, placeholders))
        if len(set(real))!=len(real):
            dup=sorted({n for n in real if real.count(n)>1})
            raise RuntimeError(f"duplicate atom names after capping: {dup}")
        write_mc(f"{a.outdir}/{a.code}.mc", ph, [cn[k] for k,_ in CAP_SPEC]+extra_omit)
        res["extra_omitted"]=extra_omit
        res["cap_names_renamed"]={k:v for k,v in cn.items() if k!=v}
        rc=sh(f"antechamber -i capped.mol2 -fi mol2 -o {a.code}.ac -fo ac -c bcc -nc {a.net_charge} "
              f"-at {a.atype} -rn {a.code} -pf y -dr no -j 4", a.outdir, log)
        res["antechamber_rc"]=rc
        if rc!=0 or not os.path.exists(f"{a.outdir}/{a.code}.ac"):
            res["status"]="antechamber_failed"; raise RuntimeError

        capmap={cn["CY"]:"C", cn["OY"]:"O", cn["NZ"]:"N", cn["HZ"]:"H"}
        res["backbone_retyped"]={k:f"{v[0]}->{v[1]}" for k,v in
                                 retype_backbone(f"{a.outdir}/{a.code}.ac", real, placeholders,
                                                 capmap, CA_TYPE).items()}
        rc=sh(f"prepgen -i {a.code}.ac -o {a.code}.prep -f prepi -m {a.code}.mc -rn {a.code}",
              a.outdir, log)
        res["prepgen_rc"]=rc
        if not os.path.exists(f"{a.outdir}/{a.code}.prep"):
            res["status"]="prepgen_failed"; raise RuntimeError
        # parmchk2's prepi reader segfaults on residues carrying Br/I/S; the .ac reader is
        # sound, and it is the .ac that holds the retyped backbone, so the CX-to-sidechain
        # cross terms (in neither ff14SB nor GAFF2) get generated from it.
        ptxt=open(f"{a.outdir}/{a.code}.prep", errors="ignore").read()
        if "nan" in ptxt.lower():
            res["status"]="prep_has_nan"
            json.dump(res, open(f"{a.outdir}/result.json","w"), indent=1); print(json.dumps(res)); return
        rc=sh(f"parmchk2 -i {a.code}.ac -f ac -o {a.code}.frcmod -s {a.atype}", a.outdir, log)
        res["parmchk2_rc"]=rc
        res["frcmod_protein_entries_dropped"]=clean_frcmod(f"{a.outdir}/{a.code}.frcmod")
        if not os.path.exists(f"{a.outdir}/{a.code}.frcmod"):
            res["status"]="parmchk2_failed"
            json.dump(res, open(f"{a.outdir}/result.json","w"), indent=1); print(json.dumps(res)); return
        restore_prep_names(f"{a.outdir}/{a.code}.prep", dict(zip(placeholders, real)))
        # how many parameters had to be guessed?
        fr=f"{a.outdir}/{a.code}.frcmod"
        txt=open(fr).read() if os.path.exists(fr) else ""
        res["frcmod_attn"]=txt.count("ATTN")
        # GAFF2 has no selenium; antechamber silently retypes Se as sulfur and parmchk2 then
        # reports nothing missing.  A silent element swap must never look like a clean build.
        # the real test: does tleap actually build a peptide with this residue in it?
        leap=(f"source {VAL_LEAPRC}\nsource leaprc.gaff2\n"
              f"loadAmberPrep {a.code}.prep\nloadamberparams {a.code}.frcmod\n"
              f"p = sequence {{ ALA {a.code} ALA }}\n"
              f"saveamberparm p val.prmtop val.inpcrd\nquit\n")
        open(f"{a.outdir}/val.leap","w").write(leap)
        sh("tleap -f val.leap > val.log 2>&1", a.outdir, log)
        vlog=open(f"{a.outdir}/val.log", errors="ignore").read() if os.path.exists(f"{a.outdir}/val.log") else ""
        res["tleap_builds_peptide"]=os.path.exists(f"{a.outdir}/val.prmtop")
        res["tleap_fatal"]=vlog.count("FATAL")
        res["tleap_missing_param"]=len(re.findall(r"[Cc]ould not find|missing", vlog))
        if not res["tleap_builds_peptide"] or res["tleap_missing_param"]:
            # tleap still writes a prmtop when a parameter is missing, substituting a default;
            # that must not count as a successful build
            res["status"]="tleap_failed" if not res["tleap_builds_peptide"] else "tleap_missing_param"
            json.dump(res, open(f"{a.outdir}/result.json","w"), indent=1); print(json.dumps(res)); return
        unsupported=sorted({x["el"].upper() for x in atoms} - GAFF2_ELEMENTS)
        res["unsupported_elements"]=unsupported
        res["status"]="ok_element_unsupported" if unsupported else "ok"
    except Exception as e:
        res.setdefault("status","error"); res["error"]=str(e)
    json.dump(res, open(f"{a.outdir}/result.json","w"), indent=1)
    print(json.dumps(res))

if __name__=="__main__": main()
