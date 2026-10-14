#!/usr/bin/env python3
"""Zhao et al. 2025 (Brief Bioinform 26:bbaf632) MM/PBSA(GBSA) protocol.

Per complex: pdb4amber clean-up -> tleap (with explicit cyclisation bonds) ->
implicit-solvent minimisation (igb=2, maxcyc 5000, ncyc 2000, cut 1000) ->
single-frame MM/PBSA and MM/GBSA over a grid of interior dielectric constants.
Entropy is excluded, as in the paper.

The released script builds filename strings only and never emits the tleap
stage, so topology construction here is ours.  Disulfides are found from the
coordinates; head-to-tail closure follows the paper's own `Chain type` column,
because two of the deposited macrocycles have a closure distance of 2.7 A that
no distance cutoff would accept.
"""
import argparse, json, os, re, subprocess, sys
import numpy as np

LEAPRC = {
    # the paper's script names "leaprc.ff03", which is the legacy leaprc: it takes
    # its terminal residues from all_amino{n,c}t94.lib and so carries NCYX/CCYX.
    # leaprc.protein.ff03.r1 uses amino{n,c}t03.lib, which has neither, and half
    # of Dataset I is a disulfide-cyclised peptide with a terminal cystine.
    "ff03": "oldff/leaprc.ff03",
    "ff14SB": "leaprc.protein.ff14SB",
    "ff14SBonlysc": "leaprc.protein.ff14SBonlysc",
    "ff99": "oldff/leaprc.ff99",
    "ff02": "oldff/leaprc.ff02",
    # off-label with GB: ff19SB's CMAP was fitted for OPC water, and CMAP applies
    # only to the 20 standard residues, so a GAFF2-hybrid ncAA never receives it
    "ff19SB": "leaprc.protein.ff19SB",
}


def run(cmd, cwd, log, timeout=None):
    with open(log, "ab") as fh:
        fh.write(b"\n$ " + cmd.encode() + b"\n")
        try:
            return subprocess.run(cmd, shell=True, cwd=cwd, stdout=fh,
                                  stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            fh.write(b"\n[TIMEOUT]\n"); return 124


STD_AA = {"ALA","ARG","ASN","ASP","CYS","CYX","GLN","GLU","GLY","HIS","HID","HIE","HIP",
          "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


def rename_to_template(lines, ncaa_dir):
    """Give each non-standard residue the atom names its Amber template uses.

    The curated library's rotamer PDB truncates four-character halogen names (CLM1 -> CLM) and
    the structures inherit that, while the .prep is built from the .params, where the names are
    unique.  tleap then cannot match the residue to its template.  Atom order is identical
    between the two, so the rename is positional."""
    import json as _json
    out=[]; buf=[]; cur=None; renamed={}
    def flush():
        if cur is None: return
        f=os.path.join(ncaa_dir, cur, "atom_names.json")
        if os.path.exists(f):
            names=_json.load(open(f))
            if len(names)==len(buf):
                for k,l in enumerate(buf):
                    nm=names[k]
                    old=l[12:16].strip()
                    if old!=nm: renamed.setdefault(cur,set()).add(f"{old}->{nm}")
                    nmf=f" {nm:<3s}" if len(nm)<4 else nm[:4]
                    buf[k]=l[:12]+nmf+l[16:]
        out.extend(buf)
    for l in lines:
        if l.startswith(("ATOM","HETATM")):
            rn=l[17:20].strip(); key=(rn, l[21], l[22:27])
            if cur is not None and key!=curkey:
                flush(); buf=[]
            if rn not in STD_AA:
                cur, curkey = rn, key; buf.append(l); continue
            flush(); buf=[]; cur=None
            out.append(l)
        else:
            flush(); buf=[]; cur=None
            out.append(l)
    flush()
    return out, {k:sorted(v) for k,v in renamed.items()}


def promote_het(src, dst):
    """Rosetta writes non-canonical residues as HETATM; pdb4amber drops them.  Promote every
    HETATM that is part of the peptide chain to ATOM, and return the residue names seen."""
    seen=set(); out=[]
    for l in open(src, errors='ignore'):
        if l.startswith('HETATM'):
            rn=l[17:20].strip()
            if rn not in ('HOH','WAT'):
                seen.add(rn); l='ATOM  '+l[6:]
        out.append(l)
    open(dst,'w').write(''.join(out))
    return sorted(seen)


def nonstandard_residues(path):
    names=[]
    for l in open(path, errors='ignore'):
        if l.startswith(('ATOM','HETATM')):
            rn=l[17:20].strip()
            if rn not in STD_AA and rn not in names: names.append(rn)
    return names


def read_pdb(path):
    out = []
    for l in open(path, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")) and l[16] in " A":
            e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            out.append(dict(nm=l[12:16].strip(), res=l[17:20].strip(), ri=l[22:27].strip(),
                            el=e, xyz=(float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out


def analyse_peptide(path, head_to_tail, charged_ref=None):
    """Cyclisation bonds and terminus state, on pdb4amber-renumbered residues."""
    a = read_pdb(path)
    order = []
    for x in a:
        if x["ri"] not in order:
            order.append(x["ri"])
    pos = {r: i + 1 for i, r in enumerate(order)}
    X = np.array([x["xyz"] for x in a])
    S = [i for i, x in enumerate(a) if x["el"] == "S" and x["nm"] == "SG"]
    ss = []
    for ii, i in enumerate(S):
        for j in S[ii + 1:]:
            if np.linalg.norm(X[i] - X[j]) < 2.5:
                ss.append((pos[a[i]["ri"]], pos[a[j]["ri"]]))
    first = {x["nm"] for x in a if x["ri"] == order[0]}
    last = {x["nm"] for x in a if x["ri"] == order[-1]}
    charged_term = bool(first & {"H1", "H2", "H3"}) or ("OXT" in last)
    if charged_ref is not None:
        charged_term = charged_ref
    amide = [(len(order), 1)] if head_to_tail else []
    return ss, amide, charged_term, len(order)


def build(prot_src, pep_src, out, ff, head_to_tail, log, ncaa_dir=None):
    # pass 1 keeps the prepared hydrogens, which is what lets pdb4amber carry the
    # histidine tautomer and cystine state over into the residue names; pass 2
    # drops the hydrogens so tleap rebuilds them to the force field's own
    # templates (the prepared files carry e.g. HXT, which ff03 has no type for)
    ncaa = []; renamed = {}
    if ncaa_dir:
        promoted = os.path.join(out, "lig_promoted.pdb")
        promote_het(pep_src, promoted)
        lines, renamed = rename_to_template(open(promoted, errors="ignore").readlines(), ncaa_dir)
        open(promoted, "w").write("".join(lines))
        pep_src = promoted
        ncaa = nonstandard_residues(pep_src)
    for src, stem in ((prot_src, "rec"), (pep_src, "lig")):
        run(f"pdb4amber -i {src} -o {stem}_h.pdb", out, log)
        run(f"pdb4amber -i {stem}_h.pdb -o {stem}_clean.pdb -y", out, log)
    pep = os.path.join(out, "lig_clean.pdb")
    if not os.path.exists(pep):
        return False, {"status": "pdb4amber_failed"}
    ref = analyse_peptide(os.path.join(out, "lig_h.pdb"), head_to_tail)[2] \
        if os.path.exists(os.path.join(out, "lig_h.pdb")) else None
    ss, amide, charged_term, nres = analyse_peptide(pep, head_to_tail, ref)
    lines = [f"source {LEAPRC[ff]}"]
    if ncaa:
        # the side chains are GAFF2-typed; the backbone was retyped to ff14SB in the prep so
        # that the peptide bond to the neighbouring standard residue has parameters
        lines.append("source leaprc.gaff2")
        for rn in ncaa:
            lines.append(f"loadAmberPrep {ncaa_dir}/{rn}/{rn}.prep")
            lines.append(f"loadamberparams {ncaa_dir}/{rn}/{rn}.frcmod")
    lines += ["set default PBRadii mbondi2", 'rec = loadpdb rec_clean.pdb']
    if amide and not charged_term:
        lines.append("clearPdbResMap")          # a macrocycle has no termini
    lines.append('lig = loadpdb lig_clean.pdb')
    # disulfides arrive through pdb4amber's CONECT records; re-issuing them here
    # is what tleap reports as "1-4: cannot add bond"
    for i, j in amide:
        lines.append(f"bond lig.{i}.C lig.{j}.N")
    lines += ["saveamberparm rec rec.prmtop rec.inpcrd",
              "saveamberparm lig lig.prmtop lig.inpcrd",
              "com = combine { rec lig }",
              "saveamberparm com com.prmtop com.inpcrd",
              "savepdb com com.pdb", "quit", ""]
    open(os.path.join(out, "leap.in"), "w").write("\n".join(lines))
    run("tleap -f leap.in", out, log)
    want = ("com.prmtop", "rec.prmtop", "lig.prmtop", "com.inpcrd")
    ok = all(os.path.exists(os.path.join(out, f)) for f in want)
    fallback = False
    if not ok and ff != "ff14SB":
        # ff03, ff99 and ff02 ship no NCYX/CCYX unit, so a peptide or chain whose
        # terminal residue is a cystine cannot be built from them at all.  Source
        # ff14SB first and let the requested force field overwrite every unit it
        # defines; only the missing terminal cystines are left coming from ff14SB.
        fallback = True
        lines[0] = "source leaprc.protein.ff14SB\n" + lines[0]
        open(os.path.join(out, "leap.in"), "w").write("\n".join(lines))
        run("tleap -f leap.in", out, log)
        ok = all(os.path.exists(os.path.join(out, f)) for f in want)
    return ok, dict(n_ss=len(ss), n_amide=len(amide), head_to_tail=bool(head_to_tail),
                    charged_term=charged_term, pep_res=nres, term_cyx_fallback=fallback,
                    ncaa=ncaa, ncaa_atoms_renamed=renamed if ncaa_dir else {})


MIN_IN = """Implicit-solvent minimisation, Zhao et al. 2025
 &cntrl
    imin=1, ntx=1, irest=0,
    maxcyc=5000, ncyc=2000,
    cut=1000.0, ntpr=5000, ntwx=0,
    igb=2, ntxo=1,
 /
"""


OMM_PYTHON = "python3"
OMM_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "min_openmm.py")


def minimise(out, log, engine="sander"):
    if engine == "openmm":
        # same force field and implicit solvent model, on the GPU; sander needs
        # hours per complex at the protocol's cut=1000 and the authors used
        # pmemd.cuda, which this conda AmberTools build does not ship
        # AmberTools sets PYTHONPATH/AMBERHOME; leaking those into the other
        # environment makes its OpenMM import the wrong numpy
        run(f"env -u PYTHONPATH -u PYTHONHOME -u AMBERHOME -u LD_LIBRARY_PATH "
            f"{OMM_PYTHON} {OMM_SCRIPT} .", out, log, timeout=7200)
        rst = os.path.join(out, "min_omm.rst")
        if not os.path.exists(rst):
            return None
        open(os.path.join(out, "traj.cpptraj"), "w").write(
            "parm com.prmtop\ntrajin min_omm.rst\ntrajout min.nc netcdf\nrun\nquit\n")
        run("cpptraj -i traj.cpptraj", out, log)
        p = os.path.join(out, "min.nc")
        return p if os.path.exists(p) else None
    open(os.path.join(out, "min.in"), "w").write(MIN_IN)
    run("sander -O -i min.in -p com.prmtop -c com.inpcrd -o min.out -r min.rst -inf min.info",
        out, log)
    if not os.path.exists(os.path.join(out, "min.rst")):
        return None
    open(os.path.join(out, "traj.cpptraj"), "w").write(
        "parm com.prmtop\ntrajin min.rst\ntrajout min.nc netcdf\nrun\nquit\n")
    run("cpptraj -i traj.cpptraj", out, log)
    p = os.path.join(out, "min.nc")
    return p if os.path.exists(p) else None


def parse_final(path):
    if not os.path.exists(path):
        return None
    m = re.findall(r"DELTA TOTAL\s+(-?\d+\.?\d*)", open(path, errors="ignore").read())
    return float(m[-1]) if m else None


def parse_components(path):
    """Keep the decomposition; an entropy-free MM/PBSA score is dominated by the
    size-extensive vdW and nonpolar terms and that has to be checkable."""
    if not os.path.exists(path):
        return {}
    txt = open(path, errors="ignore").read()
    blk = txt.split("Differences (Complex - Receptor - Ligand)")[-1]
    out = {}
    for k in ("VDWAALS", "EEL", "EPB", "EGB", "ENPOLAR", "EDISPER", "ESURF"):
        m = re.search(rf"^{k}\s+(-?\d+\.?\d*)", blk, re.M)
        if m:
            out[k] = float(m.group(1))
    return out


def _cleanup(d):
    for f in os.listdir(d):
        if f.startswith("_MM_"):
            try: os.remove(os.path.join(d, f))
            except OSError: pass


def _cleanup(d):
    for f in os.listdir(d):
        if f.startswith("_MM_"):
            try: os.remove(os.path.join(d, f))
            except OSError: pass


def score(out, method, indi, log):
    d = os.path.join(out, f"{method}_{indi}")
    os.makedirs(d, exist_ok=True)
    if method == "pb":
        # the paper's nonpolar term is LCPO with gamma 0.0072 and b 0, i.e. inp=1,
        # which was MMPBSA.py's default through Amber22 but is not AmberTools26's
        body = (f"&general\n  interval=1, verbose=1,\n/\n"
                f"&pb\n  indi={indi}, exdi=80.0, inp=1,\n"
                f"  cavity_surften=0.0072, cavity_offset=0.0,\n/\n")
    elif method == "pb2":
        body = (f"&general\n  interval=1, verbose=1,\n/\n"
                f"&pb\n  indi={indi}, exdi=80.0, inp=2,\n/\n")
    elif method.startswith("gb"):
        # AmberTools26's MMPBSA.py exposes no interior dielectric for GB: the &gb
        # namelist takes epsin but both energy engines ignore it.  The paper set it
        # by hand-writing the sander mdin, so do the same -- generate the mdin with
        # -make-mdins, inject intdiel, then run with -use-mdins.
        body = (f"&general\n  interval=1, verbose=1, use_sander=1,\n/\n"
                f"&gb\n  igb={method[-1]}, surften=0.0072,\n/\n")
        open(os.path.join(d, "mmpbsa.in"), "w").write("Zha25 replication\n" + body)
        base = ("MMPBSA.py -O -i mmpbsa.in -cp ../com.prmtop -rp ../rec.prmtop "
                "-lp ../lig.prmtop -y ../min.nc -o FINAL.dat -prefix _MM_")
        run(base + " -make-mdins", d, log)
        mdin = os.path.join(d, "_MM_gb.mdin")
        if os.path.exists(mdin):
            txt = open(mdin).read().replace(f"igb={method[-1]},",
                                            f"igb={method[-1]}, intdiel={indi}.0,")
            open(mdin, "w").write(txt)
        run(base + " -use-mdins", d, log)
        _cleanup(d)
        fin = os.path.join(d, "FINAL.dat")
        return parse_final(fin), parse_components(fin)
    else:
        body = (f"&general\n  interval=1, verbose=1,\n/\n"
                f"&gb\n  igb={method[-1]},\n/\n")
    open(os.path.join(d, "mmpbsa.in"), "w").write("Zha25 replication\n" + body)
    run("MMPBSA.py -O -i mmpbsa.in -cp ../com.prmtop -rp ../rec.prmtop -lp ../lig.prmtop "
        "-y ../min.nc -o FINAL.dat -prefix _MM_", d, log)
    _cleanup(d)
    fin = os.path.join(d, "FINAL.dat")
    return parse_final(fin), parse_components(fin)


def prmtop_natom(p):
    try:
        t = open(p, errors="ignore").read()
        return int(re.search(r"%FLAG POINTERS.*?%FORMAT\(10I8\)\s*\n\s*(\d+)", t, re.S).group(1))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdbid"); ap.add_argument("protein"); ap.add_argument("peptide")
    ap.add_argument("outdir")
    ap.add_argument("--ff", default="ff03")
    ap.add_argument("--indi", default="1,2,4,6,8")
    ap.add_argument("--methods", default="pb,gb1,gb2,gb5")
    ap.add_argument("--min-engine", default="sander", choices=["sander", "openmm"])
    ap.add_argument("--head-to-tail", action="store_true")
    ap.add_argument("--ncaa-dir", default=None,
                    help="directory of <CODE>/<CODE>.prep + .frcmod for non-standard residues")
    a = ap.parse_args()
    out = os.path.abspath(a.outdir); os.makedirs(out, exist_ok=True)
    log = os.path.join(out, "run.log"); open(log, "w").close()
    res = {"pdb": a.pdbid, "ff": a.ff, "min_engine": a.min_engine}
    ok, meta = build(os.path.abspath(a.protein), os.path.abspath(a.peptide), out,
                     a.ff, a.head_to_tail, log, a.ncaa_dir)
    res.update(meta)
    if not ok:
        res["status"] = "leap_failed"; print(json.dumps(res)); 
        json.dump(res, open(os.path.join(out, "result.json"), "w")); return 1
    res["natom"] = prmtop_natom(os.path.join(out, "com.prmtop"))
    if minimise(out, log, a.min_engine) is None:
        res["status"] = "min_failed"; print(json.dumps(res))
        json.dump(res, open(os.path.join(out, "result.json"), "w")); return 1
    res["status"] = "ok"
    for m in a.methods.split(","):
        for indi in a.indi.split(","):
            if not indi:
                continue
            tot, comp = score(out, m, indi, log)
            res[f"{m}_{indi}"] = tot
            for k, v in comp.items():
                res[f"{m}_{indi}_{k}"] = v
    json.dump(res, open(os.path.join(out, "result.json"), "w"))
    print(json.dumps(res))
    return 0


if __name__ == "__main__":
    sys.exit(main())
