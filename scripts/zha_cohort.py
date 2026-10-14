#!/usr/bin/env python3
"""Zhao et al. 2025 MM/PBSA(GBSA) protocol applied to the five-family FEP+ cohort.

Same minimisation and scoring as zha_pipeline.py.  The difference is the ligand:
the paper restricted itself to cyclic peptides built only from canonical amino
acids, so that a protein force field covers the whole complex.  This cohort is
macrocyclic peptides carrying non-canonical residues, prepared as a single
Rosetta LIG residue, so the ligand is parameterised with GAFF2 and AM1-BCC.
That is the one substantive departure from the published protocol.
"""
import argparse, json, os, re, subprocess, sys

LEAPRC = {"ff03": "oldff/leaprc.ff03", "ff14SB": "leaprc.protein.ff14SB",
          "ff14SBonlysc": "leaprc.protein.ff14SBonlysc",
          "ff99": "oldff/leaprc.ff99", "ff02": "oldff/leaprc.ff02"}

DROP_HET = {"HOH", "WAT", "GOL", "EDO", "SO4", "PO4", "ACT", "DMS", "PEG",
            "MPD", "TRS", "IMD", "FMT", "CIT", "NO3", "CL", "NA", "K", "BR", "IOD"}

# Modified amino acids deposited as HETATM but covalently inside a protein chain.
# JAK2 carries two activation-loop phosphotyrosines; dropping them opens a 7.7 A
# break that tleap then caps with a charged N- and C-terminus in the middle of the
# loop.  phosaa10 parameterises them from the ff94 charge set, which is the same
# family ff03 takes its termini from.
CHAIN_HET = {"PTR", "SEP", "TPO", "NEP", "HIP", "MSE", "CSO", "CME"}


def run(cmd, cwd, log, timeout=None):
    with open(log, "ab") as fh:
        fh.write(b"\n$ " + cmd.encode() + b"\n")
        try:
            return subprocess.run(cmd, shell=True, cwd=cwd, stdout=fh,
                                  stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            fh.write(b"\n[TIMEOUT]\n"); return 124


def split_complex(src, out):
    """Receptor = the polymer; everything heteroatomic is dropped, including the
    crystallisation additives that are not part of the binding model."""
    keep, het, kept_het = [], set(), set()
    for l in open(src, errors="ignore"):
        if l.startswith("ATOM"):
            keep.append(l)
        elif l.startswith("HETATM"):
            rn = l[17:20].strip()
            if rn in CHAIN_HET:
                keep.append("ATOM  " + l[6:]); kept_het.add(rn)
            else:
                het.add(rn)
        elif l.startswith("TER"):
            keep.append(l)
    with open(os.path.join(out, "rec_raw.pdb"), "w") as fh:
        fh.writelines(keep); fh.write("END\n")
    return sorted(het), sorted(kept_het)


def pose_drift(sdf, mol2):
    """Maximum atom displacement between the input pose and what antechamber wrote.
    Anything above a rounding tolerance means the pose was replaced."""
    try:
        import math
        n = int(open(sdf).read().splitlines()[3][:3])
        a = [tuple(float(x) for x in l.split()[:3])
             for l in open(sdf).read().splitlines()[4:4 + n]]
        lines = open(mol2).read().splitlines()
        i = lines.index("@<TRIPOS>ATOM") + 1
        b = [tuple(float(x) for x in lines[i + k].split()[2:5]) for k in range(n)]
        return round(max(math.dist(p, q) for p, q in zip(a, b)), 4)
    except Exception:
        return None


def build(prep, out, ff, charge, log):
    het, kept_het = split_complex(os.path.join(prep, "native_complex.pdb"), out)
    run("pdb4amber -i rec_raw.pdb -o rec_h.pdb", out, log)
    run("pdb4amber -i rec_h.pdb -o rec_clean.pdb -y", out, log)
    if not os.path.exists(os.path.join(out, "rec_clean.pdb")):
        return False, {"stage": "pdb4amber", "het_dropped": het, "het_kept": kept_het}
    sdf = os.path.join(prep, "ligand_v2000.sdf")
    # maxcyc=0 makes the AM1 step a single point at the input geometry.  Antechamber
    # otherwise runs a full gas-phase optimisation and writes the *optimised*
    # coordinates into the mol2, which would silently replace the bound pose the
    # whole calculation is about -- and takes tens of minutes on a macrocycle.
    ek = "maxcyc=0, qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700"
    rc = run(f"antechamber -i {sdf} -fi sdf -o lig.mol2 -fo mol2 -c bcc -nc {charge} "
             f"-at gaff2 -pf y -dr no -ek \"{ek}\"", out, log, timeout=14400)
    if not os.path.exists(os.path.join(out, "lig.mol2")):
        return False, {"stage": "antechamber", "rc": rc, "het_dropped": het, "het_kept": kept_het}
    run("parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s gaff2", out, log)
    drift = pose_drift(sdf, os.path.join(out, "lig.mol2"))
    lines = [f"source {LEAPRC[ff]}", "source leaprc.gaff2"]
    if kept_het:
        lines.append("source leaprc.phosaa10")
    lines += ["set default PBRadii mbondi2",
             "rec = loadpdb rec_clean.pdb",
             "loadamberparams lig.frcmod",
             "lig = loadmol2 lig.mol2",
             "saveamberparm rec rec.prmtop rec.inpcrd",
             "saveamberparm lig lig.prmtop lig.inpcrd",
             "com = combine { rec lig }",
             "saveamberparm com com.prmtop com.inpcrd",
             "savepdb com com.pdb", "quit", ""]
    open(os.path.join(out, "leap.in"), "w").write("\n".join(lines))
    run("tleap -f leap.in", out, log)
    ok = all(os.path.exists(os.path.join(out, f))
             for f in ("com.prmtop", "rec.prmtop", "lig.prmtop", "com.inpcrd"))
    return ok, {"stage": "leap" if not ok else "built", "het_dropped": het,
                "het_kept": kept_het, "pose_drift_A": drift}


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
        out, log, timeout=86400)
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


def score(out, method, indi, log):
    d = os.path.join(out, f"{method}_{indi}")
    os.makedirs(d, exist_ok=True)
    if method == "pb":
        body = (f"&general\n  interval=1, verbose=1,\n/\n"
                f"&pb\n  indi={indi}, exdi=80.0, inp=1,\n"
                f"  cavity_surften=0.0072, cavity_offset=0.0,\n/\n")
    elif method == "pb2":
        body = f"&general\n  interval=1, verbose=1,\n/\n&pb\n  indi={indi}, exdi=80.0, inp=2,\n/\n"
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
        "-y ../min.nc -o FINAL.dat -prefix _MM_", d, log, timeout=7200)
    for f in os.listdir(d):
        if f.startswith("_MM_"):
            try: os.remove(os.path.join(d, f))
            except OSError: pass
    fin = os.path.join(d, "FINAL.dat")
    return parse_final(fin), parse_components(fin)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slug"); ap.add_argument("prepared"); ap.add_argument("outdir")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--ff", default="ff03")
    ap.add_argument("--indi", default="1,2,4,6,8")
    ap.add_argument("--methods", default="pb,gb2")
    ap.add_argument("--min-engine", default="sander", choices=["sander", "openmm"])
    a = ap.parse_args()
    out = os.path.abspath(a.outdir); os.makedirs(out, exist_ok=True)
    log = os.path.join(out, "run.log"); open(log, "w").close()
    res = {"slug": a.slug, "ff": a.ff, "charge": a.charge, "min_engine": a.min_engine}
    ok, meta = build(a.prepared, out, a.ff, a.charge, log)
    res.update(meta)
    if not ok:
        res["status"] = "build_failed"
        json.dump(res, open(os.path.join(out, "result.json"), "w")); print(json.dumps(res)); return 1
    if minimise(out, log, a.min_engine) is None:
        res["status"] = "min_failed"
        json.dump(res, open(os.path.join(out, "result.json"), "w")); print(json.dumps(res)); return 1
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
    print(json.dumps(res)); return 0


if __name__ == "__main__":
    sys.exit(main())
