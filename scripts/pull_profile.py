#!/usr/bin/env python3
"""Score a complex by the SHAPE of its binding well, not by a single minimised geometry.

The ligand is translated rigidly away from the receptor along the outward radial direction in
a series of offsets, and MM/PBSA is evaluated at each one.  The resulting energy-versus-
distance profile gives features a single-point score cannot: how steeply the interaction
decays, how far it persists, and the area under the curve.

The motivation is that the single-point score's worst failure mode is geometric -- on KRAS the
spread caused by where the minimiser stopped was as large as the entire between-compound
signal (ratio 1.02).  A profile averages over a controlled set of geometries instead of
trusting one.

The thing to watch is that a larger ligand has a deeper and wider well for trivial reasons, so
every feature here must be checked against a heavy-atom count before it is believed.
"""
import argparse, os, subprocess, sys
import numpy as np
import parmed as pmd

MMPBSA_IN = """displacement profile, Zha25 settings
&general
  interval=1, verbose=2,
/
&pb
  indi={indi}, exdi=80.0, inp=1,
  cavity_surften=0.0072, cavity_offset=0.0,
/
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir", help="a cohort complex directory (com.prmtop, min.nc)")
    ap.add_argument("outdir")
    ap.add_argument("--offsets", default="0,0.5,1,1.5,2,2.5,3,4,5,6,8")
    ap.add_argument("--indi", type=float, default=2.0)
    ap.add_argument("--ligand-names", default="MOL,LIG")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    offs = [float(x) for x in a.offsets.split(",")]

    top = os.path.join(a.workdir, "com.prmtop")
    crd = os.path.join(a.workdir, "min.nc")
    if not (os.path.exists(top) and os.path.exists(crd)):
        print(f"SKIP {a.workdir}: missing topology or minimised coordinates")
        return 1

    p = pmd.load_file(top, xyz=crd)
    names = set(a.ligand_names.split(","))
    lig_idx = np.array([at.idx for r in p.residues if r.name in names for at in r.atoms])
    if lig_idx.size == 0:
        print(f"SKIP {a.workdir}: no ligand residue")
        return 2
    rec_idx = np.setdiff1d(np.arange(len(p.atoms)), lig_idx)

    xyz = np.array(p.coordinates, dtype=float)
    mass = np.array([at.mass for at in p.atoms], dtype=float)
    lig_com = np.average(xyz[lig_idx], axis=0, weights=mass[lig_idx])

    # Pull direction matters more than anything else here.  Ligand COM minus whole-receptor
    # COM drives a groove-bound peptide into the walls of its own pocket: VDWAALS reached
    # +2327 kcal/mol at 5 A instead of decaying to zero.  The local interface normal -- from
    # the centroid of receptor atoms actually contacting the ligand out to the ligand --
    # points out of the site.  Candidates are scored by how cleanly the ligand escapes.
    from scipy.spatial import cKDTree
    tree = cKDTree(xyz[rec_idx])
    near = set()
    for j in tree.query_ball_point(xyz[lig_idx], r=5.0):
        near.update(j)
    contact = rec_idx[sorted(near)] if near else rec_idx
    base = lig_com - xyz[contact].mean(axis=0)
    base /= max(np.linalg.norm(base), 1e-9)

    def clearance(u):
        """Closest ligand-receptor approach encountered while sliding out along u."""
        worst = 9e9
        for d in (2.0, 4.0, 6.0, 8.0):
            worst = min(worst, tree.query(xyz[lig_idx] + u * d)[0].min())
        return worst

    best, v = clearance(base), base
    rng = np.random.default_rng(0)
    for _ in range(60):
        cand = base + rng.normal(scale=0.45, size=3)
        cand /= max(np.linalg.norm(cand), 1e-9)
        if cand @ base < 0.3:
            continue
        c = clearance(cand)
        if c > best:
            best, v = c, cand
    # one frame per offset, ligand rigidly translated outward
    frames = []
    for d in offs:
        c = xyz.copy()
        c[lig_idx] += v * d
        frames.append(c)
    out_nc = os.path.join(a.outdir, "scan.nc")
    from parmed.amber import NetCDFTraj
    t = NetCDFTraj.open_new(out_nc, len(p.atoms), box=False, crds=True)
    for f in frames:
        t.add_coordinates(f)
    t.close()

    open(os.path.join(a.outdir, "mmpbsa.in"), "w").write(MMPBSA_IN.format(indi=a.indi))
    # absolute paths: the output directory need not live inside the complex directory
    W = os.path.abspath(a.workdir)
    cmd = ("MMPBSA.py -O -i mmpbsa.in -o FINAL.dat -eo per_frame.csv "
           f"-cp {W}/com.prmtop -rp {W}/rec.prmtop -lp {W}/lig.prmtop -y scan.nc")
    r = subprocess.run(cmd, shell=True, cwd=a.outdir, capture_output=True, text=True)
    ok = os.path.exists(os.path.join(a.outdir, "per_frame.csv"))
    print(f"{os.path.basename(a.workdir)}: frames={len(offs)} mmpbsa_ok={ok} "
          f"pull=({v[0]:.2f},{v[1]:.2f},{v[2]:.2f}) clearance={best:.2f}A")
    if not ok:
        sys.stderr.write((r.stdout or "")[-1200:])
        return 4
    return 0


if __name__ == "__main__":
    sys.exit(main())
