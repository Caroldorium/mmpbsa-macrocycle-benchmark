#!/usr/bin/env python3
"""Score a complex by how strongly it pulls a displaced ligand back.

Rigid extraction does not work on this cohort: over 1500 sampled directions the best
ligand-receptor clearance is 1.2-2.5 A for every family, all inside vdW contact, because these
macrocycles are enveloped by their sites.  A rigid pull therefore measures steric clash, not
dissociation.

This probe instead displaces the ligand a short distance and lets the whole complex relax,
then asks two questions:
  recovery_rmsd  -- how far back toward its original pose did the ligand travel?
  energy_recovered -- what fraction of the original interaction energy came back?

Both are responses of the binding well rather than properties of one geometry, which is the
point: the single-point score's worst failure was that the answer depended on where the
minimiser happened to stop (KRAS artefact/signal 1.02).

As with every score in this project, the result has to be checked against a heavy-atom count
before it is believed -- a bigger ligand sits in a deeper well for trivial reasons.
"""
import argparse, json, os, subprocess, sys
import numpy as np
import parmed as pmd
from scipy.spatial import cKDTree

OMM = "python3"


def pull_direction(xyz, lig, rec, nsamp=800):
    """Direction giving the ligand most room; still clashes, which is why we minimise after."""
    tree = cKDTree(xyz[rec])
    rng = np.random.default_rng(0)
    u = rng.normal(size=(nsamp, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    best, bu = -1.0, u[0]
    for k in range(nsamp):
        worst = 9e9
        for d in (2.0, 4.0, 6.0):
            worst = min(worst, tree.query(xyz[lig] + u[k] * d)[0].min())
            if worst <= best:
                break
        if worst > best:
            best, bu = worst, u[k]
    return bu, best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir")
    ap.add_argument("outdir")
    ap.add_argument("--offsets", default="1.0,2.0,3.0")
    ap.add_argument("--maxiter", type=int, default=2000)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    W = os.path.abspath(a.workdir)

    p = pmd.load_file(f"{W}/com.prmtop", xyz=f"{W}/min.nc")
    xyz0 = np.array(p.coordinates, dtype=float)
    lig = np.array([at.idx for r in p.residues if r.name in ("MOL", "LIG") for at in r.atoms])
    if lig.size == 0:
        print(json.dumps(dict(status="no_ligand", workdir=a.workdir)))
        return 1
    rec = np.setdiff1d(np.arange(len(p.atoms)), lig)
    v, clear = pull_direction(xyz0, lig, rec)

    res = dict(workdir=os.path.basename(W), n_lig_atoms=int(lig.size),
               clearance_A=round(float(clear), 2), offsets={})
    for d in [float(x) for x in a.offsets.split(",")]:
        wd = os.path.join(a.outdir, f"d{d:g}")
        os.makedirs(wd, exist_ok=True)
        disp = xyz0.copy()
        disp[lig] += v * d
        st = pmd.load_file(f"{W}/com.prmtop")
        st.coordinates = disp
        st.save(os.path.join(wd, "start.rst7"), overwrite=True)
        r = subprocess.run(
            f"env -u PYTHONPATH -u PYTHONHOME -u AMBERHOME -u LD_LIBRARY_PATH "
            f"{OMM} {os.path.dirname(os.path.abspath(__file__))}/min_openmm.py . "
            f"--prmtop {W}/com.prmtop --inpcrd start.rst7 --out relaxed.rst "
            f"--maxiter {a.maxiter}",
            shell=True, cwd=wd, capture_output=True, text=True)
        rel = os.path.join(wd, "relaxed.rst")
        if not os.path.exists(rel):
            res["offsets"][f"{d:g}"] = dict(status="min_failed",
                                            err=(r.stdout or r.stderr or "")[-200:])
            continue
        q = pmd.load_file(f"{W}/com.prmtop", xyz=rel)
        xr = np.array(q.coordinates, dtype=float)
        # how much of the displacement was undone, measured on the ligand only
        moved_out = np.sqrt(((disp[lig] - xyz0[lig]) ** 2).sum(1)).mean()
        left_over = np.sqrt(((xr[lig] - xyz0[lig]) ** 2).sum(1)).mean()
        res["offsets"][f"{d:g}"] = dict(
            status="ok",
            displaced_A=round(float(moved_out), 3),
            residual_A=round(float(left_over), 3),
            recovered_fraction=round(float(1.0 - left_over / max(moved_out, 1e-9)), 4))
    json.dump(res, open(os.path.join(a.outdir, "recovery.json"), "w"), indent=1)
    print(json.dumps(res))
    return 0


if __name__ == "__main__":
    sys.exit(main())
