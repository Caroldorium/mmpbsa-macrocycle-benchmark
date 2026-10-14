#!/usr/bin/env python3
"""Zha25's two-step workflow on one cohort complex.

Step 1 re-ranks candidate poses with the polarity-appropriate method and keeps the
top one; step 2 minimises that pose in implicit solvent under ff03 and scores it with
MM/PBSA at interior dielectric 2.  Every complex in this cohort has a polar interface
ratio below 0.6, so step 1 is the paper's low-polarity method throughout: MM/GBSA with
GB_HCT, ff03, interior dielectric 1, implicit solvent.

The cohort ships one pose per compound and its three Rosetta replicas agree to 0.2 A,
so there is nothing to re-rank as delivered.  Decoys are generated here by rigid-body
perturbation of the ligand followed by the same minimisation, with the deposited pose
kept as candidate 0 -- which reproduces the condition the paper imposes on Dataset II,
that at least one candidate be near-native.

Topologies are reused from the single-pose run; only coordinates change.
"""
import argparse, glob, json, os, re, shutil, subprocess, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zha_pipeline import run, parse_final, parse_components, _cleanup, OMM_PYTHON, OMM_SCRIPT


def prmtop_natom(p):
    t = open(p, errors="ignore").read()
    return int(re.search(r"%FLAG POINTERS.*?%FORMAT\(10I8\)\s*\n\s*(\d+)", t, re.S).group(1))


def read_rst(p):
    L = open(p, errors="ignore").read().splitlines()
    n = int(L[1].split()[0]); v = []
    for l in L[2:]:
        for i in range(0, len(l), 12):
            s = l[i:i + 12].strip()
            if s:
                v.append(float(s))
        if len(v) >= 3 * n:
            break
    return np.array(v[:3 * n]).reshape(-1, 3)


def write_rst(p, xyz, title="perturbed"):
    with open(p, "w") as fh:
        fh.write(title + "\n%5d\n" % len(xyz))
        flat = xyz.reshape(-1)
        for i in range(0, len(flat), 6):
            fh.write("".join("%12.7f" % v for v in flat[i:i + 6]) + "\n")


def rand_rotation(rng, max_deg):
    ax = rng.normal(size=3); ax /= np.linalg.norm(ax)
    th = np.deg2rad(rng.uniform(0, max_deg))
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


def make_decoys(xyz, nlig, n, rng, max_deg, max_trans, rec_heavy_idx):
    """Rigid-body perturbations of the ligand block, rejecting poses driven into the core."""
    rec, lig = xyz[:-nlig], xyz[-nlig:]
    cen = lig.mean(0)
    out, tries = [], 0
    while len(out) < n and tries < n * 20:
        tries += 1
        R = rand_rotation(rng, max_deg)
        t = rng.normal(size=3); t = t / np.linalg.norm(t) * rng.uniform(0.5, max_trans)
        new = (lig - cen) @ R.T + cen + t
        d = np.linalg.norm(rec[rec_heavy_idx][:, None, :] - new[None, :, :], axis=2)
        if (d < 1.2).sum() > 20:          # driven into the protein core
            continue
        if d.min() > 6.0:                 # floated off the surface entirely
            continue
        out.append(np.vstack([rec, new]))
    return out


def minimise(d, log):
    # com.prmtop and com.inpcrd are already staged in d
    run(f"env -u PYTHONPATH -u PYTHONHOME -u AMBERHOME -u LD_LIBRARY_PATH "
        f"{OMM_PYTHON} {OMM_SCRIPT} .", d, log, timeout=3600)
    if not os.path.exists(os.path.join(d, "min_omm.rst")):
        return None
    open(os.path.join(d, "t.cpptraj"), "w").write(
        "parm com.prmtop\ntrajin min_omm.rst\ntrajout min.nc netcdf\nrun\nquit\n")
    run("cpptraj -i t.cpptraj", d, log)
    p = os.path.join(d, "min.nc")
    return p if os.path.exists(p) else None


def score(d, method, indi, rec, lig, log):
    if method == "pb":
        body = (f"&general\n  interval=1, verbose=1,\n/\n"
                f"&pb\n  indi={indi}, exdi=80.0, inp=1,\n"
                f"  cavity_surften=0.0072, cavity_offset=0.0,\n/\n")
        open(os.path.join(d, "s.in"), "w").write("two-step\n" + body)
        run(f"MMPBSA.py -O -i s.in -cp com.prmtop -rp {rec} -lp {lig} -y min.nc "
            f"-o FINAL_{method}{indi}.dat -prefix _MM_", d, log, timeout=3600)
    else:
        body = (f"&general\n  interval=1, verbose=1, use_sander=1,\n/\n"
                f"&gb\n  igb={method[-1]}, surften=0.0072,\n/\n")
        open(os.path.join(d, "s.in"), "w").write("two-step\n" + body)
        base = (f"MMPBSA.py -O -i s.in -cp com.prmtop -rp {rec} -lp {lig} -y min.nc "
                f"-o FINAL_{method}{indi}.dat -prefix _MM_")
        run(base + " -make-mdins", d, log)
        mdin = os.path.join(d, "_MM_gb.mdin")
        if os.path.exists(mdin):
            t = open(mdin).read().replace(f"igb={method[-1]},",
                                          f"igb={method[-1]}, intdiel={indi}.0,")
            open(mdin, "w").write(t)
        run(base + " -use-mdins", d, log)
    _cleanup(d)
    return parse_final(os.path.join(d, f"FINAL_{method}{indi}.dat"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slug"); ap.add_argument("built"); ap.add_argument("outdir")
    ap.add_argument("--ndecoy", type=int, default=20)
    ap.add_argument("--max-deg", type=float, default=60.0)
    ap.add_argument("--max-trans", type=float, default=5.0)
    ap.add_argument("--rerank", default="gb1")      # low-polarity branch
    ap.add_argument("--rerank-indi", default="1")
    ap.add_argument("--final", default="pb")
    ap.add_argument("--final-indi", default="2")
    a = ap.parse_args()

    out = os.path.abspath(a.outdir); os.makedirs(out, exist_ok=True)
    log = os.path.join(out, "run.log"); open(log, "w").close()
    prm = os.path.join(a.built, "com.prmtop"); ic = os.path.join(a.built, "com.inpcrd")
    rec = os.path.abspath(os.path.join(a.built, "rec.prmtop"))
    lig = os.path.abspath(os.path.join(a.built, "lig.prmtop"))
    res = {"slug": a.slug}
    if not all(os.path.exists(p) for p in (prm, ic, rec, lig)):
        res["status"] = "no_topology"; json.dump(res, open(os.path.join(out, "result.json"), "w"))
        print(json.dumps(res)); return 1

    nlig = prmtop_natom(lig)
    xyz = read_rst(ic)
    rng = np.random.default_rng(abs(hash(a.slug)) % (2 ** 31))
    rec_xyz = xyz[:-nlig]
    step = max(1, len(rec_xyz) // 4000)
    rec_idx = np.arange(0, len(rec_xyz), step)
    decoys = [xyz] + make_decoys(xyz, nlig, a.ndecoy, rng, a.max_deg, a.max_trans, rec_idx)
    native_lig = xyz[-nlig:]

    cand = []
    for i, c in enumerate(decoys):
        d = os.path.join(out, f"pose_{i:02d}"); os.makedirs(d, exist_ok=True)
        shutil.copy(prm, os.path.join(d, "com.prmtop"))
        write_rst(os.path.join(d, "com.inpcrd"), c)
        if minimise(d, log) is None:
            continue
        s = score(d, a.rerank, a.rerank_indi, rec, lig, log)
        if s is None:
            continue
        mn = read_rst(os.path.join(d, "min_omm.rst"))
        rms_start = float(np.sqrt(((c[-nlig:] - native_lig) ** 2).sum(1).mean()))
        rms_min = float(np.sqrt(((mn[-nlig:] - native_lig) ** 2).sum(1).mean()))
        cand.append(dict(pose=i, rerank_score=s, rmsd_start=round(rms_start, 2),
                         rmsd_min=round(rms_min, 2)))

    if not cand:
        res["status"] = "no_candidate"; json.dump(res, open(os.path.join(out, "result.json"), "w"))
        print(json.dumps(res)); return 1

    cand.sort(key=lambda c: c["rerank_score"])
    top = cand[0]
    res.update(status="ok", n_candidates=len(cand), selected_pose=top["pose"],
               selected_rmsd=top["rmsd_min"], selected_rerank_score=top["rerank_score"],
               native_rank=1 + sorted(c["rerank_score"] for c in cand).index(
                   next(c["rerank_score"] for c in cand if c["pose"] == 0))
               if any(c["pose"] == 0 for c in cand) else None,
               decoy_rmsd_min=min(c["rmsd_min"] for c in cand),
               decoy_rmsd_max=max(c["rmsd_min"] for c in cand),
               candidates=cand)
    d = os.path.join(out, f"pose_{top['pose']:02d}")
    v = score(d, a.final, a.final_indi, rec, lig, log)
    res["two_step_score"] = v
    res["two_step_components"] = parse_components(
        os.path.join(d, f"FINAL_{a.final}{a.final_indi}.dat"))
    # a control: the same final scoring applied to a candidate picked at random
    ctrl = cand[len(cand) // 2]
    res["random_pose_score"] = score(os.path.join(out, f"pose_{ctrl['pose']:02d}"),
                                     a.final, a.final_indi, rec, lig, log)
    res["random_pose_rmsd"] = ctrl["rmsd_min"]
    json.dump(res, open(os.path.join(out, "result.json"), "w"))
    print(json.dumps({k: v for k, v in res.items() if k != "candidates"}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
