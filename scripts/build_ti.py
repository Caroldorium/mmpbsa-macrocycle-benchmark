#!/usr/bin/env python3
"""Turn a prepared TI edge directory into runnable Amber GTI systems and inputs.

Produces two legs -- complex (receptor + both ligands) and solvent (both ligands) -- each
solvated and neutralised, plus a minimise/heat/equilibrate/production input set per lambda
window.  dGG_bind(A->B) = dG_complex - dG_solvent.
"""
import argparse, json, os, subprocess, sys

LAMBDAS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]

LEAP = """source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
source leaprc.phosaa10
loadamberparams lga.frcmod
loadamberparams lgb.frcmod
LGA = loadmol2 lga.mol2
LGB = loadmol2 lgb.mol2
{rec_line}
sys = combine {{ {parts} }}
# no setBox here: solvateoct defines the box itself, and pre-setting one makes the buffer be
# measured from that box rather than from the solute.  tleap 25 and 26 resolved that
# differently -- 7,417 vs 12,013 waters for an identical solute -- so the build was neither
# reproducible nor the requested 12 A.
solvateoct sys TIP3PBOX {buf}
addionsrand sys Na+ 0
addionsrand sys Cl- 0
saveamberparm sys {stem}.prmtop {stem}.inpcrd
savepdb sys {stem}.pdb
quit
"""


def mdin(stage, lam, masks, nstlim, ti_extra=""):
    common = (f"  icfe=1, ifsc=1, clambda={lam:.5f},\n"
              f"  timask1='{masks['timask1']}', timask2='{masks['timask2']}',\n"
              f"  scmask1='{masks['scmask1']}', scmask2='{masks['scmask2']}',\n"
              "  gti_cut=1, gti_output=1, gti_add_sc=25, gti_scale_beta=1,\n"
              "  gti_cut_sc_on=6.0, gti_cut_sc_off=8.0, gti_lam_sch=1,\n"
              "  gti_ele_sc=1, gti_vdw_sc=1, gti_cut_sc=2,\n")
    # bar_l_min/max/incr are deprecated in Amber 22; the sampled lambdas must be listed
    lam_str = ",".join(f"{x:.5f}" for x in LAMBDAS)
    mbar = (f"  ifmbar=1, bar_intervall=5000, mbar_states={len(LAMBDAS)},\n"
            f"  mbar_lambda={lam_str},\n")
    if stage == "min":
        # softcore minimisation must use steepest descent: pmemd refuses ifsc=1 with ntmin!=2
        return ("minimisation\n &cntrl\n  imin=1, maxcyc=10000, ncyc=5000, ntmin=2,\n"
                "  ntb=1, ntp=0, cut=9.0, ntpr=500,\n"
                + common + " /\n")
    if stage == "heat":
        return ("heating 0 -> 298 K\n &cntrl\n  imin=0, ntx=1, irest=0, nstlim=%d, dt=0.001,\n"
                "  ntb=1, ntp=0, ntc=2, ntf=1, cut=9.0,\n"
                "  ntt=3, gamma_ln=5.0, tempi=0.0, temp0=298.15, ig=-1,\n"
                "  ntpr=5000, ntwx=0, ntwr=5000,\n"
                "  ntr=1, restraint_wt=2.0, restraintmask='@CA,C,N',\n"
                % nstlim + common + " /\n")
    if stage == "equil":
        return ("NPT equilibration\n &cntrl\n  imin=0, ntx=5, irest=1, nstlim=%d, dt=0.002,\n"
                "  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,\n"
                "  ntc=2, ntf=1, cut=9.0,\n"
                "  ntt=3, gamma_ln=2.0, temp0=298.15, ig=-1,\n"
                "  ntpr=5000, ntwx=0, ntwr=10000,\n" % nstlim + common + " /\n")
    return ("TI production\n &cntrl\n  imin=0, ntx=5, irest=1, nstlim=%d, dt=0.002,\n"
            "  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,\n"
            "  ntc=2, ntf=1, cut=9.0,\n"
            "  ntt=3, gamma_ln=2.0, temp0=298.15, ig=-1,\n"
            "  ntpr=5000, ntwx=0, ntwr=20000,\n" % nstlim + common + mbar + " /\n")


SBATCH = """#!/bin/bash
#SBATCH --account={account}
#SBATCH --job-name=ti_{tag}
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time={walltime}
#SBATCH --output=%x-%j.out
# amber/22.5's pmemd.cuda is compiled for sm_50-sm_80 and cannot run on the
# H100 nodes ("no kernel image available"); amber-pmemd/24.3 adds sm_86 and sm_90
module --force purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cuda/12.6 amber-pmemd/24.3
set -e
cd "$SLURM_SUBMIT_DIR"
for LEG in complex solvent; do
  for L in {lam_list}; do
    D=$LEG/lam_$L
    mkdir -p "$D"
    [ -f "$D/prod.out.gz" ] && continue
    pmemd.cuda_DPFP -O -i $D/min.in   -p $LEG.prmtop -c $LEG.inpcrd    -ref $LEG.inpcrd \\
               -o $D/min.out  -r $D/min.rst  -inf $D/min.info
    pmemd.cuda_DPFP -O -i $D/heat.in  -p $LEG.prmtop -c $D/min.rst     -ref $D/min.rst \\
               -o $D/heat.out -r $D/heat.rst -inf $D/heat.info
    pmemd.cuda_DPFP -O -i $D/equil.in -p $LEG.prmtop -c $D/heat.rst \\
               -o $D/equil.out -r $D/equil.rst -inf $D/equil.info
    pmemd.cuda_DPFP -O -i $D/prod.in  -p $LEG.prmtop -c $D/equil.rst \\
               -o $D/prod.out -r $D/prod.rst -inf $D/prod.info
    # restarts only chain the stages; every dV/dlambda and MBAR number lives in prod.out.
    # Keeping them, uncompressed, would cost 358 GB over the campaign against 34 GB free.
    rm -f $D/min.rst $D/heat.rst $D/equil.rst $D/prod.rst $D/*.info $D/min.out $D/heat.out
    gzip -f $D/prod.out $D/equil.out
  done
done
echo "TI EDGE DONE {tag}"
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("edgedir")
    ap.add_argument("--buffer", type=float, default=12.0)
    ap.add_argument("--prod-ns", type=float, default=5.0)
    ap.add_argument("--equil-ns", type=float, default=0.5)
    ap.add_argument("--heat-ps", type=float, default=100.0)
    ap.add_argument("--account", default="def-jsbrown_gpu")
    ap.add_argument("--walltime", default="23:00:00")
    ap.add_argument("--no-build", action="store_true", help="write inputs only, skip tleap")
    a = ap.parse_args()
    d = a.edgedir
    info = json.load(open(os.path.join(d, "edge.json")))
    masks = {k: info[k] for k in ("timask1", "timask2", "scmask1", "scmask2")}

    for stem, rec in (("complex", True), ("solvent", False)):
        leap = LEAP.format(rec_line="rec = loadpdb rec.pdb" if rec else "",
                           parts="rec LGA LGB" if rec else "LGA LGB",
                           buf=a.buffer, stem=stem)
        open(os.path.join(d, f"{stem}.leap"), "w").write(leap)
        if not a.no_build:
            r = subprocess.run(f"tleap -f {stem}.leap -I . > leap_{stem}.log 2>&1",
                               shell=True, cwd=d, capture_output=True, text=True)
            pt = os.path.join(d, f"{stem}.prmtop")
            ic = os.path.join(d, f"{stem}.inpcrd")
            ok = (os.path.exists(pt) and os.path.getsize(pt) > 0
                  and os.path.exists(ic) and os.path.getsize(ic) > 0)
            natom = None
            if ok:
                t = open(os.path.join(d, f"{stem}.prmtop"), errors="ignore").read(4000)
                import re
                m = re.search(r"%FLAG POINTERS.*?%FORMAT\(10I8\)\s*\n\s*(\d+)", t, re.S)
                natom = int(m.group(1)) if m else None
            print(f"  {stem}: built={ok} atoms={natom}")
            if not ok:
                lg = os.path.join(d, f"leap_{stem}.log")
                txt = open(lg, errors="ignore").read() if os.path.exists(lg) else (r.stdout or "")
                sys.stderr.write("\n".join(l for l in txt.splitlines()
                                            if "FATAL" in l or "Error" in l)[-1500:])
                for f_ in (pt, ic):
                    if os.path.exists(f_) and os.path.getsize(f_) == 0:
                        os.remove(f_)
                sys.exit(2)

    # steps = time / timestep.  heat runs at dt=0.001 ps, equil and prod at dt=0.002 ps,
    # so 1 ns of production is 500,000 steps -- dividing by 2000 instead of 2 produced runs
    # 1000x too short (5 ps rather than 5 ns), with no MBAR samples at all.
    nsteps = {"heat": int(a.heat_ps / 0.001),
              "equil": int(a.equil_ns * 1000 / 0.002),
              "prod": int(a.prod_ns * 1000 / 0.002)}
    for leg in ("complex", "solvent"):
        for lam in LAMBDAS:
            wd = os.path.join(d, leg, f"lam_{lam:.5f}")
            os.makedirs(wd, exist_ok=True)
            for stage in ("min", "heat", "equil", "prod"):
                open(os.path.join(wd, f"{stage}.in"), "w").write(
                    mdin(stage, lam, masks, nsteps.get(stage, 0)))
    tag = os.path.basename(d.rstrip("/"))
    open(os.path.join(d, "run.sbatch"), "w").write(SBATCH.format(
        account=a.account, tag=tag, walltime=a.walltime,
        lam_list=" ".join(f"{l:.5f}" for l in LAMBDAS)))
    print(json.dumps(dict(edge=tag, windows=len(LAMBDAS), legs=2,
                          prod_ns_per_window=a.prod_ns,
                          total_ns=2 * len(LAMBDAS) * (a.prod_ns + a.equil_ns))))


if __name__ == "__main__":
    main()
