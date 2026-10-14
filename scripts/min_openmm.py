#!/usr/bin/env python3
"""GPU implicit-solvent minimisation as a substitute for the published pmemd.cuda step.

The released protocol minimises with pmemd.cuda; the conda AmberTools build has no
pmemd, and sander needs hours per complex at the protocol's cut=1000.  This runs the
same force field and the same implicit solvent model (igb=2 -> OBC1, no cutoff) on
the GPU through OpenMM, and writes an Amber restart the same MMPBSA step can read.

It is only a valid substitute if it agrees with sander where both were run, which is
why every complex is scored both ways on Dataset I before this is used anywhere else.
"""
import argparse, os, sys, time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--prmtop", default="com.prmtop")
    ap.add_argument("--inpcrd", default="com.inpcrd")
    ap.add_argument("--out", default="min_omm.rst")
    ap.add_argument("--tolerance", type=float, default=0.01)   # kJ/mol/nm
    # 500 L-BFGS iterations reproduces what sander's 5000 SD+CG cycles actually do
    # to the structure: on 1SLD/1VWD/1VWN the all-atom RMSD from the starting
    # coordinates is 0.86/0.79/0.88 A against sander's 0.84/0.79/0.88, and the
    # largest single displacement 3.94/3.90/4.01 against 3.68/4.09/4.09.  Left to
    # converge, OpenMM instead drifts to a median 1.25 A, which is a different
    # protocol -- the published 5000-cycle cap is itself an implicit restraint.
    ap.add_argument("--maxiter", type=int, default=500)         # 0 = until converged
    ap.add_argument("--platform", default="CUDA")
    a = ap.parse_args()

    from openmm import app, unit, Platform, LangevinIntegrator
    import openmm

    d = os.path.abspath(a.outdir)
    prm = app.AmberPrmtopFile(os.path.join(d, a.prmtop))
    crd = app.AmberInpcrdFile(os.path.join(d, a.inpcrd))
    system = prm.createSystem(implicitSolvent=app.OBC1,          # Amber igb=2
                              implicitSolventKappa=0.0 / unit.nanometer,
                              soluteDielectric=1.0, solventDielectric=78.5,
                              nonbondedMethod=app.NoCutoff,      # Amber cut=1000
                              constraints=None, rigidWater=False)
    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               0.002 * unit.picoseconds)
    try:
        plat = Platform.getPlatformByName(a.platform)
        sim = app.Simulation(prm.topology, system, integ, plat)
    except Exception:
        sim = app.Simulation(prm.topology, system, integ)
    sim.context.setPositions(crd.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    t0 = time.time()
    sim.minimizeEnergy(tolerance=a.tolerance * unit.kilojoule_per_mole / unit.nanometer,
                       maxIterations=a.maxiter)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    dt = time.time() - t0
    pos = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

    # Amber ASCII restart: title, natom, then 6F12.7 coordinates
    with open(os.path.join(d, a.out), "w") as fh:
        fh.write("minimised with OpenMM\n")
        fh.write("%5d\n" % len(pos))
        flat = pos.reshape(-1)
        for i in range(0, len(flat), 6):
            fh.write("".join("%12.7f" % v for v in flat[i:i + 6]) + "\n")
    print(f"{os.path.basename(d)} natom={len(pos)} "
          f"E0={e0.value_in_unit(unit.kilocalorie_per_mole):.1f} "
          f"E1={st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole):.1f} "
          f"kcal/mol  {dt:.1f}s  platform={sim.context.getPlatform().getName()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
