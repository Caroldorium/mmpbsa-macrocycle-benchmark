# Running the single-structure MM/PBSA protocol

Step by step, from one protein-peptide complex to one binding score. Every command here is
what `scripts/zha_cohort.py` actually runs, pulled out so you can do it by hand once and
understand each stage before batching it.

Budget roughly **1-3 minutes per complex** with GPU minimisation. With `sander` instead, budget
hours per complex, for the reason in §6.1.

---

## 1. What you need first

### 1.1 Two conda environments, not one

```bash
conda create -y -n AmberTools26 -c conda-forge ambertools=26
conda create -y -n omm -c conda-forge openmm cudatoolkit
```

They must be separate. AmberTools sets `PYTHONPATH`, `AMBERHOME` and `LD_LIBRARY_PATH`, and if
those leak into the OpenMM process it imports the wrong numpy and dies. Step 5 unsets them
explicitly.

Check the AmberTools side has everything:

```bash
conda activate AmberTools26
which pdb4amber antechamber parmchk2 tleap sander cpptraj MMPBSA.py
```

### 1.2 Two input files per complex

Put them in a directory called `prepared/`:

| file | what it is |
|---|---|
| `native_complex.pdb` | receptor **and** ligand together, in the bound pose |
| `ligand_v2000.sdf` | the ligand alone, V2000 SDF, correct bond orders, **same coordinates as in the complex** |

You also need the ligand's **formal charge** as an integer. Get it wrong and every electrostatic
term is wrong, silently.

The pose in these files is the answer the method will give you. Minimisation will not find a
better one -- see §6.2.

### 1.3 Set up the working variables

```bash
export PREP=/path/to/prepared          # holds native_complex.pdb and ligand_v2000.sdf
export OUT=/path/to/work/my_complex    # output directory, will be created
export CHARGE=0                        # ligand formal charge
export INDI=2                          # interior dielectric; 2 is the primary estimator
mkdir -p "$OUT" && cd "$OUT"
```

---

## 2. Build the system

### Step 1 -- split the complex into receptor and ligand

The receptor is the polymer. Everything heteroatomic is dropped, **except** modified residues
that are covalently inside a chain -- phosphotyrosine and friends. Dropping those opens a gap
that `tleap` then caps with charged termini in the middle of a loop.

```bash
python3 - <<'PY'
import os
CHAIN_HET = {"PTR","SEP","TPO","NEP","HIP","MSE","CSO","CME"}
src, out = os.environ["PREP"] + "/native_complex.pdb", os.environ["OUT"]
keep, dropped, kept = [], set(), set()
for l in open(src, errors="ignore"):
    if l.startswith("ATOM"):   keep.append(l)
    elif l.startswith("TER"):  keep.append(l)
    elif l.startswith("HETATM"):
        rn = l[17:20].strip()
        if rn in CHAIN_HET: keep.append("ATOM  " + l[6:]); kept.add(rn)
        else: dropped.add(rn)
open(out + "/rec_raw.pdb", "w").writelines(keep + ["END\n"])
print("dropped HET:", sorted(dropped), "| kept in chain:", sorted(kept))
PY
```

Read that output. If something you needed is in `dropped`, stop and fix it.

### Step 2 -- clean the receptor

Twice, deliberately. The first pass keeps prepared hydrogens; the second strips them with `-y`
so `tleap` rebuilds them against its own templates.

```bash
pdb4amber -i rec_raw.pdb -o rec_h.pdb
pdb4amber -i rec_h.pdb  -o rec_clean.pdb -y
```

### Step 3 -- parameterise the ligand

```bash
antechamber -i "$PREP/ligand_v2000.sdf" -fi sdf -o lig.mol2 -fo mol2 \
  -c bcc -nc "$CHARGE" -at gaff2 -pf y -dr no \
  -ek "maxcyc=0, qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700"

parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s gaff2
```

**`maxcyc=0` is not optional.** Without it antechamber runs a full gas-phase geometry
optimisation and writes the *optimised* coordinates into the mol2, silently replacing the bound
pose the whole calculation is about. It also takes tens of minutes on a macrocycle instead of
seconds.

Verify the pose survived:

```bash
python3 - <<'PY'
import math, os
sdf = open(os.environ["PREP"] + "/ligand_v2000.sdf").read().splitlines()
n = int(sdf[3][:3])
a = [tuple(float(x) for x in l.split()[:3]) for l in sdf[4:4+n]]
m = open(os.environ["OUT"] + "/lig.mol2").read().splitlines()
i = m.index("@<TRIPOS>ATOM") + 1
b = [tuple(float(x) for x in m[i+k].split()[2:5]) for k in range(n)]
d = max(math.dist(p, q) for p, q in zip(a, b))
print(f"max atom displacement {d:.4f} A", "OK" if d < 0.01 else "<-- POSE WAS REPLACED, STOP")
PY
```

### Step 4 -- build the topologies

```bash
cat > leap.in <<'EOF'
source oldff/leaprc.ff03
source leaprc.gaff2
set default PBRadii mbondi2
rec = loadpdb rec_clean.pdb
loadamberparams lig.frcmod
lig = loadmol2 lig.mol2
saveamberparm rec rec.prmtop rec.inpcrd
saveamberparm lig lig.prmtop lig.inpcrd
com = combine { rec lig }
saveamberparm com com.prmtop com.inpcrd
savepdb com com.pdb
quit
EOF

tleap -f leap.in
ls com.prmtop rec.prmtop lig.prmtop com.inpcrd      # all four must exist
```

If you kept any `CHAIN_HET` residues in step 1, add `source leaprc.phosaa10` after the gaff2
line. For a different protein force field swap the first line: `leaprc.protein.ff14SB`,
`leaprc.protein.ff14SBonlysc`, `oldff/leaprc.ff99`, `oldff/leaprc.ff02`.

---

## 3. Minimise

One optimisation to a local minimum. No dynamics, no temperature, no time -- see §6.2.

### Step 5a -- GPU, via OpenMM (recommended)

```bash
env -u PYTHONPATH -u PYTHONHOME -u AMBERHOME -u LD_LIBRARY_PATH \
  conda run -n omm python3 /path/to/repo/scripts/min_openmm.py . \
    --prmtop com.prmtop --inpcrd com.inpcrd --out min_omm.rst \
    --maxiter 500 --platform CUDA

cat > traj.cpptraj <<'EOF'
parm com.prmtop
trajin min_omm.rst
trajout min.nc netcdf
run
quit
EOF
cpptraj -i traj.cpptraj
```

`--maxiter 500` is the calibration, not a convenience. It reproduces what sander's 5000 SD+CG
cycles actually do to the structure (0.86/0.79/0.88 A all-atom RMSD against sander's
0.84/0.79/0.88 on three test complexes). Left to converge, OpenMM drifts to a median 1.25 A,
which is a different protocol. See §6.3 for why this matters more than it sounds.

Use `--platform CPU` if you have no GPU; it still beats sander because L-BFGS beats 5000 cycles
of SD+CG.

### Step 5b -- CPU, via sander (the published path)

```bash
cat > min.in <<'EOF'
Implicit-solvent minimisation, Zhao et al. 2025
 &cntrl
    imin=1, ntx=1, irest=0,
    maxcyc=5000, ncyc=2000,
    cut=1000.0, ntpr=5000, ntwx=0,
    igb=2, ntxo=1,
 /
EOF

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
sander -O -i min.in -p com.prmtop -c com.inpcrd -o min.out -r min.rst -inf min.info

cat > traj.cpptraj <<'EOF'
parm com.prmtop
trajin min.rst
trajout min.nc netcdf
run
quit
EOF
cpptraj -i traj.cpptraj
```

Pin the thread count. sander links a threaded BLAS and will quietly take ~1.4 cores per process,
which wrecks any parallel batch.

---

## 4. Score

Single-trajectory: the receptor and ligand are split out of **the same minimised coordinates**,
never re-minimised separately. That is why all three topologies are passed against one `min.nc`.

### Step 6a -- Poisson-Boltzmann (the primary estimator)

```bash
mkdir -p "pb_$INDI" && cd "pb_$INDI"
cat > mmpbsa.in <<EOF
Zha25 replication
&general
  interval=1, verbose=1,
/
&pb
  indi=$INDI, exdi=80.0, inp=1,
  cavity_surften=0.0072, cavity_offset=0.0,
/
EOF

MMPBSA.py -O -i mmpbsa.in -cp ../com.prmtop -rp ../rec.prmtop -lp ../lig.prmtop \
  -y ../min.nc -o FINAL.dat -prefix _MM_

rm -f _MM_*
grep -A3 "DELTA TOTAL" FINAL.dat
cd ..
```

### Step 6b -- generalised Born, if you want it

GB needs a workaround. `MMPBSA.py` **silently ignores** `intdiel` passed through the namelist,
so every "different dielectric" arm returns identical numbers unless you inject it by hand.

```bash
mkdir -p "gb2_$INDI" && cd "gb2_$INDI"
cat > mmpbsa.in <<'EOF'
Zha25 replication
&general
  interval=1, verbose=1, use_sander=1,
/
&gb
  igb=2, surften=0.0072,
/
EOF

BASE="MMPBSA.py -O -i mmpbsa.in -cp ../com.prmtop -rp ../rec.prmtop -lp ../lig.prmtop \
  -y ../min.nc -o FINAL.dat -prefix _MM_"

$BASE -make-mdins
sed -i.bak "s/igb=2,/igb=2, intdiel=$INDI.0,/" _MM_gb.mdin   # -i.bak works on BSD and GNU sed
$BASE -use-mdins

rm -f _MM_*
grep -A3 "DELTA TOTAL" FINAL.dat
cd ..
```

If you skip the `-make-mdins` / `sed` / `-use-mdins` dance, your dielectric sweep is fake.

### Step 7 -- read the number

```bash
grep -E "DELTA TOTAL|VDWAALS|EEL|EPB|ENPOLAR" "pb_$INDI/FINAL.dat"
```

`DELTA TOTAL` is the binding score, in kcal/mol. `MMPBSA.py` will say
`Calculations performed using 1.0 complex frames` with zero standard error. That is correct and
expected -- one structure, scored once.

**It is not a binding free energy.** No entropy term is included, so it sits 45-75 kcal/mol from
experiment and the per-target slope is 0.10-0.21. Use it to rank compounds within one series
against one receptor. Nothing else.

---

## 5. Doing it for real

### 5.1 One complex, all stages, one command

```bash
conda activate AmberTools26
python scripts/zha_cohort.py MY_SLUG "$PREP" "$OUT" \
  --charge 0 --ff ff03 --indi 1,2,4 --methods pb,gb2 --min-engine openmm
```

Writes `$OUT/result.json` with every estimator, every energy component, the pose drift, the
dropped heteroatoms and a status field.

### 5.2 A whole panel in parallel

Build a `cohort_manifest.csv` with one row per complex
(`slug,target,structure_id,heavy_atoms,charge,dG,fep`) and run:

```bash
NPROC=10 MIN_ENGINE=openmm OUT=/path/to/work/cohort bash scripts/run_cohort.sh
```

It skips anything already marked `"status": "ok"`, so it is safe to re-run after a crash.

### 5.3 Quality control -- do not skip this

```bash
python3 - <<'PY'
import json, glob
bad = []
for f in glob.glob("work/cohort/*/result.json"):
    r = json.load(open(f))
    if r.get("status") != "ok": bad.append((f, "failed")); continue
    v = r.get("pb_2_VDWAALS")
    if v is not None and v > 0: bad.append((f, f"positive VDW {v:.0f} -- unrelieved clash"))
    d = r.get("pose_drift_A")
    if d is not None and d > 0.01: bad.append((f, f"pose drift {d} A"))
for f, why in bad: print(why, f)
print(f"\n{len(bad)} complexes need attention")
PY
```

A **positive van der Waals term** means atoms are still overlapping after minimisation -- the
build is broken, not the compound. We found one at +277,293 kcal/mol. Drop them; do not score
them. Also check your pose builder's clash count: if the score correlates with it, you are
measuring construction quality, not chemistry (see §6.2).

---

## 6. Things that will bite you

### 6.1 Why minimise in OpenMM rather than Amber

"Just do it in Amber" means `sander`, because most conda AmberTools builds ship no `pmemd.cuda`.
The protocol specifies `cut=1000.0`, which is effectively no nonbonded cutoff, making the
nonbonded term O(N^2) on a single CPU core: **hours per complex**. OpenMM with `OBC1` and
`NoCutoff` is the direct equivalent of `igb=2` with `cut=1000`, on the GPU.

It was checked, not assumed. On the 17 complexes run both ways, the primary estimator `pb_2`
gives a correlation with pKd of **-0.587 under both engines**, with score agreement r = +0.995
and a mean absolute difference of 1.44 kcal/mol. No estimator shifts by more than 0.034.

Reproduce it with `scripts/compare_minimisers.py`.

### 6.2 Minimisation does not search

It walks downhill to the **nearest local minimum** and cannot cross a barrier. The pose you put
in is the pose you get out, relieved of clashes. Its job is to make the structure scoreable, not
to find the binding mode.

So the score inherits the quality of the input pose, completely. On a sealed external panel
whose compounds were threaded onto a non-matching template, the score correlated with the pose
builder's clash count at **rho = -0.891**. It was reading how badly the pose was built. Report
pose provenance alongside any score.

### 6.3 Minimisation depth is a real free parameter

On KRAS, removing the 500-iteration cap moves the correlation from **+0.225 to +0.529**, and the
spread introduced by the stopping point alone (SD 3.60 kcal/mol) is as large as the entire
between-compound signal (SD 3.53).

It does not generalise -- uncapping improves KRAS by +0.304 and degrades PCSK9 by -0.280, net
macro change -0.009. Pin the cap to sander's drift and leave it alone. Tuning it to whichever
setting flatters your target is how you get a result that does not replicate.

### 6.4 Longer sampling does not help

Four independent 10 ns replicas per compound, 324 frames each: replicate-to-replicate SD
**6.81 kcal/mol** against a between-compound SD of **5.74**. Which replicas you happened to run
moves a compound further than real chemistry does. Ensemble averaging buys noise, expensively.
That is why this is a single-structure protocol.

### 6.5 Evaluate at matched size

Experimental affinity tracks heavy-atom count at -0.43 and the score tracks it at -0.55, so an
unmatched benchmark partly rewards reading molecular size. A heavy-atom count alone reaches
Spearman 0.487 and 60.2% pairwise accuracy. **That is the bar, not chance.**

Compare compounds within a narrow heavy-atom window, average within receptor group rather than
pooling across crystal forms, and bootstrap over **compounds, not pairs** -- pairs drawn from a
compound set are not independent, and resampling pairs understates the interval by roughly half.
`scripts/size_matched_pairwise_accuracy.py` does all three.

### 6.6 Quick reference

| symptom | cause | fix |
|---|---|---|
| ligand coordinates moved | antechamber optimised the geometry | add `-ek "maxcyc=0, ..."` |
| every dielectric gives the same number | `intdiel` silently ignored | `-make-mdins`, inject, `-use-mdins` |
| OpenMM import error / wrong numpy | AmberTools env vars leaked | `env -u PYTHONPATH -u AMBERHOME ...` |
| `tleap` caps a chain mid-loop | modified residue dropped as HETATM | keep it via `CHAIN_HET`, `source leaprc.phosaa10` |
| huge positive binding score | unrelieved clash | screen on positive VDW, drop the complex |
| parallel batch much slower than expected | sander's threaded BLAS | `export OMP_NUM_THREADS=1` |
| scores track your pose builder | minimisation inherits the input pose | fix the poses; check clash correlation |

---

## 7. Where the numbers in this file come from

| claim | source |
|---|---|
| sander vs OpenMM agreement | `scripts/compare_minimisers.py`, Dataset I, n = 17 |
| 500-iteration calibration | `scripts/min_openmm.py` docstring; methods §2.4 |
| KRAS +0.225 to +0.529 | methods §2.4; `results/capped_vs_converged.csv` |
| replica noise 6.81 vs 5.74 | `kras_all40_md_score_and_stability_metrics.csv` |
| size confound -0.43 / -0.55 | `scripts/fig_size_confound.py` |
| heavy-atom null 0.487 / 60.2% | `results/zha25_cohort_metrics.csv`; methods §5 |
| clash correlation -0.891 | `results/zha25_12ca5_built_pose_scores.csv` |

Full parameter justifications are in `methods/METHODS.md`.
