# Establishing an end-point scoring protocol for macrocyclic peptide design: size-matched benchmarking, its limits, and a registered prospective test

*Alternative title, results-forward: "Size-matched benchmarking reveals genuine but modest
ranking skill in end-point MM/PBSA for macrocyclic peptides." The methods-forward title above
reflects the paper's purpose; either works.*

**Carlos Liu¹, Joseph S. Brown¹**

¹ Leslie Dan Faculty of Pharmacy, University of Toronto, Toronto, Ontario, Canada

*Draft manuscript — prospective validation section is a registered plan, not yet executed.*

---

## Abstract

Macrocyclic peptides engage targets that small molecules cannot, but each analogue is
expensive to make, so computational ranking has unusual leverage: a method that merely sorts
candidates correctly saves real synthetic effort. Rigorous alchemical free-energy methods
achieve this but cost GPU-days per compound. End-point methods such as MM/PBSA cost seconds,
and a recent protocol reports strong performance on protein–cyclic-peptide complexes. Whether
that performance reflects physics or molecular size has not been established.

Size is not a nuisance variable for peptide binding: larger peptides bury more interface and
genuinely bind more tightly. On a five-target benchmark of 234 compounds with measured
affinities, experimental ΔG tracks ligand heavy-atom count with a mean within-target Spearman
ρ of −0.43, and the MM/PBSA score tracks it more strongly still (ρ = −0.55). Conventional
benchmarking therefore cannot separate a scoring function from a ruler.

We resolve this by evaluating only within sets of peptides matched on heavy-atom count, where a
size descriptor carries little or no information by construction, and asking how often a method
picks the tighter binder of a decisive pair. Within a three-atom window MM/PBSA succeeds in
**62.4% of pairs (95% CI 55.5–68.9)**, against 50% for chance and **83.9% (78.2–88.6)** for
published FEP+ predictions on the same compounds. At exactly equal heavy-atom count, where a
size descriptor is wholly uninformative, the estimate is 60.1% but the interval widens to
47.2–72.0 (posterior mass above chance 0.93) because few same-size decisive pairs exist.
Intervals resample compounds rather than pairs, since pairs drawn from a compound set are not
independent; resampling pairs understates the width by roughly half.

The skill is modest and unevenly distributed: only PCSK9 (75.7%) is individually above chance,
with Cyclin A/B, JAK2, KRAS and MCL-1 each spanning it.

Three further results bear on how such benchmarks should be read. The apparent accuracy depends
on the estimator: for KRAS, whose compounds span two receptor crystal forms, pooling gives
ρ = +0.225 while averaging within forms gives −0.130, the pooled value being credit for
separating two groups offset by 1.15 kcal/mol. The score is a ranking device rather than a free
energy, sitting 45–75 kcal/mol from experiment with per-target regression slopes of 0.10–0.21.
And on a sixth target held entirely outside development — 17 peptidomimetics against an
antibody, affinities sealed before prediction — every method fell below chance and the
size–affinity relationship reversed sign, showing that performance on five related targets is
not sufficient evidence of transfer.

This method was established to serve an active design campaign: 212 macrocyclic peptides
designed against the lipopolysaccharide transport protein LptA, from which candidates will be
selected for synthesis and experimental evaluation. We report the protocol applied to that
panel, and register the prospective test with its power analysis stated in advance. The purpose
of the validation reported here is to establish what the method can and cannot support
**before** committing synthetic resources, not to report a finished prediction.

**Keywords:** macrocyclic peptides, MM/PBSA, binding affinity, benchmarking, ligand size bias,
free energy perturbation

---

## 1. Introduction

### 1.1 Why ranking matters more than prediction here

Macrocyclic and cyclic peptides occupy a chemical space between small molecules and biologics,
and they can engage shallow protein–protein interfaces that conventional drugs cannot.[^1] They
are also expensive. A single non-canonical building block can cost on the order of a thousand
dollars, synthesis is laborious, and a design campaign can propose hundreds of candidates
against a budget of ten. The practical question is therefore not "what is the binding
constant?" but "of these candidates, which should we make?"

That is a ranking problem, and ranking is a weaker requirement than prediction. A method that
cannot produce a binding constant may still order candidates usefully. Conversely, a method
with an impressive correlation coefficient may be useless if that correlation comes from
something the chemist already knows.

**This work has a specific purpose.** We maintain an active design campaign against the
lipopolysaccharide transport protein LptA, comprising 212 designed macrocycles from which a
subset will be synthesised and assayed. Non-canonical building blocks cost on the order of a
thousand dollars each, so the selection is expensive and effectively irreversible. We therefore
set out to establish, in advance of that commitment, exactly what an end-point scoring protocol
can and cannot support on this chemistry — which is a methods question, and is the subject of
this paper. The design campaign is described in Section 4 and the registered prospective test
in Section 5; the validation that precedes them is what licenses either.

### 1.2 The computational options

Free energy perturbation (FEP) is the rigorous approach: it simulates one ligand transforming
into another and tracks the free-energy change. Applied to macrocycles with appropriate
enhanced sampling, FEP+ achieves a global pairwise RMSE of 1.06 kcal/mol across more than 230
analogues spanning over 10 kcal/mol of experimental range.[^2] The cost is GPU-days per
compound and a commercial licence.

End-point methods take the opposite trade. MM/PBSA evaluates the interaction energy of a
single structure and subtracts an estimate of the desolvation cost, in seconds rather than
days.[^3] Zhao and colleagues recently reported a fine-tuned MM/PBSA(GBSA) workflow for
protein–cyclic-peptide complexes achieving a Pearson correlation of −0.732 against experimental
affinity on a 50-complex dataset, roughly twice that of a conventional docking score, at about
three seconds per prediction.[^4]

The gap between these is roughly three orders of magnitude in cost. What is given up is not
well characterised for peptides.

### 1.3 The problem with the obvious benchmark

Scoring functions have a long-documented tendency to track molecular size. Jacobsson and Karlén
found "a high to intermediate linear correlation between the score and the number of heavy
atoms" for nine of ten scoring functions tested, and identified this ligand bias as a source of
false positives in virtual screening.[^5] The standard responses are normalisation by molecular
weight or heavy-atom count, ligand efficiency metrics, or matched-pair benchmark design.[^5],[^6]

For peptides the problem is sharper, for a reason that is easy to miss. Size is not an artefact
here. More residues bury more interface area, and buried area is a genuine driver of binding
affinity. A size descriptor is therefore a *mediator* of the physics, not a confounder to be
regressed away — and partialling it out over-corrects, removing exactly the signal a physics
method is supposed to capture.

This creates a genuine methodological difficulty. Comparing a scoring function against a
heavy-atom count across compounds of different sizes is unfair in both directions: the scoring
function gets credit for re-reading size, and the control gets credit for real physics.

### 1.4 What we do

We hold size fixed. Within each target we compare only peptides of **identical heavy-atom
count**, where a size descriptor has no information at all because there is no size difference
to read, and ask the question a chemist actually asks: given two candidates, which binds
tighter? This extends matched-pair benchmarking[^6] to macrocyclic peptides, where the size
range is far wider than in small-molecule series and the synthesis stakes are higher.

We report: a replication of the Zhao protocol on an independent five-target cohort; the
magnitude of the size confound; size-matched performance against a published FEP+ ceiling on
identical compounds; the sensitivity of all of this to the choice of statistical estimator; and
an application to a de novo design campaign against LptA, with a registered prospective test.

> **Box 1 — What a force field is, for readers new to this**
>
> A force field is a set of rules for computing the energy of an arrangement of atoms. It
> assigns each atom a partial electric charge and a size, and each bond a stiffness and a
> preferred length. Given the 3-D coordinates of a protein and a peptide, it returns a number:
> how favourable that arrangement is. It is classical physics — balls and springs — with
> parameters fitted to quantum calculations and experiment. It knows nothing about biology.
>
> Standard amino acids have well-tested parameters. Non-canonical residues do not, and must be
> parameterised individually before anything can be computed (Section 2.3).

---

## 2. Methods

### 2.1 Benchmark cohort

Method development and all reported statistics use a published external benchmark: the
five-target macrocycle and cyclic-peptide series of Awoonor-Williams and colleagues — KRAS,
PCSK9, MCL-1, JAK2 and Cyclin A/B — comprising over 230 peptidic and non-peptidic
analogues.[^2] We use the 241 compounds for which structures and labels could be assembled; the
attrition to the analysed sets is given below.

Two properties make this the right cohort. The affinities were measured by a third party, so
the evaluation cannot be tuned to them. And because FEP+ predictions were published alongside,
a rigorous-method ceiling is available on exactly the same compounds, which is otherwise
difficult to obtain.

**Compound attrition, stated once for reproducibility.** Of the 241 compounds for which
structures and labels could be assembled, **235** carry an experimental binding free energy;
**234** remain after excluding `MCL_1__cpd_55` for an unrelieved steric clash (Section 2.2); and
**220** remain for the pairwise analyses, which additionally exclude censored measurements
because the ordering of a censored pair is ambiguous. Correlation analyses in §3.1–3.3 use the
234; pairwise analyses in §3.4 use the 220. Every reported *n* refers to one of these two sets.

Compounds are not independent within a target: two targets span more than one receptor crystal
form (KRAS, two; PCSK9, three), giving eight receptor groups. All primary statistics are
computed **within receptor group** and averaged with equal weight per target (Section 2.6).

### 2.2 Structure preparation

Receptor and peptide are cleaned with `pdb4amber` in two passes. The first retains prepared
hydrogens so that histidine tautomer and cystine assignments carry into residue names; the
second strips hydrogens so that `tleap` rebuilds them against its own templates. A single pass
is insufficient because prepared files carry atoms (e.g. `HXT`) for which the older force
fields have no type.

Cyclisation is imposed explicitly. Head-to-tail macrocycles receive an explicit bond in
`tleap`; disulfide-cyclised peptides are detected geometrically (S–S < 2.5 Å) and left to the
CONECT records written by `pdb4amber`, because re-issuing those bonds manually produces
`1-4: cannot add bond` errors. For macrocycles with uncharged termini, `clearPdbResMap` prevents
`tleap` applying terminal templates to a molecule that has no termini.

Two preparation defects were found during validation, both of which produced plausible but
wrong energies, and both of which are corrected here:

**Spurious chain breaks.** `pdb4amber` inserted `TER` records into continuous chains, each of
which `tleap` capped, producing complexes with initial energies above 10¹⁷ kcal/mol. The problem
affects **14 of 50 reference complexes, with 9 to 25 spurious records each** (the remaining 36
are unaffected). `TER` records whose flanking C–N distance is below 2.0 Å are deleted. After the
fix, atom counts of a reference complex match published values exactly.

**Unrelieved steric clash.** One compound (`MCL_1__cpd_55`) passed the pipeline's own `ok`
status while carrying a van der Waals term of **+277,293 kcal/mol**, against ≤ −50 kcal/mol for
every other compound. A status check does not test for unrelieved overlap. Screening on
`VDWAALS > 0` is mandatory before any analysis and is applied throughout. Rank statistics are
largely insensitive to this compound (macro Spearman 0.532 → 0.545), but the pooled Pearson
correlation moves from −0.048 to **+0.265** on removing it, and MCL-1's Pearson *r* from −0.124
to +0.693. A single unrelieved clash was therefore enough to make the method look uncorrelated
with experiment.

### 2.3 Force fields and non-canonical residues

| component | assignment | reason |
|---|---|---|
| protein | `oldff/leaprc.ff03` | the leaprc named by the source protocol; it takes terminal residues from `all_amino{n,c}t94.lib` and therefore carries NCYX/CCYX. `leaprc.protein.ff03.r1` has neither and fails on 20 of 50 reference complexes whose terminal residue is a cystine |
| protein (comparison) | ff14SB,[^7] ff19SB[^8] | tested rather than assumed (Section 3.6) |
| phosphotyrosine | `leaprc.phosaa10`[^9] | JAK2's receptor carries a phosphotyrosine; without it `tleap` fails |
| ligands, ncAA side chains | GAFF2[^10] with AM1-BCC charges[^11] | standard for non-protein chemistry |
| water (alchemical legs) | TIP3P[^12] | — |

ff19SB is used off-label under implicit solvent: its CMAP correction was parameterised for OPC
water and applies only to the 20 standard residues, so a GAFF2-hybrid non-canonical residue
never receives it. This is reported rather than hidden.

#### 2.3.1 Non-canonical residue parameterisation

Designed peptides carrying non-standard residues cannot be scored until each residue has Amber
parameters. The LptA application required **170 non-standard residues** (165 unique
non-canonical codes plus five residues constituting the thanatin derivative-7 scaffold, none of
which `tleap` knows). Because this step gates everything downstream and is where most silent
failures occur, we specify it in full.

**Protocol, per residue.**

1. Take the polymer-form rotamer and the corresponding chemistry definition from the source
   library. Atom names and bond orders are read from the definition file, not perceived from
   coordinates; coordinates are matched positionally.
2. Cap with ACE at the N-terminus and NME at the C-terminus, giving a neutral dipeptide. Cap
   atoms are placed by bisector construction with a hard clash check, and cap atom names are
   chosen per residue to avoid collisions with side-chain names.
3. Derive AM1-BCC charges on the capped dipeptide with `antechamber`, under element-safe
   placeholder atom names (see defect 1 below).
4. Restore the true atom names, case-insensitively.
5. Reassign the backbone to protein atom types (N, H, CA, HA, C, O) while the side chain
   remains GAFF2. Cap atoms are retyped as protein types as well, because they stand in for
   neighbouring residues and leaving them GAFF2-typed makes `parmchk2` generate cross terms
   against the wrong types.
6. Run `prepgen` with a mainchain file (`HEAD_NAME N`, `TAIL_NAME C`, `OMIT_NAME` for the caps,
   `CHARGE 0.0`) to remove the caps and redistribute their charge, giving a residue that nets
   zero.
7. Run `parmchk2` for missing GAFF2 terms, then strip from the resulting frcmod any entry that
   would redefine a protein backbone type.

Capping matters because AM1-BCC on a free amino acid fits backbone charges in an environment
(–NH₂ / –COOH) the residue never occupies; capping places the charge fit in a peptide-like
context, which is the standard Amber construction for non-standard residues.[^13] Backbone
retyping matters because `antechamber` otherwise types the backbone as GAFF2 and `tleap` then
has no parameters across the peptide bond to a neighbouring standard residue
(`Could not find bond parameter for atom types: C - ns`). The alpha carbon is typed **CT** under
ff03/ff99 and **CX** under ff14SB/ff19SB — ff03 has no CX at all — so two complete parameter
sets are maintained and a prep typed for one force field fails on the other.

**Validation gate.** Each residue is tested by constructing ALA–X–ALA in `tleap` and requiring
**zero missing parameters**. All 170 pass under both backbone typings. This gate is not
optional: `tleap` writes a topology even when parameters are missing, substituting defaults
silently, so "the build succeeded" is not evidence of correctness.

**Failure modes.** Six defects were found, four of which produced plausible output rather than
an error. We tabulate them because every one is a property of the standard toolchain rather
than of our library, and will recur for anyone attempting this.

**Table 1.** Silent and overt failure modes in non-canonical residue parameterisation.

| # | defect | how it manifests | detection | remedy |
|---:|---|---|---|---|
| 1 | `antechamber` infers element from atom name | every alpha carbon (`CA`) read as **calcium**; 61 aromatic carbons named `CL*` read as **chlorine**; `FE`, `NI`, `CO` likewise | element census against the chemistry definition | build under element-safe placeholder names, restore afterwards |
| 2 | bond orders perceived from geometry | fails on hypervalent sulfur (sulfonamides) and iodoaryls | valence check | supply bond orders from the definition via MOL2 |
| 3 | fixed-dihedral cap placement | NME nitrogen placed 0.15 Å from the residue's own carbonyl oxygen whenever ψ ≈ 180° | inter-atomic distance check after capping | bisector placement with a clash threshold |
| 4 | cap atom names collide with side-chain names | **75 of 170** residues affected | name-set intersection | per-residue cap naming |
| 5 | library-internal name inconsistency | nine halogenated residues: the rotamer PDB truncates `CLM1` to `CLM`, and one residue names *both* chlorines `CLM`. Unmatched atoms are dropped silently by `prepgen`, surfacing later as a `nan` internal coordinate | cross-check names between definition and coordinates | take names and bonds from the definition, match coordinates positionally |
| 6 | case-sensitive name restoration | `antechamber` rewrites a two-letter placeholder to proper element case (`CL1` → `Cl1`), so a case-sensitive mapping leaves halogens unnamed | case-insensitive round-trip test | case-insensitive restoration |

**A limitation we could not remove.** GAFF2 has no selenium and `sqm` has no AM1 parameters for
it, so `antechamber` retypes selenium as sulfur and `parmchk2` reports nothing missing — the
build appears clean. Two residues are affected and are flagged `ok_element_unsupported` rather
than shipped as validated. Any score involving them is an approximation, and we say so wherever
such a residue appears.

### 2.4 Energy minimisation

Minimisation follows the published settings — `imin=1, maxcyc=5000, ncyc=2000, cut=1000.0,
igb=2` — corresponding to 2000 steps of steepest descent then conjugate gradient in the
Onufriev–Bashford–Case generalised Born model[^14] with no nonbonded cutoff.

Because the available AmberTools build does not ship `pmemd.cuda`, and `sander` requires hours
per complex at `cut=1000`, minimisation is performed in OpenMM[^15] using `OBC1` and
`NoCutoff`, the direct equivalents. The two engines are not identical, so the OpenMM L-BFGS
iteration count was **calibrated against sander by structural drift**: sander's 5000 cycles move
the complex 0.84 Å RMSD, an unrestricted OpenMM minimisation 1.25 Å. A cap of 500 iterations
reproduces sander's drift and is used throughout.

This calibration is load-bearing. On KRAS, removing the cap changes the correlation from +0.225
to +0.529, and the spread introduced by the stopping point (SD 3.60 kcal/mol) is as large as
the entire between-compound signal (SD 3.53) — a ratio of 1.02. Removing the cap does not
generalise: it improves KRAS by +0.304 and degrades PCSK9 by −0.280, for a net macro change of
−0.009 (cluster bootstrap +0.110, 95% CI [−0.066, +0.230]).

### 2.5 Free-energy evaluation

`MMPBSA.py`[^16] is run single-trajectory, so receptor internal energy cancels exactly.
Poisson–Boltzmann calculations use the internal solver with interior dielectric 2.0, exterior
80.0, and the LCPO nonpolar model (`inp=1`, `cavity_surften=0.0072`, `cavity_offset=0`).
Arms at εᵢₙ = 1 and 4 and GB arms are computed in parallel for comparison. **No entropy term is
included**, following the source protocol.

The `intdiel` keyword is silently ignored by both engines when passed through `&pb`. The
dielectric is therefore set by writing `use_sander=1`, generating inputs with `-make-mdins`,
injecting `intdiel` into the generated `_MM_gb.mdin`, and re-running with `-use-mdins`. Without
this, every "different dielectric" arm silently returns identical numbers.

The protocol is **single-structure**: one minimised geometry is scored, not an ensemble.
MMPBSA.py reports `Calculations performed using 1.0 complex frames` with zero standard error.
This is the source protocol's design, stated explicitly because "MM/PBSA" usually implies
ensemble averaging. The choice is defensible on published evidence that MM/PBSA correlation
does not generally improve with longer sampling because the limiting error is bias rather than
variance;[^3] we tested it directly (Section 3.7).

> **Box 2 — What MM/PBSA actually computes**
>
> Take the structure of the peptide sitting in the protein's binding site. Add up the
> electrostatic and van der Waals attraction between them: that is the gas-phase interaction
> energy, and it is always favourable. But in water, both surfaces are coated with ordered water
> molecules that must be stripped off before they can touch, and that costs energy. MM/PBSA
> estimates this desolvation penalty by solving the Poisson–Boltzmann equation, treating water
> as a continuous dielectric rather than individual molecules.
>
> The binding score is attraction minus desolvation. What it omits is **entropy** — the freedom
> the peptide loses when it stops tumbling in solution and locks into one conformation. That
> omission is why the numbers are not free energies (Section 3.4).

### 2.6 Statistical evaluation

**Primary estimator.** Spearman's ρ computed within receptor group and averaged across targets,
each target weighted equally. Pooling compounds across targets is rejected: it rewards a method
for knowing that one series is systematically more potent than another, which is not the
question a design campaign asks. Section 3.2 quantifies the resulting inflation within our own
data.

**Size-matched pairwise accuracy — the primary measure.** For each target we enumerate all
compound pairs whose measured affinities differ by at least a margin *m* ("decisive" pairs) and
whose ligand heavy-atom counts differ by at most a window *w*. For each pair we ask whether the
method orders them correctly. Chance is 0.5. At *w* = 0 a size descriptor is undefined, because
there is no size difference to read.

Two free parameters enter here, and neither was chosen after seeing the result.

*The affinity margin m.* Pairs whose measured affinities are nearly equal cannot discriminate
between methods and only add noise, so a margin is required. We set *m* = 0.5 kcal/mol a priori
because it is approximately the reproducibility of the assays contributing to the cohort, and
report the full sensitivity in Table 2 rather than a single value. The conclusion is stable:
MM/PBSA is above chance at every margin from 0 to 1.5 kcal/mol.

*The size window w.* We report *w* = 0, the strictest possible matching, and *w* = 3, which
admits enough pairs for the macro interval to be informative. Both are reported throughout; we
do not select between them on the basis of which is more favourable.

**Table 2.** Sensitivity to the decisive-pair margin, at *w* = 3. "Edge" is MM/PBSA minus the
heavy-atom control on the identical pair set. Intervals resample compounds.

| margin (kcal/mol) | pairs | MM/PBSA | 95% CI | heavy atoms | edge | FEP+ |
|---:|---:|---:|---|---:|---:|---:|
| 0.00 | 2384 | 0.597 | [0.542, 0.649] | 0.585 | +0.012 | 0.769 |
| 0.25 | 2002 | 0.613 | [0.554, 0.674] | 0.591 | +0.022 | 0.805 |
| **0.50** | **1645** | **0.624** | **[0.553, 0.691]** | **0.602** | **+0.022** | **0.839** |
| 0.75 | 1284 | 0.640 | [0.566, 0.718] | 0.622 | +0.018 | 0.888 |
| 1.00 | 1016 | 0.662 | [0.562, 0.750] | 0.590 | +0.072 | 0.921 |
| 1.50 | 666 | 0.697 | [0.576, 0.804] | 0.599 | +0.098 | 0.932 |
| 2.00 | 452 | 0.620 | [0.468, 0.775] | 0.420 | +0.200 | 0.956 |

Two features of Table 2 deserve comment, and one of them is unfavourable. Accuracy rises
monotonically with the margin up to 1.5 kcal/mol for both MM/PBSA and FEP+, which is the
behaviour of a real but noisy signal: pairs that are further apart in affinity are easier to
order. **But at the chosen margin the edge over a heavy-atom count is only +0.022.** It becomes
substantial (+0.07 to +0.20) only for pairs separated by 1 kcal/mol or more. The method's
advantage over the trivial descriptor within a three-atom window is therefore real but small,
and concentrated in the pairs that are easiest to tell apart. At *w* = 0, where the descriptor
has no information at all, the comparison does not arise.

**Uncertainty.** All pairwise intervals resample **compounds**, not pairs. This matters: the
*n*(*n*−1)/2 pairs formed from *n* compounds are not independent draws, because every compound
appears in many pairs. Resampling pairs treats them as independent and badly understates the
uncertainty — on the exact-match band it returns [0.530, 0.693] where the correct
compound-level interval is [0.472, 0.720]. We resample compounds within each target (2,000
draws), rebuild the pair set from each resample, and combine targets for macro values (5,000
draws). Where Spearman correlations are reported the bootstrap resamples **receptor groups**.

**Controls.** Ligand heavy-atom count is computed on the same compounds with the same
estimator throughout. Censored measurements are retained and flagged in correlation analyses
and excluded from pairwise analyses where the ordering would be ambiguous.

> **Box 3 — How to read a pairwise accuracy**
>
> Hand the method two peptides whose measured affinities genuinely differ, and ask which binds
> tighter. Do that for every such pair. The fraction it gets right is the pairwise accuracy.
>
> A coin gets 50%. A perfect method gets 100%. 60% means that in a head-to-head choice between
> two candidates, the method is right three times in five instead of one in two. That is a
> modest edge, but across a campaign of many such choices it is not nothing — and crucially it
> is measured here in a setting where the obvious shortcut, "pick the bigger one", is
> unavailable.

### 2.7 Consolidated parameters

Every setting required to reproduce the reported numbers, in one place.

| stage | parameter | value |
|---|---|---|
| preparation | `pdb4amber` passes | two; first retains prepared hydrogens, second strips with `-y` |
| | spurious `TER` removal | deleted where flanking C–N < 2.0 Å |
| | disulfide detection | geometric, S–S < 2.5 Å |
| | quality screen | exclude any complex with VDWAALS > 0 |
| force field | protein | `oldff/leaprc.ff03` (comparison arms: ff14SB, ff19SB) |
| | ligand / ncAA side chain | GAFF2, AM1-BCC charges |
| | phosphotyrosine | `leaprc.phosaa10` |
| | ncAA alpha-carbon type | CT under ff03/ff99; CX under ff14SB/ff19SB |
| minimisation | engine | OpenMM 8, CUDA, mixed precision |
| | implicit solvent | OBC1 (`igb=2` equivalent), solute dielectric 1.0, solvent 80.0 |
| | cutoff | none (`NoCutoff`, equivalent to `cut=1000`) |
| | optimiser / cap | L-BFGS, **500 iterations** (calibrated to sander drift, §2.4) |
| | reference settings | `imin=1, maxcyc=5000, ncyc=2000, cut=1000.0, igb=2` |
| scoring | engine | `MMPBSA.py`, single-trajectory |
| | solver | Poisson–Boltzmann, internal PBSA |
| | interior dielectric | **2.0** (arms at 1.0 and 4.0 also computed) |
| | exterior dielectric | 80.0 |
| | nonpolar model | LCPO, `inp=1`, `cavity_surften=0.0072`, `cavity_offset=0` |
| | entropy | **not computed** |
| | dielectric injection | `use_sander=1` → `-make-mdins` → edit `_MM_gb.mdin` → `-use-mdins` |
| statistics | primary estimator | Spearman ρ within receptor group, equal weight per target |
| | decisive-pair margin | 0.5 kcal/mol (sensitivity, Table 2) |
| | size windows | 0 and 3 heavy atoms |
| | bootstrap unit | **compounds** (pairwise) or **receptor groups** (correlations) |
| | resamples | 2,000 within target; 5,000 for macro |

### 2.8 Software, hardware and code map

AmberTools 25 and 26 (`tleap`, `antechamber`, `parmchk2`, `prepgen`, `cpptraj`, `MMPBSA.py`);
OpenMM 8; RDKit 2025.03 and 2026.03; ParmEd; NumPy, SciPy, pandas; UCSF ChimeraX 1.12 for
visualisation. Hardware: NVIDIA RTX PRO 6000 Blackwell (local) and NVIDIA H100 80GB (Digital
Research Alliance of Canada).

Each methods step maps to a specific script in the deposited repository:

| section | step | script |
|---|---|---|
| 2.2–2.5 | full protocol, one complex | `scripts/zha_pipeline.py` |
| 2.1, 3.1 | protocol applied to the five-target cohort | `scripts/zha_cohort.py`, `run_cohort.sh` |
| 2.2 | build-quality audit (TER, clashes, atom counts) | `scripts/check_build.py` |
| 2.3 | non-canonical residue parameterisation | `scripts/ncaa_param.py`, `run_ncaa_all.sh` |
| 2.4 | minimisation and the drift calibration | `scripts/min_openmm.py`, `compare_minimisers.py` |
| 2.6, 3.4 | size-matched pairwise accuracy, Tables 3 and 6, Figure 4 | `scripts/size_matched_pairwise_accuracy.py` |
| 3.3 | size-confound analysis, Figure 2 | `scripts/fig_size_confound.py` |
| 3.1–3.2 | per-target statistics and estimator comparison | `scripts/analyse_cohort.py` |
| 3.8 | external 12CA5 panel | `scripts/stage_12ca5.py`, `build_12ca5_poses.py`, `run_12ca5.sh` |
| 4 | LptA design panel | `run_lpta_canonical.sh`, `run_lpta_ncaa.sh` |

Scripts use `$PROJECT_ROOT`, `$DATA_ROOT`, `$CONDA_ENV`, `$GPU_HOST`, `$CLUSTER_HOST` and
`$CLUSTER_USER` placeholders where local infrastructure paths were removed.

---

## 3. Results

### 3.1 The protocol reproduces, with a caveat on what is being reproduced

Scoring the source protocol's 50-complex Dataset I under ff03 with PB at εᵢₙ = 2 gives a Pearson
correlation of −0.382, against the −0.545 reported for the equivalent arm.[^4] The headline
−0.732 in that work is a two-step workflow that first re-ranks docked poses and then scores the
selected pose; our evaluation scores the deposited pose directly, so the two numbers are not
comparable and we do not treat the gap as a replication failure.

Transferring the protocol unchanged to the independent five-target cohort gives a macro
Spearman ρ of **+0.545** across 234 compounds (Figure 3; per-target values in Table 3).

**Table 3.** Per-target performance on the Awoonor-Williams cohort. Pearson *r* and Spearman ρ
are computed within receptor group. "Slope" is the per-target regression slope of score on
experiment; its reciprocal gives the factor by which ΔΔG spread is inflated.

| target | n | MM/PBSA *r* | *R*² | MM/PBSA ρ | FEP+ ρ | slope | offset (kcal/mol) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Cyclin A/B | 21 | +0.609 | 0.371 | +0.591 | +0.791 | 0.202 | −59.4 |
| JAK2 | 48 | +0.588 | 0.345 | +0.595 | +0.788 | 0.103 | −44.8 |
| KRAS | 45 | +0.300 | 0.090 | +0.225 | +0.786 | 0.125 | −74.8 |
| MCL-1 | 57 | +0.693 | 0.480 | +0.728 | +0.859 | 0.099 | −69.0 |
| PCSK9 | 63 | +0.551 | 0.304 | +0.584 | +0.929 | 0.208 | −64.9 |
| **macro** | **234** | **+0.548** | **0.318** | **+0.545** | **+0.831** | — | — |

Two things in Table 3 deserve emphasis. First, the commonly quoted "R" in such analyses is a
correlation, not a coefficient of determination: *r* = 0.548 corresponds to *R*² = 0.318, so
the protocol explains under a third of the variance, against 0.719 for FEP+. Second, the
per-target slopes of 0.10–0.21 mean the score's ΔΔG spread is inflated five- to tenfold, and
the offsets of 45–75 kcal/mol are the signature of the omitted entropy term.

### 3.2 How you average changes the answer

Pooling compounds that come from different receptor structures inflates the apparent
correlation, and KRAS shows this cleanly because its 45 compounds are split across two
crystal forms.

**Table 4.** KRAS under two estimators. Pooling reverses the sign of the result.

| KRAS | n | MM/PBSA ρ | FEP+ ρ |
|---|---:|---:|---:|
| within form 7YUZ | 17 | −0.307 | +0.576 |
| within form 8JJS | 28 | +0.046 | +0.782 |
| **macro over forms** | 45 | **−0.130** | +0.679 |
| pooled over both | 45 | **+0.225** | +0.786 |

The mechanism is explicit: the MM/PBSA score is offset by −3.05 kcal/mol between the two
forms, in the same direction as a −1.15 kcal/mol offset in the experimental labels. The pooled
+0.225 is therefore credit for separating two groups, not for ranking within either. FEP+ is
far less affected (+0.679 versus +0.786), as expected of a method whose signal does not come
from a group offset.

This matters beyond KRAS because it occurs *inside* what is conventionally reported as a
within-family metric. Any benchmark of this kind must state which estimator it uses; we use
within-receptor averaging throughout.

### 3.3 Size is a real driver of binding, and the score tracks it

Figure 2 quantifies the confound. Experimental ΔG tracks ligand heavy-atom count with a mean
within-target Spearman ρ of **−0.43** (panel A; negative because more atoms give more negative
ΔG). The MM/PBSA score tracks heavy-atom count more strongly still, at **−0.55** (panel B).

The first of these is physics, not bias: more residues bury more interface area. The second
means the score and the trivial descriptor are substantially measuring the same thing. Under
the conventional estimator, a heavy-atom count reaches ρ = 0.488 pooled in-family and **0.613
within receptor** — numerically better than MM/PBSA's 0.459 within receptor.

That comparison is not evidence that the physics is useless. It is evidence that the comparison
is ill-posed. Both quantities partly measure size; the one that measures it more directly wins
on a benchmark where size varies. Panel C shows that enough compounds share a heavy-atom count
within each target to make a fair comparison possible.

### 3.4 Size-matched performance: the primary result

**Table 5.** Pairwise accuracy on decisive pairs (|ΔΔG| ≥ 0.5 kcal/mol). Chance = 0.500.
Intervals are 95% bootstrap over compounds (Section 2.6). The left block is the strictest
matching; the right admits a three-atom window and is better powered.

| | | identical heavy atoms | | within 3 heavy atoms | |
|---|---:|---:|---|---:|---|
| target | | MM/PBSA | 95% CI | MM/PBSA | 95% CI |
| Cyclin A/B | | 0.533 | [0.053, 1.000] | 0.674 | [0.454, 0.863] |
| JAK2 | | 0.637 | [0.397, 0.831] | 0.598 | [0.437, 0.732] |
| KRAS | | 0.518 | [0.289, 0.727] | 0.496 | [0.339, 0.664] |
| MCL-1 | | 0.549 | [0.328, 0.757] | 0.594 | [0.444, 0.733] |
| **PCSK9** | | **0.765** | **[0.613, 0.882]** | **0.757** | **[0.671, 0.830]** |
| **macro** | | **0.601** | [0.472, 0.720] | **0.624** | **[0.555, 0.689]** |
| FEP+ macro | | 0.841 | [0.730, 0.912] | 0.839 | [0.782, 0.886] |

**The two bands give different levels of confidence, and we report both rather than the more
favourable one.** At exact size matching the macro estimate is 0.601 but its interval includes
chance; the posterior mass above 0.5 is 0.93. There are only 337 decisive pairs in that band,
drawn from few enough compounds that the compound-level interval is wide. At a three-atom
window — still tight relative to a 25-to-127 heavy-atom cohort, and still a regime in which a
size descriptor achieves only 0.602 — the macro estimate is 0.624 with an interval that
excludes chance and posterior mass above 0.5 exceeding 0.999.

The honest summary is therefore: **the skill is established at near-matched size and suggestive
but not conclusive at exactly matched size**, the latter limited by how few same-size decisive
pairs the cohort contains rather than by the size of the effect.

The skill is also unevenly distributed, and under compound-level resampling **only PCSK9 is
individually above chance** — at both bands. JAK2's point estimate of 0.637 at exact matching
has an interval of [0.397, 0.831] and cannot be separated from chance on its own; the same is
true of Cyclin A/B, KRAS and MCL-1. The macro result at the three-atom band is therefore driven
by PCSK9 with weak support from the others, not by uniform performance.

One candidate explanation is interface polarity, since a continuum-solvent treatment should add
most where desolvation dominates. PCSK9 does have much the most polar interface in the cohort
(median buried-area polarity 42.2%, against 30.5% for the next-highest). **But the explanation
does not survive the second case:** JAK2 is the fourth most polar of the five (23.8%) and the
second-best performer, while KRAS is second most polar (30.5%) and indistinguishable from
chance. Polarity alone therefore does not account for where the method works, and with five
targets we cannot identify what does. We report the pattern without an explanation rather than
offer one the data does not support.

The behaviour of the control across matching bands is itself informative. As the permitted size
difference widens from zero to three to five heavy atoms to unrestricted, the heavy-atom count
climbs from undefined, to 0.602, to 0.629, to 0.717, while MM/PBSA moves only from 0.601 to
0.738. The control gains almost all of its apparent performance from size variation; the
scoring function does not.

### 3.5 What the score is not

With entropy omitted, the raw score is not a free energy. It sits 45–75 kcal/mol from
experimental ΔG (Table 3), and matching the experimental scale requires fitting a slope and
intercept **within a target**, which requires compounds already measured in that series. An
RMSE of 1.11 kcal/mol is achievable only after that fit. It is not a prospective accuracy, and
should not be quoted as one.

### 3.6 Force-field choice does not transfer

On the benchmark, ff03 outperforms ff19SB (0.545 versus 0.495 macro, winning four of five
targets; only JAK2 prefers ff19SB, +0.604 versus +0.595). On the LptA design panel the ordering
reverses, with ff19SB ≥ ff14SB > ff03 in every branch.

**This inversion should not be over-read.** The LptA designs have no measured affinities, so
that comparison is scored against Rosetta `dG_separated`. We tested that reference directly on
239 benchmark compounds where measured affinities exist and found a macro Spearman of **−0.149**,
wrong-signed on four of five targets, and below a heavy-atom count (+0.487). A reference
standard that does not itself track experiment cannot establish which force field is better.

We report the comparison only as a caution: force-field selections made on one dataset should
not be assumed to transfer, and we do not claim ff19SB is superior for these designs.

### 3.7 Longer sampling does not help

A natural objection to single-structure scoring is that an ensemble average would do better. We
tested this on 40 KRAS compounds for which four independent 10 ns explicit-solvent replicas were
available, scored over the 2–10 ns window at 81 frames per replica — **324 frames per compound**.

The ensemble does not converge well enough for the averaging to help. The median
replicate-to-replicate standard deviation of the ensemble score is **6.81 kcal/mol**, while the
standard deviation *between* compounds is **5.74 kcal/mol**. The noise from which replicas you
happened to run is therefore larger than the entire signal distinguishing one compound from
another, a ratio of 1.19. Splitting by panel, the 20 confirmatory compounds give 5.60 against
5.54 (ratio 1.01) and the 20 primary compounds 8.70 against 4.88 (ratio 1.78); the conclusion
does not depend on the subset.

Four 10 ns replicas is a substantial sampling budget for a system of this size, and it is not
enough to make the ensemble score reproducible at the resolution the ranking requires. This is
consistent with the general finding that MM/PBSA correlation does not improve reliably with
longer sampling because the limiting error is bias rather than variance,[^3] and it supports the
source protocol's single-structure design on pragmatic grounds: a cheaper calculation that is
reproducible is preferable to an expensive one that is not.

### 3.8 An external test on a sixth target

The five benchmark targets are congeneric medicinal-chemistry series against folded pockets,
assembled by one group for one FEP+ campaign. They may resemble each other more than they
resemble a new target. We therefore tested the pipeline on a sixth system held entirely outside
development: 17 peptidomimetics against the anti-haemagglutinin antibody 12CA5, with binding
measured by biolayer interferometry (7 binders spanning 0.21–77 nM, 10 non-binders). Affinities
were sealed before any prediction was made.

Every method failed, and several inverted.

**Table 6.** Performance on the sealed 12CA5 panel. Directions were fixed before unblinding.

| method | ρ vs pK_D (7 binders) | binder/non-binder AUC |
|---|---:|---:|
| heavy-atom count | **−0.750** | 0.257 |
| molecular weight | −0.750 | 0.271 |
| chemical model (pAffinity) | −0.107 | 0.214 |
| chemical model (rank score) | −0.306 | 0.200 |

Chance AUC is 0.5; every value is below it. The tightest binder (0.21 nM) was ranked 16th of 17
by every method, and top-5 overlap with the experimental top 5 was zero throughout.

Two findings matter for the present work. First, **the size–affinity relationship reversed
sign**: heavy-atom count correlates at −0.750 here against +0.49 to +0.61 on the benchmark.
Size is predictive in congeneric series elaborated into a pocket, and anti-predictive against a
compact epitope groove where bulk is in the way. A method whose signal is partly size therefore
carries a sign that is target-dependent and not knowable in advance without measurement — which
is a further argument for evaluating at fixed size.

Second, **MM/PBSA could not be scored reliably at all**: 3 of 17 built poses carried positive
van der Waals terms and were excluded by the screen of Section 2.2, one failed to score, and
across the remainder the score correlated with the pose builder's clash count at ρ = −0.891.
Those peptidomimetics were threaded onto a template they do not match. The lesson is that
end-point scoring inherits the quality of its input pose, and that pose provenance must be
reported alongside any score (Section 4).

We present this as a limitation rather than a general claim: one external target with 17
compounds cannot establish how the method behaves on all new systems. It does establish that
benchmark performance on five related targets is not sufficient evidence of transfer.

---

## 4. The campaign this method serves: de novo macrocycles against LptA

Lipopolysaccharide transport protein A (LptA) is the periplasmic component of the seven-protein
Lpt bridge that moves lipopolysaccharide from the inner membrane to the cell surface in
Gram-negative bacteria. The insect-derived antimicrobial peptide thanatin binds the N-terminal
β-strand of LptA, disrupting both LptA–LptA dimerisation and the LptC–LptA interaction and
causing LptA degradation.[^17] Because that β-strand edge is also the dimerisation interface,
the site is a template for designed β-hairpin binders, and thanatin analogues optimised against
carbapenem-resistant Enterobacteriaceae have been characterised structurally and
thermodynamically.[^18],[^19]

The protocol above was established in order to triage this panel. We applied it to 212 designed
macrocycles in three branches: 20 canonical
head-to-tail cyclic peptides from generative backbone design followed by sequence design; 96
single-residue non-canonical substitutions on those canonical parents; and 96 single-residue
non-canonical substitutions on a thanatin derivative-7 scaffold templated on a co-crystal
structure. A total of 213 complexes were scored — the 212 designs plus the derivative-7 parent
as a reference.

**Pose quality.** Unlike the benchmark, these poses are generated rather than crystallographic,
which raises the question of whether scores read binding or pose-construction artefact. Across
all 213 scored complexes the built poses contain **zero** interface heavy-atom contacts below
1.8 Å, and minimisation drift spans only 0.49–1.42 Å. For comparison, the 12CA5 panel of
Section 3.8, whose peptidomimetics were threaded onto a non-matching template, showed 0–44
clashes per complex,
backbone RMSD to 4.72 Å, and MM/PBSA scores correlating with clash count at ρ = −0.891 — that
is, reading pose quality rather than chemistry. The LptA poses show no such pathology, which is
a necessary condition for the scores to be interpretable, though not a sufficient one.

**Size dependence differs by branch.** The score's correlation with ligand heavy-atom count is
−0.668 for the 20 canonical de novo designs, −0.606 for their 96 non-canonical swaps, and only
**−0.291** for the 96 thanatin swaps. This follows from construction: single substitutions on a
fixed 16-residue scaffold barely change size, whereas the de novo designs vary from 14 to 18
residues. The thanatin branch is therefore the least size-confounded, and by the argument of
Section 3.4 it is the branch on which the protocol's demonstrated skill is most likely to apply.

**No measured affinities exist for any of these designs.** All statements about them are
computational, and the protocol has not demonstrated predictive skill on this system. Rankings
are reported within branch only; the branches occupy different score ranges and pooling them
produces a spurious correlation.

---

## 5. Committed prospective test (registered; synthesis pending)

*This section describes a committed experiment that has not yet been executed. No results exist
at the time of writing, and none are claimed.*

Everything above is retrospective, on compounds measured by others. The campaign of Section 4
will proceed to synthesis, and we register the forward test here, before any compound is made,
so that the prediction cannot be adjusted after the fact.

**Design.** Ten de novo LptA designs, selected and frozen before synthesis, will be made and
assayed. The ranking is registered at the point of selection and will not be revised.

**Compound selection.** The ten span a four-heavy-atom window (149–153 atoms) — essentially
size-matched, which is the regime in which Section 3.4 establishes skill — while covering
34 kcal/mol of predicted score range (−104.0 to −69.8 kcal/mol). Spanning the predicted range
rather than selecting predicted winners is deliberate: a set of predicted winners cannot
falsify the ranking.

**What ten compounds can and cannot establish.** Power analysis by simulation, at the effect
size measured here (ρ ≈ 0.30, corresponding to ~60% pairwise accuracy):

**Table 7.** Simulated power to detect a positive rank correlation at *p* < 0.05.

| n | power at ρ = 0.30 | ρ = 0.50 | ρ = 0.75 |
|---:|---:|---:|---:|
| 10 | 0.12 | 0.28 | 0.67 |
| 20 | 0.22 | 0.58 | 0.96 |
| 40 | 0.43 | 0.88 | 1.00 |

**Ten compounds cannot confirm the ranking**: at the measured effect size the power is 12%, so
a null result would be uninformative. We state this in advance rather than discovering it
afterwards. Ten compounds can, however, do two useful things: detect a catastrophic inversion
(66% power against ρ = −0.75), which is a live concern because exactly such a reversal occurred
on the external target of Section 3.8; and supply ten measured anchors that permit within-series
calibration.

This test is therefore registered as a **sanity check and anchor acquisition**, not as
confirmatory validation. Confirmatory validation would require approximately 40 compounds, and
we do not claim it here.

---

## 6. Discussion

The central finding is that end-point MM/PBSA retains genuine ranking skill on macrocyclic
peptides when size is controlled — 62.4% pairwise accuracy within a three-atom window, where
chance is 50% and a heavy-atom count reaches only 60.2% — at roughly a thousandth of the cost of
rigorous free-energy simulation, which achieves 83.9% on the same pairs. At exactly matched size
the point estimate is similar (60.1%) but the cohort does not contain enough same-size decisive
pairs to separate it from chance at 95% confidence.

The methodological contribution is the comparison itself. Ligand size bias in scoring functions
has been recognised for two decades,[^5] and matched-pair benchmarks are an established response
in small-molecule work.[^6] Neither has been applied to macrocyclic peptides, where the size
range within a series is far larger and where the usual correction — normalising by heavy-atom
count or reporting ligand efficiency — is inappropriate because size is a mediator of the
physics rather than a nuisance. Exact size matching sidesteps this: it requires no correction
factor, makes no assumption about the functional form of the size dependence, and renders the
control informative by construction.

**What this licenses us to do.** The protocol was established to decide which LptA designs to
synthesise, and the validation above determines how far that decision can lean on it. Three
constraints follow directly. Ranking is admissible within a branch, because that is what was
tested; cross-branch ranking is not, because the branches occupy different score ranges.
Size-matched selection is preferred, because that is the regime in which skill was demonstrated
— which favours the thanatin branch and, within the de novo branch, a size-matched subset.
And the selection must span the predicted range rather than take the top of it, because a set
of predicted winners cannot test the ranking that chose them. Section 5 implements all three.

The practical reading for a design campaign is that this protocol is a **triage filter, not a
ranking to defend**. A 60% edge is worth having when the alternative is arbitrary selection and
the cost of being wrong is a thousand-dollar synthesis, particularly at seconds per compound.
It is not worth defending against a rigorous calculation when the decision is expensive enough
to justify one. The natural workflow is end-point scoring to reduce hundreds of candidates to
tens, then alchemical calculation on the shortlist.

Three cautions generalise beyond this work.

**The estimator matters.** Pooling compounds drawn from different receptor structures can invert
the sign of a result, as it does for KRAS (Table 4), so benchmark numbers computed that way
should be treated sceptically — including when the pooling is hidden inside what is reported as
a within-family metric.

**Pose provenance matters.** Where poses are built rather than crystallographic, scores can
track construction artefact: on the 12CA5 panel the score correlated with the pose builder's
clash count at ρ = −0.891. Any end-point score should be reported alongside evidence about the
geometry it was computed on.

**Reference standards must themselves be validated.** We could not settle whether force-field
choice transfers between our datasets, because the only available reference on the design panel
was a computational score that we measured at ρ = −0.149 against experiment (Section 3.6).
Agreement with an unvalidated reference is not evidence.

### 6.1 Status and what follows

This paper establishes a method and reports what it supports. It does not report a validated
prediction for LptA, and no such claim should be read into it. The sequence is deliberate:
establish the method against measured data from other systems, determine its limits, apply it
to the campaign panel under those limits, register the forward test, and only then commit
synthetic resources. Synthesis and assay of the ten registered compounds is the immediate next
step; the manuscript will be updated with those results, which Section 5 anticipates in both
directions.

---

## 7. Limitations

1. **Five targets.** Macro intervals rest on five values. Skill concentrates in PCSK9 and JAK2,
   and we could not identify why: interface polarity is the obvious candidate and it fails
   (Section 3.4). With five targets we cannot distinguish a chemical explanation from chance.
2. **The strictest band is underpowered.** At exactly equal heavy-atom count the cohort yields
   337 decisive pairs in total, and the macro interval [0.472, 0.720] includes chance. Cyclin
   A/B is the extreme case, with 15 pairs and an interval spanning [0.053, 1.000]. The
   three-atom band is the better-powered analysis and the one we rely on.
3. **Only one target is individually significant.** PCSK9 excludes chance at both bands; the
   other four do not. The macro result rests substantially on PCSK9, and a cohort without a
   polar-interface target might not reproduce it.
4. **Scores are not free energies** (Section 3.5) and must not be quoted as predicted affinities.
5. **No measured affinities exist for the LptA designs.** Section 4 is an application, not a
   validation, and Section 5 is a plan.
6. **Two selenium-containing non-canonical residues** are approximations, GAFF2 having no
   selenium parameters.
7. **Single-structure scoring** inherits whatever error the input pose carries. We verify pose
   quality explicitly (Section 4) but cannot certify that generated poses are correct.
8. **Censored measurements** are excluded from pairwise analyses, which may bias the admitted
   set toward the measurable dynamic range.
9. **We evaluate the protocol as published**, not an optimised variant. We did not re-tune the
   dielectric, solvation model or minimisation depth to maximise size-matched accuracy; doing so
   on five targets would overfit. The reported 60.1% is therefore a floor for this family of
   methods rather than the best achievable.

---

## Author contributions

C.L. designed and executed the computational work, performed the analyses and wrote the
manuscript. J.S.B. supervised the project and revised the manuscript. The prospective
validation of Section 5 will be executed experimentally by the Brown laboratory.

## Competing interests

The authors declare no competing financial interests.

## Acknowledgements

Computations were performed on local GPU hardware and on the Digital Research Alliance of
Canada national systems. We thank the authors of the two benchmark studies whose data made this
evaluation possible.

## Data and code availability

All code, computed scores, per-target statistics, bootstrap intervals and figures are available
at **https://github.com/Caroldorium/mmpbsa-macrocycle-benchmark** under an MIT licence.
`scripts/size_matched_pairwise_accuracy.py` reproduces Table 5 and the primary figure.

Two third-party sources are reproduced for reproducibility and cited in place: measured
affinities and FEP+ predictions from Awoonor-Williams et al.,[^2] and the reference protocol and
dataset tables from Zhao et al.[^4] All MM/PBSA scores, statistics, figures and code were
computed in this work.

---

## References

[^1]: Vinogradov AA, Yin Y, Suga H. Macrocyclic peptides as drug candidates: recent progress and remaining challenges. *J Am Chem Soc* 2019;141(10):4167–4181. doi:10.1021/jacs.8b13178. PMID 30768253.

[^2]: Awoonor-Williams E, Beautrait A, Petridis L, Hu X. Derisking affinity optimization for macrocycles and cyclic peptides: high-precision free energy simulations across five diverse targets. *J Chem Inf Model* 2026;66(14):8545–8560. doi:10.1021/acs.jcim.6c01122. PMID 42396675.

[^3]: Genheden S, Ryde U. The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities. *Expert Opin Drug Discov* 2015;10(5):449–461.

[^4]: Zhao H, Huang J, Weng G, Jiang D, Hu R, Kang Y, Hou T. Improving the predictive performance of binding affinities and poses for protein–cyclic peptide complexes through fine-tuned MM/PBSA(GBSA)-based methods. *Brief Bioinform* 2025;26(6):bbaf632. doi:10.1093/bib/bbaf632. PMID 41319044.

[^5]: Jacobsson M, Karlén A. Ligand bias of scoring functions in structure-based virtual screening. *J Chem Inf Model* 2006;46(3):1334–1343. doi:10.1021/ci050407t. PMID 16711752.

[^6]: Kalinowsky L, Weber J, Balasupramaniam S, Baumann K, Proschak E. A diverse benchmark based on 3D matched molecular pairs for validating scoring functions. *ACS Omega* 2018;3(5):5704–5714. doi:10.1021/acsomega.7b01194. PMID 31458770.

[^7]: Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. ff14SB: improving the accuracy of protein side chain and backbone parameters from ff99SB. *J Chem Theory Comput* 2015;11(8):3696–3713.

[^8]: Tian C, Kasavajhala K, Belfon KAA, et al. ff19SB: amino-acid-specific protein backbone parameters trained against quantum mechanics energy surfaces in solution. *J Chem Theory Comput* 2020;16(1):528–552.

[^9]: Homeyer N, Horn AHC, Lanig H, Sticht H. AMBER force-field parameters for phosphorylated amino acids in different protonation states. *J Mol Model* 2006;12(3):281–289.

[^10]: Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of a general amber force field. *J Comput Chem* 2004;25(9):1157–1174.

[^11]: Jakalian A, Jack DB, Bayly CI. Fast, efficient generation of high-quality atomic charges. AM1-BCC model: II. Parameterization and validation. *J Comput Chem* 2002;23(16):1623–1641.

[^12]: Jorgensen WL, Chandrasekhar J, Madura JD, Impey RW, Klein ML. Comparison of simple potential functions for simulating liquid water. *J Chem Phys* 1983;79(2):926–935.

[^13]: Cieplak P, Cornell WD, Bayly C, Kollman PA. Application of the multimolecule and multiconformational RESP methodology to biopolymers. *J Comput Chem* 1995;16(11):1357–1377.

[^14]: Onufriev A, Bashford D, Case DA. Exploring protein native states and large-scale conformational changes with a modified generalized Born model. *Proteins* 2004;55(2):383–394.

[^15]: Eastman P, Swails J, Chodera JD, et al. OpenMM 7: rapid development of high performance algorithms for molecular dynamics. *PLoS Comput Biol* 2017;13(7):e1005659.

[^16]: Miller BR III, McGee TD, Swails JM, Homeyer N, Gohlke H, Roitberg AE. MMPBSA.py: an efficient program for end-state free energy calculations. *J Chem Theory Comput* 2012;8(9):3314–3321.

[^17]: Moura ECCM, Baeta T, Romanelli A, et al. Thanatin impairs lipopolysaccharide transport complex assembly by targeting LptC–LptA interaction and decreasing LptA stability. *Front Microbiol* 2020;11:909. doi:10.3389/fmicb.2020.00909. PMID 32477309.

[^18]: Oi KK, Moehle K, Schuster M, Zerbe O. Early molecular insights into thanatin analogues binding to *A. baumannii* LptA. *Molecules* 2023;28(11):4335. doi:10.3390/molecules28114335. PMID 37298811.

[^19]: Sinha S, Dhanabal VB, Sperandeo P, Polissi A, Bhattacharjya S. Linking dual mode of action of host defense antimicrobial peptide thanatin: structures, lipopolysaccharide and LptA binding of designed analogs. *Biochim Biophys Acta Biomembr* 2021;1864(3):183839. doi:10.1016/j.bbamem.2021.183839. PMID 34915021.


---

## Figures

**Figure 1.** Pipeline schematic and validation logic. `fig1_pipeline.png`

**Figure 2.** The size confound. (A) Experimental ΔG tracks ligand heavy-atom count (mean
within-target ρ = −0.43). (B) The MM/PBSA score tracks it more strongly (ρ = −0.55).
(C) Within-target size distributions, showing that enough compounds share a heavy-atom count to
permit direct comparison. `fig2_size_confound.png`

**Figure 3.** Per-target parity for raw MM/PBSA, calibrated MM/PBSA and published FEP+ against
experiment. `zha25_fepplus_parity_5families.png`

**Figure 4.** Size-matched pairwise accuracy, at identical heavy-atom count (left) and within
three heavy atoms (right), against the FEP+ ceiling. Chance is 0.5; a heavy-atom count is
undefined in the left panel. `size_matched_pairwise_accuracy.png`
