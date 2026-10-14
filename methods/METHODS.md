# Computational Methods

**Physics-based affinity ranking of de novo cyclic peptides and non-canonical analogues against LptA**

Prepared 27 September 2026. Every parameter choice below is stated with the reason it was made
and, where a choice was tested rather than assumed, with the result of that test.

---

## 1. Benchmark cohort and reference data

Method development was carried out on a labelled benchmark before any application to unlabelled
designs. The benchmark comprises 241 macrocyclic peptides against five protein targets
(cyclin A/B, JAK2, KRAS, MCL-1, PCSK9), of which 235 carry an experimental binding free energy
and 234 survive quality control (§2.4). Published FEP+ predictions for the same compounds
(Schrödinger, Inc.) provide an upper reference for what an alchemical free-energy method
achieves on this chemistry.

Compounds are not independent within a target. Two of the five targets are represented by more
than one receptor crystal form — KRAS by 7YUZ (n = 17) and 8JJS (n = 28), PCSK9 by 6XIC
(n = 38), 6XIE (n = 8) and 7S5G (n = 17) — giving eight receptor groups in total. All primary
statistics are therefore computed **within receptor group** and averaged (§5.1), not pooled.

---

## 2. Structure preparation and single-structure MM/PBSA

### 2.1 Protocol origin

The end-point scoring protocol reproduces Zhao *et al.* (2025), *Briefings in Bioinformatics*
26(6):bbaf632, a fine-tuned MM/PBSA and MM/GBSA workflow for protein–cyclic-peptide complexes.
The protocol is **single-structure**: one energy-minimised geometry is scored, not an ensemble
average over a molecular dynamics trajectory. This is the authors' design, and it is stated
explicitly here because "MM/PBSA" is more commonly associated with ensemble averaging. Our
minimised trajectories contain exactly one frame, and MMPBSA.py reports
`Calculations performed using 1.0 complex frames` with a standard error of zero.

The choice is defensible on published evidence: for MM/PBSA, correlation with experiment
generally does not improve with longer sampling, because the limiting error is bias in the
energy function rather than variance from short trajectories (Genheden & Ryde, *Expert Opin.
Drug Discov.* 2015, 10:449). We tested this directly rather than accepting it (§6).

### 2.2 Topology construction

Receptor and peptide are cleaned with `pdb4amber` in two passes. The first pass retains
prepared hydrogens, which is what allows histidine tautomer and cystine assignments to be
carried into the residue names; the second pass strips hydrogens with `-y` so that tleap
rebuilds them against its own templates. A single pass is insufficient because prepared files
carry atoms such as `HXT` for which the older force fields have no type.

Cyclisation is imposed explicitly. Head-to-tail macrocycles receive an explicit
`bond lig.N.C lig.1.N` in tleap; disulfide-cyclised peptides are detected geometrically
(SG–SG < 2.5 Å) and left to the CONECT records written by `pdb4amber`, since re-issuing those
bonds manually produces `1-4: cannot add bond` errors. For macrocycles with uncharged termini,
`clearPdbResMap` is issued so that tleap does not apply terminal residue templates to a
molecule that has no termini.

Two structure-preparation defects were found and corrected during validation, both of which
produced plausible-looking but wrong energies:

* **Spurious chain breaks.** `pdb4amber` inserted 9–21 `TER` records into continuous chains,
  each of which tleap capped, producing complexes with initial energies above 10^17 kcal/mol.
  `TER` records whose flanking C–N distance is below 2.0 Å are now deleted. After the fix, the
  atom count of a reference Dataset I complex matches the published value exactly.
* **Unrelieved steric clash.** One benchmark compound (`MCL_1__cpd_55`) passed the pipeline's
  `status == ok` gate carrying a van der Waals term of **+277,293 kcal/mol**, against
  ≤ −50 kcal/mol for every other compound. The `ok` status does not test for unrelieved
  overlap. Screening on `VDWAALS > 0` is now mandatory before any Pearson-based or
  absolute-error analysis. Rank statistics were unaffected, but the pooled Pearson correlation
  moved from −0.048 to +0.265 on removing this single compound.

### 2.3 Force fields

| component | assignment | reason |
|---|---|---|
| protein | `oldff/leaprc.ff03` | the legacy leaprc named by the original protocol; it takes terminal residues from `all_amino{n,c}t94.lib` and therefore carries NCYX/CCYX. `leaprc.protein.ff03.r1` uses `amino{n,c}t03.lib`, which has neither, and fails on 20 of 50 reference complexes whose terminal residue is a cystine |
| protein (comparison arms) | `leaprc.protein.ff14SB` (Maier *et al.*, *JCTC* 2015, 11:3696), `leaprc.protein.ff19SB` (Tian *et al.*, *JCTC* 2020, 16:528) | tested, not assumed — see §5.3 |
| phosphotyrosine | `leaprc.phosaa10` (Homeyer *et al.*, *J. Mol. Model.* 2006, 12:281) | JAK2's receptor carries a phosphotyrosine; without it tleap reports `Atom .R<PTR>.A<O3P> does not have a type` and fails |
| ligands and ncAA side chains | GAFF2 (Wang *et al.*, *J. Comput. Chem.* 2004, 25:1157) with AM1-BCC charges (Jakalian *et al.*, *J. Comput. Chem.* 2002, 23:1623) | standard for non-protein chemistry; charges derived per compound with `antechamber` |
| water (alchemical legs only) | TIP3P (Jorgensen *et al.*, *J. Chem. Phys.* 1983, 79:926) | — |

ff19SB is used off-label under implicit solvent: its CMAP correction was parameterised for OPC
water, and CMAP applies only to the 20 standard residues, so a GAFF2-hybrid non-canonical
residue never receives it. This is reported rather than hidden.

### 2.4 Energy minimisation

Minimisation follows the published settings — `imin=1, maxcyc=5000, ncyc=2000, cut=1000.0,
igb=2` — corresponding to 2000 steps of steepest descent followed by conjugate gradient in the
Onufriev–Bashford–Case GB model (OBC, `igb=2`; Onufriev *et al.*, *Proteins* 2004, 55:383) with
no nonbonded cutoff.

Because the AmberTools build available here does not ship `pmemd.cuda`, and `sander` requires
hours per complex at `cut=1000`, minimisation was performed with OpenMM (Eastman *et al.*,
*PLoS Comput. Biol.* 2017, 13:e1005659) using `app.OBC1` and `app.NoCutoff`, which are the
direct equivalents of `igb=2` and `cut=1000`. The two engines are not identical, so the
OpenMM L-BFGS iteration count was **calibrated against sander by structural drift**: sander's
5000 SD+CG cycles move the complex 0.84 Å RMSD, and an unrestricted OpenMM minimisation moves
it 1.25 Å. A cap of 500 L-BFGS iterations reproduces sander's drift and is used throughout.

This calibration is load-bearing rather than cosmetic. On KRAS, removing the cap changes the
correlation with experiment from +0.225 to +0.529, and the spread introduced by the stopping
point (SD 3.60 kcal/mol) is as large as the entire between-compound signal (SD 3.53), a
ratio of **1.02**. The effect does not generalise (§5.4).

### 2.5 Free-energy evaluation

`MMPBSA.py` (Miller *et al.*, *JCTC* 2012, 8:3314) is run single-trajectory: complex, receptor
and ligand energies are all extracted from one trajectory, so the receptor's internal energy
cancels exactly.

Poisson–Boltzmann calculations use the internal PBSA solver with an interior dielectric of
2.0, exterior 80.0, and the LCPO nonpolar model (`inp=1`, `cavity_surften=0.0072`,
`cavity_offset=0.0`). An interior dielectric of 2 is the value the source protocol identifies
as optimal; arms at εin = 1 and 4, and GB arms (`igb=1,2,5`), are computed in parallel for
comparison. **No entropy term is included**, following the source protocol.

Omitting entropy has a specific consequence that must be stated whenever these numbers are
used: the scores are not free energies. They sit 45–75 kcal/mol below experimental ΔG and the
regression slope against experiment is 0.10–0.21, meaning ΔΔG spread is inflated five- to
tenfold. An RMSE of 1.11 kcal/mol is achievable *only after fitting a slope and intercept per
target*, which requires measured compounds in that series. It is not a prospective accuracy.

The `intdiel` keyword is silently ignored by both engines when passed through `&pb`. The
dielectric is therefore set by writing `use_sander=1`, generating input files with
`-make-mdins`, injecting `intdiel` into the generated `_MM_gb.mdin`, and re-running with
`-use-mdins`. Without this, every "different dielectric" arm silently returns identical
numbers.

---

## 3. Non-canonical amino acid parameterisation

The LptA design panel comprises 212 molecules: 20 canonical head-to-tail cyclic peptides and
192 single-site non-canonical analogues. Scoring these required Amber parameters for **170
non-standard residues** — 165 unique ncAA codes plus five residues that constitute the thanatin
derivative-7 scaffold itself (`EU0`, `HYP`, `LE1`, `DAB`, `4FO`), none of which tleap knows.

### 3.1 Route

Chemistry is taken from a curated in-house library (Rosetta `.params` plus a rotamer PDB per
code). For each residue: the polymer-form rotamer is capped with ACE and NME, charges are
derived with AM1-BCC on the neutral capped dipeptide, `prepgen` removes the caps and
redistributes their charge to give a residue that nets zero, and `parmchk2` supplies missing
GAFF2 terms.

**Why cap.** AM1-BCC on a free amino acid fits backbone charges in an environment
(–NH₂ / –COOH) the residue never occupies. Capping places the charge fit in a peptide-like
context, which is the standard construction for Amber non-standard residues (Cieplak *et al.*,
*J. Comput. Chem.* 1995, 16:1357).

**Why the backbone is retyped.** `antechamber` types the whole molecule with GAFF2, leaving the
backbone as `ns`/`c`/`c3`. tleap then has no parameters across the peptide bond to a
neighbouring standard residue (`Could not find bond parameter for atom types: C - ns`). The
backbone is therefore reassigned to protein types (N, H, CT or CX, H1, C, O) while the side
chain remains GAFF2. The alpha carbon is typed **CT for ff03/ff99 and CX for ff14SB/ff19SB**:
ff03 has no CX at all, and a prep typed for one force field fails on the other, so two parameter
sets are maintained.

The ACE/NME cap atoms are retyped as well. They stand in for the neighbouring residues, and
leaving them GAFF2-typed makes `parmchk2` generate cross terms against the wrong types — a
proline-like residue bonded through N needs an angle `C–N–c5` and would otherwise only receive
`c–N–c5`.

### 3.2 Validation

Every residue must build a real peptide, not merely produce a file. Each parameter set is
tested by constructing ALA–X–ALA in tleap and requiring zero missing parameters. **All 170
residues pass under both backbone typings.**

This test is necessary because tleap **writes a prmtop even when parameters are missing**,
substituting defaults. A non-zero missing-parameter count is treated as failure.

### 3.3 Defects found

Six defects were identified during parameterisation, four of which produced silently wrong
output rather than an error:

1. **Atom names read as elements.** `antechamber` infers the element from the atom name. Every
   residue's alpha carbon is `CA` (calcium) and 61 aromatic carbons are named `CL*` (chlorine);
   `FE`, `NI` and `CO` also occur. One residue built "successfully" carrying two phantom
   chlorines and correspondingly wrong charges. Residues are now written under element-safe
   placeholder names and the real names restored afterwards.
2. **Silent selenium substitution.** GAFF2 has no selenium and `sqm` has no AM1 parameters for
   it, so `antechamber` retyped Se as sulfur and `parmchk2` reported nothing missing. Two
   residues (`Z0U`, `Z1P`) are affected and are flagged `ok_element_unsupported` rather than
   shipped as clean.
3. **Bond orders guessed.** `antechamber`'s perception fails on hypervalent sulfur
   (sulfonamides) and iodoaryls. Bond orders are now supplied from the library in a MOL2.
4. **Cap geometry collisions.** Fixed-dihedral NME placement puts the cap nitrogen 0.15 Å from
   the residue's own carbonyl oxygen whenever ψ is near 180°. Placement is now by bisector with
   a hard clash check, and cap names are chosen per residue to avoid collisions with side-chain
   atoms (75 of 170 required this).
5. **Library-internal name inconsistency.** For nine halogenated residues the rotamer PDB
   truncates `CLM1` to `CLM`, and one residue names *both* chlorines `CLM`. Those atoms matched
   no bond record, and `prepgen` discarded them silently, leaving an atom with a `nan` internal
   coordinate that surfaced only later as a tleap fatal error. Names and bonds are now taken
   from the `.params` file, with coordinates matched positionally.
6. **Case-sensitive name restoration.** `antechamber` rewrites a two-letter placeholder to
   proper element case (`CL1` → `Cl1`), so a case-sensitive mapping silently leaves halogens
   unnamed.

---

## 4. Alchemical relative binding free energy

### 4.1 Perturbation maps

Edges are scored LOMAP-style (Liu *et al.*, *J. Comput. Aided Mol. Des.* 2013, 27:755) from the
maximum common substructure, favouring the fewest atoms appearing or disappearing. Because MCS
on 100-atom macrocycles is combinatorially expensive, candidate pairs are pre-screened by
Morgan fingerprint similarity (Rogers & Hahn, *J. Chem. Inf. Model.* 2010, 50:742) and only the
most similar pairs go to MCS.

**Edges never cross receptor crystal forms.** A relative free energy between compounds modelled
in different receptor frames is not a well-defined quantity.

**Connectivity is enforced explicitly.** A relative map yields ΔΔG *along edges*, so a compound
in a disconnected component cannot be placed on its group's ΔG scale however well its own edges
converge. A quality filter alone (≤10 heavy atoms changed) left only 150 of 235 compounds
connected — one target had 87 edges but only 22 of 48 compounds rankable — because the filter
removes precisely the wide edges that join components. Connectivity is restored in two stages:
cheapest already-scored edge per join, then on-demand MCS on the most similar cross-component
pairs. Final coverage is 184 of 235 (78%). Edges added for connectivity are flagged and carry
lower weight, being wider perturbations expected to converge less well.

Remaining gaps are not forced. Where MCS will not converge between chemically unrelated
macrocycles, the corresponding edge would not converge in thermodynamic integration either, and
manufacturing one buys coverage at the cost of trustworthy numbers.

### 4.2 System construction

Receptor coordinates are identical for every compound within a receptor group (Cα RMSD
0.000 Å), so both ligands of an edge already occupy a common frame and no alignment is required.

Amber's softcore TI requires the common-core atoms to be the **first atoms of each TI region and
to correspond in order** between regions. Independently parameterised ligands do not satisfy
this, so both are renumbered from the MCS mapping, with hydrogens carried alongside their parent
heavy atoms and core atoms whose hydrogen counts disagree excluded from the common block.

Complexes are solvated in a truncated octahedron with a 12 Å buffer and neutralised.
`setBox` is **not** issued before `solvateoct`: pre-setting a box makes the buffer be measured
from the box rather than the solute, and AmberTools 25 and 26 resolve that ambiguity
differently — 7,417 versus 12,013 waters for an identical solute. Removing it makes the build
reproducible across versions and reduces the campaign cost by 38%.

### 4.3 Simulation

Thermodynamic integration uses Amber's GTI implementation (Lee *et al.*, *J. Chem. Inf. Model.*
2020, 60:5595) with one-step softcore (`ifsc=1`), 13 λ windows, and MBAR reweighting (Shirts &
Chodera, *J. Chem. Phys.* 2008, 129:124105) via `ifmbar=1` with all sampled λ listed in
`mbar_lambda`. Each window runs minimisation, 100 ps heating, 0.5 ns equilibration and 5 ns
production at 298.15 K under a Langevin thermostat and Monte Carlo barostat. ΔΔG_bind is the
difference between the complex and solvent legs.

Two implementation constraints are recorded because neither is documented in a way that would
be encountered before failure:

* **Minimisation with `ifsc=1` requires `ntmin=2`** (steepest descent); Amber refuses conjugate
  gradient with softcore.
* **The SPFP build cannot run a TI region of this size.** `pmemd.cuda` fails in
  `kgBuildSpecial2RestNBPreList` with `invalid argument launching kernel`. The failure was
  isolated by showing that the 7,609-atom solvent leg fails identically to the 38,928-atom
  complex (excluding system size), that `gti_cut=0` does not help, and that the same system
  runs cleanly as plain MD with TI keywords removed. The double-precision build
  `pmemd.cuda_DPFP` works. Amber TI is normally exercised on ligands of 20–60 atoms; these
  macrocycles give a 385-atom TI region across the two states.

Additionally, the `amber/22.5` module's `pmemd.cuda` is compiled for sm_50–sm_80 and cannot
execute on H100 hardware; `amber-pmemd/24.3` (CUDA 12.6) adds sm_86 and sm_90.

---

## 5. Statistical treatment

### 5.1 Estimator

The primary statistic is **Spearman's ρ computed within receptor group and averaged across
targets**, with each target weighted equally regardless of how many compounds it contributes.

Pooling compounds across targets is rejected because it rewards a method for knowing that one
series is systematically more potent than another, which is not the question a design campaign
asks. On this cohort the pooled estimator rates a randomly perturbed pose (+0.049) above the
deposited crystal pose (−0.048). The same pathology appears *within* a target when a target
spans two receptor forms: KRAS scores +0.225 pooled across its two forms but −0.130 averaged
within them, the pooled value being credit for separating two groups whose mean affinities
differ by 1.15 kcal/mol rather than for ranking inside either.

### 5.2 Controls

Every correlation is reported against two controls:

* **Molecular size.** Ligand heavy-atom count is computed on the same compounds with the same
  estimator. Partial Spearman correlations controlling heavy-atom count are reported wherever a
  method is claimed to work.
* **Cluster bootstrap.** Uncertainty is estimated by resampling **receptor groups**, not
  compounds, because compounds within a group are not independent. 4,000–10,000 resamples.

Censored measurements (assay ceilings) are retained and flagged; correlations are reported both
including and excluding them.

### 5.3 Results of force-field comparison

Force-field choice does not transfer between datasets. On the benchmark ff03 outperforms ff19SB
(0.545 versus 0.495 macro, ff03 winning four of five targets); on the LptA design panel the
ordering reverses, with ff19SB ≥ ff14SB > ff03 in every branch. The two comparisons use
different reference standards — measured affinity for the benchmark, a computational score for
the designs — so the inversion is suggestive rather than established, and is reported as such.

### 5.4 Negative results

Reported because they constrain what the method can be expected to do.

* **Minimisation depth does not generalise.** Removing the iteration cap improves KRAS by
  +0.304 and degrades PCSK9 by −0.280; averaged by target the macro correlation moves
  0.545 → 0.536. Averaged by receptor group it moves 0.459 → 0.549, but a cluster bootstrap
  over those groups puts the gain at +0.110 with a 95% CI of [−0.066, +0.230], crossing zero.
  KRAS was genuinely under-minimised and PCSK9 was not, so this is a target-specific fix
  rather than a general one.
* **Component reweighting overfits.** Scoring KRAS on the electrostatic term alone gives 0.574
  against 0.225 for the full score, but the same operation across all five targets gives 0.149
  against 0.545.
* **Longer sampling is not expected to help.** On 40 KRAS compounds with four independent 10 ns
  explicit-solvent replicas each, scored over 2–10 ns at 81 frames per replica (324 frames per
  compound), the median replicate-to-replicate SD of the ensemble score is **6.81 kcal/mol**
  against a between-compound SD of **5.74** — the noise exceeds the signal (ratio 1.19). By
  panel: confirmatory 5.60 vs 5.54 (1.01), primary 8.70 vs 4.88 (1.78). Four 10 ns replicas do
  not make the ensemble score reproducible at the resolution the ranking needs.
* **Displacement recovery carries no signal.** Rigid extraction is geometrically impossible on
  this chemistry — over 1,500 sampled directions per complex the best ligand–receptor clearance
  is 1.22–2.54 Å, inside van der Waals contact, so a rigid pull measures steric clash rather
  than dissociation. Displace-then-relax works mechanically (46/69/73% of a 1/2/3 Å displacement
  is recovered) but does not predict affinity: −0.058 by target, −0.005 by receptor group, with
  per-target signs that flip.

### 5.5 Reference methods

Rosetta `dG_separated` was evaluated as a possible reference and **rejected**. Tested against
measured affinity on 239 benchmark compounds it gives a macro Spearman of **−0.149**,
wrong-signed on four of five targets (MCL-1 −0.824) and below a heavy-atom count (+0.487). It
is therefore not used to validate any other method. Any agreement statistic computed against it
would be agreement with a non-predictor.

---

## 6. Molecular dynamics

Where dynamics are required — for the sampling-convergence test of §5.4 and for visualisation —
trajectories are generated in OpenMM with the OBC1 generalised-Born model and no cutoff.
Simulations start from the minimised complex, use a Langevin middle integrator at 298.15 K with
a 2 fs timestep and hydrogen-bond constraints, and are preceded by 100 ps of equilibration.

**These trajectories match the minimisation step, not the scoring step, and the distinction
matters.** Minimisation uses `igb=2` at the Amber default solute dielectric of 1 (§2.4), which
is what the OpenMM OBC1 setting reproduces. The reported score, however, uses the
Poisson–Boltzmann solver at an interior dielectric of 2 (§2.5). Sampling and scoring are
therefore related but **not** self-consistent, and these trajectories are not a validated
equilibrium ensemble for the scored potential. An earlier version of this section claimed they
were generated under the same solvent model as the score; that was wrong.

The construction is still closer to the scored potential than the earlier in-house KRAS
protocol, which sampled in explicit solvent and stripped water before scoring — but the
remaining mismatch (GB vs PB, and ε 1 vs 2) is a real limitation on any ensemble claim.

---

## 7. Software

AmberTools 25/26 and Amber 22/24 (`tleap`, `antechamber`, `parmchk2`, `prepgen`, `cpptraj`,
`MMPBSA.py`, `pmemd.cuda_DPFP`); OpenMM 8 (CUDA platform, mixed precision, DPFP where required);
RDKit 2025.03/2026.03 (MCS, fingerprints); ParmEd; NumPy, SciPy, pandas; UCSF ChimeraX 1.12 for
visualisation. GPU hardware: NVIDIA RTX PRO 6000 Blackwell (local) and NVIDIA H100 80GB
(Digital Research Alliance of Canada, Nibi).

---

## 8. Scope and limitations

1. **No design in the LptA panel has a measured binding affinity.** All statements about those
   molecules are computational, and no method applied to them has demonstrated predictive skill
   on that system.
2. **The scoring function does not beat molecular size on labelled data.** For scale, the
   published FEP+ predictions on the same compounds score 0.831 by target and 0.805 by
   receptor group; the MM/PBSA protocol reaches 0.545 and 0.459 respectively. Conditioned on
   receptor group, MM/PBSA scores 0.459 against 0.613 for a heavy-atom count, a deficit of
   −0.164 with a 95% CI of [−0.299, −0.026] excluding zero. Only one of five targets (PCSK9,
   the most polar interface) shows size-independent skill, at a partial ρ of 0.606
   [0.413, 0.758]; every other target's partial correlation interval includes zero.
3. **Scores are ranking devices, not free energies** (§2.5).
4. **Rank within a branch, never across.** The three design branches occupy different score
   ranges; pooling them produces a spurious correlation of −0.316.
5. **Two selenium-containing residues are approximations** (§3.3).
6. **The derivative-7 parent comparison is confounded**: the scored parent is an unrelaxed
   reference while its analogues are relaxed, so the apparent improvement of all 96 analogues
   over the parent contains a relaxation effect and is not a clean chemical comparison.

---

## References

Cieplak P, Cornell WD, Bayly C, Kollman PA. *J. Comput. Chem.* 1995, 16:1357.
Eastman P, *et al.* *PLoS Comput. Biol.* 2017, 13:e1005659.
Genheden S, Ryde U. *Expert Opin. Drug Discov.* 2015, 10:449.
Homeyer N, Horn AHC, Lanig H, Sticht H. *J. Mol. Model.* 2006, 12:281.
Jakalian A, Jack DB, Bayly CI. *J. Comput. Chem.* 2002, 23:1623.
Jorgensen WL, *et al.* *J. Chem. Phys.* 1983, 79:926.
Lee T-S, *et al.* *J. Chem. Inf. Model.* 2020, 60:5595.
Liu S, *et al.* *J. Comput. Aided Mol. Des.* 2013, 27:755.
Maier JA, *et al.* *J. Chem. Theory Comput.* 2015, 11:3696.
Miller BR III, *et al.* *J. Chem. Theory Comput.* 2012, 8:3314.
Onufriev A, Bashford D, Case DA. *Proteins* 2004, 55:383.
Rogers D, Hahn M. *J. Chem. Inf. Model.* 2010, 50:742.
Shirts MR, Chodera JD. *J. Chem. Phys.* 2008, 129:124105.
Tian C, *et al.* *J. Chem. Theory Comput.* 2020, 16:528.
Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. *J. Comput. Chem.* 2004, 25:1157.
Zhao H, Huang J, Weng G, Jiang D, Hu R, Kang Y, Hou T. *Briefings in Bioinformatics*
2025, 26(6):bbaf632. doi:10.1093/bib/bbaf632.
