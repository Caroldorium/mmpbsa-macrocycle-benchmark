# Manuscript diffs, v1 → v9

The actual text changes between versions, so the prose can be audited directly rather than
through my summary of it. Lines beginning `-` were removed, lines beginning `+` were added.

Reproduce with:

```bash
git log --format='%h' --reverse -- manuscript/MANUSCRIPT.md
git diff <prev> <next> -- manuscript/MANUSCRIPT.md
```

---

## v1 — `85b7108` — Add draft manuscript, pipeline schematic and size-confound figure

Initial version; no diff.

---

## v2 — `3a6a94a` — Manuscript v2: retract an unsupported claim, add the external-target test

`1 file changed, 137 insertions(+), 37 deletions(-)`

```diff
@@ -32,8 +32,14 @@ JAK2 (63.7%), indistinguishable from chance for Cyclin A/B, KRAS and MCL-1.
 
-We further show that the apparent accuracy of end-point scoring is sensitive to the choice of
-estimator — pooling compounds across targets rates a deliberately perturbed pose above the
-deposited crystal pose — and that the score is a ranking device rather than a free energy,
-sitting 45–75 kcal/mol from experiment with per-target regression slopes of 0.10–0.21. We
-apply the validated protocol to 308 designed macrocycles against the lipopolysaccharide
-transport protein LptA and register a prospective test.
+Three further results bear on how such benchmarks should be read. The apparent accuracy depends
+on the estimator: for KRAS, whose compounds span two receptor crystal forms, pooling gives
+ρ = +0.225 while averaging within forms gives −0.130, the pooled value being credit for
+separating two groups offset by 1.15 kcal/mol. The score is a ranking device rather than a free
+energy, sitting 45–75 kcal/mol from experiment with per-target regression slopes of 0.10–0.21.
+And on a sixth target held entirely outside development — 17 peptidomimetics against an
+antibody, affinities sealed before prediction — every method fell below chance and the
+size–affinity relationship reversed sign, showing that performance on five related targets is
+not sufficient evidence of transfer.
+
+We apply the protocol to 308 designed macrocycles against the lipopolysaccharide transport
+protein LptA and register a prospective test with its power analysis stated in advance.
 
@@ -168,4 +174,6 @@ every other compound. A status check does not test for unrelieved overlap. Scree
 `VDWAALS > 0` is mandatory before any analysis and is applied throughout. Rank statistics are
-insensitive to this compound, but the pooled Pearson correlation moves from −0.048 to +0.265 on
-removing it.
+largely insensitive to this compound (macro Spearman 0.532 → 0.545), but the pooled Pearson
+correlation moves from −0.048 to **+0.265** on removing it, and MCL-1's Pearson *r* from −0.124
+to +0.693. A single unrelieved clash was therefore enough to make the method look uncorrelated
+with experiment.
 
@@ -265,4 +273,4 @@ each target weighted equally. Pooling compounds across targets is rejected: it r
 for knowing that one series is systematically more potent than another, which is not the
-question a design campaign asks. On this cohort the pooled estimator rates a randomly perturbed
-pose (+0.049) above the deposited crystal pose (−0.048).
+question a design campaign asks. Section 3.2 quantifies the resulting inflation within our own
+data.
 
@@ -337,10 +345,25 @@ the offsets of 45–75 kcal/mol are the signature of the omitted entropy term.
 
-The pooled estimator — computing one correlation over all compounds from all targets — rates a
-deliberately perturbed pose (+0.049) above the deposited crystal pose (−0.048). The same
-pathology appears within a single target: KRAS scores +0.225 pooled across its two receptor
-crystal forms but −0.130 averaged within them, the pooled value being credit for separating two
-groups whose mean affinities differ by 1.15 kcal/mol rather than for ranking inside either.
+Pooling compounds that come from different receptor structures inflates the apparent
+correlation, and KRAS shows this cleanly because its 45 compounds are split across two
+crystal forms.
 
-Any benchmark of this kind must state which estimator it uses. We use within-receptor averaging
-throughout.
+**Table 2.** KRAS under two estimators. The pooled value is more than three times the
+within-form average, and of opposite sign.
+
+| KRAS | n | MM/PBSA ρ | FEP+ ρ |
+|---|---:|---:|---:|
+| within form 7YUZ | 17 | −0.307 | +0.576 |
+| within form 8JJS | 28 | +0.046 | +0.782 |
+| **macro over forms** | 45 | **−0.130** | +0.679 |
+| pooled over both | 45 | **+0.225** | +0.786 |
+
+The mechanism is explicit: the MM/PBSA score is offset by −3.05 kcal/mol between the two
+forms, in the same direction as a −1.15 kcal/mol offset in the experimental labels. The pooled
++0.225 is therefore credit for separating two groups, not for ranking within either. FEP+ is
+far less affected (+0.679 versus +0.786), as expected of a method whose signal does not come
+from a group offset.
+
+This matters beyond KRAS because it occurs *inside* what is conventionally reported as a
+within-family metric. Any benchmark of this kind must state which estimator it uses; we use
+within-receptor averaging throughout.
 
@@ -367,3 +390,3 @@ undefined:
 
-**Table 2.** Pairwise accuracy on decisive pairs (|ΔΔG| ≥ 0.5 kcal/mol). Chance = 0.500.
+**Table 3.** Pairwise accuracy on decisive pairs (|ΔΔG| ≥ 0.5 kcal/mol). Chance = 0.500.
 
@@ -404,7 +427,13 @@ should not be quoted as one.
 On the benchmark, ff03 outperforms ff19SB (0.545 versus 0.495 macro, winning four of five
-targets). On the LptA design panel the ordering reverses, with ff19SB ≥ ff14SB > ff03 in every
-branch. The two comparisons use different reference standards — measured affinity for the
-benchmark, a computational score for the designs — so the inversion is suggestive rather than
-established, and we report it as a caution against transferring force-field choices between
-datasets.
+targets; only JAK2 prefers ff19SB, +0.604 versus +0.595). On the LptA design panel the ordering
+reverses, with ff19SB ≥ ff14SB > ff03 in every branch.
+
+**This inversion should not be over-read.** The LptA designs have no measured affinities, so
+that comparison is scored against Rosetta `dG_separated`. We tested that reference directly on
+239 benchmark compounds where measured affinities exist and found a macro Spearman of **−0.149**,
+wrong-signed on four of five targets, and below a heavy-atom count (+0.487). A reference
+standard that does not itself track experiment cannot establish which force field is better.
+
+We report the comparison only as a caution: force-field selections made on one dataset should
+not be assumed to transfer, and we do not claim ff19SB is superior for these designs.
 
@@ -420,2 +449,43 @@ single-structure design.
 
+### 3.8 An external test on a sixth target
+
+The five benchmark targets are congeneric medicinal-chemistry series against folded pockets,
+assembled by one group for one FEP+ campaign. They may resemble each other more than they
+resemble a new target. We therefore tested the pipeline on a sixth system held entirely outside
+development: 17 peptidomimetics against the anti-haemagglutinin antibody 12CA5, with binding
+measured by biolayer interferometry (7 binders spanning 0.21–77 nM, 10 non-binders). Affinities
+were sealed before any prediction was made.
+
+Every method failed, and several inverted.
+
+**Table 4.** Performance on the sealed 12CA5 panel. Directions were fixed before unblinding.
+
+| method | ρ vs pK_D (7 binders) | binder/non-binder AUC |
+|---|---:|---:|
+| heavy-atom count | **−0.750** | 0.257 |
+| molecular weight | −0.750 | 0.271 |
+| chemical model (pAffinity) | −0.107 | 0.214 |
+| chemical model (rank score) | −0.306 | 0.200 |
+
+Chance AUC is 0.5; every value is below it. The tightest binder (0.21 nM) was ranked 16th of 17
+by every method, and top-5 overlap with the experimental top 5 was zero throughout.
+
+Two findings matter for the present work. First, **the size–affinity relationship reversed
+sign**: heavy-atom count correlates at −0.750 here against +0.49 to +0.61 on the benchmark.
+Size is predictive in congeneric series elaborated into a pocket, and anti-predictive against a
+compact epitope groove where bulk is in the way. A method whose signal is partly size therefore
+carries a sign that is target-dependent and not knowable in advance without measurement — which
+is a further argument for evaluating at fixed size.
+
+Second, **MM/PBSA could not be scored reliably at all**: 3 of 17 built poses carried positive
+van der Waals terms and were excluded by the screen of Section 2.2, one failed to score, and
+across the remainder the score correlated with the pose builder's clash count at ρ = −0.891.
+Those peptidomimetics were threaded onto a template they do not match. The lesson is that
+end-point scoring inherits the quality of its input pose, and that pose provenance must be
+reported alongside any score (Section 4).
+
+We present this as a limitation rather than a general claim: one external target with 17
+compounds cannot establish how the method behaves on all new systems. It does establish that
+benchmark performance on five related targets is not sufficient evidence of transfer.
+
 ---
@@ -424,5 +494,10 @@ single-structure design.
 
-Lipopolysaccharide transport protein A (LptA) is a periplasmic component of the LPS transport
-bridge in Gram-negative bacteria and the target of the antimicrobial peptide thanatin, which
-blocks β-strand complementation between LptA protomers.
+Lipopolysaccharide transport protein A (LptA) is the periplasmic component of the seven-protein
+Lpt bridge that moves lipopolysaccharide from the inner membrane to the cell surface in
+Gram-negative bacteria. The insect-derived antimicrobial peptide thanatin binds the N-terminal
+β-strand of LptA, disrupting both LptA–LptA dimerisation and the LptC–LptA interaction and
+causing LptA degradation.[^17] Because that β-strand edge is also the dimerisation interface,
+the site is a template for designed β-hairpin binders, and thanatin analogues optimised against
+carbapenem-resistant Enterobacteriaceae have been characterised structurally and
+thermodynamically.[^18],[^19]
 
@@ -436,6 +511,7 @@ which raises the question of whether scores read binding or pose-construction ar
 all 213 scored designs the built poses contain **zero** interface heavy-atom contacts below
-1.8 Å, and minimisation drift spans only 0.49–1.42 Å. For comparison, a separate antibody panel
-whose peptidomimetics were threaded onto a non-matching template showed 0–44 clashes per
-complex, backbone RMSD to 4.72 Å, and MM/PBSA scores correlating with clash count at ρ = −0.891
-— that is, reading pose quality rather than chemistry. The LptA poses show no such pathology.
+1.8 Å, and minimisation drift spans only 0.49–1.42 Å. For comparison, the 12CA5 panel of Section 3.8, whose
+peptidomimetics were threaded onto a non-matching template, showed 0–44 clashes per complex,
+backbone RMSD to 4.72 Å, and MM/PBSA scores correlating with clash count at ρ = −0.891 — that
+is, reading pose quality rather than chemistry. The LptA poses show no such pathology, which is
+a necessary condition for the scores to be interpretable, though not a sufficient one.
 
@@ -473,2 +549,4 @@ size measured here (ρ ≈ 0.30, corresponding to ~60% pairwise accuracy):
 
+**Table 5.** Simulated power to detect a positive rank correlation at *p* < 0.05.
+
 | n | power at ρ = 0.30 | ρ = 0.50 | ρ = 0.75 |
@@ -482,5 +560,5 @@ a null result would be uninformative. We state this in advance rather than disco
 afterwards. Ten compounds can, however, do two useful things: detect a catastrophic inversion
-(66% power against ρ = −0.75), which is a live concern given that an independent external
-antibody target showed exactly such a reversal; and supply the ten measured anchors that permit
-within-series calibration.
+(66% power against ρ = −0.75), which is a live concern because exactly such a reversal occurred
+on the external target of Section 3.8; and supply ten measured anchors that permit within-series
+calibration.
 
@@ -515,5 +593,5 @@ tens, then alchemical calculation on the shortlist.
 
-Three cautions generalise beyond this work. The estimator matters: pooling across targets can
-rate a perturbed pose above a crystal pose, so benchmark numbers computed that way should be
-treated sceptically. Pose provenance matters: where poses are built rather than crystallographic,
+Three cautions generalise beyond this work. The estimator matters: pooling compounds drawn from
+different receptor structures can invert the sign of a result, as it does for KRAS (Table 2), so
+benchmark numbers computed that way should be treated sceptically. Pose provenance matters: where poses are built rather than crystallographic,
 scores can track construction artefact — we observed ρ = −0.891 against clash count on a panel
@@ -542,2 +620,18 @@ selections made on one benchmark should not be assumed to hold on another.
 
+## Author contributions
+
+C.L. designed and executed the computational work, performed the analyses and wrote the
+manuscript. J.S.B. supervised the project and revised the manuscript. The prospective
+validation of Section 5 will be executed experimentally by the Brown laboratory.
+
+## Competing interests
+
+The authors declare no competing financial interests.
+
+## Acknowledgements
+
+Computations were performed on local GPU hardware and on the Digital Research Alliance of
+Canada national systems. We thank the authors of the two benchmark studies whose data made this
+evaluation possible.
+
 ## Data and code availability
@@ -546,3 +640,3 @@ All code, computed scores, per-target statistics, bootstrap intervals and figure
 at **https://github.com/Caroldorium/mmpbsa-macrocycle-benchmark** under an MIT licence.
-`scripts/size_matched_pairwise_accuracy.py` reproduces Table 2 and the primary figure.
+`scripts/size_matched_pairwise_accuracy.py` reproduces Table 3 and the primary figure.
 
@@ -589,2 +683,8 @@ computed in this work.
 
+[^17]: Moura ECCM, Baeta T, Romanelli A, et al. Thanatin impairs lipopolysaccharide transport complex assembly by targeting LptC–LptA interaction and decreasing LptA stability. *Front Microbiol* 2020;11:909. doi:10.3389/fmicb.2020.00909. PMID 32477309.
+
+[^18]: Oi KK, Moehle K, Schuster M, Zerbe O. Early molecular insights into thanatin analogues binding to *A. baumannii* LptA. *Molecules* 2023;28(11):4335. doi:10.3390/molecules28114335. PMID 37298811.
+
+[^19]: Sinha S, Dhanabal VB, Sperandeo P, Polissi A, Bhattacharjya S. Linking dual mode of action of host defense antimicrobial peptide thanatin: structures, lipopolysaccharide and LptA binding of designed analogs. *Biochim Biophys Acta Biomembr* 2021;1864(3):183839. doi:10.1016/j.bbamem.2021.183839. PMID 34915021.
+
 
```

---

## v3 — `15ab798` — Manuscript v3: correct three claims against source data

`1 file changed, 37 insertions(+), 20 deletions(-)`

```diff
@@ -20,4 +20,4 @@ that performance reflects physics or molecular size has not been established.
 Size is not a nuisance variable for peptide binding: larger peptides bury more interface and
-genuinely bind more tightly. On a 234-compound, five-target benchmark with measured affinities
-we find that experimental ΔG tracks ligand heavy-atom count with a mean within-target Spearman
+genuinely bind more tightly. On a five-target benchmark of 234 compounds with measured
+affinities, experimental ΔG tracks ligand heavy-atom count with a mean within-target Spearman
 ρ of −0.43, and the MM/PBSA score tracks it more strongly still (ρ = −0.55). Conventional
@@ -165,6 +165,7 @@ wrong energies, and both of which are corrected here:
 
-**Spurious chain breaks.** `pdb4amber` inserted 9–21 `TER` records into continuous chains, each
-of which `tleap` capped, producing complexes with initial energies above 10¹⁷ kcal/mol. `TER`
-records whose flanking C–N distance is below 2.0 Å are deleted. After the fix, atom counts of a
-reference complex match published values exactly.
+**Spurious chain breaks.** `pdb4amber` inserted `TER` records into continuous chains, each of
+which `tleap` capped, producing complexes with initial energies above 10¹⁷ kcal/mol. The problem
+affects **14 of 50 reference complexes, with 9 to 25 spurious records each** (the remaining 36
+are unaffected). `TER` records whose flanking C–N distance is below 2.0 Å are deleted. After the
+fix, atom counts of a reference complex match published values exactly.
 
@@ -321,4 +322,3 @@ comparable and we do not treat the gap as a replication failure.
 Transferring the protocol unchanged to the independent five-target cohort gives a macro
-Spearman ρ of **+0.545** across 234 compounds (Figure 2 of the parity analysis; per-target
-values in Table 1).
+Spearman ρ of **+0.545** across 234 compounds (Figure 3; per-target values in Table 1).
 
@@ -349,4 +349,3 @@ crystal forms.
 
-**Table 2.** KRAS under two estimators. The pooled value is more than three times the
-within-form average, and of opposite sign.
+**Table 2.** KRAS under two estimators. Pooling reverses the sign of the result.
 
@@ -406,5 +405,12 @@ substantially more pairs, confirming the result is not an artefact of the strict
 The skill is also unevenly distributed. PCSK9 (0.765) and JAK2 (0.637) are clearly above
-chance; Cyclin A/B, KRAS and MCL-1 have intervals spanning 0.5. PCSK9 has the most polar
-interface in the cohort, which is consistent with a continuum-solvent treatment adding most
-where desolvation dominates — though with five targets this is a hypothesis, not a finding.
+chance; Cyclin A/B, KRAS and MCL-1 have intervals spanning 0.5.
+
+One candidate explanation is interface polarity, since a continuum-solvent treatment should add
+most where desolvation dominates. PCSK9 does have much the most polar interface in the cohort
+(median buried-area polarity 42.2%, against 30.5% for the next-highest). **But the explanation
+does not survive the second case:** JAK2 is the fourth most polar of the five (23.8%) and the
+second-best performer, while KRAS is second most polar (30.5%) and indistinguishable from
+chance. Polarity alone therefore does not account for where the method works, and with five
+targets we cannot identify what does. We report the pattern without an explanation rather than
+offer one the data does not support.
 
@@ -441,9 +447,20 @@ not be assumed to transfer, and we do not claim ff19SB is superior for these des
 
-On 20 KRAS compounds with both available, ensemble MM/PBSA over 324 frames of 4 × 10 ns
-explicit-solvent MD scores −0.101, against −0.068 for a single minimised frame. The
-replicate-to-replicate spread (5.88 kcal/mol) equals the between-compound spread (5.69).
-Correcting for that noise (reliability 0.733) places the correlation of a perfectly converged
-ensemble at −0.118. The limiting error is bias in the energy function, not variance from short
-sampling — consistent with the general finding for MM/PBSA[^3] and with the source protocol's
-single-structure design.
+A natural objection to single-structure scoring is that an ensemble average would do better. We
+tested this on 40 KRAS compounds for which four independent 10 ns explicit-solvent replicas were
+available, scored over the 2–10 ns window at 81 frames per replica — **324 frames per compound**.
+
+The ensemble does not converge well enough for the averaging to help. The median
+replicate-to-replicate standard deviation of the ensemble score is **6.81 kcal/mol**, while the
+standard deviation *between* compounds is **5.74 kcal/mol**. The noise from which replicas you
+happened to run is therefore larger than the entire signal distinguishing one compound from
+another, a ratio of 1.19. Splitting by panel, the 20 confirmatory compounds give 5.60 against
+5.54 (ratio 1.01) and the 20 primary compounds 8.70 against 4.88 (ratio 1.78); the conclusion
+does not depend on the subset.
+
+Four 10 ns replicas is a substantial sampling budget for a system of this size, and it is not
+enough to make the ensemble score reproducible at the resolution the ranking requires. This is
+consistent with the general finding that MM/PBSA correlation does not improve reliably with
+longer sampling because the limiting error is bias rather than variance,[^3] and it supports the
+source protocol's single-structure design on pragmatic grounds: a cheaper calculation that is
+reproducible is preferable to an expensive one that is not.
 
```

---

## v4 — `78fe0a1` — Manuscript v4: correct the LptA design count, 308 -> 212

`1 file changed, 13 insertions(+), 9 deletions(-)`

```diff
@@ -42,3 +42,3 @@ not sufficient evidence of transfer.
 
-We apply the protocol to 308 designed macrocycles against the lipopolysaccharide transport
+We apply the protocol to 212 designed macrocycles against the lipopolysaccharide transport
 protein LptA and register a prospective test with its power analysis stated in advance.
@@ -520,6 +520,8 @@ thermodynamically.[^18],[^19]
 
-We applied the validated protocol to 308 designed macrocycles in two branches: 212 de novo
-designs from generative backbone design followed by sequence design, and 96 single-residue
+We applied the protocol to a panel of 212 designed macrocycles in three branches: 20 canonical
+head-to-tail cyclic peptides from generative backbone design followed by sequence design; 96
+single-residue non-canonical substitutions on those canonical parents; and 96 single-residue
 non-canonical substitutions on a thanatin derivative-7 scaffold templated on a co-crystal
-structure.
+structure. A total of 213 complexes were scored — the 212 designs plus the derivative-7 parent
+as a reference.
 
@@ -527,5 +529,6 @@ structure.
 which raises the question of whether scores read binding or pose-construction artefact. Across
-all 213 scored designs the built poses contain **zero** interface heavy-atom contacts below
-1.8 Å, and minimisation drift spans only 0.49–1.42 Å. For comparison, the 12CA5 panel of Section 3.8, whose
-peptidomimetics were threaded onto a non-matching template, showed 0–44 clashes per complex,
+all 213 scored complexes the built poses contain **zero** interface heavy-atom contacts below
+1.8 Å, and minimisation drift spans only 0.49–1.42 Å. For comparison, the 12CA5 panel of
+Section 3.8, whose peptidomimetics were threaded onto a non-matching template, showed 0–44
+clashes per complex,
 backbone RMSD to 4.72 Å, and MM/PBSA scores correlating with clash count at ρ = −0.891 — that
@@ -537,4 +540,5 @@ a necessary condition for the scores to be interpretable, though not a sufficien
 **−0.291** for the 96 thanatin swaps. This follows from construction: single substitutions on a
-fixed 16-residue scaffold barely change size, whereas de novo designs vary from 14 to 18
-residues. The thanatin branch is therefore the least size-confounded.
+fixed 16-residue scaffold barely change size, whereas the de novo designs vary from 14 to 18
+residues. The thanatin branch is therefore the least size-confounded, and by the argument of
+Section 3.4 it is the branch on which the protocol's demonstrated skill is most likely to apply.
 
```

---

## v5 — `0e055f8` — Manuscript v5: reconcile the discussion with the corrected results

`1 file changed, 23 insertions(+), 8 deletions(-)`

```diff
@@ -614,8 +614,18 @@ tens, then alchemical calculation on the shortlist.
 
-Three cautions generalise beyond this work. The estimator matters: pooling compounds drawn from
-different receptor structures can invert the sign of a result, as it does for KRAS (Table 2), so
-benchmark numbers computed that way should be treated sceptically. Pose provenance matters: where poses are built rather than crystallographic,
-scores can track construction artefact — we observed ρ = −0.891 against clash count on a panel
-where templating failed. And force-field choice did not transfer between our two datasets, so
-selections made on one benchmark should not be assumed to hold on another.
+Three cautions generalise beyond this work.
+
+**The estimator matters.** Pooling compounds drawn from different receptor structures can invert
+the sign of a result, as it does for KRAS (Table 2), so benchmark numbers computed that way
+should be treated sceptically — including when the pooling is hidden inside what is reported as
+a within-family metric.
+
+**Pose provenance matters.** Where poses are built rather than crystallographic, scores can
+track construction artefact: on the 12CA5 panel the score correlated with the pose builder's
+clash count at ρ = −0.891. Any end-point score should be reported alongside evidence about the
+geometry it was computed on.
+
+**Reference standards must themselves be validated.** We could not settle whether force-field
+choice transfers between our datasets, because the only available reference on the design panel
+was a computational score that we measured at ρ = −0.149 against experiment (Section 3.6).
+Agreement with an unvalidated reference is not evidence.
 
@@ -625,4 +635,5 @@ selections made on one benchmark should not be assumed to hold on another.
 
-1. **Five targets.** Macro intervals rest on five values. The concentration of skill in PCSK9
-   and JAK2 may be chemistry, interface polarity, or chance; we cannot distinguish these.
+1. **Five targets.** Macro intervals rest on five values. Skill concentrates in PCSK9 and JAK2,
+   and we could not identify why: interface polarity is the obvious candidate and it fails
+   (Section 3.4). With five targets we cannot distinguish a chemical explanation from chance.
 2. **Cyclin A/B is underpowered** in the strictest size-matched band (15 pairs, CI spanning
@@ -638,2 +649,6 @@ selections made on one benchmark should not be assumed to hold on another.
    set toward the measurable dynamic range.
+8. **We evaluate the protocol as published**, not an optimised variant. We did not re-tune the
+   dielectric, solvation model or minimisation depth to maximise size-matched accuracy; doing so
+   on five targets would overfit. The reported 60.1% is therefore a floor for this family of
+   methods rather than the best achievable.
 
```

---

## v6 — `7072666` — Manuscript v6: correct the bootstrap, pairs are not independent

`1 file changed, 66 insertions(+), 37 deletions(-)`

```diff
@@ -25,8 +25,14 @@ benchmarking therefore cannot separate a scoring function from a ruler.
 
-We resolve this by evaluating only within sets of peptides of **identical heavy-atom count**,
-where a size descriptor carries no information by construction. On decisive pairs, MM/PBSA
-selects the tighter binder in **60.1% of cases (95% CI 53.0–69.3)**, against **84.1%
-(75.6–91.2)** for published FEP+ predictions on the same compounds and 50% for chance. The
-skill is real but modest and unevenly distributed: clearly above chance for PCSK9 (76.5%) and
-JAK2 (63.7%), indistinguishable from chance for Cyclin A/B, KRAS and MCL-1.
+We resolve this by evaluating only within sets of peptides matched on heavy-atom count, where a
+size descriptor carries little or no information by construction, and asking how often a method
+picks the tighter binder of a decisive pair. Within a three-atom window MM/PBSA succeeds in
+**62.4% of pairs (95% CI 55.5–68.9)**, against 50% for chance and **83.9% (78.2–88.6)** for
+published FEP+ predictions on the same compounds. At exactly equal heavy-atom count, where a
+size descriptor is wholly uninformative, the estimate is 60.1% but the interval widens to
+47.2–72.0 (posterior mass above chance 0.93) because few same-size decisive pairs exist.
+Intervals resample compounds rather than pairs, since pairs drawn from a compound set are not
+independent; resampling pairs understates the width by roughly half.
+
+The skill is modest and unevenly distributed: only PCSK9 (75.7%) is individually above chance,
+with Cyclin A/B, JAK2, KRAS and MCL-1 each spanning it.
 
@@ -284,5 +290,9 @@ better-powered secondary.
 
-**Uncertainty.** Bootstrap over pairs within each target (2,000 resamples) and over targets for
-macro values (5,000 resamples). Compounds within a receptor group are not independent, so where
-Spearman correlations are reported the bootstrap resamples **receptor groups**, not compounds.
+**Uncertainty.** All pairwise intervals resample **compounds**, not pairs. This matters: the
+*n*(*n*−1)/2 pairs formed from *n* compounds are not independent draws, because every compound
+appears in many pairs. Resampling pairs treats them as independent and badly understates the
+uncertainty — on the exact-match band it returns [0.530, 0.693] where the correct
+compound-level interval is [0.472, 0.720]. We resample compounds within each target (2,000
+draws), rebuild the pair set from each resample, and combine targets for macro values (5,000
+draws). Where Spearman correlations are reported the bootstrap resamples **receptor groups**.
 
@@ -386,22 +396,34 @@ within each target to make a fair comparison possible.
 
-Restricting to pairs of peptides with identical heavy-atom count, where a size descriptor is
-undefined:
-
 **Table 3.** Pairwise accuracy on decisive pairs (|ΔΔG| ≥ 0.5 kcal/mol). Chance = 0.500.
+Intervals are 95% bootstrap over compounds (Section 2.6). The left block is the strictest
+matching; the right admits a three-atom window and is better powered.
 
-| target | pairs | MM/PBSA | 95% CI | FEP+ | 95% CI |
+| | | identical heavy atoms | | within 3 heavy atoms | |
 |---|---:|---:|---|---:|---|
-| Cyclin A/B | 15 | 0.533 | [0.267, 0.800] | 0.867 | [0.667, 1.000] |
-| JAK2 | 80 | **0.637** | [0.537, 0.738] | 0.900 | [0.838, 0.963] |
-| KRAS | 56 | 0.518 | [0.393, 0.643] | 0.818 | [0.709, 0.909] |
-| MCL-1 | 71 | 0.549 | [0.437, 0.662] | 0.671 | [0.557, 0.786] |
-| PCSK9 | 115 | **0.765** | [0.687, 0.835] | 0.947 | [0.904, 0.982] |
-| **macro** | — | **0.601** | **[0.530, 0.693]** | **0.841** | [0.756, 0.912] |
-
-The macro interval excludes chance, so **the skill is real**. Relaxing the matching to within
-three heavy atoms gives 0.624 [0.552, 0.693] for MM/PBSA and 0.839 [0.799, 0.886] for FEP+, on
-substantially more pairs, confirming the result is not an artefact of the strictest band.
-
-The skill is also unevenly distributed. PCSK9 (0.765) and JAK2 (0.637) are clearly above
-chance; Cyclin A/B, KRAS and MCL-1 have intervals spanning 0.5.
+| target | | MM/PBSA | 95% CI | MM/PBSA | 95% CI |
+| Cyclin A/B | | 0.533 | [0.053, 1.000] | 0.674 | [0.454, 0.863] |
+| JAK2 | | 0.637 | [0.397, 0.831] | 0.598 | [0.437, 0.732] |
+| KRAS | | 0.518 | [0.289, 0.727] | 0.496 | [0.339, 0.664] |
+| MCL-1 | | 0.549 | [0.328, 0.757] | 0.594 | [0.444, 0.733] |
+| **PCSK9** | | **0.765** | **[0.613, 0.882]** | **0.757** | **[0.671, 0.830]** |
+| **macro** | | **0.601** | [0.472, 0.720] | **0.624** | **[0.555, 0.689]** |
+| FEP+ macro | | 0.841 | [0.730, 0.912] | 0.839 | [0.782, 0.886] |
+
+**The two bands give different levels of confidence, and we report both rather than the more
+favourable one.** At exact size matching the macro estimate is 0.601 but its interval includes
+chance; the posterior mass above 0.5 is 0.93. There are only 337 decisive pairs in that band,
+drawn from few enough compounds that the compound-level interval is wide. At a three-atom
+window — still tight relative to a 25-to-127 heavy-atom cohort, and still a regime in which a
+size descriptor achieves only 0.602 — the macro estimate is 0.624 with an interval that
+excludes chance and posterior mass above 0.5 exceeding 0.999.
+
+The honest summary is therefore: **the skill is established at near-matched size and suggestive
+but not conclusive at exactly matched size**, the latter limited by how few same-size decisive
+pairs the cohort contains rather than by the size of the effect.
+
+The skill is also unevenly distributed, and under compound-level resampling **only PCSK9 is
+individually above chance** — at both bands. JAK2's point estimate of 0.637 at exact matching
+has an interval of [0.397, 0.831] and cannot be separated from chance on its own; the same is
+true of Cyclin A/B, KRAS and MCL-1. The macro result at the three-atom band is therefore driven
+by PCSK9 with weak support from the others, not by uniform performance.
 
@@ -594,5 +616,7 @@ we do not claim it here.
 The central finding is that end-point MM/PBSA retains genuine ranking skill on macrocyclic
-peptides when size is held fixed — 60.1% pairwise accuracy where chance is 50% and a size
-descriptor is undefined — at roughly a thousandth of the cost of rigorous free-energy
-simulation, which achieves 84.1% on the same pairs.
+peptides when size is controlled — 62.4% pairwise accuracy within a three-atom window, where
+chance is 50% and a heavy-atom count reaches only 60.2% — at roughly a thousandth of the cost of
+rigorous free-energy simulation, which achieves 83.9% on the same pairs. At exactly matched size
+the point estimate is similar (60.1%) but the cohort does not contain enough same-size decisive
+pairs to separate it from chance at 95% confidence.
 
@@ -638,14 +662,19 @@ Agreement with an unvalidated reference is not evidence.
    (Section 3.4). With five targets we cannot distinguish a chemical explanation from chance.
-2. **Cyclin A/B is underpowered** in the strictest size-matched band (15 pairs, CI spanning
-   0.267–0.800) and contributes little.
-3. **Scores are not free energies** (Section 3.5) and must not be quoted as predicted affinities.
-4. **No measured affinities exist for the LptA designs.** Section 4 is an application, not a
+2. **The strictest band is underpowered.** At exactly equal heavy-atom count the cohort yields
+   337 decisive pairs in total, and the macro interval [0.472, 0.720] includes chance. Cyclin
+   A/B is the extreme case, with 15 pairs and an interval spanning [0.053, 1.000]. The
+   three-atom band is the better-powered analysis and the one we rely on.
+3. **Only one target is individually significant.** PCSK9 excludes chance at both bands; the
+   other four do not. The macro result rests substantially on PCSK9, and a cohort without a
+   polar-interface target might not reproduce it.
+4. **Scores are not free energies** (Section 3.5) and must not be quoted as predicted affinities.
+5. **No measured affinities exist for the LptA designs.** Section 4 is an application, not a
    validation, and Section 5 is a plan.
-5. **Two selenium-containing non-canonical residues** are approximations, GAFF2 having no
+6. **Two selenium-containing non-canonical residues** are approximations, GAFF2 having no
    selenium parameters.
-6. **Single-structure scoring** inherits whatever error the input pose carries. We verify pose
+7. **Single-structure scoring** inherits whatever error the input pose carries. We verify pose
    quality explicitly (Section 4) but cannot certify that generated poses are correct.
-7. **Censored measurements** are excluded from pairwise analyses, which may bias the admitted
+8. **Censored measurements** are excluded from pairwise analyses, which may bias the admitted
    set toward the measurable dynamic range.
-8. **We evaluate the protocol as published**, not an optimised variant. We did not re-tune the
+9. **We evaluate the protocol as published**, not an optimised variant. We did not re-tune the
    dielectric, solvation model or minimisation depth to maximise size-matched accuracy; doing so
```

---

## v7 — `04461e8` — Manuscript v7: methods-first restructure, with the LptA campaign as the stated purpose

`1 file changed, 215 insertions(+), 51 deletions(-)`

```diff
@@ -1,2 +1,6 @@
-# Size-matched benchmarking reveals genuine but modest ranking skill in end-point MM/PBSA for macrocyclic peptides
+# Establishing an end-point scoring protocol for macrocyclic peptide design: size-matched benchmarking, its limits, and a registered prospective test
+
+*Alternative title, results-forward: "Size-matched benchmarking reveals genuine but modest
+ranking skill in end-point MM/PBSA for macrocyclic peptides." The methods-forward title above
+reflects the paper's purpose; either works.*
 
@@ -48,4 +52,8 @@ not sufficient evidence of transfer.
 
-We apply the protocol to 212 designed macrocycles against the lipopolysaccharide transport
-protein LptA and register a prospective test with its power analysis stated in advance.
+This method was established to serve an active design campaign: 212 macrocyclic peptides
+designed against the lipopolysaccharide transport protein LptA, from which candidates will be
+selected for synthesis and experimental evaluation. We report the protocol applied to that
+panel, and register the prospective test with its power analysis stated in advance. The purpose
+of the validation reported here is to establish what the method can and cannot support
+**before** committing synthetic resources, not to report a finished prediction.
 
@@ -72,2 +80,11 @@ something the chemist already knows.
 
+**This work has a specific purpose.** We maintain an active design campaign against the
+lipopolysaccharide transport protein LptA, comprising 212 designed macrocycles from which a
+subset will be synthesised and assayed. Non-canonical building blocks cost on the order of a
+thousand dollars each, so the selection is expensive and effectively irreversible. We therefore
+set out to establish, in advance of that commitment, exactly what an end-point scoring protocol
+can and cannot support on this chemistry — which is a methods question, and is the subject of
+this paper. The design campaign is described in Section 4 and the registered prospective test
+in Section 5; the validation that precedes them is what licenses either.
+
 ### 1.2 The computational options
@@ -200,24 +217,65 @@ never receives it. This is reported rather than hidden.
 
-**Non-canonical residues.** The LptA application required parameters for 170 non-standard
-residues. For each: the polymer-form rotamer is capped with ACE and NME, charges are derived
-with AM1-BCC on the neutral capped dipeptide, `prepgen` removes the caps and redistributes
-their charge to give a residue netting zero, and `parmchk2` supplies missing GAFF2 terms.
-Capping matters because AM1-BCC on a free amino acid fits backbone charges in an environment
-(–NH₂ / –COOH) the residue never occupies; capping places the fit in a peptide-like context, the
-standard Amber construction.[^13]
-
-The backbone is reassigned to protein atom types while the side chain remains GAFF2, because
-`antechamber` otherwise types the backbone as GAFF2 and `tleap` has no parameters across the
-peptide bond to a neighbouring standard residue. The alpha carbon is typed CT for ff03 and CX
-for ff14SB/ff19SB, so two parameter sets are maintained.
-
-Every residue is validated by constructing ALA–X–ALA and requiring zero missing parameters;
-all 170 pass under both typings. This test is necessary because `tleap` writes a topology even
-when parameters are missing, substituting defaults.
+#### 2.3.1 Non-canonical residue parameterisation
+
+Designed peptides carrying non-standard residues cannot be scored until each residue has Amber
+parameters. The LptA application required **170 non-standard residues** (165 unique
+non-canonical codes plus five residues constituting the thanatin derivative-7 scaffold, none of
+which `tleap` knows). Because this step gates everything downstream and is where most silent
+failures occur, we specify it in full.
+
+**Protocol, per residue.**
+
+1. Take the polymer-form rotamer and the corresponding chemistry definition from the source
+   library. Atom names and bond orders are read from the definition file, not perceived from
+   coordinates; coordinates are matched positionally.
+2. Cap with ACE at the N-terminus and NME at the C-terminus, giving a neutral dipeptide. Cap
+   atoms are placed by bisector construction with a hard clash check, and cap atom names are
+   chosen per residue to avoid collisions with side-chain names.
+3. Derive AM1-BCC charges on the capped dipeptide with `antechamber`, under element-safe
+   placeholder atom names (see defect 1 below).
+4. Restore the true atom names, case-insensitively.
+5. Reassign the backbone to protein atom types (N, H, CA, HA, C, O) while the side chain
+   remains GAFF2. Cap atoms are retyped as protein types as well, because they stand in for
+   neighbouring residues and leaving them GAFF2-typed makes `parmchk2` generate cross terms
+   against the wrong types.
+6. Run `prepgen` with a mainchain file (`HEAD_NAME N`, `TAIL_NAME C`, `OMIT_NAME` for the caps,
+   `CHARGE 0.0`) to remove the caps and redistribute their charge, giving a residue that nets
+   zero.
+7. Run `parmchk2` for missing GAFF2 terms, then strip from the resulting frcmod any entry that
+   would redefine a protein backbone type.
 
-Six defects were found during parameterisation, four of which produced silently wrong output.
-Most consequential: `antechamber` infers element from atom name, so every alpha carbon (`CA`)
-was read as calcium and 61 aromatic carbons named `CL*` as chlorine; and GAFF2 has no selenium,
-so `antechamber` silently retyped Se as sulfur while `parmchk2` reported nothing missing. Two
-selenium-containing residues remain approximations and are flagged as such.
+Capping matters because AM1-BCC on a free amino acid fits backbone charges in an environment
+(–NH₂ / –COOH) the residue never occupies; capping places the charge fit in a peptide-like
+context, which is the standard Amber construction for non-standard residues.[^13] Backbone
+retyping matters because `antechamber` otherwise types the backbone as GAFF2 and `tleap` then
+has no parameters across the peptide bond to a neighbouring standard residue
+(`Could not find bond parameter for atom types: C - ns`). The alpha carbon is typed **CT** under
+ff03/ff99 and **CX** under ff14SB/ff19SB — ff03 has no CX at all — so two complete parameter
+sets are maintained and a prep typed for one force field fails on the other.
+
+**Validation gate.** Each residue is tested by constructing ALA–X–ALA in `tleap` and requiring
+**zero missing parameters**. All 170 pass under both backbone typings. This gate is not
+optional: `tleap` writes a topology even when parameters are missing, substituting defaults
+silently, so "the build succeeded" is not evidence of correctness.
+
+**Failure modes.** Six defects were found, four of which produced plausible output rather than
+an error. We tabulate them because every one is a property of the standard toolchain rather
+than of our library, and will recur for anyone attempting this.
+
+**Table 1.** Silent and overt failure modes in non-canonical residue parameterisation.
+
+| # | defect | how it manifests | detection | remedy |
+|---:|---|---|---|---|
+| 1 | `antechamber` infers element from atom name | every alpha carbon (`CA`) read as **calcium**; 61 aromatic carbons named `CL*` read as **chlorine**; `FE`, `NI`, `CO` likewise | element census against the chemistry definition | build under element-safe placeholder names, restore afterwards |
+| 2 | bond orders perceived from geometry | fails on hypervalent sulfur (sulfonamides) and iodoaryls | valence check | supply bond orders from the definition via MOL2 |
+| 3 | fixed-dihedral cap placement | NME nitrogen placed 0.15 Å from the residue's own carbonyl oxygen whenever ψ ≈ 180° | inter-atomic distance check after capping | bisector placement with a clash threshold |
+| 4 | cap atom names collide with side-chain names | **75 of 170** residues affected | name-set intersection | per-residue cap naming |
+| 5 | library-internal name inconsistency | nine halogenated residues: the rotamer PDB truncates `CLM1` to `CLM`, and one residue names *both* chlorines `CLM`. Unmatched atoms are dropped silently by `prepgen`, surfacing later as a `nan` internal coordinate | cross-check names between definition and coordinates | take names and bonds from the definition, match coordinates positionally |
+| 6 | case-sensitive name restoration | `antechamber` rewrites a two-letter placeholder to proper element case (`CL1` → `Cl1`), so a case-sensitive mapping leaves halogens unnamed | case-insensitive round-trip test | case-insensitive restoration |
+
+**A limitation we could not remove.** GAFF2 has no selenium and `sqm` has no AM1 parameters for
+it, so `antechamber` retypes selenium as sulfur and `parmchk2` reports nothing missing — the
+build appears clean. Two residues are affected and are flagged `ok_element_unsupported` rather
+than shipped as validated. Any score involving them is an approximation, and we say so wherever
+such a residue appears.
 
@@ -284,7 +342,40 @@ data.
 **Size-matched pairwise accuracy — the primary measure.** For each target we enumerate all
-compound pairs whose measured affinities differ by at least 0.5 kcal/mol ("decisive" pairs) and
-whose ligand heavy-atom counts are equal. For each pair we ask whether the method orders them
-correctly. Chance is 0.5. A size descriptor is undefined on these pairs, because the sizes are
-equal. We report the band |Δheavy| = 0 as the primary analysis and |Δheavy| ≤ 3 as a
-better-powered secondary.
+compound pairs whose measured affinities differ by at least a margin *m* ("decisive" pairs) and
+whose ligand heavy-atom counts differ by at most a window *w*. For each pair we ask whether the
+method orders them correctly. Chance is 0.5. At *w* = 0 a size descriptor is undefined, because
+there is no size difference to read.
+
+Two free parameters enter here, and neither was chosen after seeing the result.
+
+*The affinity margin m.* Pairs whose measured affinities are nearly equal cannot discriminate
+between methods and only add noise, so a margin is required. We set *m* = 0.5 kcal/mol a priori
+because it is approximately the reproducibility of the assays contributing to the cohort, and
+report the full sensitivity in Table 2 rather than a single value. The conclusion is stable:
+MM/PBSA is above chance at every margin from 0 to 1.5 kcal/mol.
+
+*The size window w.* We report *w* = 0, the strictest possible matching, and *w* = 3, which
+admits enough pairs for the macro interval to be informative. Both are reported throughout; we
+do not select between them on the basis of which is more favourable.
+
+**Table 2.** Sensitivity to the decisive-pair margin, at *w* = 3. "Edge" is MM/PBSA minus the
+heavy-atom control on the identical pair set. Intervals resample compounds.
+
+| margin (kcal/mol) | pairs | MM/PBSA | 95% CI | heavy atoms | edge | FEP+ |
+|---:|---:|---:|---|---:|---:|---:|
+| 0.00 | 2384 | 0.597 | [0.542, 0.649] | 0.585 | +0.012 | 0.769 |
+| 0.25 | 2002 | 0.613 | [0.554, 0.674] | 0.591 | +0.022 | 0.805 |
+| **0.50** | **1645** | **0.624** | **[0.553, 0.691]** | **0.602** | **+0.022** | **0.839** |
+| 0.75 | 1284 | 0.640 | [0.566, 0.718] | 0.622 | +0.018 | 0.888 |
+| 1.00 | 1016 | 0.662 | [0.562, 0.750] | 0.590 | +0.072 | 0.921 |
+| 1.50 | 666 | 0.697 | [0.576, 0.804] | 0.599 | +0.098 | 0.932 |
+| 2.00 | 452 | 0.620 | [0.468, 0.775] | 0.420 | +0.200 | 0.956 |
+
+Two features of Table 2 deserve comment, and one of them is unfavourable. Accuracy rises
+monotonically with the margin up to 1.5 kcal/mol for both MM/PBSA and FEP+, which is the
+behaviour of a real but noisy signal: pairs that are further apart in affinity are easier to
+order. **But at the chosen margin the edge over a heavy-atom count is only +0.022.** It becomes
+substantial (+0.07 to +0.20) only for pairs separated by 1 kcal/mol or more. The method's
+advantage over the trivial descriptor within a three-atom window is therefore real but small,
+and concentrated in the pairs that are easiest to tell apart. At *w* = 0, where the descriptor
+has no information at all, the comparison does not arise.
 
@@ -313,7 +404,58 @@ and excluded from pairwise analyses where the ordering would be ambiguous.
 
-### 2.7 Software
+### 2.7 Consolidated parameters
 
-AmberTools 25/26 (`tleap`, `antechamber`, `parmchk2`, `prepgen`, `cpptraj`, `MMPBSA.py`);
-OpenMM 8 (CUDA, mixed precision); RDKit 2025.03/2026.03; ParmEd; NumPy, SciPy, pandas; UCSF
-ChimeraX 1.12 for visualisation. Hardware: NVIDIA RTX PRO 6000 Blackwell and NVIDIA H100 80GB.
+Every setting required to reproduce the reported numbers, in one place.
+
+| stage | parameter | value |
+|---|---|---|
+| preparation | `pdb4amber` passes | two; first retains prepared hydrogens, second strips with `-y` |
+| | spurious `TER` removal | deleted where flanking C–N < 2.0 Å |
+| | disulfide detection | geometric, S–S < 2.5 Å |
+| | quality screen | exclude any complex with VDWAALS > 0 |
+| force field | protein | `oldff/leaprc.ff03` (comparison arms: ff14SB, ff19SB) |
+| | ligand / ncAA side chain | GAFF2, AM1-BCC charges |
+| | phosphotyrosine | `leaprc.phosaa10` |
+| | ncAA alpha-carbon type | CT under ff03/ff99; CX under ff14SB/ff19SB |
+| minimisation | engine | OpenMM 8, CUDA, mixed precision |
+| | implicit solvent | OBC1 (`igb=2` equivalent), solute dielectric 1.0, solvent 80.0 |
+| | cutoff | none (`NoCutoff`, equivalent to `cut=1000`) |
+| | optimiser / cap | L-BFGS, **500 iterations** (calibrated to sander drift, §2.4) |
+| | reference settings | `imin=1, maxcyc=5000, ncyc=2000, cut=1000.0, igb=2` |
+| scoring | engine | `MMPBSA.py`, single-trajectory |
+| | solver | Poisson–Boltzmann, internal PBSA |
+| | interior dielectric | **2.0** (arms at 1.0 and 4.0 also computed) |
+| | exterior dielectric | 80.0 |
+| | nonpolar model | LCPO, `inp=1`, `cavity_surften=0.0072`, `cavity_offset=0` |
+| | entropy | **not computed** |
+| | dielectric injection | `use_sander=1` → `-make-mdins` → edit `_MM_gb.mdin` → `-use-mdins` |
+| statistics | primary estimator | Spearman ρ within receptor group, equal weight per target |
+| | decisive-pair margin | 0.5 kcal/mol (sensitivity, Table 2) |
+| | size windows | 0 and 3 heavy atoms |
+| | bootstrap unit | **compounds** (pairwise) or **receptor groups** (correlations) |
+| | resamples | 2,000 within target; 5,000 for macro |
+
+### 2.8 Software, hardware and code map
+
+AmberTools 25 and 26 (`tleap`, `antechamber`, `parmchk2`, `prepgen`, `cpptraj`, `MMPBSA.py`);
+OpenMM 8; RDKit 2025.03 and 2026.03; ParmEd; NumPy, SciPy, pandas; UCSF ChimeraX 1.12 for
+visualisation. Hardware: NVIDIA RTX PRO 6000 Blackwell (local) and NVIDIA H100 80GB (Digital
+Research Alliance of Canada).
+
+Each methods step maps to a specific script in the deposited repository:
+
+| section | step | script |
+|---|---|---|
+| 2.2–2.5 | full protocol, one complex | `scripts/zha_pipeline.py` |
+| 2.1, 3.1 | protocol applied to the five-target cohort | `scripts/zha_cohort.py`, `run_cohort.sh` |
+| 2.2 | build-quality audit (TER, clashes, atom counts) | `scripts/check_build.py` |
+| 2.3 | non-canonical residue parameterisation | `scripts/ncaa_param.py`, `run_ncaa_all.sh` |
+| 2.4 | minimisation and the drift calibration | `scripts/min_openmm.py`, `compare_minimisers.py` |
+| 2.6, 3.4 | size-matched pairwise accuracy, Tables 3 and 6, Figure 4 | `scripts/size_matched_pairwise_accuracy.py` |
+| 3.3 | size-confound analysis, Figure 2 | `scripts/fig_size_confound.py` |
+| 3.1–3.2 | per-target statistics and estimator comparison | `scripts/analyse_cohort.py` |
+| 3.8 | external 12CA5 panel | `scripts/stage_12ca5.py`, `build_12ca5_poses.py`, `run_12ca5.sh` |
+| 4 | LptA design panel | `run_lpta_canonical.sh`, `run_lpta_ncaa.sh` |
+
+Scripts use `$PROJECT_ROOT`, `$DATA_ROOT`, `$CONDA_ENV`, `$GPU_HOST`, `$CLUSTER_HOST` and
+`$CLUSTER_USER` placeholders where local infrastructure paths were removed.
 
@@ -332,5 +474,5 @@ comparable and we do not treat the gap as a replication failure.
 Transferring the protocol unchanged to the independent five-target cohort gives a macro
-Spearman ρ of **+0.545** across 234 compounds (Figure 3; per-target values in Table 1).
+Spearman ρ of **+0.545** across 234 compounds (Figure 3; per-target values in Table 3).
 
-**Table 1.** Per-target performance on the Awoonor-Williams cohort. Pearson *r* and Spearman ρ
+**Table 3.** Per-target performance on the Awoonor-Williams cohort. Pearson *r* and Spearman ρ
 are computed within receptor group. "Slope" is the per-target regression slope of score on
@@ -347,3 +489,3 @@ experiment; its reciprocal gives the factor by which ΔΔG spread is inflated.
 
-Two things in Table 1 deserve emphasis. First, the commonly quoted "R" in such analyses is a
+Two things in Table 3 deserve emphasis. First, the commonly quoted "R" in such analyses is a
 correlation, not a coefficient of determination: *r* = 0.548 corresponds to *R*² = 0.318, so
@@ -359,3 +501,3 @@ crystal forms.
 
-**Table 2.** KRAS under two estimators. Pooling reverses the sign of the result.
+**Table 4.** KRAS under two estimators. Pooling reverses the sign of the result.
 
@@ -396,3 +538,3 @@ within each target to make a fair comparison possible.
 
-**Table 3.** Pairwise accuracy on decisive pairs (|ΔΔG| ≥ 0.5 kcal/mol). Chance = 0.500.
+**Table 5.** Pairwise accuracy on decisive pairs (|ΔΔG| ≥ 0.5 kcal/mol). Chance = 0.500.
 Intervals are 95% bootstrap over compounds (Section 2.6). The left block is the strictest
@@ -447,3 +589,3 @@ scoring function does not.
 With entropy omitted, the raw score is not a free energy. It sits 45–75 kcal/mol from
-experimental ΔG (Table 1), and matching the experimental scale requires fitting a slope and
+experimental ΔG (Table 3), and matching the experimental scale requires fitting a slope and
 intercept **within a target**, which requires compounds already measured in that series. An
@@ -499,3 +641,3 @@ Every method failed, and several inverted.
 
-**Table 4.** Performance on the sealed 12CA5 panel. Directions were fixed before unblinding.
+**Table 6.** Performance on the sealed 12CA5 panel. Directions were fixed before unblinding.
 
@@ -531,3 +673,3 @@ benchmark performance on five related targets is not sufficient evidence of tran
 
-## 4. Application: a de novo design campaign against LptA
+## 4. The campaign this method serves: de novo macrocycles against LptA
 
@@ -542,3 +684,4 @@ thermodynamically.[^18],[^19]
 
-We applied the protocol to a panel of 212 designed macrocycles in three branches: 20 canonical
+The protocol above was established in order to triage this panel. We applied it to 212 designed
+macrocycles in three branches: 20 canonical
 head-to-tail cyclic peptides from generative backbone design followed by sequence design; 96
@@ -573,11 +716,13 @@ produces a spurious correlation.
 
-## 5. Prospective validation (registered plan — not yet executed)
+## 5. Committed prospective test (registered; synthesis pending)
 
-*This section describes a planned experiment. No results exist at the time of writing.*
+*This section describes a committed experiment that has not yet been executed. No results exist
+at the time of writing, and none are claimed.*
 
-Everything above is retrospective, on compounds measured by others. We therefore register a
-forward test.
+Everything above is retrospective, on compounds measured by others. The campaign of Section 4
+will proceed to synthesis, and we register the forward test here, before any compound is made,
+so that the prediction cannot be adjusted after the fact.
 
-**Design.** Ten de novo LptA designs will be selected before any measurement, synthesised, and
-assayed. The ranking is frozen at registration and will not be adjusted afterwards.
+**Design.** Ten de novo LptA designs, selected and frozen before synthesis, will be made and
+assayed. The ranking is registered at the point of selection and will not be revised.
 
@@ -592,3 +737,3 @@ size measured here (ρ ≈ 0.30, corresponding to ~60% pairwise accuracy):
 
-**Table 5.** Simulated power to detect a positive rank correlation at *p* < 0.05.
+**Table 7.** Simulated power to detect a positive rank correlation at *p* < 0.05.
 
@@ -631,2 +776,11 @@ control informative by construction.
 
+**What this licenses us to do.** The protocol was established to decide which LptA designs to
+synthesise, and the validation above determines how far that decision can lean on it. Three
+constraints follow directly. Ranking is admissible within a branch, because that is what was
+tested; cross-branch ranking is not, because the branches occupy different score ranges.
+Size-matched selection is preferred, because that is the regime in which skill was demonstrated
+— which favours the thanatin branch and, within the de novo branch, a size-matched subset.
+And the selection must span the predicted range rather than take the top of it, because a set
+of predicted winners cannot test the ranking that chose them. Section 5 implements all three.
+
 The practical reading for a design campaign is that this protocol is a **triage filter, not a
@@ -641,3 +795,3 @@ Three cautions generalise beyond this work.
 **The estimator matters.** Pooling compounds drawn from different receptor structures can invert
-the sign of a result, as it does for KRAS (Table 2), so benchmark numbers computed that way
+the sign of a result, as it does for KRAS (Table 4), so benchmark numbers computed that way
 should be treated sceptically — including when the pooling is hidden inside what is reported as
@@ -655,2 +809,12 @@ Agreement with an unvalidated reference is not evidence.
 
+### 6.1 Status and what follows
+
+This paper establishes a method and reports what it supports. It does not report a validated
+prediction for LptA, and no such claim should be read into it. The sequence is deliberate:
+establish the method against measured data from other systems, determine its limits, apply it
+to the campaign panel under those limits, register the forward test, and only then commit
+synthetic resources. Synthesis and assay of the ten registered compounds is the immediate next
+step; the manuscript will be updated with those results, which Section 5 anticipates in both
+directions.
+
 ---
@@ -705,3 +869,3 @@ All code, computed scores, per-target statistics, bootstrap intervals and figure
 at **https://github.com/Caroldorium/mmpbsa-macrocycle-benchmark** under an MIT licence.
-`scripts/size_matched_pairwise_accuracy.py` reproduces Table 3 and the primary figure.
+`scripts/size_matched_pairwise_accuracy.py` reproduces Table 5 and the primary figure.
 
```

---

## v8 — `2dd6c9a` — Manuscript v8: state the compound attrition chain, final consistency pass

`1 file changed, 9 insertions(+), 3 deletions(-)`

```diff
@@ -158,5 +158,4 @@ five-target macrocycle and cyclic-peptide series of Awoonor-Williams and colleag
 PCSK9, MCL-1, JAK2 and Cyclin A/B — comprising over 230 peptidic and non-peptidic
-analogues.[^2] We use the 241 compounds for which structures and labels could be assembled; 235
-carry an experimental binding free energy and 234 survive quality control (Section 2.4).
-Analyses requiring unambiguous labels additionally exclude censored measurements, leaving 220.
+analogues.[^2] We use the 241 compounds for which structures and labels could be assembled; the
+attrition to the analysed sets is given below.
 
@@ -167,2 +166,9 @@ difficult to obtain.
 
+**Compound attrition, stated once for reproducibility.** Of the 241 compounds for which
+structures and labels could be assembled, **235** carry an experimental binding free energy;
+**234** remain after excluding `MCL_1__cpd_55` for an unrelieved steric clash (Section 2.2); and
+**220** remain for the pairwise analyses, which additionally exclude censored measurements
+because the ordering of a censored pair is ambiguous. Correlation analyses in §3.1–3.3 use the
+234; pairwise analyses in §3.4 use the 220. Every reported *n* refers to one of these two sets.
+
 Compounds are not independent within a target: two targets span more than one receptor crystal
```

---

## v9 — `50b5c79` — Manuscript v9: repair a malformed table, verify the rendered output

`1 file changed, 12 insertions(+), 13 deletions(-)`

```diff
@@ -548,12 +548,11 @@ matching; the right admits a three-atom window and is better powered.
 
-| | | identical heavy atoms | | within 3 heavy atoms | |
-|---|---:|---:|---|---:|---|
-| target | | MM/PBSA | 95% CI | MM/PBSA | 95% CI |
-| Cyclin A/B | | 0.533 | [0.053, 1.000] | 0.674 | [0.454, 0.863] |
-| JAK2 | | 0.637 | [0.397, 0.831] | 0.598 | [0.437, 0.732] |
-| KRAS | | 0.518 | [0.289, 0.727] | 0.496 | [0.339, 0.664] |
-| MCL-1 | | 0.549 | [0.328, 0.757] | 0.594 | [0.444, 0.733] |
-| **PCSK9** | | **0.765** | **[0.613, 0.882]** | **0.757** | **[0.671, 0.830]** |
-| **macro** | | **0.601** | [0.472, 0.720] | **0.624** | **[0.555, 0.689]** |
-| FEP+ macro | | 0.841 | [0.730, 0.912] | 0.839 | [0.782, 0.886] |
+| target | MM/PBSA, identical size | 95% CI | MM/PBSA, within 3 atoms | 95% CI |
+|---|---:|---|---:|---|
+| Cyclin A/B | 0.533 | [0.053, 1.000] | 0.674 | [0.454, 0.863] |
+| JAK2 | 0.637 | [0.397, 0.831] | 0.598 | [0.437, 0.732] |
+| KRAS | 0.518 | [0.289, 0.727] | 0.496 | [0.339, 0.664] |
+| MCL-1 | 0.549 | [0.328, 0.757] | 0.594 | [0.444, 0.733] |
+| **PCSK9** | **0.765** | **[0.613, 0.882]** | **0.757** | **[0.671, 0.830]** |
+| **macro, MM/PBSA** | **0.601** | [0.472, 0.720] | **0.624** | **[0.555, 0.689]** |
+| macro, FEP+ | 0.841 | [0.730, 0.912] | 0.839 | [0.782, 0.886] |
 
@@ -580,5 +579,5 @@ most where desolvation dominates. PCSK9 does have much the most polar interface
 (median buried-area polarity 42.2%, against 30.5% for the next-highest). **But the explanation
-does not survive the second case:** JAK2 is the fourth most polar of the five (23.8%) and the
-second-best performer, while KRAS is second most polar (30.5%) and indistinguishable from
-chance. Polarity alone therefore does not account for where the method works, and with five
+does not survive the other cases:** JAK2 is the fourth most polar of the five (23.8%) yet has
+the second-highest point estimate at exact matching, while KRAS is the second most polar
+(30.5%) and sits at chance in both bands. Polarity alone therefore does not account for where the method works, and with five
 targets we cannot identify what does. We report the pattern without an explanation rather than
```
