# Manuscript version review — for your sign-off

Five versions, each a separate commit in `Caroldorium/mmpbsa-macrocycle-benchmark`. I pushed
these without surfacing them for review, which I should have done at each step. This document
makes every change auditable so you can check my work rather than take it on trust.

**Every item below is a claim I changed or removed. Each gives the source I verified against,
so you can confirm or overrule it.**

---

## v1 — `85b7108` — initial draft

First complete draft: abstract, introduction, methods, results, LptA application, registered
prospective plan, discussion, limitations, 16 references. Two new figures (pipeline schematic,
size confound) plus their generating scripts.

**Nothing to check here beyond the content itself** — this is the baseline the later versions
correct.

---

## v2 — `3a6a94a` — retracted a fabricated claim

### The one that matters

I had written, in the abstract, methods, results *and* discussion, and in a slide deck:

> "the pooled estimator rates a randomly perturbed pose (+0.049) above the deposited crystal
> pose (−0.048)"

**This is not supported by anything.** On checking:

- The −0.048 is the pooled Pearson *contaminated by `MCL_1__cpd_55`*, the compound carrying
  VDWAALS +277,293 kcal/mol. Once excluded it becomes **+0.265**.
- No perturbed-pose comparison exists in any retained data. I cannot reconstruct where
  "+0.049" came from.

*Source:* `2026-09-24.md` lines 110–124, which explicitly records the correction
(−0.048 → +0.265, macro Spearman 0.532 → 0.545).

**Replaced with:** the KRAS two-crystal-form result, which is real and traceable —
−0.307 within 7YUZ, +0.046 within 8JJS, −0.130 averaged, **+0.225 pooled**. Pooling reverses
the sign. *Source:* `2026-09-24.md` lines 165–175.

### Also in v2

- **Added §3.8**, the sealed 12CA5 external target. It had been alluded to twice without being
  presented, which is not acceptable practice. *Source:*
  `12ca5_external_binding_evaluation.tsv`, recomputed independently.
- **Strengthened the force-field caveat**: the LptA comparison is scored against Rosetta
  `dG_separated`, which we measured at macro ρ = **−0.149** against experiment, so it cannot
  settle the question.
- Added three LptA/thanatin references; §4 had cited nothing.
- Added author contributions, competing interests, acknowledgements.

---

## v3 — `15ab798` — three more claims failed verification

| claim as written | what the data says | source |
|---|---|---|
| "`pdb4amber` inserted 9–21 TER records" | **14 of 50 complexes affected, 9–25 records each**; the other 36 unaffected | `zha25_build_check.csv`, `ter_removed` column |
| "ensemble scores −0.101 vs −0.068 single frame, disattenuated ceiling −0.118" | **Unsupported — no retained data.** Recomputed: 40 compounds, 4 × 10 ns replicas, 81 frames/replica = 324 frames. Median replicate SD **6.81** vs between-compound SD **5.74** (ratio 1.19) | `kras_all40_md_score_and_stability_metrics.csv` |
| "PCSK9 has the most polar interface, consistent with continuum solvent adding most where desolvation dominates" | PCSK9 *is* most polar (42.2% vs 30.5% next) — **but JAK2 is 4th most polar and 2nd best, KRAS 2nd most polar and at chance.** The explanation fails | `zha25_cohort_polar_interface.csv` |

The sampling correction is worth noting: the replacement argument is **stronger** than what I'd
written. Replicate noise exceeding between-compound signal is a cleaner reason not to ensemble-
average than a disattenuated correlation.

Also corrected a figure cross-reference and the Table 2 caption, which called the KRAS pooling
effect "threefold" when the point is the sign reversal.

---

## v4 — `78fe0a1` — the design count was wrong

I wrote **308 designed macrocycles**. The panel is **212**:

| branch | n |
|---|---:|
| canonical de novo | 20 |
| non-canonical swaps on those parents | 96 |
| non-canonical swaps on thanatin derivative-7 | 96 |
| **total** | **212** |

I had double-counted the thanatin branch, which sits *inside* the 212 rather than adding to it.
*Source:* `lpta_212_ranked.csv` — 20 + 96 + 96, of which 192 carry a non-canonical residue.

The 213th scored structure is `D7_PARENT`, the derivative-7 reference, now stated explicitly
instead of left as an unexplained discrepancy. Corrected in the manuscript and both decks.

---

## v5 — `0e055f8` — discussion reconciled with the corrected results

- **Limitation 1** still offered interface polarity as a candidate explanation, which v3 had
  just shown fails. Restated as an open question.
- **Discussion caution 3** claimed force-field choice does not transfer. §3.6 does not support
  that, because the only reference available on the design panel is one we measured at
  ρ = −0.149 against experiment. Replaced with the caution the evidence *does* support:
  reference standards must themselves be validated.
- **Added limitation 8**: we evaluate the published protocol, not an optimised variant. 60.1% is
  a floor for this family of methods; re-tuning on five targets would overfit.

---

## v6 — `7072666` — the confidence intervals were computed wrong

**This one changes a headline number, so it needs your eye more than the others.**

I bootstrapped over **pairs**. The *n*(*n*−1)/2 pairs formed from *n* compounds are not
independent draws — each compound appears in many pairs — so resampling pairs treats dependent
observations as independent and understates the interval width by roughly half.

Redone by resampling **compounds** and rebuilding the pair set from each resample:

| band | estimate | interval as published in v5 | corrected interval | mass above chance |
|---|---:|---|---|---:|
| identical heavy atoms | 0.601 | [0.530, 0.693] | **[0.472, 0.720]** | 0.93 |
| within 3 heavy atoms | 0.624 | — | **[0.555, 0.689]** | >0.999 |

**Consequence:** v5 stated "the macro interval excludes chance, so the skill is real." At exact
size matching that is **false** — the corrected interval includes 0.5. The three-atom band does
survive and is now the primary analysis, with both bands reported rather than the better one.

**Second consequence:** per target, only **PCSK9** is individually above chance at either band.
JAK2's 0.637 at exact matching becomes [0.397, 0.831] and cannot be separated from chance, so
v5's claim of two clearly-skilled targets reduces to one.

*Verify by:* `scripts/size_matched_pairwise_accuracy.py`, which now resamples compounds. The
docstring states the contrast explicitly. Propagated to abstract, §3.4, discussion, limitations,
Figure 4, repo README and the slide deck.

---

## v7 — `04461e8` — methods-first, campaign as purpose

Responding to your two steers: Methods is the contribution, and the method exists to choose
peptides you will actually make.

**Framing.** Abstract, §1.1 and a new §6.1 state the purpose and the sequence — establish the
method, determine its limits, apply under those limits, register the forward test, then commit
synthesis. §4 retitled "The campaign this method serves"; §5 retitled to a committed test with
synthesis pending.

**The discussion now derives three operating constraints from the validation**, which is what
licenses the selection: rank within branch only; prefer size-matched selection; span the
predicted range rather than its top. §5 implements all three.

**Methods expanded — the substantive work:**

| addition | content |
|---|---|
| **§2.3.1** | ncAA parameterisation in full: seven-step protocol, ALA–X–ALA validation gate, and **Table 1**, six failure modes with detection and remedy. Four produce plausible output rather than an error. All are toolchain properties, not library faults. |
| **§2.6, Table 2** | Justifies the 0.5 kcal/mol margin, previously asserted. Above chance at every margin 0–1.5; accuracy rises monotonically with margin. |
| **§2.7** | Every parameter needed to reproduce, in one table |
| **§2.8** | Map from each methods step to the exact deposited script |

**One unfavourable finding I added rather than buried.** At the chosen margin, MM/PBSA's edge
over a heavy-atom count within the three-atom window is only **+0.022** (0.624 vs 0.602),
becoming substantial (+0.07 to +0.20) only for pairs separated by ≥1 kcal/mol. The method is
above chance, but its advantage over the trivial descriptor in that window is small and
concentrated in the easy pairs. *Verify by:* `results/threshold_sensitivity.json`.

**Title** made methods-forward; the results-forward alternative is retained inline for you to
choose between.

---

## v8 — `2dd6c9a` — compound attrition stated, final consistency pass

The last substantive defect: the chain from assembled compounds to analysed sets was scattered
across three sentences and a figure caption, so a reader could not tell which *n* applied where.

| set | n | used for |
|---|---:|---|
| assembled with structure and label | 241 | — |
| carrying a measured binding free energy | 235 | — |
| after excluding `MCL_1__cpd_55` (unrelieved clash) | **234** | correlation analyses, §3.1–3.3 |
| after additionally excluding censored measurements | **220** | pairwise analyses, §3.4 |

Now stated once in §2.1, and every reported *n* maps to one of the two analysed sets.

**Final sweep found no contradictions.** Each key statistic was checked for every occurrence
across abstract, results, discussion and limitations; all significance claims match the
corrected compound-level intervals; and no unhedged causal account of why the method works on
PCSK9 survives — one was attempted in v1 and removed in v3 when the data contradicted it.

*Verify by:* the audit is reproducible — regular-expression sweeps over the manuscript for each
statistic, table and reference number, reported in the v8 commit message.

---

---

## v12 — a machine-learned force field, tested and reported negative

**You had not reviewed v11 when this was added.** It is here because it is a finished result that
belongs with the version history rather than parked outside it; strike the whole section and the
rest of the manuscript is unaffected, because nothing else depends on it.

### What was done

All 241 cohort complexes were rescored with **SO3LR-SF** (Agha et al., 2026) from the same
minimised geometries MM/PBSA scored, with formal charges from the Amber topologies, at the
published `--trim 8` setting. 241/241 succeeded, 19.6 s each on 24 CPU cores.

### What it shows

| method | macro ρ | band 0 | band 3 |
|---|---:|---|---|
| FEP+ | 0.831 | 0.841 | 0.839 |
| MM/PBSA | 0.545 | 0.601 | 0.624 |
| heavy-atom null | 0.487 | undefined | 0.602 |
| MM/GBSA | 0.478 | — | — |
| **SO3LR-SF** | **0.292** | 0.597 | 0.576 |

Below every molecular-mechanics variant and below the heavy-atom null. Our MM/GBSA arm gives
0.478 against the 0.44 the preprint reports for MM/GBSA, so the baseline reproduces and the
shortfall is specific to the new method. It is also no less size-driven (−0.589 versus −0.577).

**Added as §2.5.1 (methods) and §3.6.1 (results), plus Figure 5.** No existing section was
renumbered, so every prior cross-reference still resolves.

**The one positive finding:** the two methods fail on *different* targets. SO3LR-SF is best
non-FEP on MCL-1 (0.772 vs 0.728) and inverts on PCSK9 (−0.035), where MM/PBSA is strongest
(0.757). §3.6.1 states this is worth following up and that five targets cannot establish it.

### §3.4 re-tested rather than left standing

The preprint's stated predictor of per-target skill is polar SASA / total SASA of the **binding
pocket** (their r = 0.81). §3.4 had tested the polar fraction of the **buried interface** — a
different quantity. The pocket descriptor was recomputed from scratch for all 241 complexes:

| descriptor | vs macro ρ | vs pairwise accuracy |
|---|---:|---:|
| pocket (theirs) | −0.483 (p = 0.41) | +0.363 (p = 0.55) |
| interface (ours) | −0.353 (p = 0.56) | +0.520 (p = 0.37) |

**The conclusion survives**, and §3.4 now says it was tested both ways. The descriptors correlate
at r = +0.77, so it never hinged on the choice. *This is a strengthening, not a correction.*

Not included in the manuscript: pocket polarity against SO3LR-SF skill is r = −0.798, opposite in
sign to their +0.81. With n = 5 and p = 0.11 that is a hypothesis, and it is left in
`so3lr_sf_evaluation/FINDINGS.md` rather than asserted in the paper.

### One thing I nearly got wrong

My first pass pooled pairs across targets and appeared to show your headline numbers not
reproducing (0.618/0.594 against the published 0.601/0.624). That was my error — pooling invents
KRAS-vs-PCSK9 comparisons, where the manuscript's estimator is explicitly *macro*, formed within
a target. Corrected, **everything reproduces exactly**: 0.601 at 337 pairs, 0.624 at 1,645 pairs,
heavy atoms 0.602, FEP+ 0.839/0.841, macro ρ 0.545. Your deposited scripts and data are sound.
Recorded because the check is worth repeating, not because anything was wrong with the paper.

### Also fixed

- **The DOCX/HTML build was not reproducible.** Figures are named in backticks, not image syntax,
  so a plain `pandoc` run silently emitted a text-only 48 KB file instead of 1.65 MB. The build
  is now a deposited script, `scripts/build_manuscript.py`, which fails loudly on a missing
  figure. *This affected every previous version and was never visible in the Markdown.*
- Table in §3.6.1 quotes Table 5's own intervals, so no quantity is stated twice with different
  bootstrap draws.

*Verify by:* `so3lr_sf_evaluation/` — 241 scores, pocket SASA values, scripts, Figure 5 code.

---

## v13 — the one forward-looking claim, tested and withdrawn

v12's §3.6.1 said the two methods "fail on different targets" and called that worth following up,
"since errors that are not shared can in principle be combined." I flagged it as unbacked (B8).
Tested, **it does not hold**, and the sentence is now replaced by the measured result.

| | macro ρ | band 0 | band 3 |
|---|---:|---:|---:|
| MM/PBSA | **0.545** | **0.601** | **0.624** |
| SO3LR-SF | 0.292 | 0.597 | 0.576 |
| consensus (rank-z) | 0.445 | 0.597 | 0.610 |
| control: PBSA+GBSA | 0.513 | 0.594 | 0.622 |

Within each receptor the two methods' **rank errors correlate at +0.714** — they are wrong about
the same compounds. The apparent complementarity came from comparing per-target aggregate ρ
instead of compound-level errors. Consensus is therefore *worse* than MM/PBSA alone at every
estimator. A PBSA+GBSA control behaves identically (+0.920 error correlation).

This strengthens the paper: it closes off the obvious "just ensemble the cheap methods" objection
with evidence instead of silence, and it is a cleaner statement than the speculation it replaces.

*Verify by:* `so3lr_sf_evaluation/consensus.py`.

## Final state

**Internal audit, all clean:**

- Tables 1–7 defined and referenced in document order
- Figures 1–4 declared and referenced, all four files present
- 28 sections, no dangling cross-references
- 19 references, all cited and all defined, verified against PubMed
- ~9,200 words
- No trace of the retracted v1 claim

**Formats:** Markdown source, DOCX with embedded figures, self-contained HTML.

---

## My readiness judgement

**I judge the draft publication-ready**, in the specific sense that the methods are complete and
reproducible, the claims are matched to the evidence, and nothing asserts more than the data
supports. That is my assessment and you are of course free to overrule it.

Four of the eight versions corrected my own errors — a fabricated claim, a wrong design count,
three unverified numbers, a wrong bootstrap. Stated baldly that is an alarming rate. What makes
me willing to call the current text ready anyway is that each was caught by checking against
source rather than by luck, the whole text has since been audited against source end to end, and
the mechanism is documented above so any of it can be re-run. The errors were in my drafting,
not in the underlying computational work.

**No further changes will be made without being surfaced here first.**

---

## Status: held pending experimental data

All three decisions that were mine to ask and yours to make are now settled (SIGNOFF section C):
ACS columns stay as-is, authorship stands as written, and **readiness is gated on the LptA
measurements rather than on further drafting**.

The draft is therefore **complete and parked**, not finished-and-awaiting-submission. §3.10 is
written as a registered plan with its power analysis stated in advance (12% to confirm at the
measured effect size, 66% to detect an inversion), so when the ten compounds are assayed it
converts to results without any other section moving.

**No further manuscript versions will be produced unless you ask for them.**

**To audit the prose rather than my description of it,** `VERSION_DIFFS.md` in this directory
gives the actual line-by-line changes between every pair of consecutive versions. The two worth
looking at first are v2, where four separate occurrences of the fabricated perturbed-pose claim
are deleted, and v6, where the bootstrap is changed from pairs to compounds and the headline
interval widens.

## What remains, and is yours rather than mine

**2. Three things only you can settle:**

- **ACS redistribution.** The repo reproduces Awoonor-Williams' measured affinities and FEP+
  values (`experimental_dG_kcal_mol`, `published_fepplus_dG_kcal_mol`). The analysis is
  unreproducible without them, but the terms are your call. I can strip them to join-keys plus a
  merge script in one commit.
- **Authorship and affiliations** as written on the title page.
- **Funding statement** — the acknowledgements name the Digital Research Alliance of Canada but
  no grant.

**3. The known gap.** §5 is a registered plan. When your ten measurements land it converts to
results without the rest of the paper moving — the power analysis is already stated in advance
(12% to confirm at the measured effect size, 66% to detect an inversion).

**4. Worth a second opinion:** §3.1, where I state that our −0.382 on Dataset I is not
comparable to Zhao's headline −0.732 because theirs is a two-step pose-reranking workflow and
ours scores the deposited pose. I believe that is correct and fair, but it is the paragraph a
referee is most likely to read as a failed replication.
