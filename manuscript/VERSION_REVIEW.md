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

## Final state

**Internal audit, all clean:**

- Tables 1–5 defined and referenced
- Figures 1–4 declared and referenced, all four files present
- 26 sections, no dangling cross-references
- 19 references, all cited and all defined, verified against PubMed
- 6,934 words
- No trace of the retracted claim

**Formats:** Markdown source, DOCX with embedded figures, self-contained HTML.

---

## What I need from you

**1. Confirm or overrule the publication-readiness judgement.** I consider the retrospective
half ready; that is my assessment, not a fact, and you have not agreed to it.

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
