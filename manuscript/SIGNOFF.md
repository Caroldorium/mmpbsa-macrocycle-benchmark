# Sign-off checklist

Reviewing nine versions is a reading task. This turns it into decisions. Strike anything you
disagree with and I will change it; silence on an item means it stands.

---

## A. Corrections I made to my own errors — confirm each is correctly resolved

These are the ones that matter. Each names the file to check it against.

| # | what I got wrong | how it stands now | verify against |
|---|---|---|---|
| A1 | Claimed the pooled estimator rates a perturbed pose (+0.049) above the crystal pose (−0.048). **Fabricated** — the −0.048 was a clash-contaminated figure, really +0.265, and no perturbed-pose comparison exists | Removed entirely; replaced with the KRAS two-crystal-form result (Table 4) | `2026-09-24.md` lines 110–124 and 165–175 |
| A2 | "308 designed macrocycles" | **212** — I double-counted the thanatin branch | `lpta_212_ranked.csv` |
| A3 | "9–21 spurious TER records" | **14 of 50 complexes, 9–25 records** | `zha25_build_check.csv` |
| A4 | Ensemble sampling figures (−0.101 / −0.068 / −0.118) — unsupported | Replaced with replicate SD 6.81 vs between-compound SD 5.74 | `kras_all40_md_score_and_stability_metrics.csv` |
| A5 | Bootstrapped over **pairs**, which aren't independent | Resamples **compounds**; exact-match CI widened to [0.472, 0.720] and **no longer excludes chance** | `scripts/size_matched_pairwise_accuracy.py` |
| A6 | Table 5 malformed (fake two-level header) | Flat five-column table; renders verified | `MANUSCRIPT.docx` |

---

## B. Judgement calls — these are arguable, tell me if you disagree

| # | call I made | the alternative |
|---|---|---|
| B1 | **Three-atom band is the primary analysis**, exact-match reported as strictest-but-underpowered | Lead with exact matching and accept the wider interval |
| B2 | **Reported that the edge over heavy atoms is only +0.022** at the chosen margin | Omit it — it's unfavourable and the band-3 result stands without it |
| B3 | **Included §3.8 (12CA5)**, where everything failed | Cut it — it's a different system and weakens the paper |
| B4 | **Stated polarity fails to explain where the method works**, rather than offering it as a plausible account | Offer polarity as a hypothesis |
| B5 | **Framed §5 as sanity-check and anchor acquisition**, not validation, given 12% power at n=10 | Call it prospective validation |
| B6 | **Methods-forward title** | Results-forward version, retained inline in the file |

---

## C. Decisions only you can make

| # | item | what I need |
|---|---|---|
| C1 | **ACS redistribution** — the repo reproduces Awoonor-Williams' measured ΔG and FEP+ columns. The analysis is unreproducible without them | Keep as-is, or one commit to strip to join-keys plus a merge script |
| C2 | **Authorship, affiliations, funding** — currently Liu & Brown, Leslie Dan; acknowledgements name the Alliance but no grant | Confirm or correct |
| C3 | **Publication readiness** | My judgement is that it is ready. Confirm or tell me what blocks it |

---

## D. Known and deliberate gaps

- **§5 is a registered plan.** Converts to results when your ten compounds are assayed; the power analysis is already stated in both directions.
- **Only PCSK9 is individually above chance.** Stated plainly in results, discussion and limitations.
- **No measured affinities exist for any LptA design.** §4 is an application under stated limits, not a validation.

---

*Nothing further will be changed without being put in front of you first.*
