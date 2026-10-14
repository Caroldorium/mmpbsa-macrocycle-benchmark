# SO3LR-SF on the macrocycle cohort

Test of the SO3LR-SF scoring function (Agha, Suárez-Dou, Kabylda, Tkatchenko & Volkamer,
ChemRxiv 10.26434/chemrxiv.15008482/v1) on the same 234/220 compound sets and under the same
two estimators as the manuscript. **No manuscript text has been changed.**

## What was run

All 241 cohort complexes were split into `receptor.pdb` + `ligand.pdb` from the **minimised**
geometry — the identical input MM/PBSA scored — with formal charges read from the prmtops, and
scored with the published CLI at its validated `--trim 8` setting. 241/241 succeeded, mean
19.6 s/complex on 24 CPU cores (JAX fell back to CPU; this matches the paper's stated CPU cost).

234 merge into the analysis set after dropping `MCL_1__cpd_55` for its unrelieved clash; 220
after additionally dropping censored measurements, exactly as in the manuscript.

## Headline

**The advantage reported on small molecules does not transfer to macrocyclic peptides.**

The preprint reports SO3LR-SF at Spearman 0.51 against MM/GBSA 0.44 and Glide 0.27 on its
8-target FEP benchmark. On our five macrocycle targets:

| method | macro Spearman | band 0 accuracy | band 3 accuracy |
|---|---:|---|---|
| FEP+ (published) | **0.831** | 0.841 [0.730, 0.910] | 0.839 [0.786, 0.885] |
| MM/PBSA (`pb_2`) | **0.545** | 0.601 [0.478, 0.722] | 0.624 [0.552, 0.689] |
| heavy-atom null | **0.487** | (undefined) | 0.602 [0.524, 0.677] |
| MM/GBSA (`gb2_1`) | **0.478** | — | — |
| **SO3LR-SF** | **0.292** | 0.597 [0.484, 0.709] | 0.576 [0.503, 0.648] |

SO3LR-SF ranks below every MM/PBSA and MM/GBSA variant tested (0.427–0.545) **and below the
heavy-atom null**. On size-matched pairs it is statistically indistinguishable from MM/PBSA at
exact matching (0.597 vs 0.601) and worse in the three-atom band (0.576 vs 0.624), where it also
falls below the trivial size descriptor (0.602).

Our MM/GBSA macro Spearman of 0.478 is close to the 0.44 the preprint reports for MM/GBSA, so
the comparator reproduces; it is SO3LR-SF that fails to.

## Where it fails, and where it wins

| target | MM/PBSA ρ | SO3LR-SF ρ |
|---|---:|---:|
| MCL-1 | 0.728 | **0.772** |
| JAK2 | 0.595 | 0.419 |
| Cyclin A/B | 0.591 | 0.409 |
| KRAS | 0.225 | **−0.106** |
| PCSK9 | 0.584 | **−0.035** |

MCL-1 is the one target where SO3LR-SF is the best non-FEP method. It inverts on KRAS and
collapses to zero on PCSK9 — the target where MM/PBSA is strongest (pairwise 0.757). The two
methods fail on different targets, which is the one genuinely encouraging result here.

SO3LR-SF is **not** less size-driven: its correlation with heavy-atom count (−0.589) is
essentially MM/PBSA's (−0.577). It does not escape the confound the manuscript is about.

## Pocket polarity — manuscript §3.4 re-tested

The preprint reports that polar SASA / total SASA of the **binding pocket** predicts per-target
performance at r = 0.81. Manuscript §3.4 tested a different quantity — polar fraction of the
**interface** that gets buried — and concluded polarity does not explain the target pattern.
The pocket descriptor was therefore recomputed from scratch (Shrake–Rupley, 241/241 complexes).

| descriptor | vs MM/PBSA ρ | vs MM/PBSA accuracy |
|---|---:|---:|
| pocket polar ratio (preprint's) | −0.483 (p = 0.41) | +0.363 (p = 0.55) |
| interface polar % (manuscript's) | −0.353 (p = 0.56) | +0.520 (p = 0.37) |

**§3.4's conclusion survives.** Neither descriptor explains where MM/PBSA works, the two sign
disagree between the two skill measures, and the descriptors themselves correlate r = +0.77 —
they order the targets similarly, so the conclusion never hinged on the choice. §3.4 can now
state it was tested both ways, which strengthens it. **No revision is required.**

One suggestive observation, reported as suggestive only: pocket polarity vs SO3LR-SF skill is
r = **−0.798**, the *opposite* sign to the preprint's +0.81 — SO3LR-SF does worse on more polar
pockets here. That is consistent with solvation/charge handling being the failure mode on
charged macrocycles, but with n = 5 targets and p = 0.11 it is a hypothesis, not a result.

## Caveat that limits all of the above

Five targets is 3 degrees of freedom. These per-target correlations cannot confirm or refute the
preprint's r = 0.81, which was measured across eight. The compound-level comparisons
(234 and 220 compounds) are adequately powered; the five-point target-level correlations are not.

## Reproduce

    python compare_so3lr.py <cohort_scores.csv> so3lr_sf_scores.csv   # pooled, diagnostic only
    python macro_pairwise.py <cohort_scores.csv> so3lr_sf_scores.csv  # the manuscript's estimator
    python pocket_polarity.py <split_dir> pocket_polar_sasa.csv
    python fig_so3lr.py

Pairs are formed **within a target** and averaged across targets. A pooled version that compares
a KRAS compound to a PCSK9 compound gives materially different numbers (0.618/0.594) and is
wrong; it is retained in `compare_so3lr.py` only as a diagnostic contrast.
