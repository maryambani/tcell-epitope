# V2 primary task specification (predeclared, before any V2 model)

Written before any V2 feature matrix, label column, split assignment, or model fit. V1 stays frozen. No labels are written by this document.

Counts are read from `audit/retrieval/full/v2_task_populations.json`, produced by `scripts/v2_task_populations.py` over `audit/retrieval/full/observations.csv` (the completed full-pool retrieval: 54,448 peptides, 157,300 observations).

V1's conflict-drop policy is **not** reused. V1 built its dataset by deleting any peptide whose historical export membership put it on both a positive and a negative list. That rule is not the V2 label rule, and nothing below inherits from it.

## Primary claim

**Given a human 8 to 11 residue peptide with annotated MHC class I restriction, predict whether it elicits a positive IFN-gamma ELISPOT response in a study not used for training and in a sequence neighborhood not used for training.** Peptides with both positive and negative eligible ELISPOT observations are left unlabeled and excluded, and their count is reported alongside every result.

Everything after this line is the audit trail for that sentence: the five decisions the reviewer demanded, in order, each with its rationale and its counts.

---

## Decision 1: Host

**Human only.** `host_is_human = 1` (`host_organism_name` is *Homo sapiens*).

HLA-transgenic mouse hosts are **not** human, as already specified in `RETRIEVAL_SPEC.md`. A transgenic mouse carrying a human allele is a mouse host and is excluded. This is the rule that produced the `VQSAKYWNL` correction recorded in `audit/retrieval/full/STATUS.md`: mouse-strain ELISPOT disagreement cannot supply a human-response label, so that peptide has no eligible evidence here rather than a label.

## Decision 2: MHC class and restriction evidence

**Annotated MHC class I required.** A peptide's observation is eligible when the record carries a class I annotation.

- Coarse annotations count. `mhc_specificity = class_only` (for example "class I" with no four-digit allele) is **eligible**. Requiring four-digit alleles would be a second, undeclared filter on curation depth rather than on biology.
- **No allele normalization is performed.** `HLA-A2` is not mapped to `HLA-A*02:01`. Serotype-level annotations (`mhc_specificity = serotype`) are **not** promoted to four-digit alleles. Raw `mhc_allele_name`, `mhc_restriction`, and `mhc_allele_evidence` strings are preserved as retrieved.
- Observations with no MHC annotation (`mhc_annotation_present = 0`) are not eligible.

Consequence to keep visible: the primary task is **not** allele-specific. It is a class I human response task. The model is not given an allele and cannot be read as predicting response in a named HLA context. Because alleles are unnormalized, any later per-allele breakdown will be a breakdown over raw strings of mixed specificity, not over a clean allele axis.

## Decision 3: Endpoint

**Primary: IFN-gamma ELISPOT** (`assay_names` = `IFNg release|ELISPOT`, within `endpoint_family = functional_response`).

**Predeclared secondary: pooled human class I functional response** (any `endpoint_family = functional_response`, which unions ELISPOT, intracellular cytokine staining, chromium release, activation, degranulation, proliferation, and other cytokine-release assays).

### Rationale

The reviewer recommended ELISPOT for **endpoint coherence, explicitly not for performance**. That reasoning is adopted as stated. ELISPOT is one measurement with one readout. The pooled functional subset is a union of assays that measure different cell behaviors at different sensitivities; a single label over that union defines the positive class by whichever assay a curator happened to run.

### What pooled functional offers, and why it is still the secondary

Recorded honestly, because it is the stronger option on every data-volume axis:

| | ELISPOT (primary) | Pooled functional (secondary) |
|---|---|---|
| Baseline labeled peptides | 14,064 | 20,634 |
| Variant B labeled peptides (primary population) | 8,544 | 15,121 |
| Variant A labeled peptides | 9,355 | 15,930 |
| Variant A positive fraction | 0.4530 | 0.4463 |
| Variant A top-study share of labeled | 8.77% (1001794) | 6.54% (1037506) |
| Variant A top 10 studies, share of labeled | 29.0% | 24.2% |

Pooled functional gives **about 70% more labeled peptides at a near-identical positive fraction and a lower top-study concentration** (70.3% more under Variant A, 77.0% more under the Variant B primary population; positive fraction differs by 0.007 under Variant A and 0.021 under Variant B).

**This was considered and rejected. Data volume does not establish construct validity.** A larger training set for a label that unions three different measurements is a larger training set for a less well defined quantity. Choosing the bigger population here would be choosing the definition likely to give the better headline score, which is the exact failure mode the reviewer has warned against repeatedly.

### One accidental consequence of the ELISPOT filter, recorded in advance

The ELISPOT filter happens to exclude reference **1037506** entirely. That study has zero rows in `per_study_human_class_I_functional_ifng_elispot.csv`; its human class I functional assays are activation, degranulation, and `IFNg release|ICS`, with no ELISPOT. In the pooled functional subset it is the single largest contributor of exclusive labeled peptides: **1,042 exclusive labeled peptides (6.54% of 15,930), 0.4% positive**.

So pooled functional retains a second near-degenerate bulk study that the primary task drops. This is a side effect of an endpoint choice made on coherence grounds, not a reason the choice was made, and it is written down here so it cannot later be presented as foresight. It does mean the pooled secondary analysis carries a confound the primary does not, and the secondary should be read with that in mind.

## Decision 4: How mixed observations become a label

**They do not. Peptides with both positive and negative eligible observations are left unlabeled and excluded from the primary analysis.**

- No any-positive-wins rule.
- No majority vote.
- No repeated-negative rule.
- No context-based conflict resolution.

The `RETRIEVAL_SPEC.md` context audit (`across_known_context` / `within_recorded_context` / `cannot_localize`) stays **descriptive**. It is not promoted into a tie-breaker.

**Mixed counts must be reported in every results table**, so that the excluded fraction is visible at every point where a metric is quoted. Under the primary population that is **1,101 mixed peptides** set aside against **8,544 labeled**, or 11.4% of eligible peptides discarded by this rule alone.

### "Not mixed" is not the same as "has a usable label"

Four states are distinct and must stay distinct in every table:

1. **positive-only**: at least one eligible positive observation, no eligible negative.
2. **negative-only**: at least one eligible negative observation, no eligible positive.
3. **no eligible evidence**: the peptide exists in the collection but has zero observations passing host, class, and endpoint filters.
4. **eligible but uninterpretable**: eligible observations present, outcome not resolvable to positive or negative.

Collapsing state 3 into state 1 or 2 is a reporting error this project has already made and corrected twice, both recorded in `audit/retrieval/full/STATUS.md`:

- **WICLLQFAY** has **no eligible evidence** under human + class I + functional. Its mixed records were tetramer and binding, filtered out by the endpoint rule. That is not a functional label.
- **VQSAKYWNL** has **no eligible evidence** under the human-host subset. Its disagreement was mouse-strain ELISPOT.

Both subsets below currently report `n_uninterpretable = 0`, so state 4 is empty as of this count. The column stays in the tables anyway; an empty category reported is auditable, an absent category is not.

## Decision 5: Confounded-study exclusion (reference 1001865)

**Reference 1001865 is excluded from the primary population.**

| Property | Value |
|---|---|
| Source | Orthopoxvirus vaccinia, Vaccinia virus NYCBH Dryvax, administration in vivo |
| PubMed identifier | none recorded |
| Eligible ELISPOT observations | 20,274, all `IFNg release\|ELISPOT` |
| Share of all eligible ELISPOT observations | 52.7% (20,274 of 38,490) |
| Peptides touched | 5,636 |
| Within-study outcome pattern | negative_only for all 5,636 (100% negative) |
| Peptides exclusive to this study | 4,804 |
| Share of baseline labeled peptides touched | 39.25% |

### Why it is excluded: entanglement, not unreliability

**The exclusion does not rest on a judgment that the study is wrong.** A large all-negative screen is a legitimate experiment and a legitimate result; most peptides in a systematic scan do not elicit a response, and reporting that is correct science. The problem is structural: with 52.7% of eligible observations and 39.25% of baseline labeled peptides coming from one experimental context, a pooled score over the baseline population is largely a measurement of how well a model reproduces that one context. That is a statement about the dataset's composition, not about the study's quality.

### Variant B is the primary population

**Variant B: drop every peptide the study touched** (all 5,636, exclusive and shared alike).

### Variant A is a sensitivity analysis only

**Variant A: drop only the 4,804 exclusive peptides, keep the 832 shared peptides with the 1001865 rows removed, and recompute labels from the remaining evidence.**

Variant A is reported as a sensitivity analysis and is **not** the primary population, because recomputing labels after removing the rows **flips 95 peptides from mixed to positive-only**. Those peptides become positive purely by deleting the evidence that called them negative. That amounts to an implicit claim that 1001865's negative calls on those peptides should be disregarded, and **this project is making no such claim about that study's reliability.** Variant B avoids the question entirely by dropping the peptides rather than re-adjudicating them.

Variant A is still worth computing and reporting: it is the more generous population (9,355 labeled against 8,544), and its study-concentration and sequence-structure diagnostics are the ones computed in the JSON. Where a diagnostic exists only for Variant A, it is labeled as such below.

---

## Populations

### Primary task: human + class I + `IFNg release|ELISPOT`

**Baseline (before any study exclusion)**

| Quantity | Value |
|---|---|
| Eligible peptides | 15,281 |
| positive-only | 4,143 |
| negative-only | 9,921 |
| mixed (left unlabeled) | 1,217 |
| eligible but uninterpretable | 0 |
| Labeled | 14,064 |
| Labeled positive fraction | 0.2946 |
| Eligible observations | 38,490 |
| Studies | 1,034 |
| Eligible rows with missing reference_id | 0 |

**Variant B, the primary population (drop every peptide 1001865 touched)**

| Quantity | Value |
|---|---|
| Labeled | 8,544 |
| Positive | 4,143 |
| Negative | 4,401 |
| Positive fraction | 0.4849 |
| mixed (left unlabeled) | 1,101 |
| Studies among labeled peptides | 760 |

The positive fraction moves from 0.2946 to 0.4849 on removing one study. That is the size of the confound, stated as a number.

**Variant A, sensitivity analysis (drop exclusive peptides only, recompute)**

| Quantity | Value |
|---|---|
| Labeled | 9,355 |
| Positive | 4,238 |
| Negative | 5,117 |
| Positive fraction | 0.4530 |
| mixed (left unlabeled) | 1,122 |
| Studies among labeled peptides | 761 |
| Peptides flipped mixed to positive-only by row removal | 95 |
| Shared peptides losing all eligible evidence | 0 |

**Variant A study concentration** (computed on the 9,355 Variant A labeled peptides)

| Quantity | Value |
|---|---|
| Largest study by exclusive labeled peptides | 1001794, 820 peptides |
| Its share of labeled population | 8.77% |
| Its positive fraction | 4.6% |
| Any single study above 10% | no |
| Top 10 studies, combined share | 29.0% |
| Studies with at least one exclusive labeled peptide | 609 |

No single remaining study dominates the way 1001865 did. Concentration is reduced, not eliminated: the top study is still 820 peptides at 4.6% positive, and 1001794 is itself a large low-positive screen.

**Variant A sequence structure** (Hamming distance <= 1 within the same length)

| Quantity | Value |
|---|---|
| Components | 8,561 |
| Singleton components | 8,173 |
| Non-singleton components | 388 |
| Largest component | 20 peptides |
| Peptides in non-singleton components | 1,182 (12.6%) |

**Variant A joint separation yield**

| Quantity | Value |
|---|---|
| Both study-exclusive and sequence-singleton | 7,646 (81.73%) |
| Positive fraction of that subset | 0.4024 |
| Study-exclusive labeled peptides (study separation alone) | 8,694 |

### Predeclared secondary: human + class I + any functional endpoint

**Baseline**

| Quantity | Value |
|---|---|
| Eligible peptides | 23,642 |
| positive-only | 7,015 |
| negative-only | 13,619 |
| mixed | 3,008 |
| Labeled | 20,634 |
| Labeled positive fraction | 0.3400 |
| Eligible observations | 71,059 |
| Studies | 2,722 |

**Variant B**

| Quantity | Value |
|---|---|
| Labeled | 15,121 |
| Positive | 7,015 |
| Negative | 8,106 |
| Positive fraction | 0.4639 |
| mixed | 2,885 |

**Variant A**

| Quantity | Value |
|---|---|
| Labeled | 15,930 |
| Positive | 7,110 |
| Negative | 8,820 |
| Positive fraction | 0.4463 |
| mixed | 2,913 |

**Variant A study concentration**

| Quantity | Value |
|---|---|
| Largest study by exclusive labeled peptides | 1037506, 1,042 peptides |
| Its share of labeled population | 6.54% |
| Its positive fraction | 0.4% |
| Any single study above 10% | no |
| Top 10 studies, combined share | 24.2% |

**Variant A sequence structure**

| Quantity | Value |
|---|---|
| Components | 14,703 |
| Singleton components | 14,025 |
| Non-singleton components | 678 |
| Largest component | 20 peptides |
| Peptides in non-singleton components | 1,905 (12.0%) |

**Variant A joint separation yield**

| Quantity | Value |
|---|---|
| Both study-exclusive and sequence-singleton | 13,151 (82.55%) |
| Positive fraction of that subset | 0.3956 |

### How settled this is

The population is defined; it is not settled in the sense of being a stable biological sample. Three things remain genuinely uncertain and are not resolved by picking a definition:

1. The positive fraction of the primary population (0.4849) is an artifact of which studies survive filtering, not a prevalence estimate. A different coherent endpoint choice moves it by 0.02 to 0.19.
2. Concentration is reduced but not low. 29.0% of Variant A labeled peptides come from 10 studies, and several of those are large low-positive screens of the same general kind as 1001865. The confound is a gradient, and 1001865 is where the gradient was steep enough to act on.
3. Why current IEDB evidence disagrees with 2025-11-03 query membership is still unresolved. These counts are current-IEDB counts and would shift on a re-pull.

---

## Evaluation design

**Generalization means both at once: a new study and a new sequence neighborhood.** The hold-out is **jointly study-separated and sequence-separated**. A peptide in the test set must satisfy both conditions.

### Study separation

Test peptides must be **exclusive to held-out studies**. A labeled peptide observed in both a held-out study and a kept study is **excluded as a leak**. It is **not** assigned to whichever study appears first in any ordering.

First-study assignment is the shortcut `RETRIEVAL_SPEC.md` already rules out, and it is the one that would quietly inflate the score: a peptide seen in training under one reference and scored in test under another is not a held-out study observation.

Under Variant A, 8,694 of 9,355 labeled peptides are study-exclusive, so study separation alone costs 661 peptides to leak exclusion.

### Sequence separation

Relatedness is **Hamming distance <= 1 within the same length**. Connected components are held out **whole**: every peptide in a component goes to the same side of the split.

Related peptides are **assigned together, not discarded.** Discarding non-singletons would throw away 1,182 peptides (12.6%) under Variant A and would bias the remaining population toward sequences with no close neighbors, which is itself a property correlated with source antigen and study design.

Note the limit of this definition: Hamming <= 1 at equal length catches point substitutions only. It does not catch length-shifted overlaps, nested peptides, or nearby windows from the same source protein. Sequence separation under this rule is a floor on relatedness control, not a guarantee of independence.

### Critical constraint: the test set cannot be drawn at random

**A balanced study-separated test set exists, but random study sampling cannot reliably reach it.**

Across **200 random study orderings** (seed 20261119) targeting a 20% test fraction on the Variant A labeled population:

| Quantity | Value |
|---|---|
| Achievable test positive fraction | 0.139 to 0.651 |
| Median test positive fraction | 0.420 |
| Studies held out, median | 136 |
| Test peptides, median | 1,897 |
| Studies available with exclusive labeled peptides | 609 |

A test positive fraction anywhere between 0.139 and 0.651 is reachable by chance alone, from the same population, at the same target size. Study positive rates are extremely heterogeneous: among the top 10 studies by exclusive labeled peptides, 1043337 is 1.8% positive and 1024519 is 100% positive. Drawing studies at random draws that heterogeneity.

Two consequences, both binding:

1. **The test set must be constructed deliberately, stratified on study positive rate.** Not sampled randomly and not reported as a random split.
2. **That construction must be fixed and written down before any model is fit.** A split chosen after seeing model scores, or reseeded until it looks reasonable, is score selection by another name, and the 0.139 to 0.651 range is exactly how much score a reseed can buy.

### Not adopted here

**This document adopts no specific fold assignment.** No study is named as held out, no component is assigned, no seed is declared final. The split rule is predeclared; the split itself is the next artifact, and it is listed below as outstanding.

---

## What this document does not claim

1. **Not a census of current IEDB.** The collection is the historical 54,448-peptide pool from `epitope_pool_8_11_preconflict.csv`, queried against current IEDB at retrieval time. Completeness is completeness for that defined collection only.
2. **No claim about peptides never in the 2025-11-03 export lists.** Nothing here measures or predicts performance on sequences outside the pool.
3. **Negative does not mean "cannot elicit a response."** A negative label means the eligible observations in this subset were negative. It is not a claim that the peptide cannot elicit a response in another host, another HLA context, another assay, or another experimental protocol.
4. **No resolution of the IEDB disagreement.** Why current IEDB evidence disagrees with 2025-11-03 query membership remains unresolved among database change, annotation change, and export difference. No choice above depends on resolving it, and no result below should be read as having resolved it.
5. **No allele-level claim.** See Decision 2: alleles are unnormalized and serotypes unmapped, so the task is class I, not per-allele.
6. **No mechanism claim.** See the anchor ablation item below.

---

## Still to be fixed before modeling

Both items are blockers in the sense that a V2 number quoted before they are done is not auditable.

1. **The deliberate stratified test-set construction.** The split rule is predeclared above; the split is not. Required: a construction procedure that is study-separated, sequence-component-separated, and stratified on study positive rate; the resulting fold assignment written to a file; and the target test positive fraction stated in advance. This must be committed before the first V2 fit. The 200-ordering sweep is the reason: without a fixed deliberate construction, the split alone moves the test positive fraction across a 0.51 range.

2. **Control masks and baselines for the notebook 06 anchor ablation.** The ablation (`notebooks/06_anchor_ablation.ipynb`) currently reports that masking P2 and PΩ costs 0.024 to 0.028 AUC on the 9-mer subset (LR 0.637 to 0.612, RF 0.673 to 0.644, CNN 0.637 to 0.613). **Still missing: equal-count non-anchor control masks, and simple baselines, to establish whether that drop exceeds what masking any two positions would cost.** Until those exist the result is information removal, not evidence of MHC anchor mechanism, and it must not be described as the latter.

Neither item is a V2 task-definition question. Both are recorded here so the predeclaration and the outstanding work sit in one place.
