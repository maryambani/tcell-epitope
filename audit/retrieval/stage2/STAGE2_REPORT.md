# Stage-2 stratified retrieval (no V2 labels)

Drawn after the pagination test passed. **48 peptides × 5 strata**, seed 43, excluding the 240 Stage-1 peptides. Retrieval time: **2026-11-07T18:11:32Z**. Current IEDB, not the 2025-11-03 export.

V1 is unchanged. No V2 labels, conflict rule, or model.

## Completeness

| Check | Result |
|-------|-------:|
| Queries | 240/240 `ok` |
| Unresolved (`fail` / `incomplete`) | 0 |
| Reported total vs rows | 0 mismatches |
| Duplicate `tcell_id` | 0 |
| Sequence mismatches | 0 |
| Queries needing a second page (page size 1,000) | 0 (max 40 assays) |
| Wall-clock | 133.6 s (0.557 s/id) |

Failed/incomplete queries would have been written to `unresolved_queries.csv` and **not** stored as “no assays.” This batch had none.

This batch still does not exercise production pagination. That path was tested separately on `NLVPMVATV` with page size 50.

## Coverage is still not a full-pool estimate

Equal n per stratum is for failure modes. Peptide-level summaries reduce domination by heavily studied peptides; stratum-specific summaries reveal heterogeneity. Neither directly estimates full-pool coverage without accounting for the sampling design.

`NLVPMVATV` is absent here (excluded with Stage 1). The top peptide is `ALDPHSGHFV` at **40 / 766 assays (5.2%)**. That does not mean concentration is resolved in the pool. Stage 1 remains 41% one peptide.

| Field | All assays (n=766) | Without top peptide (n=726) |
|-------|-------------------:|----------------------------:|
| MHC class present | 734 (95.8%) | n/a |
| Human host | 530 (69.2%) | n/a |

Assay families: 604 functional, 132 recognition/binding, 30 unresolved. Class I 694, class II 38, missing 32, non-classical 2.

## Eligibility intersections (illustrative, not V2)

| Subset | Peptides with eligible assays | Pos-only | Neg-only | Mixed | Insufficient | Assays | References |
|--------|------------------------------:|---------:|---------:|------:|-------------:|-------:|-----------:|
| All observations | 240 | 59 | 83 | 98 | 0 | 766 | 266 |
| Human host | 164 | 44 | 63 | 57 | 76 | 530 | 176 |
| Human + annotated class I | 150 | 36 | 59 | 55 | 90 | 503 | 162 |
| Human + class I + functional | 111 | 25 | 45 | 41 | 129 | 359 | 132 |

Human + class I + functional IFN-γ ELISPOT: **66 peptides** (17 pos-only, 32 neg-only, 17 mixed), 165 assays, 62 references.

## Historical conflict is not task conflict

All 96 V1-conflict peptides in this draw are mixed in the unfiltered current table. **Two V1-retained peptides are also mixed** in all observations:

- `WICLLQFAY` (retained pos): mixed in all observations and in human + class I (tetramer/binding). Under human + class I + functional it has **no eligible evidence** (0 functional assays). Losing the binding disagreement is not a usable functional label.
- `VQSAKYWNL` (retained pos): IFN-γ ELISPOT positive in C57BL/6 and negative in BALB/c. Under a human-host subset it has **no eligible evidence**. A mouse-strain disagreement cannot supply a human-response label.

Stage 1’s “every retained peptide stayed outcome-family consistent” did not hold in this draw. That is why V1 conflict, all-observation conflict, and subset conflict are stored as separate fields. Consistency after a filter is not a universally correct label.

## Context audit (unfiltered mixed peptides)

98 mixed **peptides**. Of those, 14 peptides have at least one complete recorded-context key that is itself mixed; those 14 peptides contain **21** distinct mixed complete-context keys (not 14 groups). 84 mixed peptides have no mixed complete-context key: 78 have complete context on both outcome sides that do not overlap, and 6 have incomplete context on one or both sides (3 pos-side incomplete, 2 neg-side, 1 both). Keep those 6 separate from the 78. Raw allele strings differ for 34 mixed peptides. The 14 peptides / 21 mixed keys are investigation targets, not exclusions.

Files: `stage2_sample.csv`, `completeness_log.csv`, `observations.csv`, `eligibility_intersections.csv`, `context_audit.csv`. Combined Stage-1 + Stage-2 is 480 peptides of 54,448. It is still not a candidate benchmark.
