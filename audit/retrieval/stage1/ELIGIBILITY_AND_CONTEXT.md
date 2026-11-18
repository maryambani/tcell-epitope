# Stage-1 eligibility intersections and context audit

Illustrative subsets only. **Not a V2 policy.** No conflict rule was chosen. No model was trained.

These counts are from the Stage-1 sample of **240 peptides** (48 per stratum). Equal allocation is useful for failure modes. It is not representative of the 54,448-peptide pool. Peptide-level summaries reduce domination by heavily studied peptides; stratum-specific summaries reveal heterogeneity. Neither directly estimates full-pool coverage without accounting for the sampling design.

`NLVPMVATV` is still **41.2%** of all recovered assays. Its lower share than in the 46-peptide pilot does **not** mean concentration is resolved. Tables below also show the same subsets without that peptide.

Historical V1 conflict, current all-observation conflict, and conflict **inside a candidate subset** are separate fields. A mixed peptide across all hosts and endpoints can become outcome-consistent in a subset, or have no eligible evidence left. That is not an intrinsic property of the sequence, and consistency after filtering is not a universally correct biological label.

---

## How many peptides and studies support each plausible task?

Outcome pattern is assigned only among **eligible** observations for that subset:

- `positive_only` / `negative_only` / `mixed`
- `insufficient_unknown`: no eligible assays, or none with a positive or negative outcome family

| Subset (illustrative) | Peptides with eligible assays | Pos-only | Neg-only | Mixed | Insufficient / unknown | Assays | References | Share of assays from `NLVPMVATV` |
|------------------------|------------------------------:|---------:|---------:|------:|-----------------------:|-------:|-----------:|---------------------------------:|
| All recovered observations | 240 | 67 | 77 | 96 | 0 | 2,235 | 748 | 41.2% |
| Human host | 174 | 52 | 54 | 68 | 66 | 1,927 | 650 | 45.5% |
| Human + annotated class I | 158 | 44 | 50 | 64 | 82 | 1,832 | 626 | 47.4% |
| Human + annotated class I + functional endpoint | 121 | 32 | 40 | 49 | 119 | 1,166 | 434 | 41.6% |

Without `NLVPMVATV`, the last row is **681 assays, 120 peptides, 242 references**. Concentration remains.

Among the 96 V1-conflict peptides, tightening to human + class I + functional does **not** leave them all mixed:

| Inside human + class I + functional | V1-conflict peptides |
|-------------------------------------|---------------------:|
| Still mixed | 49 |
| Positive-only | 12 |
| Negative-only | 6 |
| No eligible assays | 29 |

V1 retained peptides in this sample stay non-mixed in that subset (`n_v1_retained_now_task_mixed = 0`). That is reconciliation with historical strata, not independent biological validation.

### Same subset by sampling stratum

Human + annotated class I + functional:

| Stratum | Peptides with eligible assays | Mixed | Assays |
|---------|------------------------------:|------:|-------:|
| conflict_both_databases | 39 | 30 | 928 |
| conflict_iedb_only | 28 | 19 | 114 |
| retained_both_databases | 24 | 0 | 36 |
| retained_neg_iedb_only | 18 | 0 | 53 |
| retained_pos_iedb_only | 12 | 0 | 35 |

Most of the functional human class I assays in this draw still sit in the conflict strata.

---

## Exact functional endpoints (human + class I)

Parent subset: human host + annotated class I + functional family. Exact `assay_names` are not collapsed. Top rows by peptide count:

| Exact assay name | Peptides | Pos-only | Neg-only | Mixed | Assays | References | Share from `NLVPMVATV` |
|------------------|--------:|---------:|---------:|------:|-------:|-----------:|-----------------------:|
| IFNg release\|ELISPOT | 78 | 24 | 33 | 21 | 346 | 176 | 20.8% |
| IFNg release\|ICS | 22 | 16 | 4 | 2 | 147 | 108 | 64.0% |
| activation\|biological activity | 22 | 6 | 15 | 1 | 39 | 18 | 38.5% |
| IFNg release\|ELISA | 20 | 9 | 3 | 8 | 96 | 50 | 22.9% |
| cytotoxicity\|51 chromium | 19 | 9 | 6 | 4 | 186 | 103 | 45.7% |

IFN-γ ELISPOT is the largest exact endpoint here, and `NLVPMVATV` is less dominant in that slice than in tetramer-heavy pooled tables. That still does not choose ELISPOT as V2.

---

## Restriction-evidence category (human + class I + functional)

| Evidence category | Peptides | Pos-only | Neg-only | Mixed | Assays |
|-------------------|--------:|---------:|---------:|------:|-------:|
| prediction | 49 | 25 | 14 | 10 | 84 |
| not_determined | 35 | 13 | 19 | 3 | 101 |
| cited_reference | 31 | 6 | 9 | 16 | 721 |
| binding_assay | 29 | 12 | 5 | 12 | 107 |
| other | 18 | 8 | 8 | 2 | 52 |
| single_allele_present | 11 | 5 | 2 | 4 | 101 |

Cited-reference assays are few peptides and many rows (`NLVPMVATV` again). Prediction has more peptides and fewer rows. These are different prediction tasks if evidence quality is part of eligibility.

---

## Context-classification audit

96 peptides have both positive-family and negative observations in the **unfiltered** Stage-1 table.

| Detail | Count |
|--------|------:|
| Mixed peptides (any pos + neg observation) | 96 peptides |
| Peptides with ≥1 mixed complete recorded-context key | 14 peptides |
| Distinct mixed complete-context keys (groups) | **17 groups** |
| Across known context, both families have complete context, keys do not overlap | 78 peptides |
| Labeled across known context, but negative observations lack complete MHC class/allele | 4 peptides |
| Cannot localize (overlap only after filling missing fields) | 0 peptides |
| Raw allele **strings** differ between positive and negative groups | 48 peptides |

Do not call the 14 peptides “14 mixed groups.” Eleven of those peptides have one mixed complete-context key; three have two (14 peptides, 17 groups).

“Across known context” is **not** automatically “both families have sufficiently recorded context.” Four of the 82 unmatched peptides are incomplete on the negative side: `FGYMTFVHF`, `GQIVGGVYL`, `QLLGQTADV`, `RLVTNRAGK`. Keep them separate from the 78.

The 14 peptides / 17 mixed keys are investigation targets, not automatic exclusions. Distinct `tcell_id`s are curated records, not replicates. Exact-string allele differences (`HLA-A2` vs `HLA-A*02:01`) can block a match without proving a different experiment. Raw allele strings and `mhc_specificity` flags are stored separately. `direct_ex_vivo_bool=0` is recorded, not missing.

Files: `eligibility_intersections.csv`, `eligibility_intersection_summary.csv`, `context_audit.csv`, `eligibility_and_context.json`.
