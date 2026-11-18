# Staged metadata retrieval spec (no V2 labels)

V1 stays frozen. This retrieval produces **observations and completeness logs**, not modeling labels.

Records are **current IEDB** at request time, not the 2025-11-03 export. Keep original query membership on the peptide pool table.

## Target

All **54,448** valid 8–11 peptides in `epitope_pool_8_11_preconflict.csv`, including conflicts.

Query **every** `epitope_ids` value (pipe-separated), not only `primary_epitope_id`.

Do **not** extrapolate runtime from the 46-peptide pilot (`NLVPMVATV` was 82.9% of those assays). Measure throughput on a broader staged batch first.

Stage 1 (2026-11-07): 240 peptides, 5 strata × 48, seed 42. Results in `audit/retrieval/stage1/STAGE1_REPORT.md`.

Stage 2 (2026-11-07): another 240 peptides, seed 43. Full-pool retrieval writes to `audit/retrieval/full/` and skips IDs already `ok`/`empty`. Completing all 54,448 historical sequences is complete only for that defined collection. It is not a census of current IEDB and does not measure performance on untested peptides.

Do not treat measured seconds/id as a 54,448-peptide wall-clock.

Pagination acceptance test (2026-11-07, `NLVPMVATV`, page size 50 vs 10,000): **PASS**. See `audit/retrieval/pagination_test/PAGINATION_TEST.md`. Do not increase scale until this test has passed after client changes.

Failed and incomplete queries stay in `unresolved_queries.csv`. They are not recorded as “no assays.” Empty means HTTP 200 with zero rows.

## Completeness log (required per query)

| Field | Meaning |
|-------|---------|
| query_key | epitope_id(s) requested |
| retrieval_time_utc | |
| http_status | |
| status | `ok` / `empty` / `fail` / `incomplete` |
| content_range | API `Content-Range` if present |
| reported_total | from Content-Range `*/N` |
| n_rows_retrieved | |
| n_unique_tcell_id | |
| page_offsets | e.g. `0-999,1000-1999` |
| error | |
| raw_path | file of the JSON payload |

`ok` requires: HTTP 200, and if `reported_total` is known then `n_rows_retrieved == reported_total`. Otherwise mark `incomplete`. Matching a local observations CSV to a link-status CSV does **not** prove the API returned every assay.

## Observation identity

- `observation_id` = `IEDB_ASSAY:{tcell_id}` (stable if tcell_id is stable)
- Distinct `tcell_id` = distinct curated records, **not** independent biological replicates
- Validate `linear_sequence` vs export peptide on every row
- Keep `e_modification` / structure fields when the API returns them

## MHC annotation flags (do not normalize alleles)

Preserve raw `mhc_class`, `mhc_allele_name`, `mhc_restriction`, `mhc_allele_evidence`.

Add:

| Flag | Values |
|------|--------|
| mhc_annotation_present | 0/1 |
| mhc_specificity | `missing` / `class_only` / `serotype` / `allele` / `other` |
| mhc_modified_experiment | 1 if allele/restriction text contains mutant/substitution language |
| mhc_evidence_category | as recorded, plus a coarse bucket (`cited_reference`, `binding_assay`, `prediction`, `single_allele_present`, `not_determined`, `other`, `missing`) |

Do **not** map `HLA-A2` → `HLA-A*02:01`. Do **not** strip mutant suffixes. HLA-transgenic mouse host remains a **mouse** host (`host_organism_name` as recorded; flag `host_is_human` only for Homo sapiens).

## Endpoint dictionary (visible, not used to pick a winner)

Map `assay_names` to `endpoint_family` **without** collapsing the raw name:

| Family | Examples (substrings) |
|--------|------------------------|
| `functional_response` | ELISPOT, ICS, ELISA, cytotoxicity, chromium, proliferation, cytokine/IFNg/TNFa/IL- release, degranulation, granzyme, activation |
| `recognition_binding` | multimer, tetramer, qualitative binding, SPR, dissociation constant, KD |
| `structural` | x-ray, crystallography, 3D structure |
| `other_or_unresolved` | unmatched |

Do not choose the primary endpoint by model score.

## Conflict / disagreement (descriptive only)

For peptides with both positive-family and Negative assays, classify using recorded context  
`(host, mhc_class, mhc_allele_name, assay_names, reference_id, direct_ex_vivo_bool, immunization_description)`:

- `across_known_context`: no overlapping full context
- `within_recorded_context`: overlap on complete context keys
- `cannot_localize`: overlap only after treating missing fields as equal, or keys incomplete

Also report, separately, whether **both** outcome families have at least one complete context. Do not equate “no matching complete key” with “both sides were fully recorded.” Keep raw allele strings and `mhc_specificity` flags separate: exact-string mismatches can block a match without a meaningful experimental difference.

“Same recorded context” ≠ same experiment. Distinct `tcell_id` values are distinct curated records, not independent biological replicates. Treat `direct_ex_vivo_bool=0` as a recorded value, not missing. Do **not** apply any-positive-wins, majority vote, or repeated-negative rules.

Conflict fields are not an intrinsic property of a sequence:

| Field | Meaning |
|-------|---------|
| `v1_conflict` | Historical filename membership (both pos and neg lists) |
| `current_all_observation_conflict` | Mixed positive/negative in all recovered observations |
| `task_observation_conflict` | Mixed **inside a named eligibility subset** |

A peptide can change pattern, or lose all eligible evidence, when host/class/endpoint filters are applied. Consistency after filtering is not a universally correct biological label.

## Next summaries (before any V2 decision)

After each staged batch, and again after a complete (or explicitly incomplete) full-pool pull:

- assay-level, peptide-level, and stratum-level coverage; with and without dominant peptides
- **eligibility intersections** (illustrative, not adopted): human; human+class I; human+class I+functional; then by exact assay name and restriction-evidence category
- for each intersection: distinct peptides, references, and outcome patterns (`positive_only` / `negative_only` / `mixed` / `insufficient_unknown`)
- historical vs all-observation vs task-specific conflict
- context-audit completeness (complete vs partial keys; distinguishing fields; groups still mixed after a match)
- peptides per reference; multi-reference peptides; connectivity (no first-study assignment)
- similarity grouping remains a **separate** workstream; do not mix unverified V1 similarity files into this retrieval

Equal n-per-stratum samples do not estimate full-pool coverage without a sampling-design correction.
