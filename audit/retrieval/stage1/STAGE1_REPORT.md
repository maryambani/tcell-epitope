# Stage-1 IEDB retrieval (no V2 labels)

V1 stays frozen. This batch recovers **observations and completeness logs**, not modeling labels.

Retrieval time: **2026-11-07T10:27:24Z**. Records are **current IEDB**, not a reconstruction of the 2025-11-03 export. Original query membership stays on `audit/epitope_pool_8_11_preconflict.csv`.

No observations were removed or relabeled. No V2 filters, conflict rule, or split was selected. No model was trained.

---

## What was retrieved

Stratified sample: **48 peptides × 5 strata = 240 peptides**, seed 42.

| Stratum | Peptides |
|---------|---------:|
| retained_pos_iedb_only | 48 |
| retained_neg_iedb_only | 48 |
| retained_both_databases | 48 |
| conflict_iedb_only | 48 |
| conflict_both_databases | 48 |

Each sampled peptide had one epitope ID in this draw (240 IDs / 240 peptides). CEDAR-only 8–11 peptides remain n=0, so that stratum was not sampled.

| Check | Result |
|-------|-------:|
| HTTP 200 and status `ok` | 240/240 |
| Reported `Content-Range` total vs rows retrieved | 0 mismatches |
| Duplicate `tcell_id` | 0 |
| Assay sequence vs export peptide | 0 mismatches |
| Queried epitope ID vs returned `structure_id` | 0 mismatches |
| Queries that needed a second page | 0 (max 920 assays; page size 1,000) |

Files: `completeness_log.csv`, `observations.csv`, `raw/*.json` (sha256 per page), `stage1_sample.csv`.

Matching the observations CSV to the completeness log still cannot prove that a later full-pool pull is exhaustive. It does show that **this** batch’s reported totals matched the payloads that were stored.

---

## Throughput (measured, not extrapolated)

| Quantity | Value |
|----------|------:|
| Wall-clock | 199.43 s |
| Seconds per epitope ID | 0.831 |
| Assays per peptide | min 1, median 2, mean 9.31, max 920 |

**Do not** multiply 0.831 s by 54,448. This batch still contains `NLVPMVATV`, the median peptide has two assays, and the remaining pool still contains other high-volume peptides that will paginate. Use this number as a measured rate for typical one-page queries, not a full-pool runtime.

---

## One peptide still dominates assay-level rates

`NLVPMVATV` (stratum `conflict_both_databases`) contributes **920 / 2,235 assays (41.2%)**.

Setting it aside leaves **1,315 assays, 239 peptides, and 433 references**. That is a sensitivity check, not a recommendation to drop the peptide.

Assay-level coverage is therefore still pulled toward a heavily studied CMV peptide. Peptide-level summaries reduce domination by heavily studied peptides; stratum-specific summaries reveal heterogeneity. Neither directly estimates full-pool coverage without accounting for the sampling design. Equal allocation of 48 peptides per stratum is useful for discovering failure modes. It is not representative of the full pool’s composition. The lower share of assays from `NLVPMVATV` relative to the 46-peptide pilot does not mean concentration is resolved.

| Field | All assays (n=2,235) | Without `NLVPMVATV` (n=1,315) | ≥1 observation per peptide (n=240) |
|-------|---------------------:|------------------------------:|-----------------------------------:|
| MHC class | 98.6% | 97.6% | 230/240 (95.8%) |
| PubMed ID | 96.1% | 93.5% | 222/240 (92.5%) |
| Source-antigen name | 92.9% | 87.9% | 191/240 (79.6%) |
| MHC annotation present (class or allele text) | 98.6% | n/a | n/a |

A PubMed ID on 96.1% of assays is not “96% of peptides have a PubMed ID.”

---

## MHC annotation present is not a usable allele

Raw `mhc_allele_name` / `mhc_class` are stored unchanged. Flags only.

| `mhc_specificity` | Assays | What that flag means here |
|-------------------|-------:|---------------------------|
| allele | 1,759 | Star alleles or mouse gene-like names (`H2-Kb`). Includes **1,235 `HLA-A*02:01`** |
| serotype | 342 | e.g. `HLA-A2` (243), `HLA-A24`, `HLA-B35`. **Not** converted to `HLA-A*02:01` |
| class_only | 44 | Allele text is generic `human` |
| other | 58 | Includes `HLA-E` without a star, class II serotypes (`HLA-DR13`), species words (`rabbit`, `pig`), and haplotype-like `H2-b class I` |
| missing | 32 | No class and no allele text |

Also:

- **13** observations name an explicitly modified MHC (e.g. `HLA-A*02:01 K66A mutant`). The mutant suffix was not stripped.
- Restriction evidence is mixed: cited reference 1,296; single-allele present 259; binding assay 252; prediction 159; not determined 135; other 104; missing 30.

Class annotation in this sample: **2,087 class I, 72 class II (21 peptides), 44 non-classical, 32 missing**. Length 8–11 is still not a class I filter. Forty of the class II assays sit in `retained_pos_iedb_only`.

---

## The endpoint is heterogeneous

`T Cell Assay` in the original query is not one measurement.

| `endpoint_family` | Assays |
|-------------------|-------:|
| functional_response | 1,494 |
| recognition_binding | 644 |
| other_or_unresolved | 88 |
| structural | 9 |

Most common exact `assay_names` values:

| Assay name | Assays |
|------------|-------:|
| qualitative binding\|multimer/tetramer | 614 |
| IFNg release\|ELISPOT | 466 |
| cytotoxicity\|51 chromium | 242 |
| IFNg release\|ICS | 188 |
| IFNg release\|ELISA | 105 |

A positive tetramer, SPR, or x-ray record is not the same endpoint as cytokine release or cytotoxicity. Exact assay names are kept. The family flag is visible only. It was not used to pick a primary endpoint or a label.

---

## Hosts

| Host bucket | Assays |
|-------------|-------:|
| human (`Homo sapiens`) | 1,927 |
| transgenic mouse | 150 |
| mouse (not flagged transgenic) | 138 |
| other species | 20 |

HLA-transgenic mice are counted as **mouse hosts** (`Mus musculus HLA-A*02:01 Tg`, `Mus musculus HLA-A2 Tg`, and similar). They are not recoded as human. Other hosts in this sample include rabbit, pig, cow, rat, rhesus macaque, and sheep.

---

## By sampling stratum

| Stratum | Peptides | Assays | Median assays | Mixed-outcome peptides | Class I / II | Functional / binding / structural |
|---------|---------:|-------:|--------------:|-----------------------:|-------------:|----------------------------------:|
| conflict_both_databases | 48 | 1,633 | 4 | 48 | 1,611 / 13 | 1,022 / 542 / 5 |
| conflict_iedb_only | 48 | 326 | 4 | 48 | 256 / 18 | 264 / 52 / 4 |
| retained_pos_iedb_only | 48 | 114 | 1 | 0 | 69 / 40 | 80 / 17 / 0 |
| retained_neg_iedb_only | 48 | 92 | 1 | 0 | 85 / 1 | 84 / 7 / 0 |
| retained_both_databases | 48 | 70 | 1 | 0 | 66 / 0 | 44 / 26 / 0 |

Retained peptides in this draw usually have **one** current assay. Conflict peptides have more, and `conflict_both_databases` is still dominated by `NLVPMVATV`.

---

## Conflicts: context, not a new voting rule

All **96 V1-conflict peptides** in this sample have both positive-family and negative observations in current IEDB. **None of the 144 retained peptides** do.

That supports keeping conflicts for investigation. It does not choose any-positive-wins, majority vote, or a repeated-negative rule.

Recorded-context classification `(host, MHC class, allele text, assay name, reference_id, ex-vivo flag, immunization text)`:

| Class | Peptides | Meaning |
|-------|---------:|---------|
| across_known_context | 82 | Positive and negative observations do not share a complete context key |
| within_recorded_context | 14 | At least one complete context key is shared by a positive and a negative record |
| cannot_localize | 0 | After treating `0` / `False` as recorded, no mixed peptide was mixed only because of missing fields |

`direct_ex_vivo_bool=0` is a recorded value, not missing.

**Same recorded context is not the same experiment.** Donor, TCR/clone, dose, and timepoint are still absent. Distinct `tcell_id` values are distinct curated records, not independent biological replicates.

Example of within-recorded-context disagreement: `HYINMSLPI` has one `Positive-Low` and one `Negative` IFN-γ ELISPOT in `Mus musculus BALB/c`, `H2-Kd`, same reference, same immunization string. Those are two IEDB assay IDs. They are not enough to call an annotation error.

---

## References do not establish study-holdout

| Quantity | This sample only |
|----------|-----------------:|
| Distinct `reference_id` | 748 |
| References that contain only one of the 240 sampled peptides | 650 |
| Sampled peptides that appear in more than one reference | 71 |
| Max sampled peptides in one reference | 18 |

This is recoverable provenance among 240 peptides. It is not 748 independent evaluation cohorts, and it is not a full-pool connectivity matrix.

Held-out study outcomes must not leak into training labels through peptide-level aggregation. A peptide in several studies cannot be assigned to whichever study appears first. Those rules are **not** chosen here.

---

## What this batch does not decide

| Decision | Evidence still needed |
|----------|------------------------|
| Human-only benchmark? | Full-pool peptide counts and outcomes by host |
| Class-I benchmark? | Full-pool class annotation and evidence, including 8–11 class II |
| Functional-response benchmark? | Endpoint-specific coverage and class balance on the full pool |
| Conflict policy? | Context-resolved disagreement on the full pool, plus fields this table still lacks |
| Study-separated evaluation? | Peptide/reference overlap and usable cohort sizes on the full 54,448 |
| Similarity-separated evaluation? | Separate workstream; not mixed into this retrieval |

---

Follow-up tables for this sample: `ELIGIBILITY_AND_CONTEXT.md` (intersections, task-specific conflict, context completeness). Pagination test: `audit/retrieval/pagination_test/PAGINATION_TEST.md`.

---

## Next retrieval step

Pagination has been exercised on `NLVPMVATV` (page size 50 vs 10,000, interrupt/resume). Continue staged pulls of the **54,448** sequences, keeping **all original epitope-ID mappings**. Log every request the same way (`ok` / `empty` / `fail` / `incomplete`, `Content-Range`, raw JSON + sha256). Failed/incomplete queries stay unresolved; they are not “no assays.” Do not define V2 labels from stage 1.

Scripts: `scripts/iedb_staged_retrieve.py`, `scripts/iedb_annotate.py`, `scripts/summarize_stage1.py`, `scripts/summarize_task_candidates.py`, `scripts/test_iedb_pagination.py`. Spec: `audit/RETRIEVAL_SPEC.md`.
