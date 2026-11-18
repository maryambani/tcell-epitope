# Metadata-recovery pilot and V1 follow-ups

This is the next deliverable: **not a new model**. V1 files were not overwritten.

Retrieval date: **2026-11-07**. These assay records are the **current** IEDB, not a reconstruction of the 2025-11-03 export. Original query membership is stored separately.

---

## 1. Benchmark V1 is frozen

See `audit/v1/V1_BENCHMARK.md`.

> Peptide-level classification derived from positive/negative T-cell query membership, excluding peptides present in both lists, evaluated using a stratified random holdout.

---

## 2. Conflict counts (reconciled)

| Quantity | N | Meaning |
|----------|--:|---------|
| Unique peptides in both pos and neg lists, any sequence | **26,504** | Before amino-acid validation |
| Unique **valid-AA** peptides removed as conflicts | **26,204** | `226,064 − 199,860` |
| Rows removed at the conflict stage | **63,301** | `306,431 − 243,130` (valid-AA rows) |
| Valid 8–11 peptides **before** conflict drop | **54,448** | Target pool for metadata recovery |
| Of those, conflicts | **6,846** | |
| Of those, V1-retained (modeling set) | **47,602** | |

26,504 is not the valid-AA 8–11 conflict count. Use the table above.

---

## 3. Query wording (corrected)

The 2025-11-03 downloads used **Linear peptides, T Cell Assay, Outcome = Positive/Negative**.

> The query-level restriction is documented, but assay subtypes and experimental contexts cannot be verified from the retained epitope-level exports.

Missing assay columns do **not** prove that binding-only or elution-only records were mixed in.

---

## 4. Source field (corrected)

V1 `source=IEDB` after cleaning is **which duplicate was kept**, not exclusive membership.

In the 8–11 pre-conflict pool:

- 14,471 V1-retained peptides are in **both** IEDB and CEDAR
- **0** CEDAR-only 8–11 valid peptides
- After V1 `keep=first`, CEDAR rows disappear as a source tag

> CEDAR supplied overlapping records and participated in conflict handling; it did not contribute additional retained unique peptides after the documented cleaning stages.

Use `in_iedb` / `in_cedar` / query-membership flags (`audit/epitope_pool_8_11_preconflict.csv`).

---

## 5. Data dictionary

See `audit/DATA_DICTIONARY.md`:

- **A. Observation table**: one row per T-cell assay
- **B. Modeling table**: versioned peptide view after an explicit policy

No new conflict rule was applied.

---

## 6. Pilot answers

Sample: **46 peptides** (retained pos/neg, both databases, IEDB-only conflicts, both-database conflicts).  
Linking: epitope numeric ID → IEDB IQ-API `tcell_search?structure_id=eq.{id}`.  
CEDAR IDs matched the same IEDB `structure_id`.

| Question | Pilot result |
|----------|----------------|
| 1. Can assays be linked to original epitope IDs? | **Yes. 46/46 peptides** returned assays. |
| 2. Are required fields recoverable? | **Yes, at high coverage among 1,110 assays:** outcome, assay name, host, reference 100%; MHC class/allele ~99.5%; PubMed ~97.7%; immunization text ~98.8%. |
| 3. Do retained vs conflict peptides differ? | **Yes, in this sample.** Retained negatives were still all `Negative` in current IEDB. Conflicts had mixed Positive/Negative assays. Both-database conflicts were **heavily studied** (988 assays / 8 peptides, 385 studies) vs retained-pos IEDB-only (33 assays / 10 peptides). |
| 4. Coverage under plausible filters? | Not a full-pool count yet. Even this sample mixes **human and mouse**, **MHC I and II** (length 8–11 is not class I). Full retrieval is needed before choosing human-only or class-I-only. |
| 5. Enough studies / sequence groups? | Pilot alone has **430 distinct `reference_id`s**. Study hold-out looks **possible** after full retrieval. Similarity: see below; Hamming≤1 components max size 72, not one giant blob. |

Artifacts: `metadata_pilot_observations.csv`, `metadata_pilot_link_status.csv`, `metadata_pilot_summary.json`.

---

## 7. Similarity audit (full V1 test set, not a cutoff choice)

Quickcheck 8.1% was a **1,500-peptide subsample** (`numpy default_rng(42)`). Full test (**n = 9,521**):

| Check | Result |
|-------|--------|
| Exact overlap | 0 |
| Min same-length Hamming ≤ 1 | **761 / 9,521 = 7.99%** |
| Among those, only same-label neighbors | 586 |
| Among those, only opposite-label neighbors | 113 |
| By length, frac Hamming ≤ 1 | 8-mer 10.2%; 9-mer 8.7%; 10-mer 6.7%; 11-mer 5.6% |
| Different-length substring vs train | **2,250 / 9,521 = 23.6%** |
| Hamming≤1 connected components (all 47,602) | 44,673 components; max size **72**; 1,617 components size ≥ 2 (4,546 peptides) |

Same-length Hamming misses length variants. Substring relationships are common. Chaining Hamming≤1 does **not** produce a single giant component, but a grouping rule is **not** selected here.

Keep the random split as a **reference**. A similarity-separated evaluation would answer a different question, not “correct” V1 AUC.

Files: `similarity_summary.json`, `similarity_nearest_train.csv`.

---

## 8. Still undecided (needs full observation table)

| Decision | Why the pilot is not enough |
|----------|-----------------------------|
| Human-only vs mixed host | Host mix exists; need full counts |
| Confirmed class I vs length-only | Class II appears in 8–11mers |
| One endpoint vs pooled T-cell outcomes | Assay names vary; need distributions |
| Conflict policy | Disagreement is real in current assays; context (allele/host/study) not yet tabulated per peptide |
| Similarity grouping rule | Hamming vs substring vs antigen region still open |
| Study hold-out design | Need positives/negatives **per study** on the full pool |

**Stage 1 (done, 2026-11-07):** broader 240-peptide retrieval with completeness logs, endpoint family, MHC-resolution flags, and stratum summaries. See `audit/retrieval/stage1/STAGE1_REPORT.md`. V2 labels were not defined.

**Next science step:** continue staged retrieval of the full **54,448** pre-conflict 8–11 peptides (not only V1’s 47,602). Do not extrapolate runtime from the pilot or from stage 1. Still no transformer and no new conflict threshold.
