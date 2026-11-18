# Full-pool observation recovery (in progress)

V1 is frozen. No V2 labels, conflict rule, split, or model.

Collection: the **54,448** valid 8–11 peptides in `epitope_pool_8_11_preconflict.csv`. Finishing this pull is complete only for that historical set. It is not a census of current IEDB and does not measure performance on peptides that were never in the 2025-11-03 lists.

## Reporting corrections (Stage 1 / 2)

`WICLLQFAY` under human + class I + functional has **no eligible evidence** (0 functional assays). The mixed tetramer/binding records were filtered out. That is not a usable functional label.

`VQSAKYWNL` under a human-host subset has **no eligible evidence**. Mouse-strain ELISPOT disagreement cannot supply a human-response label.

Context counts are peptides and groups separately. Stage 1: **96 mixed peptides**, **14 peptides** with at least one mixed complete-context key, **17 distinct mixed keys**. Stage 2: **98 mixed peptides**, **14 peptides**, **21 mixed keys**. Incomplete-context unmatched peptides are kept separate from unmatched peptides that have complete context on both sides.

Why current IEDB disagrees with 2025-11-03 query membership remains unresolved (database change, annotation change, or export difference).

## Retrieval

`scripts/iedb_pool_retrieve.py` writes `audit/retrieval/full/`.

- Manifest: every peptide and epitope ID (54,448; one ID each in this pool).
- Stage 1 and Stage 2 `ok` IDs were imported and are not re-queried.
- Remaining **53,968** IDs are being requested with `order=tcell_id`, page size 1,000, `Prefer: count=exact`.
- Completeness log, raw JSON, timestamps, and unresolved queries (`fail`/`incomplete` ≠ empty) are appended as it goes.
- Restart: skip `ok`/`empty`; retry `fail`/`incomplete`; warn if page size or API changed; record `reported_total` changes instead of calling them pagination success or failure.

Progress: `progress.json`. Resume with the same command.

## Merge checks and feasibility

Run after retrieval (or anytime for a snapshot):

```
./venv/bin/python scripts/validate_merged_retrieval.py
./venv/bin/python scripts/task_feasibility.py
```

Do not describe retrieval as complete until the validator says every pool ID is `ok` or `empty`, duplicates are inventoried (not discarded), and there are no inconsistent payloads for the same assay ID.

Predeclared intersections in `task_feasibility.py` are for **feasibility**, not a chosen V2 task. Do not pick the definition that would later give the highest score.

## Intended prediction claim (parallel, not chosen)

A documented primary task still has to name:

1. Who (host).
2. Which MHC class / restriction evidence is required.
3. Which measurement counts as the endpoint (functional vs recognition vs one exact assay).
4. How mixed eligible observations become a peptide label, or whether mixed peptides are left unlabeled.
5. What “new” means at evaluation (new study, new sequence, or both), including multi-study peptides.

Until those are written down, this workstream only recovers observations and counts.

Scripts: `iedb_pool_retrieve.py`, `validate_merged_retrieval.py`, `task_feasibility.py`.
