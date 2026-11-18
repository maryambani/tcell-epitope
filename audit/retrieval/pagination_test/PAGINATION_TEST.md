# Pagination and restart acceptance test

This is a retrieval test, not a biological analysis. Peptide: `NLVPMVATV` (`structure_id=44920`). Time: **2026-11-07T18:09:00Z**. **PASS.**

V1 was not modified. No V2 labels were defined.

## What was tested

1. Multipage retrieval with **page size 50** (19 pages; last page 20 rows).
2. Fresh single-page retrieval with **page size 10,000** (920 < 10,000).
3. A second single-page fetch for order stability.
4. Stop after 2 pages, then resume from offset 100.
5. A bad query (`structure_id=not_an_integer`) must stay **fail**, not empty.

## Results

| Check | Result |
|-------|--------|
| Multipage unique assay IDs | 920 / reported total 920 |
| Internal duplicates | 0 |
| `tcell_id` order monotonic | yes |
| Content-Range continuity across pages | 0 issues |
| Page count = ceil(920/50) | 19 |
| Multipage ID set vs same-session single page | **identical, same order** |
| Two consecutive single-page fetches | identical set and order |
| Interrupted status | `incomplete` (100 rows), **not** `empty` |
| Resume vs full multipage | identical 920 IDs |
| HTTP 400 bad ID | status `fail`, not `empty` |
| Stage-1 snapshot (same morning) vs this test | same 920 IDs after coercing IDs to strings |

HTTP 206 on paged requests is expected with `Content-Range`. Same-session multipage and single-page totals were both 920. If those ever differ, record a possible database change. Do not call that a pagination success or failure.

A first draft compared Stage-1 CSV strings to JSON integers and falsely reported zero overlap. Coerced IDs match. The test script now stringifies `tcell_id` before set comparison.

## What this does not prove

- It does not prove that every future 54,448 query will paginate correctly under load or after API changes.
- It does not prove exhaustive retrieval for peptides that were not tested.
- Failed and incomplete queries must remain in `unresolved_queries.csv`. They are not “no assays.”

JSON: `pagination_test_report.json`. Script: `scripts/test_iedb_pagination.py`.
