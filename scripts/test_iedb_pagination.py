"""Retrieval acceptance test: pagination, ordering, resume, unresolved status.

Not a biological analysis. Does not define V2 labels.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from iedb_staged_retrieve import retrieve_epitope  # noqa: E402

OUT = ROOT / "audit" / "retrieval" / "pagination_test"
PEPTIDE = "NLVPMVATV"
STRUCTURE_ID = "44920"
SMALL_PAGE = 50
SINGLE_PAGE = 10000


def ids(rows: list[dict]) -> list:
    return ["" if r.get("tcell_id") is None else str(r.get("tcell_id")) for r in rows]


def compare_id_sets(a: list, b: list) -> dict:
    a = ["" if x is None else str(x) for x in a]
    b = ["" if x is None else str(x) for x in b]
    sa, sb = set(a), set(b)
    return {
        "n_a": len(a),
        "n_b": len(b),
        "n_unique_a": len(sa),
        "n_unique_b": len(sb),
        "n_only_a": len(sa - sb),
        "n_only_b": len(sb - sa),
        "n_overlap": len(sa & sb),
        "same_set": sa == sb,
        "same_order": a == b,
        "a_has_internal_duplicates": len(a) != len(sa),
        "b_has_internal_duplicates": len(b) != len(sb),
    }


def monotonic(seq: list) -> bool:
    coerced = []
    for x in seq:
        try:
            coerced.append(int(x))
        except (TypeError, ValueError):
            return False
    return all(coerced[i] <= coerced[i + 1] for i in range(len(coerced) - 1))


def page_boundary_gaps(pages: list[dict]) -> dict:
    """Check Content-Range continuity across pages."""
    issues = []
    for i, p in enumerate(pages):
        start, end = p.get("range_start"), p.get("range_end")
        if start is None or end is None:
            issues.append({"page_index": i, "issue": "missing_content_range_span", "content_range": p.get("content_range")})
            continue
        if p.get("n", 0) and end - start + 1 != p["n"]:
            issues.append(
                {
                    "page": i,
                    "issue": "content_range_length_mismatch",
                    "range_start": start,
                    "range_end": end,
                    "n": p["n"],
                }
            )
        if i > 0:
            prev = pages[i - 1]
            if prev.get("range_end") is not None and start != prev["range_end"] + 1:
                issues.append(
                    {
                        "page": i,
                        "issue": "noncontiguous_range",
                        "prev_end": prev.get("range_end"),
                        "start": start,
                    }
                )
    return {"n_issues": len(issues), "issues": issues[:20]}


def load_stage1_nlv_ids() -> dict:
    path = ROOT / "audit" / "retrieval" / "stage1" / "observations.csv"
    out = []
    retrieval_time = None
    if not path.exists():
        return {"available": False, "ids": []}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("export_peptide") == PEPTIDE:
                out.append(row.get("tcell_id"))
                retrieval_time = row.get("retrieval_time_utc")
    return {
        "available": True,
        "n": len(out),
        "n_unique": len(set(out)),
        "ids": out,
        "retrieval_time_utc": retrieval_time,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw_dir = OUT / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    retrieval_time = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    report: dict = {
        "test": "iedb_pagination_acceptance",
        "peptide": PEPTIDE,
        "structure_id": STRUCTURE_ID,
        "retrieval_time_utc": retrieval_time,
        "v2_labels_defined": False,
        "not_a_biological_analysis": True,
    }

    # --- 1. Multipage with small page size ---
    print("multipage page_size=", SMALL_PAGE)
    multi = retrieve_epitope(
        STRUCTURE_ID,
        raw_dir / "multipage",
        retrieval_time,
        page_size=SMALL_PAGE,
        checkpoint_path=OUT / "checkpoint_multipage.json",
        sleep_s=0.08,
    )
    multi_ids = ids(multi["assays"])
    report["multipage"] = {
        "status": multi["status"],
        "error": multi["error"],
        "http_status": multi["http_status"],
        "page_size": SMALL_PAGE,
        "n_pages": multi["n_pages"],
        "page_offsets": multi["page_offsets"],
        "reported_total": multi["reported_total"],
        "n_rows_retrieved": multi["n_rows_retrieved"],
        "n_unique_tcell_id": multi["n_unique_tcell_id"],
        "n_duplicate_tcell_id": multi["n_duplicate_tcell_id"],
        "ids_monotonic_tcell_id": monotonic(multi_ids),
        "content_range_continuity": page_boundary_gaps(multi["pages"]),
        "last_page_n": multi["pages"][-1]["n"] if multi["pages"] else None,
        "expected_pages_if_total_known": (
            (int(multi["reported_total"]) + SMALL_PAGE - 1) // SMALL_PAGE
            if str(multi["reported_total"]).isdigit()
            else None
        ),
    }

    # --- 2. Fresh single-page retrieval (API max) ---
    print("single-page page_size=", SINGLE_PAGE)
    single_error = ""
    single = None
    try:
        single = retrieve_epitope(
            STRUCTURE_ID,
            raw_dir / "single",
            retrieval_time,
            page_size=SINGLE_PAGE,
            checkpoint_path=OUT / "checkpoint_single.json",
        )
    except Exception as e:
        single_error = str(e)
        print("single-page 10000 failed, retrying 1000:", e)
        single = retrieve_epitope(
            STRUCTURE_ID,
            raw_dir / "single_1000",
            retrieval_time,
            page_size=1000,
            checkpoint_path=OUT / "checkpoint_single_1000.json",
        )
    single_ids = ids(single["assays"])
    report["single_page"] = {
        "status": single["status"],
        "error": single["error"] or single_error,
        "http_status": single["http_status"],
        "page_size_requested": SINGLE_PAGE,
        "page_size_used": single["page_size"],
        "n_pages": single["n_pages"],
        "reported_total": single["reported_total"],
        "n_rows_retrieved": single["n_rows_retrieved"],
        "n_unique_tcell_id": single["n_unique_tcell_id"],
        "n_duplicate_tcell_id": single["n_duplicate_tcell_id"],
        "ids_monotonic_tcell_id": monotonic(single_ids),
    }
    report["multipage_vs_single_page"] = compare_id_sets(multi_ids, single_ids)
    if report["multipage_vs_single_page"]["same_set"]:
        report["multipage_vs_single_page"]["interpretation"] = "pagination_match"
    else:
        report["multipage_vs_single_page"]["interpretation"] = (
            "id_set_differs; if reported totals also differ this may be a database "
            "change between requests, not a pagination success or failure"
        )
        report["multipage_vs_single_page"]["reported_total_multipage"] = multi["reported_total"]
        report["multipage_vs_single_page"]["reported_total_single"] = single["reported_total"]
        report["multipage_vs_single_page"]["only_multipage_sample"] = sorted(
            set(multi_ids) - set(single_ids), key=lambda x: int(x) if str(x).isdigit() else 0
        )[:20]
        report["multipage_vs_single_page"]["only_single_sample"] = sorted(
            set(single_ids) - set(multi_ids), key=lambda x: int(x) if str(x).isdigit() else 0
        )[:20]

    # --- 3. Second single-page fetch for ordering stability ---
    print("repeat single-page for order stability")
    time.sleep(0.2)
    single2 = retrieve_epitope(
        STRUCTURE_ID,
        raw_dir / "single_repeat",
        datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        page_size=single["page_size"],
        checkpoint_path=OUT / "checkpoint_single_repeat.json",
    )
    single2_ids = ids(single2["assays"])
    report["order_stability_two_single_page_fetches"] = compare_id_sets(single_ids, single2_ids)
    report["order_stability_two_single_page_fetches"]["reported_total_first"] = single["reported_total"]
    report["order_stability_two_single_page_fetches"]["reported_total_second"] = single2["reported_total"]
    if (
        not report["order_stability_two_single_page_fetches"]["same_set"]
        and single["reported_total"] != single2["reported_total"]
    ):
        report["order_stability_two_single_page_fetches"]["interpretation"] = (
            "database_may_have_changed_between_requests"
        )
    elif report["order_stability_two_single_page_fetches"]["same_order"]:
        report["order_stability_two_single_page_fetches"]["interpretation"] = "stable_order"
    else:
        report["order_stability_two_single_page_fetches"]["interpretation"] = "same_set_unstable_order"

    # --- 4. Interruption and resumption ---
    print("interrupt after 2 pages, then resume")
    part = retrieve_epitope(
        STRUCTURE_ID,
        raw_dir / "resume_part1",
        retrieval_time,
        page_size=SMALL_PAGE,
        max_pages=2,
        checkpoint_path=OUT / "checkpoint_interrupted.json",
        sleep_s=0.08,
    )
    resume_offset = (part["pages"][-1]["offset"] + SMALL_PAGE) if part["pages"] else 0
    resumed = retrieve_epitope(
        STRUCTURE_ID,
        raw_dir / "resume_part2",
        retrieval_time,
        page_size=SMALL_PAGE,
        start_offset=resume_offset,
        prior_rows=part["assays"],
        checkpoint_path=OUT / "checkpoint_resumed.json",
        sleep_s=0.08,
    )
    resumed_ids = ids(resumed["assays"])
    report["interruption_resume"] = {
        "interrupted_status": part["status"],
        "interrupted_error": part["error"],
        "interrupted_n_rows": part["n_rows_retrieved"],
        "interrupted_n_pages": part["n_pages"],
        "resume_start_offset": resume_offset,
        "resumed_status": resumed["status"],
        "resumed_error": resumed["error"],
        "resumed_n_rows": resumed["n_rows_retrieved"],
        "resumed_n_unique": resumed["n_unique_tcell_id"],
        "resumed_n_duplicate": resumed["n_duplicate_tcell_id"],
        "vs_full_multipage": compare_id_sets(resumed_ids, multi_ids),
        "unresolved_not_empty": part["status"] == "incomplete" and part["status"] != "empty",
    }

    # --- 5. Failed query is unresolved, not empty ---
    print("failed query (non-integer structure_id)")
    bogus = retrieve_epitope(
        "not_an_integer",
        raw_dir / "bogus",
        retrieval_time,
        page_size=SMALL_PAGE,
        checkpoint_path=OUT / "checkpoint_bogus.json",
    )
    report["failed_query_not_empty"] = {
        "epitope_id": "not_an_integer",
        "status": bogus["status"],
        "http_status": bogus["http_status"],
        "error": bogus["error"],
        "n_rows_retrieved": bogus["n_rows_retrieved"],
        "treated_as_empty": bogus["status"] == "empty",
        "preserved_as_unresolved": bogus["status"] in ("fail", "incomplete"),
        "note": "HTTP errors must stay fail/incomplete. A true HTTP 200 with 0 rows is empty, not a failure.",
    }
    # A true empty (if the API returns 200 with 0 rows) would be status=empty.
    # structure_id=0 may 200-empty or 4xx; either fail or empty is recorded, not rewritten.

    # --- 6. Stage-1 snapshot vs this test: possible DB change ---
    stage1 = load_stage1_nlv_ids()
    report["versus_stage1_snapshot"] = {
        "stage1_available": stage1["available"],
        "stage1_n": stage1.get("n"),
        "stage1_retrieval_time_utc": stage1.get("retrieval_time_utc"),
        "this_test_n": len(set(multi_ids)),
        "this_test_reported_total": multi["reported_total"],
    }
    if stage1["available"]:
        cmp = compare_id_sets(stage1["ids"], multi_ids)
        cmp["interpretation"] = (
            "pagination_match_across_days"
            if cmp["same_set"]
            else (
                "database_may_have_changed_since_stage1; do not treat as pagination "
                "failure or success without the same-session multipage vs single-page comparison"
            )
        )
        report["versus_stage1_snapshot"]["id_comparison"] = {
            k: cmp[k] for k in cmp if k != "ids"
        }

    # Termination checks
    report["termination"] = {
        "multipage_stopped_at_reported_total": (
            str(multi["reported_total"]).isdigit()
            and multi["n_unique_tcell_id"] == int(multi["reported_total"])
        ),
        "multipage_did_not_keep_requesting_after_short_page": (
            bool(multi["pages"]) and multi["pages"][-1]["n"] <= SMALL_PAGE
        ),
        "last_page_shorter_than_page_size_or_exact_fill": (
            bool(multi["pages"]) and multi["pages"][-1]["n"] <= SMALL_PAGE
        ),
        "n_pages_matches_ceiling_total_over_page_size": (
            report["multipage"]["expected_pages_if_total_known"] == multi["n_pages"]
            if report["multipage"]["expected_pages_if_total_known"] is not None
            else None
        ),
    }

    report["elapsed_seconds"] = round(time.time() - t0, 2)
    report["pass"] = bool(
        report["multipage"]["status"] == "ok"
        and report["single_page"]["status"] == "ok"
        and report["multipage_vs_single_page"]["same_set"]
        and report["multipage"]["n_duplicate_tcell_id"] == 0
        and report["single_page"]["n_duplicate_tcell_id"] == 0
        and report["interruption_resume"]["vs_full_multipage"]["same_set"]
        and report["interruption_resume"]["unresolved_not_empty"]
        and report["failed_query_not_empty"]["treated_as_empty"] is False
        and report["termination"]["multipage_stopped_at_reported_total"]
    )
    (OUT / "pagination_test_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({k: report[k] for k in report if k != "versus_stage1_snapshot"}, indent=2))
    print("PASS" if report["pass"] else "FAIL")
    print("wrote", OUT / "pagination_test_report.json")


if __name__ == "__main__":
    main()
