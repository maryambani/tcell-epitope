"""Reconcile the merged retrieval collection. Not a V2 decision."""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FULL = ROOT / "audit" / "retrieval" / "full"
POOL = ROOT / "audit" / "epitope_pool_8_11_preconflict.csv"

COMPARE_FIELDS = (
    "linear_sequence",
    "qualitative_measure",
    "assay_names",
    "host_organism_name",
    "mhc_class",
    "mhc_allele_name",
    "reference_id",
    "pubmed_id",
)


def read_csv(path: Path) -> list[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: list[str] = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def latest_log_by_id(logs: list[dict]) -> dict[str, dict]:
    out = {}
    for r in logs:
        eid = str(r.get("epitope_id") or "")
        if eid:
            out[eid] = r
    return out


def main() -> None:
    pool = read_csv(POOL)
    manifest = read_csv(FULL / "manifest.csv")
    logs = read_csv(FULL / "completeness_log.csv")
    obs = read_csv(FULL / "observations.csv")
    unresolved = read_csv(FULL / "unresolved_queries.csv")

    pool_ids = []
    peptide_to_ids = {}
    for r in pool:
        ids = [x for x in (r.get("epitope_ids") or "").split("|") if x]
        peptide_to_ids[r["peptide"]] = ids
        pool_ids.extend(ids)
    manifest_ids = [r["epitope_id"] for r in manifest]
    latest = latest_log_by_id(logs)

    missing_from_manifest = sorted(set(pool_ids) - set(manifest_ids))
    extra_in_manifest = sorted(set(manifest_ids) - set(pool_ids))
    never_requested = sorted(set(pool_ids) - set(latest))
    status_counts = Counter((latest[i].get("status") if i in latest else "never_requested") for i in pool_ids)

    obs_by_eid = defaultdict(list)
    for o in obs:
        obs_by_eid[str(o.get("queried_epitope_id"))].append(o)
    peptide_link_issues = []
    for pep, ids in peptide_to_ids.items():
        if not ids:
            peptide_link_issues.append({"peptide": pep, "issue": "no_epitope_ids"})
            continue
        for eid in ids:
            for o in obs_by_eid.get(str(eid), []):
                if o.get("export_peptide") != pep:
                    peptide_link_issues.append(
                        {
                            "peptide": pep,
                            "epitope_id": eid,
                            "issue": "observation_peptide_mismatch",
                            "export_peptide": o.get("export_peptide"),
                        }
                    )
                    break

    by_tcell = defaultdict(list)
    for o in obs:
        by_tcell[str(o.get("tcell_id"))].append(o)
    dup_ids = {k: v for k, v in by_tcell.items() if len(v) > 1}
    inconsistent = []
    linkage_only = []
    for tid, rows in dup_ids.items():
        payloads = [{f: r.get(f) for f in COMPARE_FIELDS} for r in rows]
        queries = {(r.get("queried_epitope_id"), r.get("export_peptide"), r.get("retrieval_batch")) for r in rows}
        same_payload = all(p == payloads[0] for p in payloads)
        if same_payload and len(queries) > 1:
            linkage_only.append(
                {
                    "tcell_id": tid,
                    "n_rows": len(rows),
                    "n_query_linkages": len(queries),
                    "linkages": "|".join(f"{a}:{b}:{c}" for a, b, c in sorted(queries, key=lambda x: str(x))),
                }
            )
        elif not same_payload:
            inconsistent.append(
                {
                    "tcell_id": tid,
                    "n_rows": len(rows),
                    "fields_that_differ": "|".join(
                        f
                        for f in COMPARE_FIELDS
                        if len({r.get(f) for r in rows}) > 1
                    ),
                }
            )

    empty = [latest[i] for i in pool_ids if i in latest and latest[i].get("status") == "empty"]
    fail = [latest[i] for i in pool_ids if i in latest and latest[i].get("status") == "fail"]
    incomplete = [latest[i] for i in pool_ids if i in latest and latest[i].get("status") == "incomplete"]
    ok = [latest[i] for i in pool_ids if i in latest and latest[i].get("status") == "ok"]

    seq_mismatch = sum(1 for r in obs if str(r.get("sequence_matches_export")) != "1")
    structure_mismatch = sum(
        1 for r in obs if str(r.get("queried_epitope_id")) != str(r.get("returned_structure_id"))
    )

    n_pool = len(pool_ids)
    n_resolved = len(ok) + len(empty)
    complete = (
        not missing_from_manifest
        and not extra_in_manifest
        and not never_requested
        and not fail
        and not incomplete
        and seq_mismatch == 0
        and not inconsistent
    )
    report = {
        "v2_labels_defined": False,
        "retrieval_described_as_complete": complete,
        "collection": "historical_8_11_preconflict_54448",
        "not_a_census_of_current_iedb": True,
        "n_pool_peptides": len(pool),
        "n_pool_epitope_ids": n_pool,
        "n_manifest_rows": len(manifest),
        "n_missing_ids_from_manifest": len(missing_from_manifest),
        "n_extra_ids_in_manifest": len(extra_in_manifest),
        "n_never_requested": len(never_requested),
        "status_latest_per_id": dict(status_counts),
        "n_ok": len(ok),
        "n_empty": len(empty),
        "n_fail": len(fail),
        "n_incomplete": len(incomplete),
        "n_historically_unresolved_file_rows": len(unresolved),
        "n_currently_unresolved": len(fail) + len(incomplete),
        "unresolved_note": (
            "unresolved_queries.csv and completeness_log.csv are append-only history. "
            "n_historically_unresolved_file_rows counts every fail/incomplete row ever written, "
            "including attempts later superseded by a successful retry. "
            "n_currently_unresolved counts pool epitope IDs whose LATEST log row is fail or "
            "incomplete; only that number reflects outstanding retrieval gaps."
        ),
        "n_observations": len(obs),
        "n_unique_tcell_id": len(by_tcell),
        "n_tcell_ids_with_duplicate_rows": len(dup_ids),
        "n_duplicate_rows_same_payload_multiple_linkages": len(linkage_only),
        "n_inconsistent_duplicate_payloads": len(inconsistent),
        "n_sequence_mismatches": seq_mismatch,
        "n_structure_id_mismatches": structure_mismatch,
        "n_peptide_link_issues": len(peptide_link_issues),
        "frac_ids_resolved_ok_or_empty": round(n_resolved / max(n_pool, 1), 4),
        "empty_is_not_fail": True,
        "duplicates_not_discarded": True,
    }
    write_csv(FULL / "duplicate_assay_linkages.csv", linkage_only)
    write_csv(FULL / "inconsistent_assay_payloads.csv", inconsistent)
    write_csv(FULL / "peptide_link_issues.csv", peptide_link_issues[:5000])
    (FULL / "merge_validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not complete:
        print("MERGED COLLECTION NOT COMPLETE")
    else:
        print("merged collection checks passed for this defined collection")


if __name__ == "__main__":
    main()
