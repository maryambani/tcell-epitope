"""Staged IEDB T-cell assay retrieval with pagination and completeness logs.

Does not overwrite V1. Does not assign V2 modeling labels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from iedb_annotate import (  # noqa: E402
    disagreement_class,
    endpoint_family,
    evidence_bucket,
    host_is_human,
    host_is_transgenic_mouse,
    mhc_modified,
    mhc_specificity,
    nonempty,
    outcome_family,
)

AUDIT = ROOT / "audit"
POOL = AUDIT / "epitope_pool_8_11_preconflict.csv"
API = "https://query-api.iedb.org/tcell_search"
PAGE = 1000
SELECT = ",".join(
    [
        "tcell_id",
        "tcell_iri",
        "structure_id",
        "structure_iri",
        "linear_sequence",
        "e_modification",
        "qualitative_measure",
        "assay_names",
        "assay_description",
        "host_organism_name",
        "host_organism_iri",
        "mhc_class",
        "mhc_allele_name",
        "mhc_restriction",
        "mhc_allele_evidence",
        "reference_id",
        "reference_iri",
        "pubmed_id",
        "immunization_description",
        "direct_ex_vivo_bool",
        "parent_source_antigen_name",
        "source_organism_name",
    ]
)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def peptide_stratum(row: dict) -> str:
    if row.get("v1_retained") == "1" and row.get("v1_query_label") == "1" and row.get("in_cedar") == "0":
        return "retained_pos_iedb_only"
    if row.get("v1_retained") == "1" and row.get("v1_query_label") == "0" and row.get("in_cedar") == "0":
        return "retained_neg_iedb_only"
    if row.get("v1_retained") == "1" and row.get("in_iedb") == "1" and row.get("in_cedar") == "1":
        return "retained_both_databases"
    if row.get("v1_conflict") == "1" and row.get("in_cedar") == "0":
        return "conflict_iedb_only"
    if row.get("v1_conflict") == "1" and row.get("in_iedb") == "1" and row.get("in_cedar") == "1":
        return "conflict_both_databases"
    if row.get("in_iedb") == "0" and row.get("in_cedar") == "1":
        return "cedar_only"
    return "other"


def draw_stage_sample(
    pool: list[dict],
    n_per: int,
    seed: int,
    exclude_peptides: set[str] | None = None,
) -> list[dict]:
    rng = random.Random(seed)
    skip = exclude_peptides or set()
    by = defaultdict(list)
    for r in pool:
        if r.get("peptide") in skip:
            continue
        by[peptide_stratum(r)].append(r)
    out = []
    wanted = [
        "retained_pos_iedb_only",
        "retained_neg_iedb_only",
        "retained_both_databases",
        "conflict_iedb_only",
        "conflict_both_databases",
    ]
    for tag in wanted:
        cands = list(by[tag])
        rng.shuffle(cands)
        for row in cands[:n_per]:
            item = dict(row)
            item["stage_stratum"] = tag
            out.append(item)
    return out


def parse_content_range(header: str | None) -> tuple[int | None, int | None, int | None]:
    if not header:
        return None, None, None
    # e.g. 0-999/1234 or 0-9/*
    try:
        span, total_s = header.split("/")
        start_s, end_s = span.split("-")
        total = None if total_s.strip() == "*" else int(total_s)
        return int(start_s), int(end_s), total
    except Exception:
        return None, None, None


def fetch_page(
    structure_id: str,
    offset: int,
    page_size: int = PAGE,
    order: str = "tcell_id",
    timeout: int = 90,
):
    params = {
        "structure_id": f"eq.{structure_id}",
        "select": SELECT,
        "limit": str(page_size),
        "offset": str(offset),
    }
    if order:
        params["order"] = order
    url = API + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "Prefer": "count=exact",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()
        headers = {k.lower(): v for k, v in resp.headers.items()}
        data = json.loads(raw.decode("utf-8"))
        return resp.status, headers, data if isinstance(data, list) else [], raw


def _write_checkpoint(path: Path | None, payload: dict) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def retrieve_epitope(
    structure_id: str,
    raw_dir: Path,
    retrieval_time: str,
    page_size: int = PAGE,
    start_offset: int = 0,
    prior_rows: list | None = None,
    max_pages: int | None = None,
    checkpoint_path: Path | None = None,
    order: str = "tcell_id",
    sleep_s: float = 0.12,
) -> dict:
    """Fetch all pages for one structure_id.

    fail/incomplete are unresolved. They are never rewritten as empty.
    0/False-length successful payloads with HTTP 200 are status=empty.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = list(prior_rows or [])
    seen = {r.get("tcell_id") for r in rows}
    pages = []
    reported_total = None
    error = ""
    status_name = "ok"
    http_status = None
    offset = start_offset
    n_new_pages = 0
    interrupted = False
    repeated_page = False
    try:
        while True:
            http_status, headers, data, raw = fetch_page(
                structure_id, offset, page_size=page_size, order=order
            )
            start, end, tot = parse_content_range(headers.get("content-range"))
            if tot is not None:
                reported_total = tot
            digest = hashlib.sha256(raw).hexdigest()
            raw_path = raw_dir / f"structure_{structure_id}_limit_{page_size}_offset_{offset}.json"
            raw_path.write_bytes(raw)
            page_ids = [r.get("tcell_id") for r in data]
            n_new_ids = sum(1 for i in page_ids if i not in seen)
            page_info = {
                "offset": offset,
                "n": len(data),
                "n_new_ids": n_new_ids,
                "sha256": digest,
                "content_range": headers.get("content-range", ""),
                "range_start": start,
                "range_end": end,
                "reported_total": tot,
                "first_tcell_id": page_ids[0] if page_ids else None,
                "last_tcell_id": page_ids[-1] if page_ids else None,
                "raw_path": str(raw_path.relative_to(ROOT)),
                "retrieval_time_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            }
            pages.append(page_info)
            rows.extend(data)
            n_new_pages += 1
            if data and n_new_ids == 0:
                repeated_page = True
            seen.update(page_ids)
            _write_checkpoint(
                checkpoint_path,
                {
                    "structure_id": structure_id,
                    "page_size": page_size,
                    "order": order,
                    "next_offset": offset + page_size,
                    "n_rows": len(rows),
                    "n_unique": len(seen),
                    "reported_total": reported_total,
                    "pages": pages,
                    "status": "in_progress",
                },
            )
            if repeated_page:
                break
            if reported_total is not None:
                if len(seen) >= reported_total or len(data) == 0:
                    break
            elif len(data) < page_size:
                break
            if max_pages is not None and n_new_pages >= max_pages:
                interrupted = True
                break
            offset += page_size
            time.sleep(sleep_s)
        if interrupted:
            status_name = "incomplete"
            error = f"stopped after max_pages={max_pages}"
        elif repeated_page:
            status_name = "incomplete"
            error = "page contained only already-seen tcell_ids"
        elif reported_total is not None and len(seen) < reported_total:
            status_name = "incomplete"
            error = f"retrieved unique {len(seen)} < reported_total {reported_total}"
        elif not rows:
            status_name = "empty"
    except urllib.error.HTTPError as e:
        status_name = "fail"
        http_status = e.code
        error = str(e)
    except Exception as e:
        status_name = "fail"
        error = str(e)

    unique_ids = [r.get("tcell_id") for r in rows]
    result = {
        "epitope_id": structure_id,
        "retrieval_time_utc": retrieval_time,
        "http_status": http_status,
        "status": status_name,
        "content_range": pages[0]["content_range"] if pages else "",
        "reported_total": reported_total if reported_total is not None else "",
        "n_rows_retrieved": len(rows),
        "n_unique_tcell_id": len(set(unique_ids)),
        "n_duplicate_tcell_id": len(unique_ids) - len(set(unique_ids)),
        "page_size": page_size,
        "page_offsets": ",".join(str(p["offset"]) for p in pages),
        "n_pages": len(pages),
        "error": error,
        "raw_paths": "|".join(p["raw_path"] for p in pages),
        "page_sha256": "|".join(p["sha256"] for p in pages),
        "start_offset": start_offset,
        "assays": rows,
        "pages": pages,
    }
    _write_checkpoint(
        checkpoint_path,
        {
            "structure_id": structure_id,
            "page_size": page_size,
            "order": order,
            "status": status_name,
            "error": error,
            "reported_total": reported_total,
            "n_rows": len(rows),
            "n_unique_tcell_id": result["n_unique_tcell_id"],
            "pages": pages,
        },
    )
    return result


def annotate_obs(assay: dict, pep_row: dict, queried_id: str, retrieval_time: str) -> dict:
    allele = assay.get("mhc_allele_name")
    klass = assay.get("mhc_class")
    seq = assay.get("linear_sequence") or ""
    return {
        "observation_id": assay.get("tcell_iri") or f"IEDB_ASSAY:{assay.get('tcell_id')}",
        "tcell_id": assay.get("tcell_id"),
        "tcell_iri": assay.get("tcell_iri"),
        "queried_epitope_id": queried_id,
        "returned_structure_id": assay.get("structure_id"),
        "structure_iri": assay.get("structure_iri"),
        "export_peptide": pep_row["peptide"],
        "linear_sequence": seq,
        "sequence_matches_export": int(seq == pep_row["peptide"]),
        "e_modification": assay.get("e_modification"),
        "qualitative_measure": assay.get("qualitative_measure"),
        "outcome_family": outcome_family(assay.get("qualitative_measure") or ""),
        "assay_names": assay.get("assay_names"),
        "assay_description": assay.get("assay_description"),
        "endpoint_family": endpoint_family(assay.get("assay_names") or ""),
        "host_organism_name": assay.get("host_organism_name"),
        "host_organism_iri": assay.get("host_organism_iri"),
        "host_is_human": host_is_human(assay.get("host_organism_name") or ""),
        "host_is_transgenic_mouse": host_is_transgenic_mouse(assay.get("host_organism_name") or ""),
        "mhc_class": klass,
        "mhc_allele_name": allele,
        "mhc_restriction": assay.get("mhc_restriction"),
        "mhc_allele_evidence": assay.get("mhc_allele_evidence"),
        "mhc_annotation_present": int(nonempty(klass) or nonempty(allele)),
        "mhc_specificity": mhc_specificity(allele or "", klass or ""),
        "mhc_modified_experiment": mhc_modified(allele or "", assay.get("mhc_restriction") or ""),
        "mhc_evidence_category": evidence_bucket(assay.get("mhc_allele_evidence") or ""),
        "reference_id": assay.get("reference_id"),
        "reference_iri": assay.get("reference_iri"),
        "pubmed_id": assay.get("pubmed_id"),
        "immunization_description": assay.get("immunization_description"),
        "direct_ex_vivo_bool": assay.get("direct_ex_vivo_bool"),
        "parent_source_antigen_name": assay.get("parent_source_antigen_name"),
        "source_organism_name": assay.get("source_organism_name"),
        "stage_stratum": pep_row.get("stage_stratum"),
        "v1_conflict": pep_row.get("v1_conflict"),
        "v1_query_label": pep_row.get("v1_query_label"),
        "retrieval_time_utc": retrieval_time,
        "retrieval_source": "iedb_query_api.tcell_search",
    }


def summarize(sample: list[dict], logs: list[dict], obs: list[dict], elapsed_s: float) -> dict:
    by_pep = defaultdict(list)
    for r in obs:
        by_pep[r["export_peptide"]].append(r)
    mixed = []
    for pep, rs in by_pep.items():
        if any(x["outcome_family"] == "positive" for x in rs) and any(x["outcome_family"] == "negative" for x in rs):
            mixed.append({"peptide": pep, "disagreement_class": disagreement_class(rs), "n": len(rs)})
    ok_peps = {r["peptide"] for r in logs if r.get("status") in ("ok", "empty")}
    unresolved_peps = {r["peptide"] for r in logs if r.get("status") in ("fail", "incomplete")}
    n_assays = [len(by_pep[p]) for p in ok_peps]
    return {
        "n_peptides_queried": len(sample),
        "n_epitope_ids_queried": sum(len([x for x in (r.get("epitope_ids") or "").split("|") if x]) for r in sample),
        "n_observations": len(obs),
        "n_unique_tcell_id": len({r["tcell_id"] for r in obs}),
        "log_status": dict(Counter(r["status"] for r in logs)),
        "elapsed_seconds": round(elapsed_s, 2),
        "seconds_per_epitope_id": round(elapsed_s / max(len(logs), 1), 3),
        "assays_per_peptide_min": min(n_assays) if n_assays else None,
        "assays_per_peptide_median": sorted(n_assays)[len(n_assays) // 2] if n_assays else None,
        "assays_per_peptide_max": max(n_assays) if n_assays else None,
        "top_peptide": max(by_pep, key=lambda p: len(by_pep[p])) if by_pep else None,
        "top_peptide_n": max((len(v) for v in by_pep.values()), default=0),
        "frac_obs_from_top_peptide": round(max((len(v) for v in by_pep.values()), default=0) / max(len(obs), 1), 4),
        "endpoint_family": dict(Counter(r["endpoint_family"] for r in obs)),
        "outcome_family": dict(Counter(r["outcome_family"] for r in obs)),
        "mhc_specificity": dict(Counter(r["mhc_specificity"] for r in obs)),
        "mhc_class": dict(Counter((r.get("mhc_class") or "MISSING") for r in obs)),
        "host_is_human_assays": sum(int(r["host_is_human"]) for r in obs),
        "n_mixed_outcome_peptides": len(mixed),
        "disagreement_classes": dict(Counter(m["disagreement_class"] for m in mixed)),
        "n_distinct_reference_id": len({r["reference_id"] for r in obs if nonempty(r.get("reference_id"))}),
        "n_unresolved_peptides": len(unresolved_peps),
        "do_not_extrapolate_from_pilot": True,
        "v2_labels_defined": False,
    }


def parse_args(argv: list[str] | None = None):
    p = argparse.ArgumentParser(description="Staged IEDB T-cell observation retrieval")
    p.add_argument("--stage", default="stage1")
    p.add_argument("--n-per", type=int, default=48)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--page-size", type=int, default=PAGE)
    p.add_argument(
        "--exclude-from",
        action="append",
        default=[],
        help="CSV with a peptide column; may be passed more than once",
    )
    p.add_argument(
        "--redraw",
        action="store_true",
        help="Draw a new sample even if the stage sample CSV already exists",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    out_dir = AUDIT / "retrieval" / args.stage
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    sample_path = out_dir / f"{args.stage}_sample.csv"

    pool = read_csv(POOL)
    exclude: set[str] = set()
    for path_s in args.exclude_from:
        for row in read_csv(Path(path_s)):
            if row.get("peptide"):
                exclude.add(row["peptide"])

    if sample_path.exists() and not args.redraw:
        sample = read_csv(sample_path)
        print("Reusing", sample_path, "n=", len(sample))
    else:
        sample = draw_stage_sample(pool, n_per=args.n_per, seed=args.seed, exclude_peptides=exclude)
        write_csv(sample_path, sample)
        print("Wrote", sample_path, "n=", len(sample), "excluded", len(exclude))

    t0 = time.time()
    retrieval_time = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    logs = []
    observations = []
    unresolved = []
    n_ids = 0
    for i, pep in enumerate(sample, 1):
        ids = [x for x in (pep.get("epitope_ids") or "").split("|") if x]
        print(f"[{i}/{len(sample)}] {pep['peptide']} ids={ids}")
        if not ids:
            unresolved.append(
                {
                    "peptide": pep["peptide"],
                    "epitope_id": "",
                    "stage_stratum": pep.get("stage_stratum"),
                    "status": "fail",
                    "error": "no epitope_ids on pool row",
                    "n_rows_retrieved": 0,
                }
            )
            continue
        for eid in ids:
            n_ids += 1
            result = retrieve_epitope(
                eid,
                raw_dir,
                retrieval_time,
                page_size=args.page_size,
                checkpoint_path=out_dir / "checkpoints" / f"structure_{eid}.json",
            )
            assays = result.pop("assays")
            result.pop("pages", None)
            result["peptide"] = pep["peptide"]
            result["stage_stratum"] = pep.get("stage_stratum")
            logs.append(result)
            if result["status"] in ("fail", "incomplete"):
                unresolved.append(
                    {
                        "peptide": pep["peptide"],
                        "epitope_id": eid,
                        "stage_stratum": pep.get("stage_stratum"),
                        "status": result["status"],
                        "error": result.get("error"),
                        "http_status": result.get("http_status"),
                        "n_rows_retrieved": result.get("n_rows_retrieved"),
                        "reported_total": result.get("reported_total"),
                    }
                )
                # Partial rows stay in raw JSON / completeness log, not the observation table.
                continue
            for a in assays:
                observations.append(annotate_obs(a, pep, eid, retrieval_time))
            time.sleep(0.12)

    write_csv(out_dir / "completeness_log.csv", logs)
    write_csv(out_dir / "observations.csv", observations)
    write_csv(out_dir / "unresolved_queries.csv", unresolved)
    summary = summarize(sample, logs, observations, time.time() - t0)
    summary["retrieval_time_utc"] = retrieval_time
    summary["n_epitope_id_requests"] = n_ids
    summary["n_unresolved_queries"] = len(unresolved)
    summary["n_peptides_with_ok_or_empty"] = sum(1 for r in logs if r["status"] in ("ok", "empty"))
    summary["n_peptides_queried_unresolved"] = sum(
        1 for pep in sample if pep["peptide"] in {u["peptide"] for u in unresolved}
    )
    (out_dir / f"{args.stage}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
