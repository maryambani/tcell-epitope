"""Retrieve T-cell observations for the full 8–11 pre-conflict pool.

Skips epitope IDs already ok/empty in prior batches. Appends completeness
and observations. fail/incomplete stay unresolved, not empty.

Does not assign V2 labels.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from iedb_staged_retrieve import (  # noqa: E402
    API,
    AUDIT,
    PAGE,
    POOL,
    SELECT,
    annotate_obs,
    peptide_stratum,
    read_csv,
    retrieve_epitope,
    write_csv,
)

OUT = AUDIT / "retrieval" / "full"
CONDITIONS = {
    "api": API,
    "select": SELECT,
    "page_size": PAGE,
    "order": "tcell_id",
    "prefer": "count=exact",
    "sleep_s": 0.12,
}


def load_done_ids(log_path: Path) -> dict[str, dict]:
    if not log_path.exists() or log_path.stat().st_size == 0:
        return {}
    out = {}
    for r in read_csv(log_path):
        eid = str(r.get("epitope_id") or "")
        if eid:
            out[eid] = r
    return out


def append_rows(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if new_file:
            w.writeheader()
        w.writerows(rows)
        f.flush()


def import_prior_batch(batch: str, fieldnames_holder: dict) -> tuple[list[dict], list[dict]]:
    bdir = AUDIT / "retrieval" / batch
    logs = read_csv(bdir / "completeness_log.csv")
    obs = read_csv(bdir / "observations.csv")
    for r in logs:
        r["retrieval_batch"] = batch
        r.setdefault("page_size", "" if batch == "stage1" else r.get("page_size", PAGE))
    for r in obs:
        r["retrieval_batch"] = batch
    if obs:
        for k in obs[0]:
            if k not in fieldnames_holder["obs"]:
                fieldnames_holder["obs"].append(k)
    if logs:
        for k in logs[0]:
            if k not in fieldnames_holder["log"]:
                fieldnames_holder["log"].append(k)
    return logs, obs


def write_progress(payload: dict) -> None:
    (OUT / "progress.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--max-ids", type=int, default=0, help="0 = remaining pool")
    p.add_argument("--sleep", type=float, default=CONDITIONS["sleep_s"])
    p.add_argument("--page-size", type=int, default=PAGE)
    args = p.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    raw_dir = OUT / "raw"
    raw_dir.mkdir(exist_ok=True)
    cond_path = OUT / "retrieval_conditions.json"
    if cond_path.exists():
        prev = json.loads(cond_path.read_text(encoding="utf-8"))
        if prev.get("page_size") != args.page_size or prev.get("api") != API:
            print("WARNING: retrieval conditions differ from", cond_path)
            print(" previous", prev)
            print(" current page_size", args.page_size, "api", API)
    conditions = dict(CONDITIONS)
    conditions["page_size"] = args.page_size
    conditions["sleep_s"] = args.sleep
    cond_path.write_text(json.dumps(conditions, indent=2), encoding="utf-8")

    pool = read_csv(POOL)
    manifest = []
    for r in pool:
        ids = [x for x in (r.get("epitope_ids") or "").split("|") if x]
        if not ids and r.get("primary_epitope_id"):
            ids = [r["primary_epitope_id"]]
        for eid in ids:
            manifest.append(
                {
                    "peptide": r["peptide"],
                    "epitope_id": eid,
                    "peptide_stratum": peptide_stratum(r),
                    "v1_conflict": r.get("v1_conflict"),
                    "v1_retained": r.get("v1_retained"),
                    "v1_query_label": r.get("v1_query_label"),
                    "in_iedb": r.get("in_iedb"),
                    "in_cedar": r.get("in_cedar"),
                }
            )
    write_csv(OUT / "manifest.csv", manifest)
    print("manifest", len(manifest), "peptides", len(pool))

    log_path = OUT / "completeness_log.csv"
    obs_path = OUT / "observations.csv"
    unresolved_path = OUT / "unresolved_queries.csv"
    fields = {"log": [], "obs": []}

    if not log_path.exists() or log_path.stat().st_size == 0:
        imported_logs = []
        imported_obs = []
        for batch in ("stage1", "stage2"):
            logs, obs = import_prior_batch(batch, fields)
            imported_logs.extend(logs)
            imported_obs.extend(obs)
            print("imported", batch, "logs", len(logs), "obs", len(obs))
        if "retrieval_batch" not in fields["log"]:
            fields["log"].append("retrieval_batch")
        if "retrieval_batch" not in fields["obs"]:
            fields["obs"].append("retrieval_batch")
        write_csv(log_path, imported_logs)
        write_csv(obs_path, imported_obs)
        write_csv(unresolved_path, [])

    done = load_done_ids(log_path)
    skip = {eid for eid, r in done.items() if r.get("status") in ("ok", "empty")}
    retry = {eid for eid, r in done.items() if r.get("status") in ("fail", "incomplete")}
    todo = [m for m in manifest if m["epitope_id"] not in skip]
    if args.max_ids:
        todo = todo[: args.max_ids]
    print("skip ok/empty", len(skip), "retry", len(retry), "todo this run", len(todo))

    # Fieldnames from existing files if present
    if obs_path.exists() and obs_path.stat().st_size:
        with obs_path.open(encoding="utf-8") as f:
            obs_fields = next(csv.reader(f))
    else:
        obs_fields = fields["obs"] or ["observation_id"]
    if log_path.exists() and log_path.stat().st_size:
        with log_path.open(encoding="utf-8") as f:
            log_fields = next(csv.reader(f))
    else:
        log_fields = fields["log"]
    for extra in (
        "retrieval_batch",
        "page_size",
        "order",
        "n_duplicate_tcell_id",
        "start_offset",
    ):
        if extra not in log_fields:
            log_fields.append(extra)
    if "retrieval_batch" not in obs_fields:
        obs_fields.append("retrieval_batch")

    unresolved_fields = [
        "peptide",
        "epitope_id",
        "peptide_stratum",
        "status",
        "error",
        "http_status",
        "n_rows_retrieved",
        "reported_total",
        "retrieval_time_utc",
        "retrieval_batch",
    ]

    t0 = time.time()
    batch_name = "full"
    n_ok = n_empty = n_fail = n_inc = n_obs = 0
    for i, item in enumerate(todo, 1):
        eid = item["epitope_id"]
        pep_row = {
            "peptide": item["peptide"],
            "stage_stratum": item["peptide_stratum"],
            "v1_conflict": item["v1_conflict"],
            "v1_query_label": item["v1_query_label"],
        }
        retrieval_time = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        print(f"[{i}/{len(todo)}] {item['peptide']} id={eid}", flush=True)
        prev = done.get(eid)
        result = retrieve_epitope(
            eid,
            raw_dir,
            retrieval_time,
            page_size=args.page_size,
            sleep_s=args.sleep,
        )
        assays = result.pop("assays")
        result.pop("pages", None)
        result["peptide"] = item["peptide"]
        result["stage_stratum"] = item["peptide_stratum"]
        result["retrieval_batch"] = batch_name
        result["order"] = "tcell_id"
        if prev and str(prev.get("reported_total") or "") and str(result.get("reported_total") or ""):
            if str(prev["reported_total"]) != str(result["reported_total"]):
                result["error"] = (
                    (result.get("error") or "")
                    + f" reported_total_changed prev={prev['reported_total']} now={result['reported_total']}"
                ).strip()
                print("  reported_total changed", prev["reported_total"], "->", result["reported_total"])
        append_rows(log_path, [result], log_fields)
        status = result["status"]
        if status == "ok":
            n_ok += 1
            rows = [annotate_obs(a, pep_row, eid, retrieval_time) for a in assays]
            for r in rows:
                r["retrieval_batch"] = batch_name
                r["stage_stratum"] = item["peptide_stratum"]
            for r in rows:
                for k in r:
                    if k not in obs_fields:
                        obs_fields.append(k)
            append_rows(obs_path, rows, obs_fields)
            n_obs += len(rows)
        elif status == "empty":
            n_empty += 1
        else:
            if status == "fail":
                n_fail += 1
            else:
                n_inc += 1
            append_rows(
                unresolved_path,
                [
                    {
                        "peptide": item["peptide"],
                        "epitope_id": eid,
                        "peptide_stratum": item["peptide_stratum"],
                        "status": status,
                        "error": result.get("error"),
                        "http_status": result.get("http_status"),
                        "n_rows_retrieved": result.get("n_rows_retrieved"),
                        "reported_total": result.get("reported_total"),
                        "retrieval_time_utc": retrieval_time,
                        "retrieval_batch": batch_name,
                    }
                ],
                unresolved_fields,
            )
        if i % 10 == 0 or i == len(todo):
            write_progress(
                {
                    "retrieval_time_utc": retrieval_time,
                    "n_todo_this_run": len(todo),
                    "n_finished_this_run": i,
                    "n_ok": n_ok,
                    "n_empty": n_empty,
                    "n_fail": n_fail,
                    "n_incomplete": n_inc,
                    "n_observations_appended": n_obs,
                    "elapsed_seconds": round(time.time() - t0, 1),
                    "last_peptide": item["peptide"],
                    "last_epitope_id": eid,
                    "v2_labels_defined": False,
                    "collection": "historical_8_11_preconflict_54448",
                    "not_a_census_of_current_iedb": True,
                }
            )
        time.sleep(args.sleep)

    write_progress(
        {
            "done_run": True,
            "n_todo_this_run": len(todo),
            "n_ok": n_ok,
            "n_empty": n_empty,
            "n_fail": n_fail,
            "n_incomplete": n_inc,
            "n_observations_appended": n_obs,
            "elapsed_seconds": round(time.time() - t0, 1),
            "v2_labels_defined": False,
        }
    )
    print("finished run ok", n_ok, "empty", n_empty, "fail", n_fail, "incomplete", n_inc)


if __name__ == "__main__":
    main()
