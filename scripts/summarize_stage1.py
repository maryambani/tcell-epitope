"""Peptide- and stratum-level summaries for stage-1 retrieval. No V2 labels."""
from __future__ import annotations

import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from iedb_annotate import (  # noqa: E402
    disagreement_class,
    mhc_specificity,
    nonempty,
)

STAGE = ROOT / "audit" / "retrieval" / "stage1"


def read_csv(path: Path) -> list[dict]:
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


def pct(n: int, d: int) -> float:
    return round(100.0 * n / d, 1) if d else 0.0


def median_int(xs: list[int]) -> float:
    return float(statistics.median(xs)) if xs else 0.0


def field_presence(rows: list[dict], field: str) -> dict:
    n = sum(1 for r in rows if nonempty(r.get(field)))
    return {"n": n, "pct_assays": pct(n, len(rows)), "n_assays": len(rows)}


def peptide_any(by_pep: dict[str, list[dict]], field: str) -> dict:
    n = sum(1 for rs in by_pep.values() if any(nonempty(r.get(field)) for r in rs))
    return {"n_peptides": n, "of_peptides": len(by_pep), "pct_peptides": pct(n, len(by_pep))}


def host_bucket(host: str) -> str:
    h = (host or "").lower()
    if "homo sapiens" in h:
        return "human"
    if "mus musculus" in h and "tg" in h:
        return "transgenic_mouse"
    if "mus musculus" in h:
        return "mouse"
    if not h.strip():
        return "missing"
    return "other"


def main() -> None:
    sample = read_csv(STAGE / "stage1_sample.csv")
    logs = read_csv(STAGE / "completeness_log.csv")
    obs = read_csv(STAGE / "observations.csv")
    summary = json.loads((STAGE / "stage1_summary.json").read_text(encoding="utf-8"))
    # Re-apply specificity flags from raw allele strings (no allele normalization).
    for r in obs:
        r["mhc_specificity"] = mhc_specificity(r.get("mhc_allele_name") or "", r.get("mhc_class") or "")
    write_csv(STAGE / "observations.csv", obs)

    by_pep: dict[str, list[dict]] = defaultdict(list)
    for r in obs:
        by_pep[r["export_peptide"]].append(r)

    tcell_ids = [r["tcell_id"] for r in obs]
    seq_mismatch = [r for r in obs if r.get("sequence_matches_export") != "1"]
    structure_mismatch = [
        r for r in obs if str(r.get("queried_epitope_id")) != str(r.get("returned_structure_id"))
    ]

    completeness = []
    for log in logs:
        reported = log.get("reported_total")
        retrieved = int(log["n_rows_retrieved"])
        unique = int(log["n_unique_tcell_id"])
        reported_i = int(reported) if str(reported).strip().isdigit() else None
        ok_count = reported_i is None or reported_i == retrieved
        completeness.append(
            {
                "peptide": log["peptide"],
                "status": log["status"],
                "http_status": log["http_status"],
                "reported_total": reported_i,
                "retrieved": retrieved,
                "unique": unique,
                "n_pages": int(log.get("n_pages") or 0),
                "match_reported": int(ok_count and unique == retrieved),
            }
        )

    famous = "NLVPMVATV"
    obs_wo = [r for r in obs if r["export_peptide"] != famous]
    by_pep_wo = {k: v for k, v in by_pep.items() if k != famous}

    stratum_of = {r["peptide"]: r.get("stage_stratum") for r in sample}
    by_stratum_pep: dict[str, list[str]] = defaultdict(list)
    for pep, tag in stratum_of.items():
        by_stratum_pep[tag].append(pep)

    stratum_rows = []
    for tag, peps in sorted(by_stratum_pep.items()):
        assays = [a for p in peps for a in by_pep.get(p, [])]
        n_per = [len(by_pep.get(p, [])) for p in peps]
        mixed = sum(
            1
            for p in peps
            if any(x["outcome_family"] == "positive" for x in by_pep.get(p, []))
            and any(x["outcome_family"] == "negative" for x in by_pep.get(p, []))
        )
        stratum_rows.append(
            {
                "stage_stratum": tag,
                "n_peptides": len(peps),
                "n_assays": len(assays),
                "assays_min": min(n_per) if n_per else 0,
                "assays_median": median_int(n_per),
                "assays_max": max(n_per) if n_per else 0,
                "n_mixed_outcome_peptides": mixed,
                "n_human_assays": sum(int(r["host_is_human"]) for r in assays),
                "n_class_I": sum(1 for r in assays if r.get("mhc_class") == "I"),
                "n_class_II": sum(1 for r in assays if r.get("mhc_class") == "II"),
                "n_functional": sum(1 for r in assays if r["endpoint_family"] == "functional_response"),
                "n_binding": sum(1 for r in assays if r["endpoint_family"] == "recognition_binding"),
                "n_structural": sum(1 for r in assays if r["endpoint_family"] == "structural"),
                "n_distinct_refs": len({r["reference_id"] for r in assays if nonempty(r.get("reference_id"))}),
            }
        )

    peptide_rows = []
    for pep in [r["peptide"] for r in sample]:
        rs = by_pep.get(pep, [])
        peptide_rows.append(
            {
                "peptide": pep,
                "stage_stratum": stratum_of.get(pep),
                "v1_conflict": next((r.get("v1_conflict") for r in sample if r["peptide"] == pep), ""),
                "v1_query_label": next((r.get("v1_query_label") for r in sample if r["peptide"] == pep), ""),
                "n_assays": len(rs),
                "n_positive": sum(1 for r in rs if r["outcome_family"] == "positive"),
                "n_negative": sum(1 for r in rs if r["outcome_family"] == "negative"),
                "mixed_outcome": int(
                    any(r["outcome_family"] == "positive" for r in rs)
                    and any(r["outcome_family"] == "negative" for r in rs)
                ),
                "n_references": len({r["reference_id"] for r in rs if nonempty(r.get("reference_id"))}),
                "any_pubmed": int(any(nonempty(r.get("pubmed_id")) for r in rs)),
                "any_mhc_class": int(any(nonempty(r.get("mhc_class")) for r in rs)),
                "any_allele_name": int(any(nonempty(r.get("mhc_allele_name")) for r in rs)),
                "any_source_antigen": int(any(nonempty(r.get("parent_source_antigen_name")) for r in rs)),
                "n_human": sum(int(r["host_is_human"]) for r in rs),
                "n_transgenic_mouse": sum(int(r["host_is_transgenic_mouse"]) for r in rs),
                "n_class_I": sum(1 for r in rs if r.get("mhc_class") == "I"),
                "n_class_II": sum(1 for r in rs if r.get("mhc_class") == "II"),
                "n_allele_level": sum(1 for r in rs if r.get("mhc_specificity") == "allele"),
                "n_serotype": sum(1 for r in rs if r.get("mhc_specificity") == "serotype"),
                "n_functional": sum(1 for r in rs if r["endpoint_family"] == "functional_response"),
                "n_binding": sum(1 for r in rs if r["endpoint_family"] == "recognition_binding"),
                "n_structural": sum(1 for r in rs if r["endpoint_family"] == "structural"),
            }
        )

    pep_to_refs: dict[str, set] = defaultdict(set)
    ref_to_peps: dict[str, set] = defaultdict(set)
    for r in obs:
        rid = r.get("reference_id")
        if nonempty(rid):
            pep_to_refs[r["export_peptide"]].add(rid)
            ref_to_peps[rid].add(r["export_peptide"])

    n_refs_single_pep = sum(1 for peps in ref_to_peps.values() if len(peps) == 1)
    n_peps_multi_ref = sum(1 for refs in pep_to_refs.values() if len(refs) > 1)

    assay_names = Counter(r.get("assay_names") or "MISSING" for r in obs)
    hosts = Counter(host_bucket(r.get("host_organism_name") or "") for r in obs)
    evidence = Counter(r.get("mhc_evidence_category") or "missing" for r in obs)
    modified = sum(int(r.get("mhc_modified_experiment") or 0) for r in obs)
    disagreement_rows = sorted(
        [
            {
                "peptide": p,
                "n_assays": len(rs),
                "disagreement_class": disagreement_class(rs),
                "stage_stratum": stratum_of.get(p),
            }
            for p, rs in by_pep.items()
            if any(x["outcome_family"] == "positive" for x in rs)
            and any(x["outcome_family"] == "negative" for x in rs)
        ],
        key=lambda x: (x["disagreement_class"], -x["n_assays"], x["peptide"]),
    )

    coverage = {
        "all_assays": {
            "n": len(obs),
            "mhc_class": field_presence(obs, "mhc_class"),
            "mhc_allele_name": field_presence(obs, "mhc_allele_name"),
            "pubmed_id": field_presence(obs, "pubmed_id"),
            "parent_source_antigen_name": field_presence(obs, "parent_source_antigen_name"),
            "mhc_annotation_present": {
                "n": sum(int(r["mhc_annotation_present"]) for r in obs),
                "pct_assays": pct(sum(int(r["mhc_annotation_present"]) for r in obs), len(obs)),
            },
        },
        "without_NLVPMVATV": {
            "n": len(obs_wo),
            "n_peptides": len(by_pep_wo),
            "n_references": len({r["reference_id"] for r in obs_wo if nonempty(r.get("reference_id"))}),
            "mhc_class": field_presence(obs_wo, "mhc_class"),
            "pubmed_id": field_presence(obs_wo, "pubmed_id"),
            "parent_source_antigen_name": field_presence(obs_wo, "parent_source_antigen_name"),
        },
        "peptide_level_any_observation": {
            "pubmed_id": peptide_any(by_pep, "pubmed_id"),
            "mhc_class": peptide_any(by_pep, "mhc_class"),
            "mhc_allele_name": peptide_any(by_pep, "mhc_allele_name"),
            "parent_source_antigen_name": peptide_any(by_pep, "parent_source_antigen_name"),
        },
    }

    out = {
        "retrieval_time_utc": summary.get("retrieval_time_utc"),
        "n_peptides_queried": len(sample),
        "n_epitope_ids_queried": len(logs),
        "n_observations": len(obs),
        "n_unique_tcell_id": len(set(tcell_ids)),
        "n_duplicate_tcell_id": len(tcell_ids) - len(set(tcell_ids)),
        "n_sequence_mismatches": len(seq_mismatch),
        "n_structure_id_mismatches": len(structure_mismatch),
        "completeness_status": dict(Counter(c["status"] for c in completeness)),
        "n_reported_total_mismatch": sum(1 for c in completeness if not c["match_reported"]),
        "n_queries_gt_one_page": sum(1 for c in completeness if c["n_pages"] > 1),
        "max_pages": max(c["n_pages"] for c in completeness) if completeness else 0,
        "elapsed_seconds": summary.get("elapsed_seconds"),
        "seconds_per_epitope_id": summary.get("seconds_per_epitope_id"),
        "assays_per_peptide": {
            "min": min(r["n_assays"] for r in peptide_rows),
            "median": median_int([r["n_assays"] for r in peptide_rows]),
            "mean": round(sum(r["n_assays"] for r in peptide_rows) / len(peptide_rows), 2),
            "max": max(r["n_assays"] for r in peptide_rows),
        },
        "top_peptide": famous if famous in by_pep else summary.get("top_peptide"),
        "top_peptide_n": len(by_pep.get(famous, [])),
        "frac_obs_from_top_peptide": round(len(by_pep.get(famous, [])) / max(len(obs), 1), 4),
        "do_not_extrapolate_runtime_from_this_batch": True,
        "reason_not_to_extrapolate": (
            "Median assays/peptide is 2, but the remaining 54k pool still contains other "
            "high-volume peptides; wall-clock is not 54448 times seconds_per_epitope_id."
        ),
        "coverage": coverage,
        "mhc_specificity": dict(Counter(r["mhc_specificity"] for r in obs)),
        "mhc_specificity_notes": {
            "n_star_allele": sum(1 for r in obs if "*" in (r.get("mhc_allele_name") or "")),
            "n_HLA_A0201": sum(1 for r in obs if (r.get("mhc_allele_name") or "") == "HLA-A*02:01"),
            "n_HLA_A2_serotype": sum(1 for r in obs if (r.get("mhc_allele_name") or "") == "HLA-A2"),
            "n_human_generic_allele_name": sum(
                1 for r in obs if (r.get("mhc_allele_name") or "").lower() == "human"
            ),
            "n_H2_haplotype_like": sum(
                1
                for r in obs
                if (r.get("mhc_allele_name") or "").startswith("H2-")
                and r.get("mhc_specificity") != "allele"
            ),
            "do_not_treat_as_usable_allele_rate": True,
        },
        "mhc_class": dict(Counter((r.get("mhc_class") or "MISSING") for r in obs)),
        "mhc_modified_experiment_n": modified,
        "mhc_evidence_category": dict(evidence),
        "endpoint_family": dict(Counter(r["endpoint_family"] for r in obs)),
        "top_assay_names": assay_names.most_common(15),
        "host_bucket": dict(hosts),
        "n_class_II_peptides": len({r["export_peptide"] for r in obs if r.get("mhc_class") == "II"}),
        "n_mixed_outcome_peptides": sum(r["mixed_outcome"] for r in peptide_rows),
        "disagreement_classes": dict(Counter(r["disagreement_class"] for r in disagreement_rows)),
        "n_disagreement_peptides": len(disagreement_rows),
        "mixed_vs_v1_conflict": {
            "v1_conflict_and_mixed": sum(
                1 for r in peptide_rows if r["v1_conflict"] == "1" and r["mixed_outcome"]
            ),
            "v1_conflict_not_mixed": sum(
                1 for r in peptide_rows if r["v1_conflict"] == "1" and not r["mixed_outcome"]
            ),
            "v1_retained_and_mixed": sum(
                1 for r in peptide_rows if r["v1_conflict"] != "1" and r["mixed_outcome"]
            ),
        },
        "references": {
            "n_distinct_reference_id": len(ref_to_peps),
            "n_refs_with_only_one_sampled_peptide": n_refs_single_pep,
            "n_sampled_peptides_in_multiple_refs": n_peps_multi_ref,
            "max_sampled_peptides_per_ref": max((len(v) for v in ref_to_peps.values()), default=0),
            "note": (
                "Counts are among the 240 sampled peptides only. They do not establish "
                "study-holdout feasibility for the full 54,448 pool."
            ),
        },
        "v2_labels_defined": False,
        "strata": stratum_rows,
    }

    write_csv(STAGE / "peptide_summary.csv", peptide_rows)
    write_csv(STAGE / "stratum_summary.csv", stratum_rows)
    write_csv(STAGE / "disagreement_peptides.csv", disagreement_rows)
    (STAGE / "stage1_diagnostics.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
