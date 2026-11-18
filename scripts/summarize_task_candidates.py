"""Peptide-level eligibility intersections and context-classification audit.

Illustrative subsets only. Does not choose a V2 policy, conflict rule, or label.
"""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from iedb_annotate import (  # noqa: E402
    CONTEXT_KEYS,
    disagreement_class,
    mixed_complete_context_groups,
    nonempty,
    outcome_pattern,
    recorded_value,
)

FAMOUS = "NLVPMVATV"


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


def as_int(v) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


def context_missing(row: dict) -> list[str]:
    return [k for k in CONTEXT_KEYS if recorded_value(row.get(k)) is None]


def eligible(row: dict, spec: dict) -> bool:
    if spec.get("human") and as_int(row.get("host_is_human")) != 1:
        return False
    if spec.get("class_i") and (row.get("mhc_class") or "") != "I":
        return False
    if spec.get("endpoint_family") and row.get("endpoint_family") != spec["endpoint_family"]:
        return False
    if spec.get("assay_names") is not None and (row.get("assay_names") or "") != spec["assay_names"]:
        return False
    if spec.get("evidence") is not None and (row.get("mhc_evidence_category") or "") != spec["evidence"]:
        return False
    return True


def peptide_task_row(
    peptide: str,
    all_rows: list[dict],
    subset_rows: list[dict],
    spec_id: str,
    spec_label: str,
    extra: dict | None = None,
) -> dict:
    refs = {r.get("reference_id") for r in subset_rows if nonempty(r.get("reference_id"))}
    pos_refs = {
        r.get("reference_id")
        for r in subset_rows
        if r.get("outcome_family") == "positive" and nonempty(r.get("reference_id"))
    }
    neg_refs = {
        r.get("reference_id")
        for r in subset_rows
        if r.get("outcome_family") == "negative" and nonempty(r.get("reference_id"))
    }
    v1 = next((r.get("v1_conflict") for r in all_rows), "")
    stratum = next((r.get("stage_stratum") for r in all_rows), "")
    out = {
        "subset_id": spec_id,
        "subset_label": spec_label,
        "illustrative_not_v2_policy": 1,
        "peptide": peptide,
        "stage_stratum": stratum,
        "v1_conflict": v1,
        "current_all_observation_conflict": int(outcome_pattern(all_rows) == "mixed"),
        "task_observation_conflict": int(outcome_pattern(subset_rows) == "mixed"),
        "n_eligible_assays": len(subset_rows),
        "n_positive": sum(1 for r in subset_rows if r.get("outcome_family") == "positive"),
        "n_negative": sum(1 for r in subset_rows if r.get("outcome_family") == "negative"),
        "n_other_or_missing_outcome": sum(
            1 for r in subset_rows if r.get("outcome_family") not in ("positive", "negative")
        ),
        "outcome_pattern": outcome_pattern(subset_rows),
        "n_references": len(refs),
        "n_references_with_positive": len(pos_refs),
        "n_references_with_negative": len(neg_refs),
        "n_references_with_both_outcomes": len(pos_refs & neg_refs),
    }
    if extra:
        out.update(extra)
    return out


def summarize_subset(rows: list[dict], queried_peptides: list[str]) -> dict:
    by = defaultdict(list)
    for r in rows:
        by[r["peptide"]].append(r)
    patterns = Counter(r["outcome_pattern"] for r in rows)
    n_peps_with_obs = len({r["peptide"] for r in rows if r["n_eligible_assays"] > 0})
    return {
        "n_queried_peptides_in_stage": len(queried_peptides),
        "n_peptides_with_eligible_assays": n_peps_with_obs,
        "n_peptides_no_eligible_evidence": patterns.get("no_eligible_evidence", 0)
        + patterns.get("insufficient_unknown", 0),
        "n_peptides_eligible_uninterpretable_outcome": patterns.get(
            "eligible_uninterpretable_outcome", 0
        ),
        "n_peptides_positive_only": patterns.get("positive_only", 0),
        "n_peptides_negative_only": patterns.get("negative_only", 0),
        "n_peptides_mixed": patterns.get("mixed", 0),
        "n_assays": sum(r["n_eligible_assays"] for r in rows),
        "n_references": len({}),  # filled by caller if needed
        "n_v1_conflict_still_task_mixed": sum(
            1 for r in rows if r["v1_conflict"] == "1" and r["outcome_pattern"] == "mixed"
        ),
        "n_v1_conflict_now_positive_only": sum(
            1 for r in rows if r["v1_conflict"] == "1" and r["outcome_pattern"] == "positive_only"
        ),
        "n_v1_conflict_now_negative_only": sum(
            1 for r in rows if r["v1_conflict"] == "1" and r["outcome_pattern"] == "negative_only"
        ),
        "n_v1_conflict_now_no_eligible_evidence": sum(
            1 for r in rows if r["v1_conflict"] == "1" and r["outcome_pattern"] == "no_eligible_evidence"
        ),
        "n_v1_retained_now_task_mixed": sum(
            1 for r in rows if r["v1_conflict"] != "1" and r["outcome_pattern"] == "mixed"
        ),
    }


def distinguishing_fields(pos: list[dict], neg: list[dict]) -> list[dict]:
    out = []
    for k in CONTEXT_KEYS:
        pv = {recorded_value(r.get(k)) for r in pos}
        nv = {recorded_value(r.get(k)) for r in neg}
        out.append(
            {
                "field": k,
                "n_distinct_pos": len(pv),
                "n_distinct_neg": len(nv),
                "n_overlap": len(pv & nv),
                "values_equal": pv == nv,
                "pos_has_missing": None in pv,
                "neg_has_missing": None in nv,
            }
        )
    return out


def context_audit_row(peptide: str, rows: list[dict]) -> dict:
    pos = [r for r in rows if r.get("outcome_family") == "positive"]
    neg = [r for r in rows if r.get("outcome_family") == "negative"]
    pos_complete = [r for r in pos if not context_missing(r)]
    neg_complete = [r for r in neg if not context_missing(r)]
    pos_partial = [r for r in pos if context_missing(r)]
    neg_partial = [r for r in neg if context_missing(r)]
    missing_counts = Counter()
    for r in pos + neg:
        for k in context_missing(r):
            missing_counts[k] += 1
    dist = distinguishing_fields(pos, neg)
    differing = [d["field"] for d in dist if not d["values_equal"]]
    allele_raw_pos = {r.get("mhc_allele_name") or "" for r in pos}
    allele_raw_neg = {r.get("mhc_allele_name") or "" for r in neg}
    spec_pos = {r.get("mhc_specificity") or "" for r in pos}
    spec_neg = {r.get("mhc_specificity") or "" for r in neg}

    both_complete = bool(pos_complete) and bool(neg_complete)
    coarse = disagreement_class(rows)
    if coarse == "within_recorded_context":
        detail = "mixed_after_matching_complete_context"
    elif coarse == "cannot_localize":
        detail = "overlap_only_after_filling_missing"
    elif both_complete:
        detail = "across_known_context_both_families_have_complete_context"
    elif pos_complete and not neg_complete:
        detail = "across_or_unmatched_neg_context_incomplete"
    elif neg_complete and not pos_complete:
        detail = "across_or_unmatched_pos_context_incomplete"
    else:
        detail = "across_or_unmatched_both_families_incomplete"

    return {
        "peptide": peptide,
        "n_assays": len(rows),
        "n_positive": len(pos),
        "n_negative": len(neg),
        "n_pos_complete_context": len(pos_complete),
        "n_neg_complete_context": len(neg_complete),
        "n_pos_partial_context": len(pos_partial),
        "n_neg_partial_context": len(neg_partial),
        "both_outcome_families_have_complete_context": int(both_complete),
        "disagreement_class": coarse,
        "disagreement_detail": detail,
        "fields_that_differ_between_pos_and_neg": "|".join(differing),
        "n_fields_that_differ": len(differing),
        "missing_context_fields": "|".join(f"{k}:{n}" for k, n in missing_counts.most_common()),
        "raw_allele_strings_equal": int(allele_raw_pos == allele_raw_neg),
        "raw_allele_n_pos": len(allele_raw_pos),
        "raw_allele_n_neg": len(allele_raw_neg),
        "raw_allele_overlap": len(allele_raw_pos & allele_raw_neg),
        "specificity_flag_equal": int(spec_pos == spec_neg),
        "v1_conflict": next((r.get("v1_conflict") for r in rows), ""),
        "stage_stratum": next((r.get("stage_stratum") for r in rows), ""),
        "still_mixed_after_matching_recorded_context": int(coarse == "within_recorded_context"),
        "n_mixed_complete_context_groups": len(mixed_complete_context_groups(rows)),
    }


def main(argv: list[str] | None = None) -> None:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--stage", default="stage1")
    args = p.parse_args(argv)
    STAGE = ROOT / "audit" / "retrieval" / args.stage
    sample_path = STAGE / f"{args.stage}_sample.csv"
    sample = read_csv(sample_path)
    obs = read_csv(STAGE / "observations.csv")
    queried = [r["peptide"] for r in sample]
    by_pep: dict[str, list[dict]] = defaultdict(list)
    for r in obs:
        by_pep[r["export_peptide"]].append(r)
        r["peptide"] = r["export_peptide"]

    base_specs = [
        ("all_observations", "All recovered observations", {}),
        ("human_host", "Human host", {"human": True}),
        ("human_class_I", "Human host + annotated class I", {"human": True, "class_i": True}),
        (
            "human_class_I_functional",
            "Human host + annotated class I + functional endpoint",
            {"human": True, "class_i": True, "endpoint_family": "functional_response"},
        ),
    ]

    peptide_rows = []
    for spec_id, label, spec in base_specs:
        for pep in queried:
            all_rows = by_pep.get(pep, [])
            sub = [r for r in all_rows if eligible(r, spec)]
            peptide_rows.append(peptide_task_row(pep, all_rows, sub, spec_id, label))

    # Exact endpoints under human + class I + functional
    parent = [
        r
        for r in obs
        if eligible(r, {"human": True, "class_i": True, "endpoint_family": "functional_response"})
    ]
    assay_names = sorted({r.get("assay_names") or "" for r in parent})
    for name in assay_names:
        spec = {
            "human": True,
            "class_i": True,
            "endpoint_family": "functional_response",
            "assay_names": name,
        }
        sid = "human_class_I_functional_by_assay"
        label = f"Human + class I + functional; assay_names={name or 'MISSING'}"
        peps_in = {r["export_peptide"] for r in parent if (r.get("assay_names") or "") == name}
        for pep in sorted(peps_in):
            all_rows = by_pep.get(pep, [])
            sub = [r for r in all_rows if eligible(r, spec)]
            peptide_rows.append(
                peptide_task_row(
                    pep,
                    all_rows,
                    sub,
                    sid,
                    label,
                    extra={"assay_names": name},
                )
            )

    # Restriction-evidence under human + class I + functional
    evidence_vals = sorted({r.get("mhc_evidence_category") or "missing" for r in parent})
    for ev in evidence_vals:
        spec = {
            "human": True,
            "class_i": True,
            "endpoint_family": "functional_response",
            "evidence": ev,
        }
        sid = "human_class_I_functional_by_evidence"
        label = f"Human + class I + functional; mhc_evidence_category={ev}"
        peps_in = {r["export_peptide"] for r in parent if (r.get("mhc_evidence_category") or "missing") == ev}
        for pep in sorted(peps_in):
            all_rows = by_pep.get(pep, [])
            sub = [r for r in all_rows if eligible(r, spec)]
            peptide_rows.append(
                peptide_task_row(
                    pep,
                    all_rows,
                    sub,
                    sid,
                    label,
                    extra={"mhc_evidence_category": ev},
                )
            )

    write_csv(STAGE / "eligibility_intersections.csv", peptide_rows)

    rollup = []
    for spec_id, label, spec in base_specs:
        sub_rows = [r for r in peptide_rows if r["subset_id"] == spec_id]
        assays = [r for r in obs if eligible(r, spec)]
        assays_wo = [r for r in assays if r["export_peptide"] != FAMOUS]
        stats = summarize_subset(sub_rows, queried)
        stats.update(
            {
                "subset_id": spec_id,
                "subset_label": label,
                "illustrative_not_v2_policy": 1,
                "n_assays": len(assays),
                "n_references": len({r.get("reference_id") for r in assays if nonempty(r.get("reference_id"))}),
                "n_assays_without_NLVPMVATV": len(assays_wo),
                "n_peptides_without_NLVPMVATV": len({r["export_peptide"] for r in assays_wo}),
                "n_references_without_NLVPMVATV": len(
                    {r.get("reference_id") for r in assays_wo if nonempty(r.get("reference_id"))}
                ),
                "frac_assays_from_NLVPMVATV": round(
                    sum(1 for r in assays if r["export_peptide"] == FAMOUS) / max(len(assays), 1), 4
                ),
                "sampling_note": (
                    "Equal n per stratum. Peptide-level summaries reduce domination by "
                    "heavily studied peptides; stratum-specific summaries reveal heterogeneity. "
                    "Neither directly estimates full-pool coverage without accounting for the sampling design."
                ),
            }
        )
        # stratum breakdown
        for tag in sorted({r.get("stage_stratum") or "" for r in sample}):
            tag_peps = {r["peptide"] for r in sample if r.get("stage_stratum") == tag}
            tag_rows = [r for r in sub_rows if r["peptide"] in tag_peps]
            stats[f"stratum__{tag}__n_peptides_with_eligible"] = sum(
                1 for r in tag_rows if r["n_eligible_assays"] > 0
            )
            stats[f"stratum__{tag}__n_mixed"] = sum(1 for r in tag_rows if r["outcome_pattern"] == "mixed")
            stats[f"stratum__{tag}__n_assays"] = sum(r["n_eligible_assays"] for r in tag_rows)
        rollup.append(stats)

    # Exact endpoint rollup
    for name in assay_names:
        sub_rows = [
            r
            for r in peptide_rows
            if r["subset_id"] == "human_class_I_functional_by_assay" and r.get("assay_names") == name
        ]
        assays = [
            r
            for r in parent
            if (r.get("assay_names") or "") == name
        ]
        stats = summarize_subset(sub_rows, queried)
        stats.update(
            {
                "subset_id": "human_class_I_functional_by_assay",
                "subset_label": f"Human + class I + functional; assay_names={name or 'MISSING'}",
                "assay_names": name,
                "illustrative_not_v2_policy": 1,
                "n_assays": len(assays),
                "n_references": len({r.get("reference_id") for r in assays if nonempty(r.get("reference_id"))}),
                "n_assays_without_NLVPMVATV": sum(1 for r in assays if r["export_peptide"] != FAMOUS),
                "frac_assays_from_NLVPMVATV": round(
                    sum(1 for r in assays if r["export_peptide"] == FAMOUS) / max(len(assays), 1), 4
                ),
            }
        )
        rollup.append(stats)

    for ev in evidence_vals:
        sub_rows = [
            r
            for r in peptide_rows
            if r["subset_id"] == "human_class_I_functional_by_evidence"
            and r.get("mhc_evidence_category") == ev
        ]
        assays = [r for r in parent if (r.get("mhc_evidence_category") or "missing") == ev]
        stats = summarize_subset(sub_rows, queried)
        stats.update(
            {
                "subset_id": "human_class_I_functional_by_evidence",
                "subset_label": f"Human + class I + functional; mhc_evidence_category={ev}",
                "mhc_evidence_category": ev,
                "illustrative_not_v2_policy": 1,
                "n_assays": len(assays),
                "n_references": len({r.get("reference_id") for r in assays if nonempty(r.get("reference_id"))}),
                "n_assays_without_NLVPMVATV": sum(1 for r in assays if r["export_peptide"] != FAMOUS),
                "frac_assays_from_NLVPMVATV": round(
                    sum(1 for r in assays if r["export_peptide"] == FAMOUS) / max(len(assays), 1), 4
                ),
            }
        )
        rollup.append(stats)

    write_csv(STAGE / "eligibility_intersection_summary.csv", rollup)

    mixed_peps = [
        pep
        for pep, rs in by_pep.items()
        if any(r.get("outcome_family") == "positive" for r in rs)
        and any(r.get("outcome_family") == "negative" for r in rs)
    ]
    audit_rows = [context_audit_row(pep, by_pep[pep]) for pep in sorted(mixed_peps)]
    write_csv(STAGE / "context_audit.csv", audit_rows)

    detail_counts = Counter(r["disagreement_detail"] for r in audit_rows)
    class_counts = Counter(r["disagreement_class"] for r in audit_rows)
    both_complete = sum(r["both_outcome_families_have_complete_context"] for r in audit_rows)
    still_mixed = sum(r["still_mixed_after_matching_recorded_context"] for r in audit_rows)
    allele_string_diff = sum(1 for r in audit_rows if not r["raw_allele_strings_equal"])

    headline = []
    for spec_id, label, _ in base_specs:
        row = next(r for r in rollup if r["subset_id"] == spec_id)
        headline.append(
            {
                "subset_id": spec_id,
                "subset_label": label,
                "n_peptides_with_eligible_assays": row["n_peptides_with_eligible_assays"],
                "n_peptides_positive_only": row["n_peptides_positive_only"],
                "n_peptides_negative_only": row["n_peptides_negative_only"],
                "n_peptides_mixed": row["n_peptides_mixed"],
                "n_peptides_no_eligible_evidence": row["n_peptides_no_eligible_evidence"],
                "n_peptides_eligible_uninterpretable_outcome": row[
                    "n_peptides_eligible_uninterpretable_outcome"
                ],
                "n_assays": row["n_assays"],
                "n_references": row["n_references"],
                "frac_assays_from_NLVPMVATV": row["frac_assays_from_NLVPMVATV"],
                "n_v1_conflict_now_positive_only": row["n_v1_conflict_now_positive_only"],
                "n_v1_conflict_now_negative_only": row["n_v1_conflict_now_negative_only"],
                "n_v1_conflict_still_task_mixed": row["n_v1_conflict_still_task_mixed"],
                "n_v1_conflict_now_no_eligible_evidence": row[
                    "n_v1_conflict_now_no_eligible_evidence"
                ],
            }
        )

    out = {
        "v2_labels_defined": False,
        "illustrative_not_v2_policy": True,
        "sampling_note": (
            "Equal n per stratum. Peptide-level summaries reduce domination by heavily "
            "studied peptides; stratum-specific summaries reveal heterogeneity. Neither "
            "directly estimates full-pool coverage without accounting for the sampling design."
        ),
        "concentration_not_resolved": True,
        "eligibility_headline": headline,
        "context_audit": {
            "n_mixed_peptides": len(audit_rows),
            "n_peptides_with_mixed_complete_context_key": sum(
                1 for r in audit_rows if int(r["n_mixed_complete_context_groups"]) > 0
            ),
            "n_mixed_complete_context_groups": sum(
                int(r["n_mixed_complete_context_groups"]) for r in audit_rows
            ),
            "disagreement_class": dict(class_counts),
            "disagreement_detail": dict(detail_counts),
            "n_both_families_have_complete_context": both_complete,
            "n_peptides_with_at_least_one_mixed_complete_key": still_mixed,
            "n_raw_allele_strings_differ": allele_string_diff,
            "note": (
                "across_known_context is not automatically 'both families have complete "
                "context'. See disagreement_detail. Exact allele-string differences can "
                "prevent a match without a meaningful experimental difference. Raw allele "
                "strings and specificity flags are stored separately."
            ),
        },
    }
    (STAGE / "eligibility_and_context.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
