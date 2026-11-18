"""Task-feasibility inventory for predeclared endpoint intersections.

Does not choose a V2 policy or a split. Counts what remains under candidate
study- and sequence-separated evaluations.
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
    mixed_complete_context_groups,
    nonempty,
    outcome_pattern,
)
from iedb_staged_retrieve import peptide_stratum  # noqa: E402
from summarize_task_candidates import eligible  # noqa: E402

FULL = ROOT / "audit" / "retrieval" / "full"
POOL = ROOT / "audit" / "epitope_pool_8_11_preconflict.csv"

CANDIDATES = [
    ("all_observations", "All recovered observations", {}),
    ("human_host", "Human host", {"human": True}),
    ("human_class_I", "Human host + annotated class I", {"human": True, "class_i": True}),
    (
        "human_class_I_functional",
        "Human host + annotated class I + functional endpoint",
        {"human": True, "class_i": True, "endpoint_family": "functional_response"},
    ),
    (
        "human_class_I_functional_ifng_elispot",
        "Human + class I + IFNg ELISPOT",
        {
            "human": True,
            "class_i": True,
            "endpoint_family": "functional_response",
            "assay_names": "IFNg release|ELISPOT",
        },
    ),
    (
        "human_class_I_functional_ifng_ics",
        "Human + class I + IFNg ICS",
        {
            "human": True,
            "class_i": True,
            "endpoint_family": "functional_response",
            "assay_names": "IFNg release|ICS",
        },
    ),
    (
        "human_class_I_functional_chromium",
        "Human + class I + chromium cytotoxicity",
        {
            "human": True,
            "class_i": True,
            "endpoint_family": "functional_response",
            "assay_names": "cytotoxicity|51 chromium",
        },
    ),
    (
        "human_class_I_recognition_binding",
        "Human + class I + recognition/binding (contrast, not adopted)",
        {"human": True, "class_i": True, "endpoint_family": "recognition_binding"},
    ),
]


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


def union_find(items: list[str], edges: list[tuple[str, str]]) -> dict[str, int]:
    parent = {x: x for x in items}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for a, b in edges:
        if a in parent and b in parent:
            union(a, b)
    size: Counter = Counter()
    for x in items:
        size[find(x)] += 1
    return size


def hamming_components(seqs: list[str], max_d: int = 1) -> list[int]:
    seqs = list(dict.fromkeys(seqs))
    if len(seqs) >= 1500:
        from similarity_audit import union_find_components

        return union_find_components(seqs, max_d)
    n = len(seqs)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    by_len = defaultdict(list)
    for i, s in enumerate(seqs):
        by_len[len(s)].append(i)
    for idxs in by_len.values():
        for a in range(len(idxs)):
            sa = seqs[idxs[a]]
            for b in range(a + 1, len(idxs)):
                sb = seqs[idxs[b]]
                if sum(x != y for x, y in zip(sa, sb)) <= max_d:
                    union(idxs[a], idxs[b])
    size = Counter(find(i) for i in range(n))
    return sorted(size.values(), reverse=True)


def study_eval_stats(peps: list[str], by_pep_obs: dict[str, list[dict]]) -> dict:
    pep_to_refs = {}
    ref_to_peps = defaultdict(set)
    for pep in peps:
        refs = {r.get("reference_id") for r in by_pep_obs[pep] if nonempty(r.get("reference_id"))}
        pep_to_refs[pep] = refs
        for ref in refs:
            ref_to_peps[ref].add(pep)
    n_multi_study = sum(1 for refs in pep_to_refs.values() if len(refs) > 1)
    n_single_study = sum(1 for refs in pep_to_refs.values() if len(refs) == 1)
    n_no_study = sum(1 for refs in pep_to_refs.values() if not refs)
    exclusive = [pep for pep, refs in pep_to_refs.items() if len(refs) == 1]
    studies_with_exclusive = defaultdict(list)
    for pep in exclusive:
        ref = next(iter(pep_to_refs[pep]))
        studies_with_exclusive[ref].append(pep)
    # Holding out a study: peptides exclusive to it can be test-only; multi-study peptides leak if kept.
    holdout_sizes = sorted((len(v) for v in studies_with_exclusive.values()), reverse=True)
    edges = []
    items = list(peps)
    for ref, group in ref_to_peps.items():
        g = list(group)
        for i in range(len(g) - 1):
            edges.append((g[i], g[i + 1]))
    # Study graph: nodes = peptides, edge whenever two peptides share a reference_id.
    # Peptides without a reference are their own singleton components.
    comp = union_find(items, edges) if items else Counter()
    comp_sizes = sorted(comp.values(), reverse=True)
    return {
        "n_peptides": len(peps),
        "n_references": len(ref_to_peps),
        "n_peptides_single_study": n_single_study,
        "n_peptides_multi_study": n_multi_study,
        "n_peptides_missing_reference": n_no_study,
        "n_studies_with_any_exclusive_peptide": len(studies_with_exclusive),
        "max_exclusive_peptides_in_one_study": holdout_sizes[0] if holdout_sizes else 0,
        "median_exclusive_peptides_in_one_study": (
            holdout_sizes[len(holdout_sizes) // 2] if holdout_sizes else 0
        ),
        "n_peptides_if_purge_multi_study": n_single_study,
        "n_study_connected_components": len(comp_sizes),
        "n_peptides_in_multi_study_components": sum(s for s in comp_sizes if s > 1),
        "largest_study_component_n_peptides": comp_sizes[0] if comp_sizes else 0,
        "n_peptides_if_collapse_to_component_representatives": len(comp_sizes),
        "note": (
            "A peptide in several studies cannot be assigned to whichever study appears first. "
            "Hold-out of exclusive-study peptides is one candidate; purging multi-study peptides "
            "is another; grouping connected studies is a third. None is selected here. "
            "Study graph: peptides are nodes, linked when they share a reference_id; components "
            "are its connected components. n_peptides_in_multi_study_components counts peptides "
            "in components of size >1 (i.e. sharing at least one study with another peptide). "
            "n_peptides_if_collapse_to_component_representatives equals the component count by "
            "construction (one representative per component); it is an illustrative lower bound "
            "on retained peptides, not an adopted split."
        ),
    }


def sequence_eval_stats(peps: list[str]) -> dict:
    if not peps:
        return {
            "n_peptides": 0,
            "hamming_le1_n_components": 0,
            "hamming_le1_largest_component": 0,
            "hamming_le1_n_peptides_in_components_gt1": 0,
            "n_if_drop_non_singletons": 0,
        }
    sizes = hamming_components(peps, max_d=1)
    n_in_gt1 = sum(s for s in sizes if s > 1)
    n_singletons = sum(1 for s in sizes if s == 1)
    return {
        "n_peptides": len(peps),
        "hamming_le1_n_components": len(sizes),
        "hamming_le1_largest_component": sizes[0] if sizes else 0,
        "hamming_le1_n_peptides_in_components_gt1": n_in_gt1,
        "n_if_drop_non_singletons": n_singletons,
        "note": (
            "Hamming ≤1, same length only. Substring grouping is not applied here. "
            "No cutoff is adopted. Counts describe remaining peptides if non-singleton "
            "components were dropped; that is not a chosen split."
        ),
    }


def main() -> None:
    pool = read_csv(POOL)
    pool_by_pep = {r["peptide"]: r for r in pool}
    obs = read_csv(FULL / "observations.csv")
    logs = read_csv(FULL / "completeness_log.csv")
    latest = {}
    for r in logs:
        latest[str(r.get("epitope_id") or "")] = r
    n_ok = sum(1 for r in latest.values() if r.get("status") in ("ok", "empty"))
    retrieved_peps = {r["export_peptide"] for r in obs}

    by_pep = defaultdict(list)
    for r in obs:
        by_pep[r["export_peptide"]].append(r)

    inventory = []
    for spec_id, label, spec in CANDIDATES:
        pep_rows = []
        for pep, all_rows in by_pep.items():
            sub = [r for r in all_rows if eligible(r, spec)]
            pool_row = pool_by_pep.get(pep, {})
            pat = outcome_pattern(sub)
            pep_rows.append(
                {
                    "subset_id": spec_id,
                    "peptide": pep,
                    "peptide_stratum": peptide_stratum(pool_row) if pool_row else "",
                    "v1_conflict": pool_row.get("v1_conflict", ""),
                    "outcome_pattern": pat,
                    "n_eligible_assays": len(sub),
                    "n_positive": sum(1 for r in sub if r.get("outcome_family") == "positive"),
                    "n_negative": sum(1 for r in sub if r.get("outcome_family") == "negative"),
                    "n_references": len({r.get("reference_id") for r in sub if nonempty(r.get("reference_id"))}),
                    "n_mixed_complete_context_groups": (
                        len(mixed_complete_context_groups(sub)) if pat == "mixed" else 0
                    ),
                }
            )
        patterns = Counter(r["outcome_pattern"] for r in pep_rows)
        labeled = [r["peptide"] for r in pep_rows if r["outcome_pattern"] in ("positive_only", "negative_only")]
        mixed = [r["peptide"] for r in pep_rows if r["outcome_pattern"] == "mixed"]
        eligible_peps = [r["peptide"] for r in pep_rows if r["n_eligible_assays"] > 0]
        labeled_obs = {
            pep: [x for x in by_pep[pep] if eligible(x, spec)]
            for pep in labeled + mixed
        }
        # Study/sequence stats on peptides with an interpretable exclusive outcome (pos-only/neg-only),
        # and separately on all peptides with any eligible assay.
        study_labeled = study_eval_stats(labeled, labeled_obs)
        seq_labeled = sequence_eval_stats(labeled)
        study_eligible = study_eval_stats(
            eligible_peps, {pep: [x for x in by_pep[pep] if eligible(x, spec)] for pep in eligible_peps}
        )
        seq_eligible = sequence_eval_stats(eligible_peps)
        stratum_counts = Counter(r["peptide_stratum"] for r in pep_rows if r["n_eligible_assays"] > 0)
        stratum_patterns = defaultdict(Counter)
        for r in pep_rows:
            stratum_patterns[r["peptide_stratum"]][r["outcome_pattern"]] += 1
        inventory.append(
            {
                "subset_id": spec_id,
                "subset_label": label,
                "illustrative_not_v2_policy": 1,
                "n_peptides_with_eligible_assays": len(eligible_peps),
                "n_positive_only": patterns.get("positive_only", 0),
                "n_negative_only": patterns.get("negative_only", 0),
                "n_mixed": patterns.get("mixed", 0),
                "n_no_eligible_evidence_among_retrieved": patterns.get("no_eligible_evidence", 0),
                "n_not_yet_retrieved": len(pool) - len(retrieved_peps),
                "n_eligible_uninterpretable_outcome": patterns.get(
                    "eligible_uninterpretable_outcome", 0
                ),
                "n_assays": sum(r["n_eligible_assays"] for r in pep_rows),
                "n_mixed_peptides_with_mixed_complete_context_key": sum(
                    1 for r in pep_rows if r["n_mixed_complete_context_groups"] > 0
                ),
                "n_mixed_complete_context_groups": sum(r["n_mixed_complete_context_groups"] for r in pep_rows),
                "stratum_peptides_with_eligible": dict(stratum_counts),
                "study_eval_positive_or_negative_only": study_labeled,
                "sequence_eval_positive_or_negative_only": seq_labeled,
                "study_eval_any_eligible_assay": study_eligible,
                "sequence_eval_any_eligible_assay": seq_eligible,
            }
        )
        write_csv(FULL / f"feasibility_peptides_{spec_id}.csv", pep_rows)

    merge = json.loads((FULL / "merge_validation.json").read_text()) if (FULL / "merge_validation.json").exists() else {}
    out = {
        "v2_labels_defined": False,
        "primary_task_chosen": False,
        "collection": "historical_8_11_preconflict_54448",
        "not_a_census_of_current_iedb": True,
        "not_performance_on_untested_peptides": True,
        "n_pool_peptides": len(pool),
        "n_peptides_with_any_observation_row": len(retrieved_peps),
        "n_epitope_ids_ok_or_empty": n_ok,
        "retrieval_complete_for_defined_collection": bool(merge.get("retrieval_described_as_complete")),
        "candidates": inventory,
        "note": (
            "These intersections are predeclared for feasibility. Do not pick the definition "
            "that would later give the highest model score. Labels should be derived from "
            "eligible observations for a named task, not from a global conflicting-sequence flag."
        ),
    }
    (FULL / "task_feasibility.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    write_csv(
        FULL / "task_feasibility_summary.csv",
        [
            {
                "subset_id": c["subset_id"],
                "subset_label": c["subset_label"],
                "n_peptides_with_eligible_assays": c["n_peptides_with_eligible_assays"],
                "n_positive_only": c["n_positive_only"],
                "n_negative_only": c["n_negative_only"],
                "n_mixed": c["n_mixed"],
                "n_no_eligible_evidence_among_retrieved": c["n_no_eligible_evidence_among_retrieved"],
                "n_not_yet_retrieved": c["n_not_yet_retrieved"],
                "n_eligible_uninterpretable_outcome": c["n_eligible_uninterpretable_outcome"],
                "n_assays": c["n_assays"],
                "n_mixed_peptides_with_mixed_complete_context_key": c[
                    "n_mixed_peptides_with_mixed_complete_context_key"
                ],
                "n_mixed_complete_context_groups": c["n_mixed_complete_context_groups"],
                "n_single_study_posneg_only": c["study_eval_positive_or_negative_only"][
                    "n_peptides_single_study"
                ],
                "n_multi_study_posneg_only": c["study_eval_positive_or_negative_only"][
                    "n_peptides_multi_study"
                ],
                "n_study_connected_components_posneg_only": c["study_eval_positive_or_negative_only"][
                    "n_study_connected_components"
                ],
                "n_peptides_in_multi_study_components_posneg_only": c[
                    "study_eval_positive_or_negative_only"
                ]["n_peptides_in_multi_study_components"],
                "largest_study_component_n_peptides_posneg_only": c[
                    "study_eval_positive_or_negative_only"
                ]["largest_study_component_n_peptides"],
                "n_posneg_only_if_collapse_to_component_representatives": c[
                    "study_eval_positive_or_negative_only"
                ]["n_peptides_if_collapse_to_component_representatives"],
                "largest_hamming1_component_posneg_only": c["sequence_eval_positive_or_negative_only"][
                    "hamming_le1_largest_component"
                ],
                "n_posneg_only_if_drop_hamming1_nonsingletons": c[
                    "sequence_eval_positive_or_negative_only"
                ]["n_if_drop_non_singletons"],
            }
            for c in inventory
        ],
    )
    print(json.dumps({k: out[k] for k in out if k != "candidates"}, indent=2))
    for c in inventory:
        print(
            c["subset_id"],
            "eligible",
            c["n_peptides_with_eligible_assays"],
            "P",
            c["n_positive_only"],
            "N",
            c["n_negative_only"],
            "M",
            c["n_mixed"],
            "none",
            c["n_no_eligible_evidence_among_retrieved"],
            "unint",
            c["n_eligible_uninterpretable_outcome"],
        )


if __name__ == "__main__":
    main()
