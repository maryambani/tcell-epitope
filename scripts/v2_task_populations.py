"""Exact population counts for two predeclared V2 task candidates.

Counting exercise only. Does not train a model, does not write a label column,
does not choose a task or a split. Supports a human decision about whether the
IFNg ELISPOT intersection is usable once reference 1001865 (a single
all-negative vaccinia Dryvax study) is accounted for.

Reuses, without reimplementing:
  summarize_task_candidates.eligible  (host / class / endpoint / assay filter)
  iedb_annotate.outcome_pattern, nonempty
  task_feasibility.hamming_components, read_csv, CANDIDATES
  per_study_table.MISSING_REF convention for rows without a reference_id

The only new graph code is a degree-zero test used to identify Hamming<=1
singleton peptides (membership, not sizes); it is cross-checked against
hamming_components() and the check is reported in the JSON.
"""
from __future__ import annotations

import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from iedb_annotate import nonempty, outcome_pattern  # noqa: E402
from per_study_table import MISSING_REF, ref_of  # noqa: E402
from summarize_task_candidates import eligible  # noqa: E402
from task_feasibility import CANDIDATES, hamming_components, read_csv  # noqa: E402

FULL = ROOT / "audit" / "retrieval" / "full"
OUT = FULL / "v2_task_populations.json"

SUBSETS = ("human_class_I_functional_ifng_elispot", "human_class_I_functional")
DOMINANT_REF = "1001865"
LABELED = ("positive_only", "negative_only")
TEST_FRACTION_TARGET = 0.20
N_ORDERINGS = 200
SEED = 20261119

KEEP_FIELDS = (
    "export_peptide",
    "outcome_family",
    "assay_names",
    "endpoint_family",
    "host_is_human",
    "mhc_class",
    "mhc_evidence_category",
    "reference_id",
)


def frac(num: int, den: int) -> float | None:
    return round(num / den, 6) if den else None


def pos_fraction(peps, pattern: dict[str, str]) -> float | None:
    lab = [p for p in peps if pattern[p] in LABELED]
    if not lab:
        return None
    return frac(sum(1 for p in lab if pattern[p] == "positive_only"), len(lab))


def singleton_flags(peps: list[str], max_d: int = 1) -> dict[str, bool]:
    """True when a peptide has no same-length neighbour within Hamming<=max_d.

    Degree zero in the same graph hamming_components() builds, so these are
    exactly its singleton components.
    """
    flags = {p: True for p in peps}
    by_len: dict[int, list[str]] = defaultdict(list)
    for p in peps:
        by_len[len(p)].append(p)
    for seqs in by_len.values():
        if len(seqs) < 2:
            continue
        arr = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(
            len(seqs), len(seqs[0])
        )
        for i in range(len(seqs)):
            d = (arr != arr[i]).sum(axis=1)
            d[i] = max_d + 1
            if int((d <= max_d).sum()):
                flags[seqs[i]] = False
    return flags


def study_holdout_sweep(
    exclusive_by_study: dict[str, list[str]],
    pattern: dict[str, str],
    n_labeled_total: int,
    rng: random.Random,
) -> dict:
    """Random study orderings held out until ~TEST_FRACTION_TARGET is reached.

    Only peptides exclusive to a held-out study can enter the test set; a
    labeled peptide shared with a kept study would leak and is not counted.
    """
    studies = sorted(exclusive_by_study)
    target = TEST_FRACTION_TARGET * n_labeled_total
    pos_by_study = {
        s: sum(1 for p in peps if pattern[p] == "positive_only")
        for s, peps in exclusive_by_study.items()
    }
    n_by_study = {s: len(peps) for s, peps in exclusive_by_study.items()}
    draws = []
    for _ in range(N_ORDERINGS):
        order = studies[:]
        rng.shuffle(order)
        n = pos = 0
        used = 0
        for s in order:
            n += n_by_study[s]
            pos += pos_by_study[s]
            used += 1
            if n >= target:
                break
        draws.append(
            {
                "n_studies_held_out": used,
                "n_test_peptides": n,
                "test_fraction_of_labeled": frac(n, n_labeled_total),
                "test_positive_fraction": frac(pos, n),
            }
        )
    pf = sorted(d["test_positive_fraction"] for d in draws)
    ns = sorted(d["n_studies_held_out"] for d in draws)
    npep = sorted(d["n_test_peptides"] for d in draws)
    mid = len(pf) // 2
    representative = min(draws, key=lambda d: abs(d["test_positive_fraction"] - pf[mid]))
    return {
        "n_random_orderings": N_ORDERINGS,
        "seed": SEED,
        "target_test_fraction": TEST_FRACTION_TARGET,
        "n_studies_available_with_exclusive_labeled_peptides": len(studies),
        "test_positive_fraction_median": pf[mid],
        "test_positive_fraction_min": pf[0],
        "test_positive_fraction_max": pf[-1],
        "n_studies_held_out_median": ns[mid],
        "n_studies_held_out_min": ns[0],
        "n_studies_held_out_max": ns[-1],
        "n_test_peptides_median": npep[mid],
        "n_test_peptides_min": npep[0],
        "n_test_peptides_max": npep[-1],
        "representative_draw": representative,
        "note": (
            "Test set = union of the held-out studies' exclusive labeled peptides. "
            "Labeled peptides shared between a held-out and a kept study are excluded "
            "(they would leak), so the achieved test fraction can fall short or overshoot "
            "the 20% target in a single step when a large study is drawn early. "
            "Sequence separation is not imposed here."
        ),
    }


def analyse(spec_id: str, label: str, spec: dict, by_pep: dict[str, list[dict]]) -> dict:
    # ---- 1. baseline population -------------------------------------------
    elig: dict[str, list[dict]] = {}
    for pep, rows in by_pep.items():
        sub = [r for r in rows if eligible(r, spec)]
        if sub:
            elig[pep] = sub
    pattern = {pep: outcome_pattern(rows) for pep, rows in elig.items()}
    pats = Counter(pattern.values())
    labeled = [p for p in elig if pattern[p] in LABELED]
    n_lab = len(labeled)
    pep_refs = {pep: {ref_of(r) for r in rows} for pep, rows in elig.items()}
    all_refs = set().union(*pep_refs.values()) if pep_refs else set()
    baseline = {
        "n_eligible_peptides": len(elig),
        "n_positive_only": pats.get("positive_only", 0),
        "n_negative_only": pats.get("negative_only", 0),
        "n_mixed": pats.get("mixed", 0),
        "n_uninterpretable": pats.get("eligible_uninterpretable_outcome", 0),
        "n_labeled": n_lab,
        "labeled_positive_fraction": frac(pats.get("positive_only", 0), n_lab),
        "n_eligible_observations": sum(len(v) for v in elig.values()),
        "n_studies": len(all_refs - {MISSING_REF}),
        "n_eligible_rows_missing_reference_id": sum(
            1 for rows in elig.values() for r in rows if not nonempty(r.get("reference_id"))
        ),
    }

    # ---- 2. reference 1001865 contribution --------------------------------
    touched = [p for p in elig if DOMINANT_REF in pep_refs[p]]
    exclusive = [p for p in touched if pep_refs[p] == {DOMINANT_REF}]
    within = Counter(
        outcome_pattern([r for r in elig[p] if ref_of(r) == DOMINANT_REF]) for p in touched
    )
    within_excl = Counter(
        outcome_pattern([r for r in elig[p] if ref_of(r) == DOMINANT_REF]) for p in exclusive
    )
    global_excl = Counter(pattern[p] for p in exclusive)
    global_touched = Counter(pattern[p] for p in touched)
    dominant = {
        "reference_id": DOMINANT_REF,
        "n_peptides_touched": len(touched),
        "n_peptides_exclusive": len(exclusive),
        "n_peptides_shared_with_other_studies": len(touched) - len(exclusive),
        "n_eligible_observations_from_this_study": sum(
            1 for p in touched for r in elig[p] if ref_of(r) == DOMINANT_REF
        ),
        "within_study_pattern_all_touched": dict(within),
        "within_study_pattern_exclusive_only": dict(within_excl),
        "entirely_negative_within_study": within.get("positive_only", 0) == 0
        and within.get("mixed", 0) == 0,
        "global_pattern_of_exclusive_peptides": dict(global_excl),
        "global_pattern_of_all_touched_peptides": dict(global_touched),
        "exclusive_global_matches_within_study": dict(global_excl) == dict(within_excl),
        "n_labeled_touched": sum(1 for p in touched if pattern[p] in LABELED),
        "n_labeled_exclusive": sum(1 for p in exclusive if pattern[p] in LABELED),
        "share_of_baseline_labeled_touched": frac(
            sum(1 for p in touched if pattern[p] in LABELED), n_lab
        ),
        "share_of_baseline_labeled_exclusive": frac(
            sum(1 for p in exclusive if pattern[p] in LABELED), n_lab
        ),
    }

    # ---- 3. removal variants ----------------------------------------------
    excl_set = set(exclusive)
    touched_set = set(touched)

    # Variant A: drop exclusives only; recompute pattern for shared peptides
    # from the remaining studies' observations.
    a_obs = {
        p: ([r for r in elig[p] if ref_of(r) != DOMINANT_REF] if p in touched_set else elig[p])
        for p in elig
        if p not in excl_set
    }
    a_obs = {p: rows for p, rows in a_obs.items() if rows}
    a_pattern = {p: outcome_pattern(rows) for p, rows in a_obs.items()}
    a_pats = Counter(a_pattern.values())
    a_labeled = [p for p in a_obs if a_pattern[p] in LABELED]
    a_refs = {ref_of(r) for rows in a_obs.values() for r in rows}
    a_refs_labeled = {ref_of(r) for p in a_labeled for r in a_obs[p]}
    changes = Counter(
        f"{pattern[p]}->{a_pattern[p]}"
        for p in (touched_set - excl_set)
        if p in a_pattern and a_pattern[p] != pattern[p]
    )
    variant_a = {
        "definition": "drop peptides exclusive to 1001865; keep shared peptides with 1001865 rows removed",
        "n_eligible_peptides": len(a_obs),
        "n_labeled": len(a_labeled),
        "n_positive_only": a_pats.get("positive_only", 0),
        "n_negative_only": a_pats.get("negative_only", 0),
        "n_mixed": a_pats.get("mixed", 0),
        "n_uninterpretable": a_pats.get("eligible_uninterpretable_outcome", 0),
        "labeled_positive_fraction": frac(a_pats.get("positive_only", 0), len(a_labeled)),
        "n_studies_remaining_any_eligible": len(a_refs - {MISSING_REF}),
        "n_studies_remaining_among_labeled": len(a_refs_labeled - {MISSING_REF}),
        "n_shared_peptides_kept": len(touched_set - excl_set),
        "n_shared_peptides_relabeled_by_removal": sum(changes.values()),
        "shared_peptide_pattern_change": dict(changes),
        "n_shared_peptides_losing_all_eligible_evidence": len(
            (touched_set - excl_set) - set(a_obs)
        ),
    }

    # Variant B: drop every peptide 1001865 touches at all.
    b_obs = {p: elig[p] for p in elig if p not in touched_set}
    b_labeled = [p for p in b_obs if pattern[p] in LABELED]
    b_refs = {ref_of(r) for rows in b_obs.values() for r in rows}
    b_refs_labeled = {ref_of(r) for p in b_labeled for r in b_obs[p]}
    b_pats = Counter(pattern[p] for p in b_obs)
    variant_b = {
        "definition": "drop every peptide with any eligible 1001865 observation",
        "n_eligible_peptides": len(b_obs),
        "n_labeled": len(b_labeled),
        "n_positive_only": b_pats.get("positive_only", 0),
        "n_negative_only": b_pats.get("negative_only", 0),
        "n_mixed": b_pats.get("mixed", 0),
        "n_uninterpretable": b_pats.get("eligible_uninterpretable_outcome", 0),
        "labeled_positive_fraction": frac(b_pats.get("positive_only", 0), len(b_labeled)),
        "n_studies_remaining_any_eligible": len(b_refs - {MISSING_REF}),
        "n_studies_remaining_among_labeled": len(b_refs_labeled - {MISSING_REF}),
    }

    # ---- 4. concentration after removal (variant A) -----------------------
    a_pep_refs = {p: {ref_of(r) for r in rows} for p, rows in a_obs.items()}
    a_excl_by_study: dict[str, list[str]] = defaultdict(list)
    for p in a_labeled:
        refs = a_pep_refs[p]
        if len(refs) == 1:
            a_excl_by_study[next(iter(refs))].append(p)
    a_any_by_study: dict[str, list[str]] = defaultdict(list)
    for p in a_labeled:
        for ref in a_pep_refs[p]:
            a_any_by_study[ref].append(p)
    n_a_lab = len(a_labeled)
    top = sorted(a_excl_by_study.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:10]
    top_rows = []
    for ref, peps in top:
        npos = sum(1 for p in peps if a_pattern[p] == "positive_only")
        any_peps = a_any_by_study.get(ref, [])
        top_rows.append(
            {
                "reference_id": ref,
                "n_exclusive_labeled_peptides": len(peps),
                "n_positive_only": npos,
                "n_negative_only": len(peps) - npos,
                "positive_fraction": frac(npos, len(peps)),
                "share_of_reduced_labeled_population": frac(len(peps), n_a_lab),
                "n_labeled_peptides_touched_any": len(set(any_peps)),
                "share_touched_any_of_reduced_labeled": frac(len(set(any_peps)), n_a_lab),
            }
        )
    max_share = max((r["share_of_reduced_labeled_population"] for r in top_rows), default=None)
    max_share_any = max((r["share_touched_any_of_reduced_labeled"] for r in top_rows), default=None)
    concentration = {
        "population": "variant_A labeled peptides",
        "n_labeled_reduced_population": n_a_lab,
        "top10_studies_by_exclusive_labeled_peptides": top_rows,
        "largest_single_study_share_exclusive": max_share,
        "largest_single_study_share_touched_any": max_share_any,
        "any_single_study_exceeds_10pct_of_reduced_population_exclusive": bool(
            max_share is not None and max_share > 0.10
        ),
        "any_single_study_exceeds_10pct_of_reduced_population_touched_any": bool(
            max_share_any is not None and max_share_any > 0.10
        ),
        "n_studies_with_any_exclusive_labeled_peptide": len(a_excl_by_study),
    }

    # ---- 5. sequence structure after removal (variant A labeled) ----------
    sizes = hamming_components(a_labeled, max_d=1)
    flags = singleton_flags(a_labeled)
    n_singletons_from_sizes = sum(1 for s in sizes if s == 1)
    n_singletons_from_degree = sum(1 for v in flags.values() if v)
    sequence = {
        "population": "variant_A labeled peptides",
        "n_peptides": len(a_labeled),
        "hamming_le1_same_length_n_components": len(sizes),
        "n_singleton_components": n_singletons_from_sizes,
        "n_non_singleton_components": sum(1 for s in sizes if s > 1),
        "largest_component_size": sizes[0] if sizes else 0,
        "n_peptides_in_non_singleton_components": sum(s for s in sizes if s > 1),
        "singleton_count_cross_check_degree_zero": n_singletons_from_degree,
        "singleton_counts_agree": n_singletons_from_sizes == n_singletons_from_degree,
    }

    # ---- 6. joint split feasibility (variant A) ---------------------------
    study_exclusive = {p for peps in a_excl_by_study.values() for p in peps}
    joint = sorted(p for p in a_labeled if p in study_exclusive and flags[p])
    joint_pos = sum(1 for p in joint if a_pattern[p] == "positive_only")
    rng = random.Random(SEED)
    sweep = study_holdout_sweep(dict(a_excl_by_study), a_pattern, n_a_lab, rng)
    joint_block = {
        "population": "variant_A labeled peptides",
        "n_study_exclusive_and_sequence_singleton": len(joint),
        "n_positive_only": joint_pos,
        "n_negative_only": len(joint) - joint_pos,
        "positive_fraction": frac(joint_pos, len(joint)),
        "share_of_reduced_labeled_population": frac(len(joint), n_a_lab),
        "n_study_exclusive_labeled": len(study_exclusive),
        "study_holdout_sweep": sweep,
    }

    return {
        "subset_id": spec_id,
        "subset_label": label,
        "spec": spec,
        "baseline_population": baseline,
        "dominant_study_contribution": dominant,
        "variant_A_drop_exclusive_only": variant_a,
        "variant_B_drop_all_touched": variant_b,
        "concentration_after_removal_variant_A": concentration,
        "sequence_structure_variant_A": sequence,
        "joint_split_feasibility_variant_A": joint_block,
    }


def print_subset(res: dict) -> None:
    b = res["baseline_population"]
    d = res["dominant_study_contribution"]
    a = res["variant_A_drop_exclusive_only"]
    vb = res["variant_B_drop_all_touched"]
    c = res["concentration_after_removal_variant_A"]
    s = res["sequence_structure_variant_A"]
    j = res["joint_split_feasibility_variant_A"]
    print(f"\n{'=' * 78}\n{res['subset_id']}  ({res['subset_label']})\n{'=' * 78}")
    print("[1] baseline population")
    print(
        f"    eligible peptides {b['n_eligible_peptides']}  "
        f"pos-only {b['n_positive_only']}  neg-only {b['n_negative_only']}  "
        f"mixed {b['n_mixed']}  uninterpretable {b['n_uninterpretable']}"
    )
    print(
        f"    labeled {b['n_labeled']}  positive fraction {b['labeled_positive_fraction']}  "
        f"studies {b['n_studies']}  eligible obs {b['n_eligible_observations']}  "
        f"rows missing reference_id {b['n_eligible_rows_missing_reference_id']}"
    )
    print(f"[2] reference {d['reference_id']} contribution")
    print(
        f"    peptides touched {d['n_peptides_touched']}  exclusive {d['n_peptides_exclusive']}  "
        f"shared {d['n_peptides_shared_with_other_studies']}  "
        f"eligible obs {d['n_eligible_observations_from_this_study']}"
    )
    print(f"    within-study pattern (touched): {d['within_study_pattern_all_touched']}")
    print(f"    entirely negative within study: {d['entirely_negative_within_study']}")
    print(
        f"    global pattern of exclusives: {d['global_pattern_of_exclusive_peptides']}  "
        f"matches within-study: {d['exclusive_global_matches_within_study']}"
    )
    print(
        f"    labeled touched {d['n_labeled_touched']} "
        f"({d['share_of_baseline_labeled_touched']} of baseline labeled); "
        f"labeled exclusive {d['n_labeled_exclusive']} "
        f"({d['share_of_baseline_labeled_exclusive']})"
    )
    print("[3] removal variants")
    print(
        f"    A labeled {a['n_labeled']}  P {a['n_positive_only']}  N {a['n_negative_only']}  "
        f"posfrac {a['labeled_positive_fraction']}  studies {a['n_studies_remaining_any_eligible']} "
        f"(labeled-only {a['n_studies_remaining_among_labeled']})"
    )
    print(
        f"      shared kept {a['n_shared_peptides_kept']}  relabeled by removal "
        f"{a['n_shared_peptides_relabeled_by_removal']} {a['shared_peptide_pattern_change']}  "
        f"lost all evidence {a['n_shared_peptides_losing_all_eligible_evidence']}"
    )
    print(
        f"    B labeled {vb['n_labeled']}  P {vb['n_positive_only']}  N {vb['n_negative_only']}  "
        f"posfrac {vb['labeled_positive_fraction']}  studies {vb['n_studies_remaining_any_eligible']} "
        f"(labeled-only {vb['n_studies_remaining_among_labeled']})"
    )
    print(f"[4] concentration after removal (variant A, {c['n_labeled_reduced_population']} labeled)")
    print("    ref | excl labeled | P | N | posfrac | share | touched-any | share-any")
    for r in c["top10_studies_by_exclusive_labeled_peptides"]:
        print(
            f"    {r['reference_id']} | {r['n_exclusive_labeled_peptides']} | "
            f"{r['n_positive_only']} | {r['n_negative_only']} | {r['positive_fraction']} | "
            f"{r['share_of_reduced_labeled_population']} | {r['n_labeled_peptides_touched_any']} | "
            f"{r['share_touched_any_of_reduced_labeled']}"
        )
    print(
        f"    any single study >10% of reduced labeled population? exclusive: "
        f"{c['any_single_study_exceeds_10pct_of_reduced_population_exclusive']} "
        f"(max {c['largest_single_study_share_exclusive']}); touched-any: "
        f"{c['any_single_study_exceeds_10pct_of_reduced_population_touched_any']} "
        f"(max {c['largest_single_study_share_touched_any']})"
    )
    print("[5] sequence structure (variant A labeled)")
    print(
        f"    Hamming<=1 same-length components {s['hamming_le1_same_length_n_components']}  "
        f"singletons {s['n_singleton_components']}  non-singletons {s['n_non_singleton_components']}  "
        f"largest {s['largest_component_size']}  peptides in non-singletons "
        f"{s['n_peptides_in_non_singleton_components']}  (singleton cross-check "
        f"{s['singleton_counts_agree']})"
    )
    print("[6] joint split feasibility (variant A)")
    print(
        f"    study-exclusive AND sequence-singleton: {j['n_study_exclusive_and_sequence_singleton']}  "
        f"P {j['n_positive_only']}  N {j['n_negative_only']}  posfrac {j['positive_fraction']}  "
        f"share {j['share_of_reduced_labeled_population']}"
    )
    w = j["study_holdout_sweep"]
    print(
        f"    study holdout to ~20%: studies {w['n_studies_held_out_median']} "
        f"[{w['n_studies_held_out_min']}-{w['n_studies_held_out_max']}]  peptides "
        f"{w['n_test_peptides_median']} [{w['n_test_peptides_min']}-{w['n_test_peptides_max']}]"
    )
    print(
        f"    test positive fraction median {w['test_positive_fraction_median']}  "
        f"min {w['test_positive_fraction_min']}  max {w['test_positive_fraction_max']}  "
        f"({w['n_random_orderings']} orderings, {w['n_studies_available_with_exclusive_labeled_peptides']} "
        f"studies available)"
    )
    print(f"    representative draw: {w['representative_draw']}")


def main() -> None:
    specs = {sid: (label, spec) for sid, label, spec in CANDIDATES}
    obs = read_csv(FULL / "observations.csv")
    by_pep: dict[str, list[dict]] = defaultdict(list)
    for r in obs:
        by_pep[r["export_peptide"]].append({k: r.get(k) for k in KEEP_FIELDS})
    del obs

    results = []
    for sid in SUBSETS:
        label, spec = specs[sid]
        results.append(analyse(sid, label, spec, by_pep))
        print_subset(results[-1])

    out = {
        "purpose": (
            "Exact population counts for two predeclared V2 task candidates, with and "
            "without reference 1001865. Input to a human task-definition decision."
        ),
        "v2_task_chosen": False,
        "counts_only_no_labels_written": True,
        "v1_unchanged": True,
        "split_chosen": False,
        "collection": "historical_8_11_preconflict_54448",
        "source_observations": str((FULL / "observations.csv").relative_to(ROOT)),
        "dominant_reference_id": DOMINANT_REF,
        "missing_reference_convention": (
            f"rows without a reference_id are grouped under {MISSING_REF} as a pseudo-study, "
            "as in per_study_table.py; counts of such rows are reported per subset"
        ),
        "subsets": results,
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    print("v2_task_chosen=false  counts_only_no_labels_written=true")


if __name__ == "__main__":
    main()
