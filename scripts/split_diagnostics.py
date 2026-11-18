"""Candidate-partition feasibility for human + class I + IFNg ELISPOT.

Exploratory only. No model is fitted, no split is adopted, and no label is
assigned to mixed peptides. "pos"/"neg" below are the descriptive
positive_only / negative_only outcome patterns from outcome_pattern() over
eligible observations in this intersection; they are not V2 labels.
"""
from __future__ import annotations

import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from iedb_annotate import nonempty, outcome_pattern  # noqa: E402
from similarity_audit import encode  # noqa: E402
from summarize_task_candidates import eligible  # noqa: E402
from task_feasibility import CANDIDATES, hamming_components, read_csv  # noqa: E402

FULL = ROOT / "audit" / "retrieval" / "full"
SUBSET = "human_class_I_functional_ifng_elispot"
TEST_FRAC = 0.2
N_RANDOM_ORDERS = 200
SEED = 0


def dist(values: list[int]) -> dict:
    if not values:
        return {"n": 0}
    v = sorted(values)
    q = lambda p: v[min(len(v) - 1, int(p * (len(v) - 1) + 0.5))]  # noqa: E731
    return {"n": len(v), "min": v[0], "p25": q(0.25), "median": q(0.5), "p75": q(0.75), "max": v[-1],
            "mean": round(sum(v) / len(v), 2)}


def hamming_component_ids(seqs: list[str], max_d: int = 1) -> dict[str, int]:
    """Same rule as task_feasibility.hamming_components (same length, Hamming <= max_d),
    but returns membership. Sizes are cross-checked against that function in main()."""
    seqs = list(dict.fromkeys(seqs))
    parent = list(range(len(seqs)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    by_len = defaultdict(list)
    for i, s in enumerate(seqs):
        by_len[len(s)].append(i)
    for idxs in by_len.values():
        arr = np.array([encode(seqs[i]) for i in idxs], dtype=np.int8)
        for a in range(len(idxs)):
            hits = np.where((arr[a + 1:] != arr[a]).sum(axis=1) <= max_d)[0]
            for h in hits:
                ra, rb = find(idxs[a]), find(idxs[a + 1 + int(h)])
                if ra != rb:
                    parent[rb] = ra
    return {s: find(i) for i, s in enumerate(seqs)}


def components_to_reach(sizes: list[int], target: int) -> int:
    total = 0
    for k, s in enumerate(sizes, 1):
        total += s
        if total >= target:
            return k
    return len(sizes)


def main() -> None:
    spec = {sid: sp for sid, _, sp in CANDIDATES}[SUBSET]
    obs = read_csv(FULL / "observations.csv")
    elig_by_pep = defaultdict(list)
    for r in obs:
        if eligible(r, spec):
            elig_by_pep[r["export_peptide"]].append(r)
    pattern = {p: outcome_pattern(rows) for p, rows in elig_by_pep.items()}
    pos = {p for p, pat in pattern.items() if pat == "positive_only"}
    neg = {p for p, pat in pattern.items() if pat == "negative_only"}
    labeled = sorted(pos | neg)
    n_lab = len(labeled)
    pep_refs = {
        p: {r["reference_id"] for r in elig_by_pep[p] if nonempty(r.get("reference_id"))} for p in labeled
    }

    # ---- A. study hold-out (exclusive peptides of one study) ----
    excl_by_ref = defaultdict(list)
    shared_by_ref = Counter()
    for p, refs in pep_refs.items():
        if len(refs) == 1:
            excl_by_ref[next(iter(refs))].append(p)
        else:
            for ref in refs:
                shared_by_ref[ref] += 1
    per_study = []
    for ref, test in excl_by_ref.items():
        t_pos = sum(1 for p in test if p in pos)
        t_neg = len(test) - t_pos
        frac = t_pos / len(test)
        per_study.append({
            "reference_id": ref,
            "n_train": n_lab - len(test),
            "n_test": len(test),
            "n_train_pos": len(pos) - t_pos,
            "n_train_neg": len(neg) - t_neg,
            "n_test_pos": t_pos,
            "n_test_neg": t_neg,
            "test_pos_frac": round(frac, 4),
            "flag_test_pos_frac_lt0.1_or_gt0.9": frac < 0.1 or frac > 0.9,
            "n_train_peptides_also_in_this_study": shared_by_ref.get(ref, 0),
        })
    per_study.sort(key=lambda r: -r["n_test"])
    flagged = [r for r in per_study if r["flag_test_pos_frac_lt0.1_or_gt0.9"]]
    split_a = {
        "description": "Test = peptides whose only eligible study is this one; train = all other "
                       "pos-only/neg-only peptides. Multi-study peptides stay in train.",
        "n_labeled_peptides": n_lab,
        "n_studies_with_exclusive_peptides": len(per_study),
        "n_multi_study_labeled_peptides": sum(1 for r in pep_refs.values() if len(r) > 1),
        "test_size_distribution": dist([r["n_test"] for r in per_study]),
        "train_size_distribution": dist([r["n_train"] for r in per_study]),
        "n_studies_flagged_test_pos_frac_lt0.1_or_gt0.9": len(flagged),
        "n_flagged_with_n_test_ge_20": sum(1 for r in flagged if r["n_test"] >= 20),
        "n_studies_with_n_test_ge_20": sum(1 for r in per_study if r["n_test"] >= 20),
        "leakage_note": "n_train_peptides_also_in_this_study counts multi-study peptides kept in train "
                        "that were also assayed by the held-out study (study-level leakage).",
        "per_study": per_study,
    }

    # ---- B. sequence hold-out (Hamming <= 1, same length) ----
    comp_of = hamming_component_ids(labeled, 1)
    members = defaultdict(list)
    for p, c in comp_of.items():
        members[c].append(p)
    sizes = sorted((len(m) for m in members.values()), reverse=True)
    assert sizes == hamming_components(labeled, 1), "component sizes disagree with task_feasibility"
    target = math.ceil(TEST_FRAC * n_lab)
    rng = random.Random(SEED)
    rand_k, rand_test = [], []
    shuffled = list(sizes)
    for _ in range(N_RANDOM_ORDERS):
        rng.shuffle(shuffled)
        k = components_to_reach(shuffled, target)
        rand_k.append(k)
        rand_test.append(sum(shuffled[:k]))
    comp_label_mix = Counter(
        "pos_only" if all(p in pos for p in m) else "neg_only" if all(p in neg for p in m) else "both"
        for m in members.values() if len(m) > 1
    )
    k_large = components_to_reach(sizes, target)
    k_small = components_to_reach(sorted(sizes), target)
    split_b = {
        "description": "Each Hamming<=1 same-length component is a group; held-out components form test.",
        "graph_built_on": "pos-only/neg-only peptides only (mixed peptides excluded; they could bridge components)",
        "n_components": len(sizes),
        "component_size_distribution": dist(sizes),
        "component_size_counts": {str(s): n for s, n in sorted(Counter(sizes).items())},
        "n_singletons": sum(1 for s in sizes if s == 1),
        "n_non_singletons": sum(1 for s in sizes if s > 1),
        "n_peptides_in_non_singletons": sum(s for s in sizes if s > 1),
        "non_singleton_component_outcome_mix": dict(comp_label_mix),
        "test_fraction_target": TEST_FRAC,
        "target_n_test": target,
        "n_components_to_hold_out": {
            "largest_first": k_large,
            "smallest_first": k_small,
            f"random_order_median_over_{N_RANDOM_ORDERS}": sorted(rand_k)[len(rand_k) // 2],
        },
        "estimated_n_test": sorted(rand_test)[len(rand_test) // 2],
        "estimated_n_train": n_lab - sorted(rand_test)[len(rand_test) // 2],
        "note": "Component sizes are small, so the 20% target is met to within one component "
                "(n_test overshoots target by at most the largest component size).",
    }

    # ---- C. combined: exclusive to one study AND sequence singleton ----
    comp_size = {c: len(m) for c, m in members.items()}
    both = [p for p in labeled if len(pep_refs[p]) == 1 and comp_size[comp_of[p]] == 1]
    c_pos = sum(1 for p in both if p in pos)
    by_ref_c = Counter(next(iter(pep_refs[p])) for p in both)
    top_ref, top_n = by_ref_c.most_common(1)[0] if by_ref_c else ("", 0)
    split_c = {
        "description": "Peptides exclusive to one study AND in a singleton Hamming<=1 component.",
        "n_peptides": len(both),
        "n_pos": c_pos,
        "n_neg": len(both) - c_pos,
        "frac_of_labeled": round(len(both) / max(n_lab, 1), 4),
        "n_studies_represented": len(by_ref_c),
        "largest_study_reference_id": top_ref,
        "largest_study_n_peptides": top_n,
    }

    out = {
        "v2_labels_defined": False,
        "split_adopted": False,
        "model_fitted": False,
        "mixed_peptides_labeled": False,
        "subset_id": SUBSET,
        "n_eligible_peptides": len(elig_by_pep),
        "n_positive_only": len(pos),
        "n_negative_only": len(neg),
        "n_mixed_excluded": sum(1 for p in pattern.values() if p == "mixed"),
        "A_study_holdout_exclusive": split_a,
        "B_sequence_holdout_hamming1": split_b,
        "C_study_and_sequence": split_c,
        "note": "Exploratory partition feasibility. pos/neg are descriptive outcome patterns, not V2 labels.",
    }
    (FULL / "split_diagnostics.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    print(f"Split diagnostics: {SUBSET} (exploratory; no split adopted; no model fitted)")
    print(f"labeled (pos-only/neg-only): {n_lab}  pos {len(pos)}  neg {len(neg)}  mixed excluded {out['n_mixed_excluded']}")
    print("\nA. Study hold-out (exclusive peptides)")
    print(f"  studies with exclusive peptides: {len(per_study)}; multi-study labeled peptides: {split_a['n_multi_study_labeled_peptides']}")
    print(f"  test size: {split_a['test_size_distribution']}")
    print(f"  studies with n_test>=20: {split_a['n_studies_with_n_test_ge_20']}; "
          f"flagged pos_frac<0.1 or >0.9: {len(flagged)} ({split_a['n_flagged_with_n_test_ge_20']} with n_test>=20)")
    for r in per_study[:5]:
        print(f"  ref {r['reference_id']}: test {r['n_test']} (pos {r['n_test_pos']}, neg {r['n_test_neg']}, "
              f"frac {r['test_pos_frac']}), train {r['n_train']}, shared-in-train {r['n_train_peptides_also_in_this_study']}")
    print("\nB. Sequence hold-out (Hamming<=1, same length)")
    print(f"  components {len(sizes)}; singletons {split_b['n_singletons']}; non-singletons {split_b['n_non_singletons']} "
          f"({split_b['n_peptides_in_non_singletons']} peptides); largest {sizes[0] if sizes else 0}")
    print(f"  non-singleton outcome mix: {dict(comp_label_mix)}")
    print(f"  20% target n_test={target}: components to hold out {split_b['n_components_to_hold_out']}; "
          f"est n_train {split_b['estimated_n_train']}, n_test {split_b['estimated_n_test']}")
    print("\nC. Study-exclusive AND sequence-singleton")
    print(f"  {len(both)} peptides ({split_c['frac_of_labeled']:.1%} of labeled): pos {c_pos}, neg {len(both) - c_pos}; "
          f"{len(by_ref_c)} studies, largest {top_ref} with {top_n}")
    print(f"\nwrote {(FULL / 'split_diagnostics.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
