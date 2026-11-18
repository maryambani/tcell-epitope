"""Fuller V1 similarity audit. Does not choose a cutoff or retune models."""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audit"

AA = "ACDEFGHIKLMNPQRSTVWY"


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def encode(seq: str) -> np.ndarray:
    table = {c: i for i, c in enumerate(AA)}
    x = np.array([table.get(c, -1) for c in seq], dtype=np.int8)
    return x


def min_hamming_to_train(test_seq: str, train_arr: np.ndarray, train_seqs: list[str]) -> tuple[int, int]:
    """Return (min_distance, n_neighbors_at_that_distance) for same-length train set."""
    t = encode(test_seq)
    diffs = (train_arr != t).sum(axis=1)
    dmin = int(diffs.min())
    n = int((diffs == dmin).sum())
    return dmin, n


def union_find_components(seqs: list[str], max_d: int) -> list[int]:
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
    for L, idxs in by_len.items():
        arr = np.array([encode(seqs[i]) for i in idxs], dtype=np.int8)
        m = len(idxs)
        for a in range(m):
            diffs = (arr[a + 1 :] != arr[a]).sum(axis=1)
            hits = np.where(diffs <= max_d)[0]
            for h in hits:
                union(idxs[a], idxs[a + 1 + int(h)])
    size = Counter(find(i) for i in range(n))
    return sorted(size.values(), reverse=True)


def main() -> None:
    split_path = AUDIT / "split_assignments.csv"
    rows = list(csv.DictReader(split_path.open()))
    train = [r for r in rows if r["partition"] == "train"]
    test = [r for r in rows if r["partition"] == "test"]
    train_by_len = defaultdict(list)
    train_lab_by_len = defaultdict(list)
    for r in train:
        train_by_len[int(r["length"])].append(r["peptide"])
        train_lab_by_len[int(r["length"])].append(int(r["observed_label"]))
    train_arr = {L: np.array([encode(s) for s in seqs], dtype=np.int8) for L, seqs in train_by_len.items()}

    nearest_rows = []
    dist_hist = Counter()
    dist_hist_by_len = defaultdict(Counter)
    n_d_le1 = 0
    n_d_le1_same_label = 0
    n_d_le1_diff_label = 0
    n_containment = 0

    train_set = {r["peptide"] for r in train}
    train_list = [r["peptide"] for r in train]

    for r in test:
        p = r["peptide"]
        L = int(r["length"])
        y = int(r["observed_label"])
        arr = train_arr[L]
        dmin, n_at = min_hamming_to_train(p, arr, train_by_len[L])
        # labels among Hamming<=1 neighbors
        diffs = (arr != encode(p)).sum(axis=1)
        neigh = np.where(diffs <= 1)[0]
        neigh_labs = [train_lab_by_len[L][int(i)] for i in neigh]
        n_same = sum(lab == y for lab in neigh_labs)
        n_diff = len(neigh_labs) - n_same
        # different-length containment
        contained_in_train = any((p in t or t in p) and t != p for t in train_list if abs(len(t) - L) >= 1)
        if dmin <= 1:
            n_d_le1 += 1
            if n_same and not n_diff:
                n_d_le1_same_label += 1
            elif n_diff and not n_same:
                n_d_le1_diff_label += 1
        if contained_in_train:
            n_containment += 1
        dist_hist[dmin] += 1
        dist_hist_by_len[L][dmin] += 1
        nearest_rows.append(
            {
                "peptide": p,
                "length": L,
                "observed_label": y,
                "min_hamming_same_length_train": dmin,
                "n_train_at_min_distance": n_at,
                "n_train_hamming_le_1": int(len(neigh)),
                "n_hamming_le_1_same_label": n_same,
                "n_hamming_le_1_diff_label": n_diff,
                "different_length_containment_in_train": int(contained_in_train),
            }
        )

    # Component sizes under Hamming<=1 (same length only), on full V1 modeling set
    all_seqs = [r["peptide"] for r in rows]
    sizes_d1 = union_find_components(all_seqs, max_d=1)
    sizes_d0 = union_find_components(all_seqs, max_d=0)  # should be all 1s if unique

    def size_summary(sizes):
        arr = np.array(sizes)
        return {
            "n_components": int(len(sizes)),
            "max_size": int(arr.max()) if len(arr) else 0,
            "median_size": float(np.median(arr)) if len(arr) else None,
            "n_size_ge_2": int((arr >= 2).sum()),
            "n_size_ge_10": int((arr >= 10).sum()),
            "n_peptides_in_size_ge_2": int(arr[arr >= 2].sum()) if (arr >= 2).any() else 0,
            "top10_sizes": [int(x) for x in sizes[:10]],
        }

    n_test = len(test)
    summary = {
        "split": "V1 reconstructed stratified random holdout (random_state=42)",
        "n_train": len(train),
        "n_test": n_test,
        "exact_overlap": len({r["peptide"] for r in train} & {r["peptide"] for r in test}),
        "full_test_same_length_hamming": {
            "n_test_with_min_d_le_1": n_d_le1,
            "frac_test_with_min_d_le_1": round(n_d_le1 / n_test, 4),
            "n_those_only_same_label_neighbors": n_d_le1_same_label,
            "n_those_only_diff_label_neighbors": n_d_le1_diff_label,
            "distance_histogram": {str(k): dist_hist[k] for k in sorted(dist_hist)},
            "by_length_frac_d_le_1": {
                str(L): round(dist_hist_by_len[L][0] + dist_hist_by_len[L][1], 4) / Counter(int(r["length"]) for r in test)[L]
                for L in sorted(dist_hist_by_len)
            },
        },
        "different_length_containment": {
            "n_test_contained_or_containing_train": n_containment,
            "frac": round(n_containment / n_test, 4),
            "note": "Substring relationship across lengths; not a grouping rule.",
        },
        "connected_components_same_length_hamming_le_1": size_summary(sizes_d1),
        "connected_components_exact_sequence": size_summary(sizes_d0),
        "prior_quickcheck": {
            "n_sampled": 1500,
            "frac_d_le_1": 0.0807,
            "selection": "numpy default_rng(42) choice without replacement from test peptides",
            "note": "Quickcheck is a subsample. Use full-test numbers above for design, not score chasing.",
        },
        "do_not": (
            "Do not pick a similarity cutoff to maximize or minimize V1 AUC. "
            "Keep the random split as a reference; a similarity-separated evaluation answers a different question."
        ),
    }
    # fix by_length to be cleaner
    n_by_len = Counter(int(r["length"]) for r in test)
    summary["full_test_same_length_hamming"]["by_length"] = {}
    for L in sorted(n_by_len):
        nL = n_by_len[L]
        n_close = dist_hist_by_len[L][0] + dist_hist_by_len[L][1]
        summary["full_test_same_length_hamming"]["by_length"][str(L)] = {
            "n_test": nL,
            "n_min_d_le_1": int(n_close),
            "frac": round(n_close / nL, 4),
            "distance_histogram": {str(k): int(dist_hist_by_len[L][k]) for k in sorted(dist_hist_by_len[L])},
        }
    del summary["full_test_same_length_hamming"]["by_length_frac_d_le_1"]

    with (AUDIT / "similarity_nearest_train.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(nearest_rows[0].keys()))
        w.writeheader()
        w.writerows(nearest_rows)
    (AUDIT / "similarity_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2)[:2500])
    print("Wrote similarity_nearest_train.csv and similarity_summary.json")


if __name__ == "__main__":
    main()
