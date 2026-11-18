"""Build reviewer audit package artifacts under audit/.

Uses csv + numpy only (avoids pandas/TF sandbox crashes). Run predictions
separately if models need TensorFlow.
"""
from __future__ import annotations

import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audit"
DATA = ROOT / "data"
PROC = DATA / "processed"
FIXED = PROC / "fixed"
MODEL_DIR = ROOT / "outputs" / "models"

AUDIT.mkdir(exist_ok=True)
random.seed(42)
np.random.seed(42)

EXPORTS = [
    ("iedb_tcell_pos_fixed.csv", "IEDB", 1),
    ("iedb_tcell_neg_fixed.csv", "IEDB", 0),
    ("cedar_tcell_pos_fixed.csv", "CEDAR", 1),
    ("cedar_tcell_neg_fixed.csv", "CEDAR", 0),
]

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def is_valid(seq: str) -> bool:
    if not isinstance(seq, str):
        return False
    seq = seq.strip().upper()
    return bool(seq) and all(ch in STANDARD_AA for ch in seq)


def read_csv_dicts(path: Path) -> tuple[list[str], list[dict]]:
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def peptide_col(fieldnames: list[str]) -> str:
    for c in fieldnames:
        if c.endswith("|Name") and c.startswith("Epitope"):
            # Prefer Epitope.1|Name over Related Object names
            if "Related" not in c:
                return c
    raise KeyError(f"No peptide name column in {fieldnames[:10]}...")


def epitope_id_col(fieldnames: list[str]) -> str:
    for c in fieldnames:
        if c.startswith("Epitope ID|"):
            return c
    raise KeyError("No epitope id column")


def col_endswith(fieldnames: list[str], suffix: str, prefix: str | None = None) -> str | None:
    for c in fieldnames:
        if c.endswith(suffix) and (prefix is None or c.startswith(prefix)):
            if prefix == "Epitope" and c.startswith("Related"):
                continue
            return c
    return None


def main() -> None:
    inventory_rows = []
    field_rows = []
    sample_rows = []
    label_traces = []

    for fname, source, label in EXPORTS:
        path = FIXED / fname
        print(f"Loading {fname}...")
        fields, rows = read_csv_dicts(path)
        pcol = peptide_col(fields)
        idcol = epitope_id_col(fields)
        modcol = col_endswith(fields, "|Modifications", "Epitope")
        orgcol = col_endswith(fields, "|Source Organism", "Epitope")
        spcol = col_endswith(fields, "|Species", "Epitope")
        molcol = col_endswith(fields, "|Source Molecule", "Epitope")
        otcol = col_endswith(fields, "|Object Type", "Epitope")
        relname = next((c for c in fields if c.startswith("Related Object") and c.endswith("|Name")), None)
        relsp = next((c for c in fields if c.startswith("Related Object") and c.endswith("|Species")), None)

        inventory_rows.append(
            {
                "file": fname,
                "source": source,
                "file_label": label,
                "n_rows": len(rows),
                "n_cols": len(fields),
                "peptide_col_used": pcol,
                "n_unique_peptide_name": len({r[pcol] for r in rows if r.get(pcol)}),
                "n_unique_epitope_id": len({r[idcol] for r in rows if r.get(idcol)}),
            }
        )

        for c in fields:
            non_null = sum(1 for r in rows if r.get(c) not in (None, "", "nan"))
            field_rows.append(
                {
                    "file": fname,
                    "column": c,
                    "non_null_frac": round(non_null / max(len(rows), 1), 4),
                    "n_unique": len({r.get(c) for r in rows if r.get(c) not in (None, "", "nan")}),
                }
            )

        # buckets for sampling
        buckets = {"short_nomod": [], "short_mod": [], "other": []}
        for r in rows:
            pep = (r.get(pcol) or "").strip()
            L = len(pep)
            has_mod = bool((r.get(modcol) or "").strip()) if modcol else False
            if 8 <= L <= 11 and not has_mod:
                buckets["short_nomod"].append(r)
            elif 8 <= L <= 11 and has_mod:
                buckets["short_mod"].append(r)
            else:
                buckets["other"].append(r)

        picks = []
        for key in ("short_nomod", "short_mod", "other"):
            pool = buckets[key]
            random.Random(42).shuffle(pool)
            picks.extend(pool[:4])
        for r in picks[:12]:
            sample_rows.append(
                {
                    "export_file": fname,
                    "source": source,
                    "assigned_label": label,
                    "epitope_id": r.get(idcol),
                    "peptide": r.get(pcol),
                    "object_type": r.get(otcol) if otcol else None,
                    "modifications": r.get(modcol) if modcol else None,
                    "source_molecule": r.get(molcol) if molcol else None,
                    "source_organism": r.get(orgcol) if orgcol else None,
                    "epitope_species": r.get(spcol) if spcol else None,
                    "related_object_name": r.get(relname) if relname else None,
                    "related_object_species": r.get(relsp) if relsp else None,
                }
            )

        # label traces: 5 valid 8-11mers
        valid_short = [
            r
            for r in rows
            if is_valid((r.get(pcol) or "")) and 8 <= len((r.get(pcol) or "")) <= 11
        ]
        random.Random(42).shuffle(valid_short)
        for r in valid_short[:5]:
            label_traces.append(
                {
                    "export_file": fname,
                    "source": source,
                    "assigned_label_from_filename": label,
                    "epitope_id": r.get(idcol),
                    "peptide": r.get(pcol),
                    "object_type": r.get(otcol) if otcol else None,
                    "modifications": r.get(modcol) if modcol else None,
                    "source_molecule": r.get(molcol) if molcol else None,
                    "source_organism": r.get(orgcol) if orgcol else None,
                    "epitope_species": r.get(spcol) if spcol else None,
                    "related_object_name": r.get(relname) if relname else None,
                    "related_object_species": r.get(relsp) if relsp else None,
                    "label_rule": "label taken from export filename (_pos_=1, _neg_=0); no assay-outcome column in file",
                }
            )

    write_csv(AUDIT / "export_inventory.csv", inventory_rows)
    write_csv(AUDIT / "field_missingness.csv", field_rows)
    write_csv(AUDIT / "sample_export_rows.csv", sample_rows)
    write_csv(AUDIT / "label_trace_examples.csv", label_traces)
    print("Wrote export inventory / missingness / samples / traces")

    # ------------------------------------------------------------------
    # Attrition from combined_raw / combined_cleaned
    # ------------------------------------------------------------------
    _, raw_rows = read_csv_dicts(PROC / "combined_raw.csv")
    _, clean_rows = read_csv_dicts(PROC / "combined_cleaned.csv")

    # label conflicts
    labels_by_pep: dict[str, set[int]] = defaultdict(set)
    rows_by_pep: dict[str, list[dict]] = defaultdict(list)
    for r in raw_rows:
        pep = r["peptide"]
        labels_by_pep[pep].add(int(r["label"]))
        rows_by_pep[pep].append(r)
    conflict_peps = {p for p, labs in labels_by_pep.items() if len(labs) > 1}

    def stage_stats(rows: list[dict], stage: str, **extra):
        return {
            "stage": stage,
            "n_rows": len(rows),
            "n_unique_peptides": len({r["peptide"] for r in rows}),
            "n_pos_rows": sum(1 for r in rows if int(r["label"]) == 1),
            "n_neg_rows": sum(1 for r in rows if int(r["label"]) == 0),
            "n_IEDB": sum(1 for r in rows if r["source"] == "IEDB"),
            "n_CEDAR": sum(1 for r in rows if r["source"] == "CEDAR"),
            **extra,
        }

    attrition = [stage_stats(raw_rows, "0_raw_combined_rows")]

    valid_rows = [r for r in raw_rows if is_valid(r["peptide"])]
    attrition.append(
        stage_stats(
            valid_rows,
            "1_valid_standard_AA_only",
            removed_reason="non-standard AA / empty / non-string",
            n_removed_rows=len(raw_rows) - len(valid_rows),
        )
    )

    no_conf = [r for r in valid_rows if r["peptide"] not in conflict_peps]
    attrition.append(
        stage_stats(
            no_conf,
            "2_drop_conflicting_peptides",
            removed_reason="same peptide has both labels in combined_raw",
            n_removed_peptides=len(conflict_peps),
        )
    )

    seen = set()
    dedup = []
    for r in no_conf:
        if r["peptide"] in seen:
            continue
        seen.add(r["peptide"])
        dedup.append(r)
    attrition.append(
        stage_stats(
            dedup,
            "3_dedupe_to_one_row_per_peptide",
            removed_reason="duplicate peptide rows after conflict drop",
            n_removed_rows=len(no_conf) - len(dedup),
        )
    )

    len_filt = [r for r in dedup if 8 <= len(r["peptide"]) <= 11]
    attrition.append(
        stage_stats(
            len_filt,
            "4_length_filter_8_to_11",
            removed_reason="peptide length not in [8,11]",
            n_removed_rows=len(dedup) - len(len_filt),
        )
    )
    # cleaned file check
    attrition.append(
        stage_stats(
            clean_rows,
            "3b_combined_cleaned_csv_on_disk",
            note="output of clean_data.py (should match stage 3)",
        )
    )
    write_csv(AUDIT / "attrition_table.csv", attrition)

    # conflict examples
    cross = [
        p
        for p, rs in rows_by_pep.items()
        if p in conflict_peps and {r["source"] for r in rs} == {"IEDB", "CEDAR"}
    ]
    picks = cross[:5]
    other = [p for p in conflict_peps if p not in picks]
    random.Random(0).shuffle(other)
    picks += other[:7]
    conflict_trace = []
    for pep in picks:
        for r in rows_by_pep[pep]:
            conflict_trace.append(
                {
                    "peptide": pep,
                    "source": r["source"],
                    "label": r["label"],
                    "note": "CONFLICT in combined_raw; peptide dropped entirely in clean_data.py",
                }
            )
    write_csv(AUDIT / "conflict_trace_examples.csv", conflict_trace)
    print("Wrote attrition + conflict traces")

    # Mark which label traces survive cleaning
    clean_peps = {r["peptide"] for r in clean_rows}
    for t in label_traces:
        t["survives_to_cleaned"] = t["peptide"] in clean_peps and t["peptide"] not in conflict_peps
    write_csv(AUDIT / "label_trace_examples.csv", label_traces)

    # ------------------------------------------------------------------
    # Reconstruct modeling split (matches features.npz)
    # ------------------------------------------------------------------
    print("Reconstructing split...")
    model_rows = [r for r in clean_rows if 8 <= len(r["peptide"]) <= 11]
    # preserve cleaned file order
    MAX_LEN = 11
    AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
    aa_to_idx = {aa: i for i, aa in enumerate(AA_LIST)}

    def one_hot(seq: str) -> np.ndarray:
        padded = seq.ljust(MAX_LEN, "-")
        enc = np.zeros((MAX_LEN, 21))
        for i, aa in enumerate(padded):
            if aa in aa_to_idx:
                enc[i, aa_to_idx[aa]] = 1
            else:
                enc[i, 20] = 1
        return enc.flatten()

    peptides = np.array([r["peptide"] for r in model_rows])
    y = np.array([int(r["label"]) for r in model_rows])
    sources = np.array([r["source"] for r in model_rows])
    lengths = np.array([len(p) for p in peptides])
    ids = np.arange(len(peptides))
    X = np.vstack([one_hot(p) for p in peptides])

    X_tr, X_te, y_tr, y_te, pep_tr, pep_te, src_tr, src_te, len_tr, len_te, id_tr, id_te = train_test_split(
        X, y, peptides, sources, lengths, ids, test_size=0.2, random_state=42, stratify=y
    )

    feat = np.load(PROC / "features.npz")
    assert np.allclose(X_te, feat["X_test"]), "X_test mismatch"
    assert np.array_equal(y_te, feat["y_test"]), "y_test mismatch"
    print("Split reconstruction MATCHES features.npz")

    split_rows = []
    for i in range(len(y_tr)):
        split_rows.append(
            {
                "peptide_id": int(id_tr[i]),
                "peptide": pep_tr[i],
                "observed_label": int(y_tr[i]),
                "source": src_tr[i],
                "length": int(len_tr[i]),
                "partition": "train",
            }
        )
    for i in range(len(y_te)):
        split_rows.append(
            {
                "peptide_id": int(id_te[i]),
                "peptide": pep_te[i],
                "observed_label": int(y_te[i]),
                "source": src_te[i],
                "length": int(len_te[i]),
                "partition": "test",
            }
        )
    write_csv(AUDIT / "split_assignments.csv", split_rows)

    # leakage quickcheck
    train_by_len: dict[int, list[str]] = defaultdict(list)
    for p in pep_tr:
        train_by_len[len(p)].append(p)

    def hamming(a: str, b: str) -> int:
        return sum(x != y for x, y in zip(a, b))

    rng = np.random.default_rng(42)
    sample_idx = rng.choice(len(pep_te), size=min(1500, len(pep_te)), replace=False)
    near = 0
    for i in sample_idx:
        p = pep_te[i]
        if any(hamming(p, q) <= 1 for q in train_by_len[len(p)]):
            near += 1
    leakage = {
        "exact_peptide_overlap_train_test": int(len(set(pep_tr) & set(pep_te))),
        "test_sample_size_for_hamming": int(len(sample_idx)),
        "test_frac_with_hamming_le_1_same_length_train_neighbor": float(near / len(sample_idx)),
        "note": (
            "Approximate; subsampled test peptides; same-length Hamming only. "
            "Study/HLA leakage cannot be checked: those fields were absent from exports."
        ),
    }
    (AUDIT / "split_leakage_quickcheck.json").write_text(json.dumps(leakage, indent=2), encoding="utf-8")
    print("Leakage:", leakage)

    # Save arrays needed for prediction step
    np.savez(
        AUDIT / "test_arrays_for_scoring.npz",
        X_test=X_te,
        y_test=y_te,
        peptide_id=id_te,
        peptide=pep_te,
        source=src_te,
        length=len_te,
    )
    print("Wrote split assignments + test arrays")

    # Try scoring without TF first for LR/RF; CNN optional
    try:
        import joblib

        lr = joblib.load(MODEL_DIR / "logistic_regression.joblib")
        rf = joblib.load(MODEL_DIR / "random_forest.joblib")
        scores = {
            "LogisticRegression": lr.predict_proba(X_te)[:, 1],
            "RandomForest": rf.predict_proba(X_te)[:, 1],
        }
        try:
            import tensorflow as tf

            cnn = tf.keras.models.load_model(MODEL_DIR / "cnn_epitope.keras")
            scores["CNN"] = cnn.predict(X_te.reshape(-1, 11, 21), verbose=0).flatten()
        except Exception as e:
            print("CNN scoring skipped:", e)

        pred_rows = []
        metric_rows = []
        for name, sc in scores.items():
            metric_rows.append(
                {
                    "model": name,
                    "partition": "test",
                    "n": len(y_te),
                    "roc_auc": float(roc_auc_score(y_te, sc)),
                    "average_precision": float(average_precision_score(y_te, sc)),
                    "ap_definition": "sklearn.metrics.average_precision_score(y_true, continuous scores)",
                }
            )
            for i in range(len(y_te)):
                pred_rows.append(
                    {
                        "peptide_id": int(id_te[i]),
                        "peptide": pep_te[i],
                        "observed_label": int(y_te[i]),
                        "model_score": float(sc[i]),
                        "model_name": name,
                        "partition": "test",
                        "source": src_te[i],
                        "length": int(len_te[i]),
                    }
                )
        write_csv(AUDIT / "test_predictions.csv", pred_rows)
        write_csv(AUDIT / "test_metrics_from_scores.csv", metric_rows)

        rf_sc = scores["RandomForest"]
        breakdown = []
        for L in sorted(set(len_te.tolist())):
            m = len_te == L
            breakdown.append(
                {
                    "by": "length",
                    "value": int(L),
                    "n": int(m.sum()),
                    "pos_rate": float(y_te[m].mean()),
                    "rf_auc": float(roc_auc_score(y_te[m], rf_sc[m])) if len(set(y_te[m])) > 1 else None,
                    "rf_ap": float(average_precision_score(y_te[m], rf_sc[m])),
                }
            )
        for S in sorted(set(src_te.tolist())):
            m = src_te == S
            breakdown.append(
                {
                    "by": "source",
                    "value": S,
                    "n": int(m.sum()),
                    "pos_rate": float(y_te[m].mean()),
                    "rf_auc": float(roc_auc_score(y_te[m], rf_sc[m])) if len(set(y_te[m])) > 1 else None,
                    "rf_ap": float(average_precision_score(y_te[m], rf_sc[m])),
                }
            )
        write_csv(AUDIT / "rf_test_breakdown_by_length_source.csv", breakdown)
        print("Wrote predictions")
    except Exception as e:
        print("Model scoring failed (run later with TF available):", e)

    print("DONE. Artifacts in", AUDIT)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    # union of keys preserving order from first row then extras
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
