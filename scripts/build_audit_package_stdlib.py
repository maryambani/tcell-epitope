"""Stdlib-only audit artifacts (no numpy/pandas). Safe in restricted environments."""
from __future__ import annotations

import csv
import json
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audit"
PROC = ROOT / "data" / "processed"
FIXED = PROC / "fixed"
AUDIT.mkdir(exist_ok=True)
random.seed(42)

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


def read_csv_dicts(path: Path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def peptide_col(fieldnames):
    for c in fieldnames:
        if c.endswith("|Name") and c.startswith("Epitope") and "Related" not in c:
            return c
    raise KeyError("peptide col")


def epitope_id_col(fieldnames):
    for c in fieldnames:
        if c.startswith("Epitope ID|"):
            return c
    raise KeyError("id col")


def col_endswith(fieldnames, suffix, must_start=None, must_not_start=None):
    for c in fieldnames:
        if not c.endswith(suffix):
            continue
        if must_start and not c.startswith(must_start):
            continue
        if must_not_start and c.startswith(must_not_start):
            continue
        return c
    return None


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    inventory_rows, field_rows, sample_rows, label_traces = [], [], [], []

    for fname, source, label in EXPORTS:
        print("Loading", fname)
        fields, rows = read_csv_dicts(FIXED / fname)
        pcol = peptide_col(fields)
        idcol = epitope_id_col(fields)
        modcol = col_endswith(fields, "|Modifications", "Epitope", "Related")
        orgcol = col_endswith(fields, "|Source Organism", "Epitope", "Related")
        # Species columns: prefer Epitope.*.Species not Related
        spcol = None
        for c in fields:
            if c.endswith("|Species") and c.startswith("Epitope") and "Related" not in c:
                spcol = c
                break
        molcol = col_endswith(fields, "|Source Molecule", "Epitope", "Related")
        otcol = col_endswith(fields, "|Object Type", "Epitope", "Related")
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
            pool = list(buckets[key])
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
                    "object_type": r.get(otcol) if otcol else "",
                    "modifications": r.get(modcol) if modcol else "",
                    "source_molecule": r.get(molcol) if molcol else "",
                    "source_organism": r.get(orgcol) if orgcol else "",
                    "epitope_species": r.get(spcol) if spcol else "",
                    "related_object_name": r.get(relname) if relname else "",
                    "related_object_species": r.get(relsp) if relsp else "",
                }
            )

        valid_short = [
            r for r in rows if is_valid(r.get(pcol) or "") and 8 <= len(r.get(pcol) or "") <= 11
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
                    "object_type": r.get(otcol) if otcol else "",
                    "modifications": r.get(modcol) if modcol else "",
                    "source_molecule": r.get(molcol) if molcol else "",
                    "source_organism": r.get(orgcol) if orgcol else "",
                    "epitope_species": r.get(spcol) if spcol else "",
                    "related_object_name": r.get(relname) if relname else "",
                    "related_object_species": r.get(relsp) if relsp else "",
                    "label_rule": "label from filename (_pos_=1/_neg_=0); no assay-outcome column in export",
                }
            )

    write_csv(AUDIT / "export_inventory.csv", inventory_rows)
    write_csv(AUDIT / "field_missingness.csv", field_rows)
    write_csv(AUDIT / "sample_export_rows.csv", sample_rows)

    _, raw_rows = read_csv_dicts(PROC / "combined_raw.csv")
    _, clean_rows = read_csv_dicts(PROC / "combined_cleaned.csv")
    labels_by_pep = defaultdict(set)
    rows_by_pep = defaultdict(list)
    for r in raw_rows:
        labels_by_pep[r["peptide"]].add(int(r["label"]))
        rows_by_pep[r["peptide"]].append(r)
    conflict_peps = {p for p, labs in labels_by_pep.items() if len(labs) > 1}
    clean_peps = {r["peptide"] for r in clean_rows}

    for t in label_traces:
        t["survives_to_cleaned"] = t["peptide"] in clean_peps and t["peptide"] not in conflict_peps
    write_csv(AUDIT / "label_trace_examples.csv", label_traces)

    def stage_stats(rows, stage, **extra):
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
            removed_reason="same peptide has both labels",
            n_removed_peptides=len(conflict_peps),
        )
    )
    seen, dedup = set(), []
    for r in no_conf:
        if r["peptide"] in seen:
            continue
        seen.add(r["peptide"])
        dedup.append(r)
    attrition.append(
        stage_stats(
            dedup,
            "3_dedupe_one_row_per_peptide",
            removed_reason="duplicate peptide rows",
            n_removed_rows=len(no_conf) - len(dedup),
        )
    )
    len_filt = [r for r in dedup if 8 <= len(r["peptide"]) <= 11]
    attrition.append(
        stage_stats(
            len_filt,
            "4_length_filter_8_to_11",
            removed_reason="length not in [8,11]",
            n_removed_rows=len(dedup) - len(len_filt),
        )
    )
    attrition.append(stage_stats(clean_rows, "3b_combined_cleaned_on_disk", note="clean_data.py output"))
    write_csv(AUDIT / "attrition_table.csv", attrition)

    cross = [
        p
        for p, rs in rows_by_pep.items()
        if p in conflict_peps and {r["source"] for r in rs} >= {"IEDB", "CEDAR"}
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
                    "note": "CONFLICT in combined_raw; dropped entirely in clean_data.py",
                }
            )
    write_csv(AUDIT / "conflict_trace_examples.csv", conflict_trace)

    # Metadata availability summary (critical for reviewer)
    meta = {
        "download_date_documented": "2025-11-03",
        "download_filters_documented": "Linear Peptides, T Cell Assay, Outcome = Positive/Negative (data/README.md)",
        "export_grain": "epitope-level list rows, NOT assay-level records",
        "fields_present_in_exports": [
            "Epitope ID",
            "Epitope Name (sequence)",
            "Modifications / modified residues",
            "Source molecule / parent / organism",
            "Epitope Species",
            "Related Object fields",
        ],
        "fields_absent_from_exports": [
            "Assay ID",
            "Qualitative measure / assay type detail",
            "Host organism (assay host)",
            "MHC allele / MHC class field",
            "Publication / reference ID",
            "Donor / sample ID",
            "Stimulation / processing context",
        ],
        "how_labels_assigned": (
            "File membership only: rows in *_tcell_pos_* get label=1; "
            "rows in *_tcell_neg_* get label=0 (src/prepare_data.py). "
            "No outcome column is read from the CSV."
        ),
        "n_conflicting_peptides": len(conflict_peps),
        "n_raw_rows": len(raw_rows),
        "n_cleaned_peptides": len(clean_rows),
        "n_modeling_8_11mers": len(len_filt),
    }
    (AUDIT / "metadata_inventory.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("Wrote stdlib audit artifacts to", AUDIT)


if __name__ == "__main__":
    main()
