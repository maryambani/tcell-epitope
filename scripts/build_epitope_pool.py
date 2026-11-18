"""Build the pre-conflict 8-11mer epitope pool with provenance flags.

Does not overwrite V1 cleaned data, splits, or models.
"""
from __future__ import annotations

import csv
import json
import random
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXED = ROOT / "data" / "processed" / "fixed"
AUDIT = ROOT / "audit"
AUDIT.mkdir(exist_ok=True)

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
ID_RE = re.compile(r"/epitope/(\d+)")

EXPORTS = [
    ("iedb_tcell_pos_fixed.csv", "IEDB", "positive"),
    ("iedb_tcell_neg_fixed.csv", "IEDB", "negative"),
    ("cedar_tcell_pos_fixed.csv", "CEDAR", "positive"),
    ("cedar_tcell_neg_fixed.csv", "CEDAR", "negative"),
]


def is_valid(seq: str) -> bool:
    if not isinstance(seq, str):
        return False
    seq = seq.strip().upper()
    return bool(seq) and all(ch in STANDARD_AA for ch in seq)


def peptide_col(fields):
    for c in fields:
        if c.endswith("|Name") and c.startswith("Epitope") and "Related" not in c:
            return c
    raise KeyError("peptide column")


def id_col(fields):
    for c in fields:
        if c.startswith("Epitope ID|"):
            return c
    raise KeyError("epitope id column")


def parse_epitope_id(iri: str) -> str:
    if not iri:
        return ""
    m = ID_RE.search(iri)
    return m.group(1) if m else iri.strip()


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    # peptide -> accumulators
    recs = defaultdict(
        lambda: {
            "in_iedb_pos": False,
            "in_iedb_neg": False,
            "in_cedar_pos": False,
            "in_cedar_neg": False,
            "iedb_ids": set(),
            "cedar_ids": set(),
            "n_export_rows": 0,
        }
    )
    n_export_rows = 0
    n_invalid_rows = 0

    for fname, database, membership in EXPORTS:
        path = FIXED / fname
        print("Reading", fname)
        with path.open(newline="", encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f)
            fields = list(reader.fieldnames or [])
            pcol = peptide_col(fields)
            icol = id_col(fields)
            for row in reader:
                n_export_rows += 1
                pep = (row.get(pcol) or "").strip()
                eid = parse_epitope_id(row.get(icol) or "")
                recs[pep]["n_export_rows"] += 1
                flag = f"in_{database.lower()}_{membership[:3]}"
                recs[pep][flag] = True
                if database == "IEDB" and eid:
                    recs[pep]["iedb_ids"].add(eid)
                if database == "CEDAR" and eid:
                    recs[pep]["cedar_ids"].add(eid)
                if not is_valid(pep):
                    n_invalid_rows += 1

    all_peptides = list(recs)
    valid = [p for p in all_peptides if is_valid(p)]
    conflict_all = [
        p
        for p in all_peptides
        if (recs[p]["in_iedb_pos"] or recs[p]["in_cedar_pos"])
        and (recs[p]["in_iedb_neg"] or recs[p]["in_cedar_neg"])
    ]
    conflict_valid = [p for p in valid if p in set(conflict_all)]
    retained_valid = [p for p in valid if p not in set(conflict_valid)]
    pool_811 = [p for p in valid if 8 <= len(p) <= 11]
    pool_811_conflict = [p for p in pool_811 if p in set(conflict_valid)]
    pool_811_retained = [p for p in pool_811 if p not in set(conflict_valid)]

    counts = {
        "n_export_rows": n_export_rows,
        "n_unique_peptides_any": len(all_peptides),
        "n_unique_valid_aa": len(valid),
        "n_conflicting_peptides_any_including_invalid_aa": len(conflict_all),
        "n_conflicting_peptides_valid_aa": len(conflict_valid),
        "n_retained_valid_aa_any_length": len(retained_valid),
        "n_valid_8_11": len(pool_811),
        "n_valid_8_11_conflict": len(pool_811_conflict),
        "n_valid_8_11_retained": len(pool_811_retained),
        "note": (
            "26,504-style count is unique peptides with both pos and neg membership "
            "before requiring standard AA. The valid-AA conflict count is smaller. "
            "V1 modeling set is the length 8-11 slice of retained valid peptides."
        ),
    }
    (AUDIT / "conflict_count_reconciliation.json").write_text(
        json.dumps(counts, indent=2), encoding="utf-8"
    )
    print(json.dumps(counts, indent=2))

    pool_rows = []
    for p in pool_811:
        r = recs[p]
        pos = r["in_iedb_pos"] or r["in_cedar_pos"]
        neg = r["in_iedb_neg"] or r["in_cedar_neg"]
        conflict = pos and neg
        v1_label = ""
        if not conflict:
            v1_label = "1" if pos else "0"
        ids = sorted(r["iedb_ids"] | r["cedar_ids"], key=lambda x: int(x) if x.isdigit() else x)
        pool_rows.append(
            {
                "peptide": p,
                "length": len(p),
                "in_iedb": int(r["in_iedb_pos"] or r["in_iedb_neg"]),
                "in_cedar": int(r["in_cedar_pos"] or r["in_cedar_neg"]),
                "in_iedb_pos": int(r["in_iedb_pos"]),
                "in_iedb_neg": int(r["in_iedb_neg"]),
                "in_cedar_pos": int(r["in_cedar_pos"]),
                "in_cedar_neg": int(r["in_cedar_neg"]),
                "v1_conflict": int(conflict),
                "v1_retained": int(not conflict),
                "v1_query_label": v1_label,
                "n_export_rows": r["n_export_rows"],
                "epitope_ids": "|".join(ids),
                "primary_epitope_id": ids[0] if ids else "",
            }
        )
    write_csv(AUDIT / "epitope_pool_8_11_preconflict.csv", pool_rows)
    print("Wrote", AUDIT / "epitope_pool_8_11_preconflict.csv", "n=", len(pool_rows))

    # Stratified pilot sample (epitope IDs, not only peptides)
    rng = random.Random(42)

    def pick(pred, n, tag):
        cands = [row for row in pool_rows if pred(row)]
        rng.shuffle(cands)
        out = []
        for row in cands[:n]:
            item = dict(row)
            item["pilot_stratum"] = tag
            out.append(item)
        return out, len(cands)

    sample, sizes = [], {}
    chunks = [
        (lambda r: r["v1_retained"] == 1 and r["v1_query_label"] == "1" and r["in_cedar"] == 0, 10, "retained_pos_iedb_only"),
        (lambda r: r["v1_retained"] == 1 and r["v1_query_label"] == "0" and r["in_cedar"] == 0, 10, "retained_neg_iedb_only"),
        (lambda r: r["v1_retained"] == 1 and r["in_iedb"] == 1 and r["in_cedar"] == 1, 10, "retained_both_databases"),
        (lambda r: r["v1_conflict"] == 1 and r["in_iedb"] == 1 and r["in_cedar"] == 0, 8, "conflict_iedb_only"),
        (lambda r: r["v1_conflict"] == 1 and r["in_iedb"] == 1 and r["in_cedar"] == 1, 8, "conflict_both_databases"),
        (lambda r: r["in_iedb"] == 0 and r["in_cedar"] == 1, 4, "cedar_only"),
    ]
    for pred, n, tag in chunks:
        part, n_avail = pick(pred, n, tag)
        sample.extend(part)
        sizes[tag] = {"n_available": n_avail, "n_sampled": len(part)}

    write_csv(AUDIT / "metadata_pilot_sample.csv", sample)
    (AUDIT / "metadata_pilot_sample_sizes.json").write_text(json.dumps(sizes, indent=2), encoding="utf-8")
    print("Pilot sample", len(sample), "strata", sizes)


if __name__ == "__main__":
    main()
