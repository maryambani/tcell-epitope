"""Per-study (reference_id) breakdown for the two shortlisted task intersections.

Exploratory inspection of dominant studies only. Does not define V2 labels,
change the conflict rule, or choose a split. Outcome patterns per study are
computed from that study's eligible observations only and are descriptive.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from iedb_annotate import nonempty, outcome_pattern  # noqa: E402
from summarize_task_candidates import eligible  # noqa: E402
from task_feasibility import CANDIDATES, read_csv, write_csv  # noqa: E402

FULL = ROOT / "audit" / "retrieval" / "full"
SHORTLIST = ("human_class_I_functional", "human_class_I_functional_ifng_elispot")
MISSING_REF = "(missing_reference_id)"
LABELED = ("positive_only", "negative_only")


def top_k(counter: Counter, k: int = 3) -> str:
    return "; ".join(f"{v} ({n})" for v, n in counter.most_common(k))


def ref_of(row: dict) -> str:
    ref = row.get("reference_id")
    return str(ref).strip() if nonempty(ref) else MISSING_REF


def per_study_rows(obs_by_pep: dict[str, list[dict]], spec: dict) -> tuple[list[dict], dict]:
    elig_by_pep = {}
    for pep, rows in obs_by_pep.items():
        sub = [r for r in rows if eligible(r, spec)]
        if sub:
            elig_by_pep[pep] = sub
    global_pattern = {pep: outcome_pattern(rows) for pep, rows in elig_by_pep.items()}

    pep_refs: dict[str, set[str]] = {}
    study_obs: dict[str, list[dict]] = defaultdict(list)
    for pep, rows in elig_by_pep.items():
        pep_refs[pep] = {ref_of(r) for r in rows}
        for r in rows:
            study_obs[ref_of(r)].append(r)

    out = []
    for ref, rows in study_obs.items():
        by_pep = defaultdict(list)
        for r in rows:
            by_pep[r["export_peptide"]].append(r)
        pats = Counter(outcome_pattern(prs) for prs in by_pep.values())
        exclusive = [p for p in by_pep if pep_refs[p] == {ref}]
        out.append(
            {
                "reference_id": ref,
                "n_peptides": len(by_pep),
                "n_peptides_exclusive": len(exclusive),
                "n_positive_only": pats.get("positive_only", 0),
                "n_negative_only": pats.get("negative_only", 0),
                "n_mixed": pats.get("mixed", 0),
                "n_uninterpretable": pats.get("eligible_uninterpretable_outcome", 0),
                "n_peptides_globally_posneg_only": sum(1 for p in by_pep if global_pattern[p] in LABELED),
                "n_exclusive_globally_posneg_only": sum(1 for p in exclusive if global_pattern[p] in LABELED),
                "n_assays": len(rows),
                "top3_host_organism_name": top_k(Counter(r.get("host_organism_name") or "" for r in rows)),
                "top3_mhc_allele_name": top_k(Counter(r.get("mhc_allele_name") or "" for r in rows)),
                "top3_assay_names": top_k(Counter(r.get("assay_names") or "" for r in rows)),
                "pubmed_id": next((r.get("pubmed_id") for r in rows if nonempty(r.get("pubmed_id"))), ""),
            }
        )
    out.sort(key=lambda r: (-r["n_peptides_exclusive"], -r["n_peptides"], r["reference_id"]))
    stats = {
        "n_eligible_peptides": len(elig_by_pep),
        "n_labeled_peptides": sum(1 for p in global_pattern.values() if p in LABELED),
        "n_studies": len(out),
        "n_missing_reference_rows": len(study_obs.get(MISSING_REF, [])),
    }
    return out, stats


def main() -> None:
    obs = read_csv(FULL / "observations.csv")
    obs_by_pep = defaultdict(list)
    for r in obs:
        obs_by_pep[r["export_peptide"]].append(r)
    specs = {sid: (label, spec) for sid, label, spec in CANDIDATES}

    print("Per-study tables (exploratory; v2_labels_defined=false; no split chosen)")
    for sid in SHORTLIST:
        label, spec = specs[sid]
        rows, stats = per_study_rows(obs_by_pep, spec)
        path = FULL / f"per_study_{sid}.csv"
        write_csv(path, rows)
        n_lab = stats["n_labeled_peptides"]
        print(f"\n== {sid} ({label}) -> {path.relative_to(ROOT)}")
        print(
            f"studies: {stats['n_studies']}  eligible peptides: {stats['n_eligible_peptides']}  "
            f"globally pos-only/neg-only peptides: {n_lab}  "
            f"eligible rows missing reference_id: {stats['n_missing_reference_rows']}"
        )
        print("top 5 by exclusive peptides (ref | pmid | n_pep | n_excl | P/N/M within study | top allele):")
        for r in rows[:5]:
            print(
                f"  {r['reference_id']} | {r['pubmed_id']} | {r['n_peptides']} | {r['n_peptides_exclusive']} | "
                f"{r['n_positive_only']}/{r['n_negative_only']}/{r['n_mixed']} | "
                f"{r['top3_mhc_allele_name'].split(';')[0]}"
            )
        if rows and n_lab:
            top = rows[0]
            print(
                f"top study {top['reference_id']}: {top['n_peptides_globally_posneg_only']}/{n_lab} "
                f"({top['n_peptides_globally_posneg_only'] / n_lab:.1%}) of labeled peptides appear in it; "
                f"{top['n_exclusive_globally_posneg_only']}/{n_lab} "
                f"({top['n_exclusive_globally_posneg_only'] / n_lab:.1%}) are exclusive to it"
            )


if __name__ == "__main__":
    main()
