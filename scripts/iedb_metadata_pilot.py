"""Pilot: retrieve current IEDB T-cell assay records for sampled epitope IDs.

Stores retrieval_date. This is NOT a reconstruction of the 2025-11-03 snapshot.
Does not overwrite V1 artifacts.
"""
from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audit"

API = "https://query-api.iedb.org/tcell_search"
SELECT = ",".join(
    [
        "tcell_id",
        "tcell_iri",
        "structure_id",
        "structure_iri",
        "linear_sequence",
        "qualitative_measure",
        "assay_names",
        "assay_description",
        "host_organism_name",
        "host_organism_iri",
        "mhc_class",
        "mhc_allele_name",
        "mhc_restriction",
        "mhc_allele_evidence",
        "reference_id",
        "reference_iri",
        "pubmed_id",
        "immunization_description",
        "direct_ex_vivo_bool",
        "parent_source_antigen_name",
        "source_organism_name",
    ]
)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
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


def fetch_assays(structure_id: str, timeout: int = 60) -> list[dict]:
    params = {
        "structure_id": f"eq.{structure_id}",
        "select": SELECT,
        "limit": "1000",
        "order": "tcell_id",
    }
    url = API + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    if not isinstance(data, list):
        raise RuntimeError(f"Unexpected response for {structure_id}: {data}")
    return data


def nonempty(v) -> bool:
    if v is None:
        return False
    if isinstance(v, str) and v.strip() == "":
        return False
    return True


def main() -> None:
    sample = read_csv(AUDIT / "metadata_pilot_sample.csv")
    retrieval_date = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    observations = []
    link_rows = []

    for i, pep_row in enumerate(sample, 1):
        ids = [x for x in (pep_row.get("epitope_ids") or "").split("|") if x]
        if not ids and pep_row.get("primary_epitope_id"):
            ids = [pep_row["primary_epitope_id"]]
        print(f"[{i}/{len(sample)}] {pep_row['peptide']} ids={ids} stratum={pep_row['pilot_stratum']}")
        assays_for_peptide = []
        errors = []
        for eid in ids:
            try:
                assays = fetch_assays(eid)
            except Exception as e:
                errors.append(f"{eid}:{e}")
                assays = []
            time.sleep(0.15)
            for a in assays:
                seq = a.get("linear_sequence") or ""
                observations.append(
                    {
                        "observation_id": a.get("tcell_iri") or f"IEDB_ASSAY:{a.get('tcell_id')}",
                        "tcell_id": a.get("tcell_id"),
                        "epitope_id": a.get("structure_id"),
                        "export_peptide": pep_row["peptide"],
                        "linear_sequence": seq,
                        "sequence_matches_export": int(seq == pep_row["peptide"]),
                        "qualitative_measure": a.get("qualitative_measure"),
                        "assay_names": a.get("assay_names"),
                        "host_organism_name": a.get("host_organism_name"),
                        "mhc_class": a.get("mhc_class"),
                        "mhc_allele_name": a.get("mhc_allele_name"),
                        "mhc_restriction": a.get("mhc_restriction"),
                        "mhc_allele_evidence": a.get("mhc_allele_evidence"),
                        "reference_id": a.get("reference_id"),
                        "pubmed_id": a.get("pubmed_id"),
                        "immunization_description": a.get("immunization_description"),
                        "direct_ex_vivo_bool": a.get("direct_ex_vivo_bool"),
                        "parent_source_antigen_name": a.get("parent_source_antigen_name"),
                        "source_organism_name": a.get("source_organism_name"),
                        "pilot_stratum": pep_row["pilot_stratum"],
                        "v1_conflict": pep_row["v1_conflict"],
                        "v1_query_label": pep_row["v1_query_label"],
                        "retrieval_date": retrieval_date,
                        "retrieval_source": "iedb_query_api.tcell_search",
                    }
                )
                assays_for_peptide.append(a)
        link_rows.append(
            {
                "peptide": pep_row["peptide"],
                "pilot_stratum": pep_row["pilot_stratum"],
                "v1_conflict": pep_row["v1_conflict"],
                "v1_query_label": pep_row["v1_query_label"],
                "in_iedb": pep_row["in_iedb"],
                "in_cedar": pep_row["in_cedar"],
                "n_epitope_ids_queried": len(ids),
                "n_assays": len(assays_for_peptide),
                "linked": int(len(assays_for_peptide) > 0),
                "errors": ";".join(errors),
            }
        )

    write_csv(AUDIT / "metadata_pilot_observations.csv", observations)
    write_csv(AUDIT / "metadata_pilot_link_status.csv", link_rows)

    # Coverage summary
    fields = [
        "qualitative_measure",
        "assay_names",
        "host_organism_name",
        "mhc_class",
        "mhc_allele_name",
        "reference_id",
        "pubmed_id",
        "immunization_description",
        "parent_source_antigen_name",
        "source_organism_name",
    ]
    field_cov = {}
    n_obs = len(observations)
    for fld in fields:
        n = sum(1 for r in observations if nonempty(r.get(fld)))
        field_cov[fld] = {"n_non_null": n, "frac": round(n / n_obs, 4) if n_obs else None}

    def by_stratum(rows, key="pilot_stratum"):
        g = defaultdict(list)
        for r in rows:
            g[r[key]].append(r)
        return g

    strata = {}
    for tag, rows in by_stratum(link_rows).items():
        assays = [r for r in observations if r["pilot_stratum"] == tag]
        hosts = Counter(r["host_organism_name"] or "MISSING" for r in assays)
        classes = Counter(r["mhc_class"] or "MISSING" for r in assays)
        outcomes = Counter(r["qualitative_measure"] or "MISSING" for r in assays)
        n_studies = len({r["reference_id"] for r in assays if nonempty(r.get("reference_id"))})
        strata[tag] = {
            "n_peptides_sampled": len(rows),
            "n_peptides_with_assays": sum(int(r["linked"]) for r in rows),
            "n_assays": len(assays),
            "n_studies": n_studies,
            "host_top": hosts.most_common(8),
            "mhc_class": dict(classes),
            "qualitative_measure": dict(outcomes),
        }

    summary = {
        "retrieval_date": retrieval_date,
        "snapshot_warning": (
            "These assay records reflect IEDB at retrieval_date, not the 2025-11-03 export. "
            "Keep original query membership separate from newly retrieved outcomes."
        ),
        "n_peptides_sampled": len(sample),
        "n_peptides_linked": sum(int(r["linked"]) for r in link_rows),
        "n_observations": n_obs,
        "n_distinct_studies": len({r["reference_id"] for r in observations if nonempty(r.get("reference_id"))}),
        "field_coverage_among_observations": field_cov,
        "by_stratum": strata,
        "wording_for_exports": (
            "The query-level restriction is documented (Linear peptides, T Cell Assay, "
            "Outcome Positive/Negative), but assay subtypes and experimental contexts "
            "cannot be verified from the retained epitope-level exports."
        ),
    }
    (AUDIT / "metadata_pilot_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in summary if k != "by_stratum"}, indent=2))
    print("Wrote observations", n_obs, "linked", summary["n_peptides_linked"], "/", len(sample))


if __name__ == "__main__":
    main()
