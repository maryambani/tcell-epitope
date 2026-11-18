"""Independent validation of the four metadata-pilot files. Does not modify V1."""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

AUDIT = Path(__file__).resolve().parents[1] / "audit"


def read_csv(name: str) -> list[dict]:
    with (AUDIT / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def nonempty(v) -> bool:
    return v is not None and str(v).strip() not in ("", "nan", "None")


def is_specific_host(h: str) -> bool:
    if not nonempty(h):
        return False
    h = h.strip()
    vague = {"mammal", "vertebrate", "unknown"}
    return h.lower() not in vague


def mhc_specific(v: str) -> bool:
    if not nonempty(v):
        return False
    v = v.strip()
    # too coarse to use as restriction
    coarse = {"human", "HLA class I", "HLA class II", "MHC class I", "MHC class II", "class I", "class II"}
    return v not in coarse


def main() -> dict:
    sizes = json.loads((AUDIT / "metadata_pilot_sample_sizes.json").read_text())
    sample = read_csv("metadata_pilot_sample.csv")
    links = read_csv("metadata_pilot_link_status.csv")
    obs = read_csv("metadata_pilot_observations.csv")

    out: dict = {}
    out["n_sample"] = len(sample)
    out["n_unique_peptides_sample"] = len({r["peptide"] for r in sample})
    out["n_links"] = len(links)
    out["n_obs"] = len(obs)
    out["sample_cols"] = list(sample[0].keys()) if sample else []
    out["link_cols"] = list(links[0].keys()) if links else []
    out["obs_cols"] = list(obs[0].keys()) if obs else []

    # stratum sizes
    sample_strata = Counter(r["pilot_stratum"] for r in sample)
    size_sampled = {k: v.get("n_sampled") for k, v in sizes.items()}
    out["stratum_sample_counts"] = dict(sample_strata)
    out["stratum_json_n_sampled"] = size_sampled
    out["stratum_match_json"] = dict(sample_strata) == {k: v for k, v in size_sampled.items() if v}

    # sample vs link peptide identity
    s_peps = [r["peptide"] for r in sample]
    l_peps = [r["peptide"] for r in links]
    out["sample_link_peptides_equal_order"] = s_peps == l_peps
    out["duplicate_peptides_in_sample"] = [p for p, n in Counter(s_peps).items() if n > 1]

    # epitope ids
    n_multi_id = 0
    id_mismatch_primary = 0
    for r in sample:
        ids = [x for x in (r.get("epitope_ids") or "").split("|") if x]
        if len(ids) != 1:
            n_multi_id += 1
        if r.get("primary_epitope_id") and ids and r["primary_epitope_id"] != ids[0]:
            id_mismatch_primary += 1
    out["n_peptides_with_multiple_epitope_ids"] = n_multi_id
    out["n_primary_id_mismatch"] = id_mismatch_primary

    # stratum flag consistency
    flag_issues = []
    for r in sample:
        st = r["pilot_stratum"]
        if st == "retained_pos_iedb_only":
            ok = r["v1_retained"] == "1" and r["v1_query_label"] == "1" and r["in_cedar"] == "0" and r["in_iedb"] == "1"
        elif st == "retained_neg_iedb_only":
            ok = r["v1_retained"] == "1" and r["v1_query_label"] == "0" and r["in_cedar"] == "0"
        elif st == "retained_both_databases":
            ok = r["v1_retained"] == "1" and r["in_iedb"] == "1" and r["in_cedar"] == "1"
        elif st == "conflict_iedb_only":
            ok = r["v1_conflict"] == "1" and r["in_cedar"] == "0"
        elif st == "conflict_both_databases":
            ok = r["v1_conflict"] == "1" and r["in_iedb"] == "1" and r["in_cedar"] == "1"
        else:
            ok = True
        if not ok:
            flag_issues.append(r["peptide"] + ":" + st)
    out["n_stratum_flag_inconsistencies"] = len(flag_issues)
    out["stratum_flag_issues"] = flag_issues[:10]

    # observations
    out["n_unique_observation_id"] = len({r["observation_id"] for r in obs})
    out["n_unique_tcell_id"] = len({r["tcell_id"] for r in obs})
    out["n_dup_observation_id"] = out["n_obs"] - out["n_unique_observation_id"]
    tcell_counts = Counter(r["tcell_id"] for r in obs)
    out["n_dup_tcell_id_values"] = sum(1 for n in tcell_counts.values() if n > 1)
    out["max_tcell_id_repeats"] = max(tcell_counts.values()) if tcell_counts else 0

    # reconcile n_assays
    link_n = {r["peptide"]: int(r["n_assays"]) for r in links}
    obs_n = Counter(r["export_peptide"] for r in obs)
    mismatches = []
    for pep, n in link_n.items():
        if obs_n[pep] != n:
            mismatches.append({"peptide": pep, "link_n_assays": n, "obs_rows": obs_n[pep]})
    out["n_link_obs_count_mismatches"] = len(mismatches)
    out["link_obs_mismatches"] = mismatches[:15]
    out["sum_link_n_assays"] = sum(link_n.values())
    out["n_obs_rows"] = len(obs)
    out["n_peptides_linked"] = sum(int(r["linked"]) for r in links)
    out["n_peptides_with_errors"] = sum(1 for r in links if nonempty(r.get("errors")))

    # sequence agreement
    seq_mismatch = sum(1 for r in obs if r.get("sequence_matches_export") != "1")
    true_mismatch = [
        {"peptide": r["export_peptide"], "linear_sequence": r["linear_sequence"], "epitope_id": r["epitope_id"]}
        for r in obs
        if (r.get("linear_sequence") or "") != (r.get("export_peptide") or "")
    ]
    out["n_obs_flag_sequence_mismatch"] = seq_mismatch
    out["n_obs_actual_sequence_mismatch"] = len(true_mismatch)
    out["sequence_mismatch_examples"] = true_mismatch[:5]

    # peptides in obs vs sample
    out["peptides_in_obs_not_sample"] = sorted(set(obs_n) - set(s_peps))
    out["peptides_in_sample_not_obs"] = sorted(set(s_peps) - set(obs_n))

    # epitope_id vs sample ids
    sample_ids = {}
    for r in sample:
        sample_ids[r["peptide"]] = set(x for x in (r.get("epitope_ids") or "").split("|") if x)
    id_not_in_sample = 0
    for r in obs:
        if str(r.get("epitope_id") or "") not in sample_ids.get(r["export_peptide"], set()):
            id_not_in_sample += 1
    out["n_obs_epitope_id_not_in_sample_ids"] = id_not_in_sample

    # retrieval provenance
    out["retrieval_dates"] = sorted({r.get("retrieval_date") for r in obs})
    out["retrieval_sources"] = sorted({r.get("retrieval_source") for r in obs})

    # field coverage assay-level
    fields = [
        "qualitative_measure",
        "assay_names",
        "host_organism_name",
        "mhc_class",
        "mhc_allele_name",
        "mhc_restriction",
        "mhc_allele_evidence",
        "reference_id",
        "pubmed_id",
        "immunization_description",
        "direct_ex_vivo_bool",
        "parent_source_antigen_name",
        "source_organism_name",
    ]
    cov = {}
    for fld in fields:
        n = sum(1 for r in obs if nonempty(r.get(fld)))
        cov[fld] = {"n": n, "frac": round(n / len(obs), 4)}
    out["assay_field_nonempty"] = cov

    # peptide-level recoverability: at least one nonempty assay field
    by_pep = defaultdict(list)
    for r in obs:
        by_pep[r["export_peptide"]].append(r)

    def pep_any(fld):
        return sum(1 for pep, rs in by_pep.items() if any(nonempty(x.get(fld)) for x in rs))

    out["peptide_level_any_nonempty"] = {fld: pep_any(fld) for fld in fields}
    out["n_peptides_with_obs"] = len(by_pep)

    # interpretable specificity
    hosts = Counter(r["host_organism_name"] for r in obs if nonempty(r.get("host_organism_name")))
    out["n_distinct_hosts"] = len(hosts)
    out["host_top15"] = hosts.most_common(15)
    humanish = sum(n for h, n in hosts.items() if "Homo sapiens" in h or h.startswith("human"))
    mouseish = sum(n for h, n in hosts.items() if "Mus musculus" in h)
    out["n_assays_human_host"] = humanish
    out["n_assays_mouse_host"] = mouseish

    classes = Counter(r["mhc_class"] if nonempty(r.get("mhc_class")) else "MISSING" for r in obs)
    out["mhc_class_counts"] = dict(classes)
    alleles = Counter(r["mhc_allele_name"] if nonempty(r.get("mhc_allele_name")) else "MISSING" for r in obs)
    out["n_distinct_allele_name"] = len(alleles)
    out["n_allele_coarse"] = sum(1 for r in obs if nonempty(r.get("mhc_allele_name")) and not mhc_specific(r["mhc_allele_name"]))
    out["n_allele_specific"] = sum(1 for r in obs if mhc_specific(r.get("mhc_allele_name") or ""))
    out["allele_top15"] = alleles.most_common(15)
    evidence = Counter(r["mhc_allele_evidence"] if nonempty(r.get("mhc_allele_evidence")) else "MISSING" for r in obs)
    out["mhc_allele_evidence"] = dict(evidence)

    outcomes = Counter(r["qualitative_measure"] if nonempty(r.get("qualitative_measure")) else "MISSING" for r in obs)
    out["qualitative_measure"] = dict(outcomes)

    assays = Counter(r["assay_names"] if nonempty(r.get("assay_names")) else "MISSING" for r in obs)
    out["n_distinct_assay_names"] = len(assays)
    out["assay_names_top15"] = assays.most_common(15)

    # stimulation present vs informative
    stim_empty = sum(1 for r in obs if not nonempty(r.get("immunization_description")))
    out["n_immunization_empty"] = stim_empty

    refs = Counter(r["reference_id"] for r in obs if nonempty(r.get("reference_id")))
    out["n_distinct_reference_id"] = len(refs)
    out["max_assays_one_reference"] = max(refs.values()) if refs else 0

    # concentration by peptide
    pep_n = Counter(r["export_peptide"] for r in obs)
    out["assays_per_peptide_min"] = min(pep_n.values())
    out["assays_per_peptide_median"] = sorted(pep_n.values())[len(pep_n) // 2]
    out["assays_per_peptide_max"] = max(pep_n.values())
    top_pep, top_n = pep_n.most_common(1)[0]
    out["top_peptide"] = top_pep
    out["top_peptide_n_assays"] = top_n
    out["frac_obs_from_top_peptide"] = round(top_n / len(obs), 4)
    out["top5_peptides"] = pep_n.most_common(5)

    # stratum assay counts with/without top peptide
    stratum_n = Counter(r["pilot_stratum"] for r in obs)
    out["assays_by_stratum"] = dict(stratum_n)
    obs_wo_top = [r for r in obs if r["export_peptide"] != top_pep]
    out["assays_by_stratum_excluding_top_peptide"] = dict(Counter(r["pilot_stratum"] for r in obs_wo_top))
    out["n_obs_excluding_top"] = len(obs_wo_top)
    out["n_studies_excluding_top"] = len({r["reference_id"] for r in obs_wo_top if nonempty(r.get("reference_id"))})

    # peptide-level study counts
    pep_refs = defaultdict(set)
    for r in obs:
        if nonempty(r.get("reference_id")):
            pep_refs[r["export_peptide"]].add(r["reference_id"])
    nrefs = [len(pep_refs[p]) for p in by_pep]
    out["studies_per_peptide_min"] = min(nrefs) if nrefs else None
    out["studies_per_peptide_median"] = sorted(nrefs)[len(nrefs) // 2] if nrefs else None
    out["studies_per_peptide_max"] = max(nrefs) if nrefs else None

    # V1 label vs current assay outcomes
    def is_pos(m: str) -> bool:
        return nonempty(m) and str(m).lower().startswith("positive")

    def is_neg(m: str) -> bool:
        return str(m).strip() == "Negative"

    label_vs = []
    for pep, rs in by_pep.items():
        v1 = rs[0]["v1_query_label"]
        conflict = rs[0]["v1_conflict"]
        stratum = rs[0]["pilot_stratum"]
        n_pos = sum(1 for x in rs if is_pos(x.get("qualitative_measure") or ""))
        n_neg = sum(1 for x in rs if is_neg(x.get("qualitative_measure") or ""))
        other = len(rs) - n_pos - n_neg
        label_vs.append(
            {
                "peptide": pep,
                "stratum": stratum,
                "v1_conflict": conflict,
                "v1_query_label": v1,
                "n_pos": n_pos,
                "n_neg": n_neg,
                "n_other": other,
                "n_assays": len(rs),
            }
        )
    out["n_retained_pos_with_any_current_neg"] = sum(
        1 for r in label_vs if r["v1_conflict"] == "0" and r["v1_query_label"] == "1" and r["n_neg"] > 0
    )
    out["n_retained_neg_with_any_current_pos"] = sum(
        1 for r in label_vs if r["v1_conflict"] == "0" and r["v1_query_label"] == "0" and r["n_pos"] > 0
    )
    out["n_conflict_peptides_with_both_pos_neg_now"] = sum(
        1 for r in label_vs if r["v1_conflict"] == "1" and r["n_pos"] > 0 and r["n_neg"] > 0
    )
    out["n_conflict_peptides_only_pos_now"] = sum(
        1 for r in label_vs if r["v1_conflict"] == "1" and r["n_pos"] > 0 and r["n_neg"] == 0
    )
    out["n_conflict_peptides_only_neg_now"] = sum(
        1 for r in label_vs if r["v1_conflict"] == "1" and r["n_neg"] > 0 and r["n_pos"] == 0
    )

    # disagreement context: same host+mhc_allele+assay_names+reference vs across
    def ctx(r):
        return (
            r.get("host_organism_name") or "",
            r.get("mhc_allele_name") or "",
            r.get("assay_names") or "",
            r.get("reference_id") or "",
        )

    within_ctx = 0
    across_ctx = 0
    mixed_unknown = 0
    for pep, rs in by_pep.items():
        pos_ctx = {ctx(x) for x in rs if is_pos(x.get("qualitative_measure") or "")}
        neg_ctx = {ctx(x) for x in rs if is_neg(x.get("qualitative_measure") or "")}
        if not pos_ctx or not neg_ctx:
            continue
        overlap = pos_ctx & neg_ctx
        if overlap:
            # if any overlapping context has empty key fields, bound as unknown
            if any(not all(nonempty(part) for part in c) for c in overlap):
                mixed_unknown += 1
            else:
                within_ctx += 1
        else:
            across_ctx += 1
    out["n_peptides_pos_and_neg_assays"] = sum(1 for pep, rs in by_pep.items() if any(is_pos(x.get("qualitative_measure") or "") for x in rs) and any(is_neg(x.get("qualitative_measure") or "") for x in rs))
    out["n_disagreement_within_observed_context"] = within_ctx
    out["n_disagreement_across_context"] = across_ctx
    out["n_disagreement_overlap_but_incomplete_context"] = mixed_unknown

    # pagination risk: any peptide with n_assays == 1000
    out["n_peptides_at_limit_1000"] = sum(1 for n in pep_n.values() if n >= 1000)
    out["n_peptides_at_limit_500"] = sum(1 for n in pep_n.values() if n >= 500)

    (AUDIT / "pilot_file_validation.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2)[:8000])
    print("\n... wrote audit/pilot_file_validation.json")
    return out


if __name__ == "__main__":
    main()
