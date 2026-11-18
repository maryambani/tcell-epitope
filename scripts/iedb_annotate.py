"""Shared annotation flags for IEDB T-cell observations. No allele normalization."""
from __future__ import annotations

import re

POSITIVE_PREFIX = "positive"
COARSE_ALLELE = {
    "human",
    "hla class i",
    "hla class ii",
    "mhc class i",
    "mhc class ii",
    "class i",
    "class ii",
    "hla-i",
    "hla-ii",
}

SEROTYPE_RE = re.compile(
    r"^HLA-[ABC](\d{1,2})(\:\d{2})?$",
    re.I,
)
# HLA-A2, HLA-A24 kept as serotype if no *allele. Do not map to four-digit alleles.
STAR_ALLELE_RE = re.compile(r"HLA-[A-Z]+\*\d+", re.I)
# Mouse gene-like names (H2-Kb, H2-Db). Haplotypes such as "H2-b class I" are not alleles.
H2_GENE_RE = re.compile(r"^H2-[KDLQA]", re.I)
MUTANT_RE = re.compile(r"mutant|substitution|\b[A-Z]\d+[A-Z]\b", re.I)

CONTEXT_KEYS = (
    "host_organism_name",
    "mhc_class",
    "mhc_allele_name",
    "assay_names",
    "reference_id",
    "direct_ex_vivo_bool",
    "immunization_description",
)

FUNCTIONAL_TOKENS = (
    "elispots",
    "elispot",
    "ics",
    "elisa",
    "cytotoxicity",
    "chromium",
    "proliferation",
    "degranulation",
    "granzyme",
    "activation",
    "ifng",
    "ifn-g",
    "ifnγ",
    "tnfa",
    "tnf-a",
    "il-2",
    "il-4",
    "il-5",
    "il-10",
    "il-13",
    "cytokine",
    "release",
)
BINDING_TOKENS = (
    "multimer",
    "tetramer",
    "qualitative binding",
    "surface plasmon",
    "spr",
    "dissociation constant",
    "kd",
)
STRUCTURAL_TOKENS = ("x-ray", "crystall", "3d structure")


def nonempty(v) -> bool:
    if v is None:
        return False
    if isinstance(v, bool):
        return True
    if isinstance(v, (int, float)) and v == v:
        return True
    return str(v).strip() not in ("", "nan", "None")


def recorded_value(v) -> str | None:
    """Presence-preserving stringify. 0/False are recorded values, not missing."""
    if v is None:
        return None
    if isinstance(v, bool):
        return "1" if v else "0"
    s = str(v).strip()
    if s.lower() in ("", "nan", "none"):
        return None
    if s.lower() == "false":
        return "0"
    if s.lower() == "true":
        return "1"
    return s


def outcome_family(measure: str) -> str:
    if not nonempty(measure):
        return "missing"
    m = str(measure).strip().lower()
    if m.startswith(POSITIVE_PREFIX):
        return "positive"
    if m == "negative":
        return "negative"
    return "other"


def endpoint_family(assay_names: str) -> str:
    if not nonempty(assay_names):
        return "other_or_unresolved"
    s = assay_names.lower().replace("γ", "g").replace("α", "a")
    if any(t in s for t in STRUCTURAL_TOKENS):
        return "structural"
    if any(t in s for t in BINDING_TOKENS):
        return "recognition_binding"
    if any(t in s for t in FUNCTIONAL_TOKENS):
        return "functional_response"
    return "other_or_unresolved"


def mhc_specificity(allele_name: str, mhc_class: str) -> str:
    if not nonempty(allele_name):
        return "missing" if not nonempty(mhc_class) else "class_only"
    a = allele_name.strip()
    al = a.lower()
    if al in COARSE_ALLELE:
        return "class_only"
    if STAR_ALLELE_RE.search(a) or H2_GENE_RE.search(a) or "*" in a:
        return "allele"
    # HLA-A2, HLA-A24, HLA-B7
    if re.match(r"^HLA-[ABC]\d", a, re.I) and "*" not in a:
        return "serotype"
    if re.match(r"^HLA-[ABC]$", a, re.I):
        return "class_only"
    return "other"


def mhc_modified(allele_name: str, restriction: str) -> int:
    text = f"{allele_name or ''} {restriction or ''}"
    return int(bool(MUTANT_RE.search(text)))


def evidence_bucket(evidence: str) -> str:
    if not nonempty(evidence):
        return "missing"
    e = evidence.strip().lower()
    if "not determined" in e:
        return "not_determined"
    if "cited reference" in e:
        return "cited_reference"
    if "binding assay" in e:
        return "binding_assay"
    if "prediction" in e:
        return "prediction"
    if "single allele" in e:
        return "single_allele_present"
    return "other"


def host_is_human(host: str) -> int:
    if not nonempty(host):
        return 0
    return int("homo sapiens" in host.lower())


def host_is_transgenic_mouse(host: str) -> int:
    if not nonempty(host):
        return 0
    h = host.lower()
    return int("mus musculus" in h and "tg" in h)


def outcome_pattern(rows: list[dict]) -> str:
    """Distinguish consistency from absence. Not a modeling label."""
    if not rows:
        return "no_eligible_evidence"
    fam = {r.get("outcome_family") for r in rows}
    has_p = "positive" in fam
    has_n = "negative" in fam
    if has_p and has_n:
        return "mixed"
    if has_p:
        return "positive_only"
    if has_n:
        return "negative_only"
    return "eligible_uninterpretable_outcome"


def complete_context_key(row: dict):
    parts = []
    for k in CONTEXT_KEYS:
        v = recorded_value(row.get(k))
        if v is None:
            return None
        parts.append(v)
    return tuple(parts)


def mixed_complete_context_groups(rows: list[dict]) -> set[tuple]:
    """Distinct complete recorded-context keys that contain both outcome families."""
    pos = {complete_context_key(r) for r in rows if r.get("outcome_family") == "positive"}
    neg = {complete_context_key(r) for r in rows if r.get("outcome_family") == "negative"}
    pos.discard(None)
    neg.discard(None)
    return pos & neg


def disagreement_class(rows: list[dict]) -> str:
    """Descriptive only. Same recorded context is not the same experiment."""
    pos = [r for r in rows if r.get("outcome_family") == "positive"]
    neg = [r for r in rows if r.get("outcome_family") == "negative"]
    if not pos or not neg:
        return "no_mixed_outcome"

    def ctx(r, require_complete: bool):
        parts = []
        for k in CONTEXT_KEYS:
            v = recorded_value(r.get(k))
            if require_complete and v is None:
                return None
            parts.append(v if v is not None else "")
        return tuple(parts)

    pos_full = {c for r in pos if (c := ctx(r, True))}
    neg_full = {c for r in neg if (c := ctx(r, True))}
    if pos_full & neg_full:
        return "within_recorded_context"
    pos_any = {ctx(r, False) for r in pos}
    neg_any = {ctx(r, False) for r in neg}
    if pos_any & neg_any:
        return "cannot_localize"
    return "across_known_context"
