import pandas as pd
from pathlib import Path

DATA = Path("data")
OUT = DATA / "processed"
FIXED = OUT / "fixed"
NORMAL = OUT / "normalized"
OUT.mkdir(parents=True, exist_ok=True)
FIXED.mkdir(parents=True, exist_ok=True)
NORMAL.mkdir(parents=True, exist_ok=True)

FILES = [
    ("iedb_tcell_pos_fixed.csv",  "IEDB",  1),
    ("iedb_tcell_neg_fixed.csv",  "IEDB",  0),
    ("cedar_tcell_pos_fixed.csv", "CEDAR", 1),
    ("cedar_tcell_neg_fixed.csv", "CEDAR", 0),
]

PEPTIDE_COL_OVERRIDES = {
    "iedb_tcell_pos_fixed.csv":  "Epitope.1|Name",
    "iedb_tcell_neg_fixed.csv":  "Epitope.1|Name",
    "cedar_tcell_pos_fixed.csv": "Epitope.1|Name",
    "cedar_tcell_neg_fixed.csv": "Epitope.1|Name",
}

def pick_peptide_col(cols):
    cols_lower = [c.lower() for c in cols]
    # Best: Epitope|Linear sequence
    for i,c in enumerate(cols_lower):
        if c == "epitope|linear sequence":
            return i
    # Next: any Epitope|*sequence*
    for i,c in enumerate(cols_lower):
        if c.startswith("epitope|") and "sequence" in c:
            return i
    # Fallback: any *sequence*
    for i,c in enumerate(cols_lower):
        if "sequence" in c:
            return i
    return None

def main():
    parts = []
    for fname, source, label in FILES:
        path = FIXED / fname
        if not path.exists():
            print(f"[WARN] Missing {fname}, skipping.")
            continue

        df = pd.read_csv(path, low_memory=False)
        
        # Overriding first
        override = PEPTIDE_COL_OVERRIDES.get(fname)
        if override and override in df.columns:
            pep_col = override
            print(f"[INFO] Using override peptide column '{pep_col}' for {fname}")
        else:
            # Falling back to heuristic if no override
            idx = pick_peptide_col(df.columns)
            if idx is None:
                print(f"[WARN] No peptide column found in {fname}")
                continue
            pep_col = df.columns[idx]

        out = df[[pep_col]].copy()
        out.columns = ["peptide"]
        out["source"] = source
        out["label"] = int(label)
        out_path = NORMAL / f"{fname.replace('_fixed.csv','')}_norm.csv"
        out.to_csv(out_path, index=False)
        parts.append(out)
        print(f"[OK] {fname}: using '{pep_col}' → {len(out)} rows")

    if parts:
        combined = pd.concat(parts, ignore_index=True)
        combined.to_csv(OUT / "combined_raw.csv", index=False)
        print(f"[OK] Wrote {OUT/'combined_raw.csv'} with {len(combined)} rows")

if __name__ == "__main__":
    main()