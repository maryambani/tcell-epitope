# src/prepare_data.py  (Milestone A)
import pandas as pd
from pathlib import Path

DATA = Path("data")
OUT = DATA / "processed"
FIXED = OUT / "fixed"
OUT.mkdir(parents=True, exist_ok=True)
FIXED.mkdir(parents=True, exist_ok=True)

FILES = [
    ("iedb_tcell_pos.csv",  "IEDB",  1),
    ("iedb_tcell_neg.csv",  "IEDB",  0),
    ("cedar_tcell_pos.csv", "CEDAR", 1),
    ("cedar_tcell_neg.csv", "CEDAR", 0),
]

def fix_two_row_header(df: pd.DataFrame) -> pd.DataFrame:
    """Combine top header and first row into 'Group|Subheader' column names, drop first row."""
    if df.empty:
        return df
    sub = df.iloc[0].astype(str).fillna("")
    new_cols = []
    for col in df.columns:
        group = str(col).strip()
        subname = sub.get(col, "").strip()
        new_cols.append(f"{group}|{subname}" if subname else group)
    df.columns = new_cols
    return df.iloc[1:].reset_index(drop=True)

def main():
    for fname, source, label in FILES:
        path = DATA / fname
        if not path.exists():
            print(f"[WARN] Missing {fname}, skipping.")
            continue
        print(f"[INFO] Loading {fname} …")
        df = pd.read_csv(path, low_memory=False)
        df_fixed = fix_two_row_header(df)
        out_path = FIXED / f"{fname.replace('.csv','')}_fixed.csv"
        df_fixed.to_csv(out_path, index=False)
        print(f"[OK] Wrote {out_path} with {len(df_fixed)} rows and {len(df_fixed.columns)} cols")

if __name__ == "__main__":
    main()
