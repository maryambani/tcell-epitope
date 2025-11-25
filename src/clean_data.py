import pandas as pd
from pathlib import Path

DATA_DIR = Path("data") / "processed"
INPUT_FILE = DATA_DIR / "combined_raw.csv"
OUTPUT_FILE = DATA_DIR / "combined_cleaned.csv"

# Only use the 20 standard amino acids
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if not {"peptide", "label", "source"}.issubset(df.columns):
        raise ValueError("Input file must contain 'peptide', 'label', 'source'")
    return df

def drop_invalid_chars(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    def is_valid(seq: str) -> bool:
        if not isinstance(seq, str):
            return False
        seq = seq.strip().upper()
        if not seq:
            return False
        return all(ch in STANDARD_AA for ch in seq)

    valid_mask = df["peptide"].apply(is_valid)
    removed_count = (~valid_mask).sum()
    cleaned = df[valid_mask].copy()
    return cleaned, removed_count

def drop_duplicates(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    before = len(df)
    deduped = df.drop_duplicates(subset="peptide", keep="first")
    removed_count = before - len(deduped)
    return deduped, removed_count

def main():
    print(f"[INFO] Loading {INPUT_FILE}")
    df = load_data(INPUT_FILE)
    print(f"[INFO] Starting rows: {len(df):,}")

    df, invalid_removed = drop_invalid_chars(df)
    print(f"[INFO] Removed invalid sequences: {invalid_removed:,}")

    df, dup_removed = drop_duplicates(df)
    print(f"[INFO] Removed duplicate sequences: {dup_removed:,}")

    df.to_csv(OUTPUT_FILE, index=False)
    print(f"[OK] Wrote cleaned data to {OUTPUT_FILE} with {len(df):,} rows")

if __name__ == "__main__":
    main()