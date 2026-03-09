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


def resolve_conflicts(df: pd.DataFrame) -> tuple[pd.DataFrame, int, int]:
    """Handle peptides that appear with both label 0 and label 1.

    - If all entries for a peptide agree, keep one copy with that label.
    - If entries disagree (conflicting labels), drop the peptide entirely.
    """
    label_counts = df.groupby("peptide")["label"].nunique()
    conflicting_peptides = set(label_counts[label_counts > 1].index)
    num_conflicts = len(conflicting_peptides)

    # Remove conflicting peptides
    df_clean = df[~df["peptide"].isin(conflicting_peptides)].copy()

    # Deduplicate the remaining (non-conflicting) peptides
    before = len(df_clean)
    df_clean = df_clean.drop_duplicates(subset="peptide", keep="first")
    num_deduped = before - len(df_clean)

    return df_clean, num_conflicts, num_deduped

def main():
    print(f"[INFO] Loading {INPUT_FILE}")
    df = load_data(INPUT_FILE)
    print(f"[INFO] Starting rows: {len(df):,}")

    df, invalid_removed = drop_invalid_chars(df)
    print(f"[INFO] Removed invalid sequences: {invalid_removed:,}")

    df, num_conflicts, num_deduped = resolve_conflicts(df)
    print(f"[INFO] Dropped conflicting peptides: {num_conflicts:,}")
    print(f"[INFO] Removed remaining duplicates: {num_deduped:,}")

    df.to_csv(OUTPUT_FILE, index=False)
    print(f"[OK] Wrote cleaned data to {OUTPUT_FILE} with {len(df):,} rows")

if __name__ == "__main__":
    main()