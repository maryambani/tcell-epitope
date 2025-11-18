import pandas as pd
from pathlib import Path

FIXED = Path("data/processed/fixed")

files = [
    "iedb_tcell_pos_fixed.csv",
    "iedb_tcell_neg_fixed.csv",
    "cedar_tcell_pos_fixed.csv",
    "cedar_tcell_neg_fixed.csv",
]

for fname in files:
    path = FIXED / fname
    print(f"\n=== {fname} ===")
    try:
        df = pd.read_csv(path, low_memory=False)
        print("Columns:")
        for c in df.columns:
            print("  ", c)
        print("\nFirst 2 rows:")
        print(df.head(2).T)  # transpose so you can see values under each col
    except Exception as e:
        print("  ERROR:", e)
