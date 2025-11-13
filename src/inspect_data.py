import pandas as pd
from pathlib import Path

data_dir = Path("data")

files = [
    "iedb_tcell_pos.csv",
    "iedb_tcell_neg.csv",
    "cedar_tcell_pos.csv",
    "cedar_tcell_neg.csv"
]

for fname in files:
    path = data_dir / fname
    print(f"\n=== {fname} ===")
    try:
        df = pd.read_csv(path)
        print("Columns:", list(df.columns))
        print(df.head(3))
    except Exception as e:
        print(f"Could not read {fname}: {e}")