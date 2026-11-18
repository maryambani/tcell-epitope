# Benchmark V1 (frozen)

Do **not** overwrite these artifacts when building later benchmarks.

## Definition

> Peptide-level classification derived from positive/negative T-cell query membership, excluding peptides present in both lists, evaluated using a stratified random holdout.

This is **not** a claim of a single experimental endpoint. It is query-membership after dropping disagreements.

## Snapshot

| Item | Value |
|------|--------|
| Original export date | 2025-11-03 |
| Query filters | Linear peptides, T Cell Assay, Outcome = Positive/Negative |
| Sources | IEDB + CEDAR epitope-level lists |
| Label rule | Filename membership (`*_pos_*` → 1, `*_neg_*` → 0) |
| Conflict rule | Drop peptide if it appears with both labels |
| Dedup | One row per peptide, `keep="first"` (IEDB listed first) |
| Length | 8–11 residues after cleaning |
| Encoding | Right-pad to 11, 21-way one-hot (231 features) |
| Split | `train_test_split(test_size=0.2, random_state=42, stratify=y)` |
| Evaluation status | Exploratory (test set already inspected) |

## Frozen files (do not replace)

| Path | Role |
|------|------|
| `data/iedb_tcell_*.csv`, `data/cedar_tcell_*.csv` | Original 2025-11-03 exports |
| `data/processed/combined_raw.csv` | Peptide/source/label rows before cleaning |
| `data/processed/combined_cleaned.csv` | After invalid-AA + conflict drop + dedupe |
| `data/processed/features.npz` | V1 train/test arrays |
| `outputs/models/` | V1 LR / RF / CNN |
| `audit/split_assignments.csv` | Reconstructed V1 peptide partitions |
| `audit/test_predictions.csv` | V1 per-peptide test scores |
| `audit/test_metrics_from_scores.csv` | V1 headline metrics from scores |

## Implementation points already resolved

- Class weights were converted to a dict for Keras from **training** labels.
- CNN early stopping used `validation_split=0.1` on the **training** partition.
- Average precision used **continuous** scores.
- Exact peptide overlap train/test = **0**.
- LR/RF had no hyperparameter search (limited baseline, not evidence of leakage).

## What V1 does not measure

Transfer across related sequences, studies, hosts, alleles, or assay subtypes. Those questions need recovered assay metadata and a later versioned modeling table.
