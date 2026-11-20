## T-Cell Epitope Immunogenicity Classifier

Personal learning / research-portfolio project: **benchmarking peptide-only prediction of curated T-cell assay outcomes**, and investigating what sequence-only models learn under a simple encoding.

### Prediction target (operational definition)

This project predicts **curated positive versus negative T-cell assay outcomes from peptide sequence**.

A database “negative” means **no detected response under particular experimental conditions**. It is **not** proof that the peptide can never elicit a response in another host, HLA context, or assay. Positive labels come from MHC-dependent experimental contexts, but omitting HLA does **not** control for HLA-associated bias or isolate TCR-recognition signal.

You can still call this *peptide immunogenicity prediction*, with that definition in mind.

### Biological background

T cells recognize infected or abnormal cells when their receptors bind a small peptide presented on the cell surface by MHC molecules. Peptide–MHC recognition is an early step toward an immune response, but it is not sufficient on its own. [1]

Predicting which peptides are associated with T-cell responses matters for vaccine design, autoimmunity research, drug safety, and neoantigen discovery. [2]

### Scope

| Item | Choice |
|------|--------|
| Task | Supervised binary classification |
| Input | Peptide sequence only (no explicit MHC allele) |
| Length filter | 8–11 amino acids (typical MHC class I range) [3] |
| Models | Logistic regression, random forest, small 1D CNN |
| Current evaluation | Stratified random holdout (exploratory) |

**What this is:** an end-to-end prototype plus analysis of ranking signal, class imbalance, positional encoding, and planned biological ablations.

**What this is not yet:** a claim of a generalizable immunogenicity predictor across related peptides, studies, alleles, or experimental contexts.

### Suggested portfolio framing

> Benchmarking peptide-only prediction of curated T-cell assay outcomes: investigating generalization, dataset bias, and positional sequence signal.

---

### Environment setup

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Use the project `venv` kernel for notebooks (Python 3.11).

### Data sources

- [IEDB](https://www.iedb.org/): T-cell epitope exports (positive and negative files)
- [CEDAR](https://cedar.iedb.org): cancer epitope exports (positive and negative files)

Exports were prepared as T-cell assay–oriented files (`*_tcell_pos_*` / `*_tcell_neg_*`). The modeling table currently retains only `peptide`, `source`, and `label`. Host species, HLA allele, and assay subtype are **not** carried into training features yet (documented limitation).

### Data preparation

1. `src/fix_headers.py`: merge two-row IEDB/CEDAR headers  
2. `src/prepare_data.py`: combine four exports into `data/processed/combined_raw.csv`  
3. `src/clean_data.py`: drop invalid amino-acid characters; **drop peptides with conflicting 0/1 labels**; deduplicate

**Attrition (from EDA / cleaning, notebook 01 + `clean_data.py`):**

| Step | Count |
|------|------:|
| Raw combined rows | 308,790 |
| Invalid / non-standard sequences removed | 2,359 |
| Unique peptides in both pos and neg lists (any AA) | 26,504 |
| Valid-AA peptides dropped as conflicts | 26,204 |
| Rows removed at conflict stage | 63,301 |
| Remaining duplicates collapsed | (see cleaning log) |
| Cleaned unique peptides | ~199,860 |
| Class ratio (neg:pos) after cleaning | ~2.49:1 |

**Conflict policy:** if the same peptide appears with both labels, it is dropped rather than majority-voted. That makes a cleaner binary benchmark, but the same sequence can genuinely differ across experimental contexts. Assay-level provenance should be preserved in future work even if the modeling table stays one row per peptide.

**Length filter (notebook 02):** keep 8–11mers only (~47.6k sequences after filtering), right-pad to length 11, one-hot encode to 231 features (11 × 21), stratified 80/20 split.

### Modeling

All models use `class_weight='balanced'` (sklearn string; for Keras this is converted to a class→weight dict from **training** labels only).

| Model | Notes |
|-------|--------|
| Logistic regression | Linear baseline; interpretable coefficients |
| Random forest | 100 trees; captures nonlinear interactions; can overfit |
| 1D CNN | Input shaped `(11, 21)`; Conv1D → GlobalMaxPooling → Dense |

### Evaluation (notebook 04): exploratory holdout

On the stratified random test set (n = 9,521; positive rate ≈ 0.297):

| Model | Accuracy | ROC-AUC | Average precision |
|-------|----------|---------|-------------------|
| Logistic Regression | 0.596 | 0.634 | 0.423 |
| Random Forest | 0.753 | **0.702** | **0.571** |
| CNN | 0.577 | 0.649 | 0.451 |

At threshold 0.5, RF positive-class **recall ≈ 0.24** (misses most labeled positives) despite high positive precision. CNN has higher positive recall/F1 at that threshold.

**How to read these numbers**

> On a stratified random holdout, random forest achieved the highest observed ROC-AUC (0.702) and average precision (0.571). RF AP is substantially above the prevalence baseline (0.297), so there is within-dataset ranking signal. These results do **not** establish generalization across related peptides, studies, or experimental contexts.

The current test set has already been inspected (metrics, plots, threshold sweep). Treat results as **exploratory**. Do not treat a reshuffled random split as a fresh unconsulted test. Prefer development folds for future model / threshold / encoding choices.

### V2 observation recovery (in progress, no new model)

V1 is frozen. The V2 workstream recovers full assay-level provenance before choosing a new prediction task.

**Completed:**
- Staged IEDB T-cell observation retrieval for all 54,448 valid 8|11-mer peptides in `audit/epitope_pool_8_11_preconflict.csv`: 157,300 observations, 0 unresolved IDs, merge validation passed (`audit/retrieval/full/merge_validation.json`)
- Two reporting corrections: `WICLLQFAY` and `VQSAKYWNL` have no eligible human evidence under functional subsets; mixed records are tetramer/binding and cannot supply functional labels
- Full task-feasibility inventory across predeclared endpoint intersections (`audit/retrieval/full/task_feasibility.json`)
- Per-study tables and split diagnostics for ELISPOT and functional candidates (`per_study_*.csv`, `split_diagnostics.json`)

**Not yet decided:** primary prediction claim, label policy for mixed observations, evaluation split design, V2 model.

The V1 conflict-drop policy (drop any peptide appearing in both positive and negative exports) is documented in V1 but should not be applied as the V2 label rule. Labels should be derived from eligible observations for a named task.

---

### Known limitations (priority for next work)

1. **Label definition & provenance.** *In progress (V2 workstream).* Full-pool IEDB T-cell observation recovery is complete (54,448 peptides, 157,300 assays). Staged retrieval with per-query completeness checks, raw JSON, and append-only history is documented in `audit/RETRIEVAL_SPEC.md`. Feasibility inventory across candidate endpoint intersections is in `audit/retrieval/full/task_feasibility.json`. Primary prediction claim (host, MHC class, endpoint, mixed-label rule) is **not yet chosen**.
2. **Evaluation design.** *In progress (V2 workstream).* Study concentration analysis and sequence-separation diagnostics are done (`audit/retrieval/full/split_diagnostics.json`, `per_study_*.csv`). One all-negative vaccinia study (reference 1001865) accounts for ~39% of ELISPOT labeled peptides; a predeclared stress test holding it out is required before committing to a split. Hamming ≤1 sequence groups are computed; largest group is 20 peptides. No split has been adopted.
3. **Anchor ablation (notebook 06).** *Partially done.* Ablation ran on the 9-mer subset only (20,690 train / 5,166 test), avoiding the right-padding ambiguity by fixing PΩ at index 8. Masking P2 and PΩ costs 0.024 to 0.028 AUC across all three models (LR 0.637 to 0.612, RF 0.673 to 0.644, CNN 0.637 to 0.613). The drop is consistent but modest, so most predictive signal survives anchor removal. **Still missing:** control masks (equal-count non-anchor positions) and simple baselines to establish whether the drop exceeds what masking any two positions would cost. Interpret as information removal, not proof of MHC mechanism.  
4. **Encoding:** compare right-padding vs terminally aligned padding; report metrics by length.  
5. **RF regularization / weighting / calibration:** tune on validation folds; treat impurity importances as training diagnostics, not biology.

A larger neural net / transformer is **not** the priority.

### Project structure

```
tcell-epitope/
├── data/processed/          # combined_raw, combined_cleaned, features.npz
├── src/                     # fix_headers, prepare_data, clean_data (V1)
├── scripts/                 # iedb_annotate, iedb_pool_retrieve, task_feasibility,
│                            #   validate_merged_retrieval, per_study_table,
│                            #   split_diagnostics, similarity_audit, …
├── audit/                   # staged retrieval outputs, feasibility inventory,
│                            #   epitope pool, RETRIEVAL_SPEC.md
├── notebooks/               # 01 EDA → 06 anchor ablation
├── outputs/                 # figures (+ models/, regeneratable)
└── requirements.txt
```

### References

[1] Cano RLE, Lopera HDE. Introduction to T and B lymphocytes. In: Anaya JM et al., eds. *Autoimmunity: From Bench to Bedside*. El Rosario University Press; 2013. https://www.ncbi.nlm.nih.gov/books/NBK459471/

[2] Xu Z et al. Applying artificial intelligence for cancer immunotherapy. *Acta Pharm Sin B.* 2021;11(11):3393-3405. https://doi.org/10.1016/j.apsb.2021.02.007

[3] Alberts B et al. *Molecular Biology of the Cell.* 4th ed. Garland Science; 2002. T Cells and MHC Proteins. https://www.ncbi.nlm.nih.gov/books/NBK26926/

**Methods literature informing next steps (not direct AUC comparators):** Calis et al. 2013 (PLOS Comp Biol); PRIME / PRIME2.0; TRAP (Genome Medicine 2023); leakage-aware IEDB evaluation work; MUNIS (Nat Mach Intell 2025). Endpoints and splits differ across papers, so do not compare headline AUCs directly.
