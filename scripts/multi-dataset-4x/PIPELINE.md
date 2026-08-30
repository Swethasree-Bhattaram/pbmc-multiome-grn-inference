# Multi-dataset (3k + 6k) Experiment Workflow

Reproduce the four-dataset experiments A / B1 / B2 from raw 10x data.

Set `PBSC4K_ROOT` to the repo root:
```
export PBSC4K_ROOT=/path/to/repo
```
Set `SCAGA_REPO` to your AluruLab scSAGA checkout (contains `scmint/scsaga_saveH.py`).

## 0. Input data

```
data/3k_rna/counts.mtx barcodes.txt features.txt pca_50.txt   # original 3k RNA  (2,711 cells)
data/3k_atac/...                                                # original 3k ATAC (2,711 cells)
data/6k_rna/...                                                 # 6k RNA  (2,711 cells, subsampled 10k)
data/6k_atac/...                                                # 6k ATAC (2,711 cells, same cells as 6k RNA)
```

The 3k files are produced by `scripts/common/extract_data.py` from the 10x 3k h5.
The 6k files are produced by `scripts/multi-dataset-4x/preprocess_10k.py` from the
10x 10k h5 (splits RNA/ATAC and deterministically subsamples the same 2,711 cells
in both modalities).

```
python preprocess_10k.py
```

## 1. scSAGA integration (4 datasets, anchor rna3k)

```
python scripts/multi-dataset-4x/run_integration.py
```
Writes `results/integration/joint_embedding_H.npy` (10,844 x 30: rna3k, atac3k,
rna6k, atac6k) and per-dataset `aligned_*.npy`.

## 2. Build the all-cells matrices for each experiment

```
python scripts/multi-dataset-4x/build_imputation.py expA combined
python scripts/multi-dataset-4x/build_imputation.py expB1 rna3k
python scripts/multi-dataset-4x/build_imputation.py expB2 rna6k
```
- `expA` needs the SCEMENT combined reference first (step 2a).
- Each writes `results/<exp>/all_cells_gene_expression.npy` (10,844 x 36,601) and
  `imputed_expression_genes_x_atac.npy`.

### 2a. SCEMENT combined RNA reference (Experiment A only)

```
python scripts/multi-dataset-4x/run_scement_combine.py
```
Writes `results/expA/combined_ref_expression.npy` (36,601 x 5,422) and
`combined_ref_embedding.npy` / `combined_ref_barcodes.txt`.

## 3. GRN inference (Arboreto GRNBoost2, TRRUST regulators)

Use the old-stack arboreto env (see ENVIRONMENT.md). Each experiment:
```
python scripts/multi-dataset-4x/run_arboreto.py <exp>       # exp in {expA,expB1,expB2}
```
Writes `results/<exp>/grn/grnboost2_network.tsv`.

## 4. Evaluation vs PBMC ground truth

```
python scripts/multi-dataset-4x/evaluate_grn.py <exp>
```
Writes `results/<exp>/grn/evaluation/evaluation_summary.txt` and `top_edges.csv`.

## 5. Reports

```
python scripts/multi-dataset-4x/write_reports.py    # per-experiment reports/
python scripts/multi-dataset-4x/build_report.py     # consolidated comparison + figures
```

## Environment

Three Python stacks are needed; see `ENVIRONMENT.md`.
