# Unpaired 10k workflow — multiome RNA × external 10x ATAC v1.1

Reproduce the **unpaired** experiment: the RNA-seq of the 10x PBMC 10k multiome
paired against a *separate* 10x ATAC dataset (no shared nuclei), then integrate →
reverse-imputeKNN → GRNBoost2.

Unlike the other workflows in this repo, the two modalities come from **different
cells**. There is no cell-level pairing, so the alignment is partial by construction.

## 0. Input data

```
raw/atac10k_v1.1_peaks.h5          # downloaded, see below
data/10k_rna/                      # 10k multiome RNA (11,898 cells, GRCh38) — symlink/regenerate
trrust_tf.txt                      # TARGET genes
data/tf_only.txt                   # REGULATOR TFs
data/ground_truth/{PBMC-TRRUST,PBMC-Blood}.csv
```

Download the external ATAC "cells by peaks" matrix:

```sh
curl -L -o raw/atac10k_v1.1_peaks.h5 \
  https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_v1/atac_pbmc_10k_v1_filtered_peak_bc_matrix.h5
```

Dataset page: <https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-atac-v-1-1-chromium-controller-1-1-standard-2-0-0>

Note the naming: the page says **v1.1**, but the 10x sample slug is
`atac_pbmc_10k_v1` and the pipeline version directory is `1.1.0`. This is 8,161 cells ×
87,863 peaks, **hg19** (the multiome is GRCh38) — *not* the v2 `nextgem` build.

## 1. Preprocess the external ATAC into scSAGA format

```sh
python scripts/preprocess_atac10k_ext.py
```
Writes `data/atac10k_ext/{counts.mtx, barcodes.txt, features.txt, pca_50.txt}` using the
same PCA recipe as the in-repo modalities (log1p CPM/1e4, top-2000 variable features,
50 PCs).

## 2. scSAGA integration (unpaired, rna10k anchor)

```sh
python scripts/run_integration_unpaired.py
```
Writes `results/integration_unpaired/joint_embedding_H.npy` — **20,059 × 30**
(11,898 rna10k + 8,161 atac10k_ext), plus the transport plan and `aligned_*.npy`.
`s_shared_cells = 8,161 = min(n_rna, n_atac)`.

## 3. Reverse-imputeKNN + all-cells matrix

```sh
python scripts/build_imputation_unpaired.py
```
Reference = real RNA10k expression; query = the 8,161 external ATAC cells.
Writes `results/unpaired_grn/`:

```
imputed_expression_genes_x_atac.npy   36601 x 8161
all_cells_gene_expression.npy         20059 x 36601   (rows: rna10k real, atac-ext imputed)
genes.npy / genes.txt / all_cells_barcodes.txt / n_cells.txt
```

## 4. GRN inference (Arboreto GRNBoost2)

```sh
bash scripts/run_all_grn_unpaired.sh          # caffeinate + ulimit wrapper
```
- **Targets** = genes present in `trrust_tf.txt` (2,827)
- **Regulators** (`tf_names`) = TFs present in `data/tf_only.txt` (816)
- Union columns fed to GRNBoost2 = 2,852

Writes `results/unpaired_grn/grn_tfonly_reg/grnboost2_network.tsv`.
Runs in the old-stack env (`~/.venv39`: python 3.9, numpy 1.21, dask 2021.10,
arboreto 0.1.6); `GRN_WORKERS` env var selects worker count (default 4).

## 5. Evaluation

```sh
python scripts/evaluate_grn_unpaired.py
```
Writes `results/unpaired_grn/grn_tfonly_reg/evaluation/{evaluation_summary.txt,
top_edges.csv, evaluation.json}` scored against deduplicated PBMC-TRRUST / PBMC-Blood.

## 6. QC

```sh
python scripts/qc_unpaired.py     # distributional sanity, marker recovery
python scripts/qc_mixing.py       # block-mixing benchmark vs the paired 4-dataset run
python scripts/qc_biology.py      # independent ATAC clustering → imputed marker coherence
```

`qc_biology.py` is the decisive one for an unpaired run: mixing statistics are
block-segregated in the paired runs too, so they cannot diagnose alignment quality.
Cluster the ATAC cells on their own PCA and check that imputed marker profiles match
independently-derived RNA clusters against a shuffled control.

## 7. Report

```sh
python scripts/build_report_unpaired.py     # -> reports/unpaired-10k/REPORT.md
```

## Environments

| Step | Interpreter |
|---|---|
| preprocessing, integration, imputation, QC, evaluation, report | scSAGA `.venv` (py3.12; h5py, torch, pot, geosketch, pydantic) |
| Arboreto GRNBoost2 | `.venv39` (py3.9, numpy 1.21, dask 2021.10, arboreto 0.1.6) |

The py3.9 env has no h5py — do the h5 split with the py3.12 venv.

## Interpretive caveats

1. **No shared cells** — partial alignment only; the two cell sets cannot be matched.
2. **hg19 vs GRCh38** — peak coordinates never enter scSAGA (it integrates on PCA), but
   biological comparability of the ATAC modality differs.
3. **40.7% of GRN rows are imputed**, and the imputation is concentrated: the nearest
   reference of all 8,161 queries comes from only 2,831 of 11,898 RNA cells.
4. Cell-type structure (B, mono, NK, T) survives the alignment; rare populations
   (platelet PPBP, basophil/DC FCER1A) do not. See `reports/unpaired-10k/REPORT.md`.
