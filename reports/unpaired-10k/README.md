# Unpaired 10k — multiome RNA × external 10x ATAC v1.1 (TF-only regulators)

Variant of the full-10k TF-only-regulator experiments in which the ATAC modality is a
**separate 10x dataset** instead of the multiome's own ATAC. Same regulator/target
definition as `reports/full_10k_3k_onlyTf/` (regulators = `data/tf_only.txt`,
targets = `trrust_tf.txt` genes present).

## Setup (differs from the paired full-10k runs)

| Dataset | Cells | Source |
|---|---|---|
| rna10k (anchor) | 11,898 | 10x PBMC 10k multiome (granulocyte-sorted), GRCh38 |
| atac10k_ext | 8,161 | 10x **10k Human PBMCs ATAC v1.1** "cells by peaks", hg19 |

The two modalities share **no nuclei** (barcode overlap 0): different donor/protocol and
different peak genome build. The alignment is therefore partial by construction.

- External ATAC download:
  `https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_v1/atac_pbmc_10k_v1_filtered_peak_bc_matrix.h5`
  (the dataset page says "v1.1"; the 10x sample slug is `atac_pbmc_10k_v1`, Cell Ranger
  ATAC 1.1.0 — not the v2 `nextgem` build)
- scSAGA anchor = rna10k; joint embedding **H = 20,059 × 30**
- All-cells matrix: **20,059 × 36,601** (11,898 real RNA + 8,161 imputed ATAC)
- Reference for reverse-imputeKNN = the single 10k RNA reference (expB2 analogue)

## GRN inference

Arboreto GRNBoost2, regulators = 816 `tf_only` TFs, targets = 2,827 `trrust_tf.txt` genes
(union columns = 2,852). 20,059 cells, **4 dask workers, 2h10m** wall clock.

## Results summary

| Ground truth | Deduplicated | Recovered |
|---|---|---|
| PBMC-TRRUST | 8,751 | 3,110 (35.5%) |
| PBMC-Blood | 96,846 | 4,486 (4.6%) |

Inferred edges: 659,314.

For comparison, the paired full-10k TF-only runs (29,218 cells) recovered 39.5–41.9%
(PBMC-TRRUST) and 5.3–5.5% (PBMC-Blood). See
`comparison_report_tfonly_unpaired.md` for the full table, Precision@K tables and figures.

## QC

`qc_imputation.txt`, `qc_mixing.txt`, `qc_biology.txt` — plus `REPORT.md` §6.

The decisive check for an unpaired run is `qc_biology.txt`: clustering the ATAC cells on
their **own peak PCA** and matching imputed marker profiles to independently-derived RNA
clusters gives mean matched r = +0.741 vs +0.485 for a shuffled control (B, mono, NK and
T lineages assigned correctly). Block-mixing statistics are **not** diagnostic: they are
segregated (0.988–1.000) in the repo's paired runs too.

Known degradation: attenuated marker amplitude, rare populations lost (PPBP, FCER1A),
concentrated imputation (nearest reference of all queries from only 2,831 of 11,898 RNA
cells), 40.7% of GRN rows imputed.

## Reproducing

Scripts: `scripts/unpaired-10k/` (see `PIPELINE.md` there).
Environments: scSAGA `.venv` (py3.12) for preprocessing/integration/imputation/QC;
`.venv39` (py3.9, numpy 1.21, dask 2021.10, arboreto 0.1.6) for GRNBoost2.

Working checkout: `/Volumes/samsung_ssd/tmp/pbmc-10k-unpaired-grn/`
