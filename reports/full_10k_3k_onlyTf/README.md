# Full-10k scSAGA Integration & GRN — TF-only Regulators (tf_only.txt)

This folder documents the **full-10k** GRN runs where the Arboreto GRNBoost2
**regulators are restricted to the tf_only.txt transcription-factor list (816 TFs
present in the expression matrix)** instead of the full TRRUST mixed list (2,827).
Target genes remain the TRRUST TFs present (2,827). Regulator list curated from the
TF columns of the PBMC ground-truth files (see the tf-only-list curation note).

## Setup (identical integration to the TRRUST-only full-10k runs)

| Dataset | Cells | Source |
|---|---|---|
| rna3k  | 2,711 | 10x PBMC multiome 3k (granulocyte-sorted) |
| atac3k | 2,711 | 10x PBMC multiome 3k (granulocyte-sorted) |
| rna10k | 11,898 | 10x PBMC multiome 10k — **all cells** |
| atac10k| 11,898 | 10x PBMC multiome 10k — **all cells** |

Anchor for scSAGA integration = rna3k. All-cells matrix for every experiment:
**29,218 x 36,601** (rna3k, rna10k, atac3k-imputed, atac10k-imputed).

## Experiments (reverse-imputeKNN reference strategy)

| Experiment | Reference | Reference cells |
|---|---|---|
| A  | SCEMENT-integrated 3k+10k RNA | 14,609 |
| B1 | 3k RNA only | 2,711 |
| B2 | 10k RNA only | 11,898 |

## GRN inference

Arboreto GRNBoost2, regulators = 816 tf_only TFs, targets = 2,827 TRRUST genes
(union columns = 2,852 so every regulator is a matrix column). 29,218 cells,
4 dask workers, ~4 h CPU per experiment (run sequentially).

## Results summary

| Experiment | Inferred edges | PBMC-TRRUST recovered | PBMC-Blood recovered |
|---|---|---|---|
| A  | 829,982 | 3,660 / 8,751 (41.8%) | 5,269 / 96,846 (5.4%) |
| B1 | 779,701 | 3,454 / 8,751 (39.5%) | 5,132 / 96,846 (5.3%) |
| B2 | 832,664 | 3,664 / 8,751 (41.9%) | 5,362 / 96,846 (5.5%) |

See `comparison_report_tfonly.md` for full Precision@K tables and per-experiment detail.

## Reproducing

Scripts (in this working checkout under `scripts/`):
```
run_arboreto_tfonly_reg_10k.py   # arboreto per experiment (writes results/<exp>_10k/grn_tfonly_reg/)
run_all_grn_tfonly_reg_10k.sh    # sequential runner (expA, expB1, expB2)
```
Environments: `.venv39` for Arboreto (python 3.9).
