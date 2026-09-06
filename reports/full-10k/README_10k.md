# Full-10k Multi-Dataset scSAGA Integration & GRN (TRRUST-only targets)

This folder documents the **full-10k** variant of the multi-dataset GRN experiments:
the same four-dataset scSAGA integration and reverse-imputeKNN pipeline as the
`multi-dataset-4x` (3k + 6k) runs, but using the **entire 10k multiome (11,898 cells)**
instead of a 2,711-cell subsample.

## Datasets

| Dataset | Cells | Source |
|---|---|---|
| rna3k  | 2,711 | 10x PBMC multiome 3k (granulocyte-sorted) |
| atac3k | 2,711 | 10x PBMC multiome 3k (granulocyte-sorted) |
| rna10k | 11,898 | 10x PBMC multiome 10k — **all cells** |
| atac10k| 11,898 | 10x PBMC multiome 10k — **all cells** |

The 10k RNA/ATAC are the same 11,898 physical cells (paired multiome). Anchor for
scSAGA integration = rna3k.

## Experiments (reverse-imputeKNN reference strategy)

| Experiment | Reference | Reference cells |
|---|---|---|
| A  | SCEMENT-integrated 3k+10k RNA | 14,609 |
| B1 | 3k RNA only | 2,711 |
| B2 | 10k RNA only | 11,898 |

All-cells matrix for every experiment: **29,218 x 36,601** (rows: rna3k, rna10k,
atac3k-imputed, atac10k-imputed). Only the imputed ATAC values vary by experiment.

## GRN inference

Arboreto GRNBoost2 with **TRRUST-only target genes** (2,827 TFs present in the
expression matrix). Each run: 29,218 cells x 2,827 targets, 2 dask workers,
~18-19 h CPU per experiment (run sequentially).

## Results summary

| Experiment | Inferred edges | PBMC-TRRUST recovered | PBMC-Blood recovered |
|---|---|---|---|
| A  | 1,177,622 | 1,957 / 8,751 (22.4%) | 2,642 / 96,846 (2.7%) |
| B1 | 1,095,408 | 1,908 / 8,751 (21.8%) | 2,528 / 96,846 (2.6%) |
| B2 | 1,189,398 | 1,939 / 8,751 (22.2%) | 2,694 / 96,846 (2.8%) |

See `comparison_report_trrust_10k.md` for the full Precision@K tables and
per-experiment detail.

## Comparison vs the 6k (subsampled) run

Using the full 10k (11,898 cells) instead of a 2,711-cell subsample improves
ground-truth edge recovery:

| Experiment | PBMC-TRRUST recovered (6k) | PBMC-TRRUST recovered (full-10k) |
|---|---|---|
| A  | 1,603 / 8,751 (18.3%) | 1,957 / 8,751 (22.4%) |
| B1 | 1,657 / 8,751 (18.9%) | 1,908 / 8,751 (21.8%) |
| B2 | 1,652 / 8,751 (18.9%) | 1,939 / 8,751 (22.2%) |

## Reproducing

Scripts live in `scripts/multi-dataset-4x/` (the `*_10k` variants). Set
`PBSC4K_ROOT` to the repo root, then run in order:

```
python scripts/multi-dataset-4x/preprocess_10k_full.py      # -> data/10k_rna, data/10k_atac
python scripts/multi-dataset-4x/run_integration_10k.py     # -> results/integration_10k/
python scripts/multi-dataset-4x/run_scement_combine_10k.py  # -> results/expA_10k/ (SCEMENT ref)
python scripts/multi-dataset-4x/build_imputation_10k.py expA combined
python scripts/multi-dataset-4x/build_imputation_10k.py expB1 rna3k
python scripts/multi-dataset-4x/build_imputation_10k.py expB2 rna10k
bash scripts/multi-dataset-4x/run_all_grn_trrust_10k.sh     # sequential GRN + eval
python scripts/multi-dataset-4x/build_report_trrust_10k.py  # -> reports/
python scripts/multi-dataset-4x/make_figures_10k.py          # -> reports/figures_10k/
```

Environments: see `ENVIRONMENT.md` (scSAGA .venv for integration/imputation/eval,
.venv39 for Arboreto, .venv-scement for SCEMENT).
