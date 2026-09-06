# PBMC Multiome GRN inference — single- and multi-dataset

End-to-end gene regulatory network (GRN) inference on 10x Genomics PBMC multiome
(paired RNA + ATAC) data using **scSAGA** integration, **reverse-imputeKNN**
imputation, and **Arboreto GRNBoost2**, evaluated against PBMC ground-truth
regulatory edges.

This repository covers three workflows:

1. **Single-dataset (3k)** — the original pbmc_3k granulocyte-sorted multiome
   (2,711 shared cells): scSAGA integrates RNA+ATAC, reverse-imputeKNN fills ATAC
   expression from RNA, GRNBoost2 infers the network, evaluated against
   PBMC-TRRUST / PBMC-Blood.
2. **Multi-dataset (3k + 6k)** — extends integration to four datasets
   (2,711 cells each: 3k-RNA, 3k-ATAC, and a 2,711-cell subset of the 10k dataset
   as 6k-RNA, 6k-ATAC) and compares **three reverse-imputeKNN reference
   strategies**:
   - **Experiment A**: a single SCEMENT-integrated combined RNA reference
     (3k + 6k RNA, 5,422 cells).
   - **Experiment B1**: the single 3k RNA reference.
   - **Experiment B2**: the single 6k RNA reference.
3. **Multi-dataset full-10k** — the same four-dataset integration and three
   reference strategies, but using the **entire 10k multiome (11,898 cells)**
   instead of a 2,711-cell subsample. All-cells matrix is 29,218 x 36,601.
   See `reports/full-10k/` and the `*_10k` scripts in `scripts/multi-dataset-4x/`.

The only variable across Experiments A/B1/B2 is which RNA reference feeds
reverse-imputeKNN (the all-cells matrix is always the same real RNA cells plus
the imputed ATAC cells), so the experiments isolate the effect of the
imputation reference strategy.

## Pipeline

```
preprocess → scSAGA integration (joint embedding H) → reverse-imputeKNN
  → all-cells gene-expression matrix → Arboreto GRNBoost2 → evaluation
```

## Repository layout

```
├── scripts/
│   ├── common/               # shared steps (h5 split, scSAGA runner, single-3k)
│   ├── single-dataset-3k/    # original 3k downstream (reverse-impute, arboreto, eval)
│   └── multi-dataset-4x/     # 4-dataset experiments A/B1/B2 (+ vendored SCEMENT)
├── reports/                  # per-experiment + consolidated comparison reports
│   ├── single-dataset-3k/    # original 3k workflow reports (REPORT.md, REPORT_DETAILED.md)
│   └── full-10k/             # full-10k (11,898-cell) multi-dataset reports + figures
├── data/
│   ├── trrust_tf.txt         # TRRUST transcription factors (regulators)
│   └── ground_truth/         # PBMC-TRRUST.csv, PBMC-Blood.csv
├── results/                  # integration, per-experiment GRN + evaluations
├── ENVIRONMENT.md            # per-step Python environment recipes
├── README.md
```

Raw 10x `.h5` files and large regenerable matrices are **not** committed
(see `.gitignore`); they are reproduced by the scripts.

## Datasets

- **3k**: 10x "PBMC from a Healthy Donor — Granulocytes Removed Through Cell
  Sorting (3k)", Multiome ATAC+GEX (Cell Ranger ARC 2.0.0). 134,920 features
  (36,601 RNA + 98,319 peaks) × 2,711 cells, shared nuclei.
- **10k → 6k**: the same protocol at 10k nominal (11,898 cells); we deterministically
  subsample the **same 2,711 cells** in both modalities to form the "6k" datasets.
- **10k (full)**: the same 10k multiome using **all 11,898 cells** (no subsampling),
  used in the full-10k workflow.

## Results summary

See `reports/00_COMPARISON.md` for the consolidated table,
`reports/exp{EXP}_REPORT.md` for per-experiment detail, and
`reports/full-10k/` for the full-10k (11,898-cell) runs.

| Experiment | Imputation reference | GRN evaluation |
|---|---|---|
| A  | SCEMENT-combined 3k+6k RNA (5,422) | PBMC-TRRUST / PBMC-Blood |
| B1 | 3k RNA only (2,711)                 | PBMC-TRRUST / PBMC-Blood |
| B2 | 6k RNA only (2,711)                 | PBMC-TRRUST / PBMC-Blood |

Full-10k variants (A: SCEMENT 3k+10k RNA, 14,609 cells; B1: 3k RNA; B2: 10k RNA)
are in `reports/full-10k/`.

## Reproducing

See `ENVIRONMENT.md` for the exact Python stacks. Set `PBSC4K_ROOT` to the repo
root, then run the scripts in the order given in `scripts/multi-dataset-4x/`
and `scripts/single-dataset-3k/`.
