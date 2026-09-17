# RNA-only baseline #2 (PBMC multiome 3k RNA + PBMC multiome 10k RNA → GRNBoost2)

A second RNA-only baseline control for the integration/imputation workflows:
**no integration and no imputation**. Two multiome Gene Expression matrices are
stacked and handed directly to Arboreto GRNBoost2.

This run stacks the two *multiome* RNA matrices (3k + 10k, same reference,
GRCh38-2020-A). The original RNA-only baseline
(`scripts/rna-baseline/`) stacked the multiome 3k RNA with the *external*
10x v2 4k RNA (3k+4k, two different references).

## Datasets

| Dataset | Source | Cells | Genes | Reference |
|---|---|---|---|---|
| 3k RNA | PBMC multiome 3k granulocyte-sorted, Gene Expression | 2,711 | 36,601 | GRCh38-2020-A |
| 10k RNA | PBMC multiome 10k granulocyte-sorted, Gene Expression | 11,898 | 36,601 | GRCh38-2020-A |

Both files come from `scripts/common/extract_data.py` (10x multiome h5) and the
10k full-multiome preprocessor. Point the data dirs via `RNA_3K_DIR` /
`RNA_10K_DIR` (the matrices are too large to live under the repo `data/`).

## Stacking

1. Both references are the same (GRCh38-2020-A) with the *identical* feature
   list, so no reordering is required.
2. Duplicate gene symbols within a reference (10 per reference) are collapsed by
   summing their counts → keep the **36,591** common symbols (shared order).
3. log1p(CPM)-normalize each matrix independently (per-cell library-size scaling).
4. Stack cells: **3k rows first, then 10k** → all-cells matrix **14,609 x 36,591**.

## GRN

GRNBoost2 with **tf_only.txt as regulators** and **trrust_tf.txt genes as
targets** (2,827 targets present; 816 regulators present; 2,852 union columns —
every regulator must be a matrix column, so the expression subset is the union).
Input to GRNBoost2 is 14,609 cells x 2,852 genes.

```bash
export PBSC4K_ROOT=/path/to/repo
RNA_3K_DIR=/Volumes/samsung_ssd/tmp/pbsc4k-multiome-experiment/data/3k_rna \
RNA_10K_DIR=/Volumes/samsung_ssd/tmp/pbmc-full10k-grn/data/10k_rna \
  python scripts/rna-baseline-3k10k/prepare_rna_baseline_3k10k.py    # scSAGA .venv
PBSC4K_ROOT=$PWD python scripts/rna-baseline-3k10k/run_arboreto_tfonly_3k10k.py  # .venv39 (arboreto)
python scripts/rna-baseline-3k10k/evaluate_grn_tfonly_3k10k.py           # scSAGA .venv
python scripts/rna-baseline-3k10k/write_report_tfonly_3k10k.py            # scSAGA .venv
```

`run_arboreto_tfonly_3k10k.py` must run in the old-stack `.venv39` env (python
3.9, dask 2021.10, arboreto 0.1.6); see `ENVIRONMENT.md`. The network `.tsv` is
gitignored.

## Outputs

```
results/rna_baseline_3k10k/grn_tfonly_reg/grnboost2_network.tsv   (not committed)
results/rna_baseline_3k10k/grn_tfonly_reg/{tf_regulators,target_genes,all_columns}.txt
results/rna_baseline_3k10k/grn_tfonly_reg/n_cells.txt
results/rna_baseline_3k10k/grn_tfonly_reg/evaluation/evaluation_summary.txt
results/rna_baseline_3k10k/grn_tfonly_reg/evaluation/top_edges.csv
reports/rna-baseline-3k10k/REPORT_tfonly.md
```
