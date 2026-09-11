# RNA-only baseline (PBMC 3k RNA + PBMC 4k RNA → GRNBoost2)

Baseline control for the integration/imputation workflows: **no integration and
no imputation**. Two RNA gene-expression matrices are stacked and handed directly
to Arboreto GRNBoost2, so the resulting network reflects only the raw
concatenation of two RNA datasets.

## Datasets

| Dataset | Source | Cells | Genes | Reference |
|---|---|---|---|---|
| 3k RNA | PBMC multiome 3k granulocyte-sorted, Gene Expression | 2,711 | 36,601 | GRCh38-2020-A |
| 4k RNA | 10x "4k PBMCs from a Healthy Donor" (v2, Cell Ranger 2.1.0) | 4,340 | 33,694 | GRCh38-3.0.0 |

4k download (filtered matrices, ~18 MB):
```
https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz
```
The 3k files come from `scripts/common/extract_data.py` (10x multiome h5).

## Stacking

The two references share neither gene count nor gene order, so:

1. Gene symbols that repeat within a reference (multiple Ensembl IDs sharing a
   symbol) are collapsed by summing their counts.
2. Keep the **gene symbols common to both** references — 21,932 genes
   (3k 36,601 -> 36,591 unique; 4k 33,694 -> 33,660 unique; intersection 21,932).
3. log1p(CPM)-normalize each matrix independently (per-cell library-size scaling).
4. Stack cells: **3k rows first, then 4k** → all-cells matrix 7,051 x 21,932.

## GRN

GRNBoost2 with **TRRUST TFs as regulators and TRRUST-only targets** (2,808 TFs
present in the common gene set; no top-2000 HVG). Input to GRNBoost2 is
7,051 cells x 2,808 genes.

```
export PBSC4K_ROOT=/path/to/repo
python scripts/rna-baseline/prepare_rna_baseline.py     # scSAGA .venv
N_WORKERS=5 python scripts/rna-baseline/run_arboreto_trrust.py   # .venv39 (arboreto)
python scripts/rna-baseline/evaluate_grn_trrust.py      # scSAGA .venv
```

`run_arboreto_trrust.py` reads the TRRUST list from `data/trrust_tf.txt` and must
run in the old-stack `.venv39` env (python 3.9, dask 2021.10, arboreto 0.1.6);
see `ENVIRONMENT.md`.

## Outputs

```
results/rna_baseline/grn_trrust/grnboost2_network.tsv   (not committed; ~751k edges)
results/rna_baseline/grn_trrust/{tf_regulators,target_genes}.txt
results/rna_baseline/grn_trrust/n_cells.txt
results/rna_baseline/grn_trrust/evaluation/evaluation_summary.txt
results/rna_baseline/grn_trrust/evaluation/top_edges.csv
reports/rna-baseline/REPORT.md
```

## Result summary

- Edges inferred: **751,233**; regulators/targets: 2,447 each.
- PBMC-TRRUST: recovered **1433/8751 (16.4%)**.
- PBMC-Blood: recovered **1841/96846 (1.9%)**.
- Top edges are lineage-coherent (CD8A↔CD8B, HLA-B↔HLA-C, CD79A↔MS4A1, LYZ→S100A9,
  GNLY↔PRF1), i.e. the baseline recovers real PBMC co-expression structure even
  without integration.
