# RNA-only baseline (3k + 10k multiome) — tf_only.txt regulators — GRNBoost2

**Baseline:** no integration, no imputation — the PBMC multiome 3k Gene Expression and the PBMC multiome 10k Gene Expression matrices are stacked and handed directly to Arboreto.
**GRN:** Arboreto GRNBoost2; targets = `trrust_tf.txt` genes present in the matrix, **regulators = `tf_only.txt` TFs present in the matrix**.
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated).

---

## 1. Inputs

- **3k RNA:** PBMC multiome 3k granulocyte-sorted, Gene Expression (36601 genes, GRCh38-2020-A), 2711 cells.
- **10k RNA:** PBMC multiome 10k granulocyte-sorted, Gene Expression (36601 genes, GRCh38-2020-A), 11898 cells.
- **Stacking:** both matrices use the same Cell Ranger reference (GRCh38-2020-A) and the identical feature list, so no reordering is needed. Duplicate symbols (10 per reference) summed per reference; keep the **36591** common symbols. Each matrix log1p(CPM)-normalized independently, then cells stacked (3k rows first, then 10k).
- **All-cells matrix:** 14609 cells x 36591 genes.

## 2. Regulator / target specification

| Regulator list | `data/tf_only.txt` |
|---|---|
| Regulator list size | 829 |
| Regulators present in matrix | 816 |
| Regulators appearing in the network | 751 |
| Target genes (`trrust_tf.txt` present) | 2827 |
| Expression subset handed to GRNBoost2 | 2852 columns (union of targets + tf_only regulators) |
| Cells | 14609 |
| Inferred edges | 616,865 |

## 3. Evaluation vs ground truth (deduplicated)

| Ground truth | Dedup. GT edges | Recovered |
|---|---|---|
| PBMC-TRRUST | 8751 | 2920/8751 (33.4%) |
| PBMC-Blood | 96846 | 4210/96846 (4.3%) |

### Precision@K / Recall@K — PBMC-TRRUST

| Top-K | Prec@K | Recall@K |
|---|---|---|
| 100 | 0.0500 | 0.0006 |
| 500 | 0.0380 | 0.0022 |
| 1000 | 0.0280 | 0.0032 |
| 5000 | 0.0188 | 0.0107 |
| 10000 | 0.0161 | 0.0184 |
| 616865 | 0.0047 | 0.3337 |

### Precision@K / Recall@K — PBMC-Blood

| Top-K | Prec@K | Recall@K |
|---|---|---|
| 100 | 0.0300 | 0.0000 |
| 500 | 0.0300 | 0.0002 |
| 1000 | 0.0270 | 0.0003 |
| 5000 | 0.0190 | 0.0010 |
| 10000 | 0.0180 | 0.0019 |
| 616865 | 0.0068 | 0.0435 |

## 4. Top 25 edges

| TF -> target | importance |
|---|---|
| PAX5 -> EBF1 | 218.5450 |
| EBF1 -> PAX5 | 190.2514 |
| PAX5 -> MS4A1 | 188.4677 |
| FOSB -> JUN | 172.5060 |
| TCF7L2 -> CDKN1C | 141.8560 |
| LEF1 -> NELL2 | 138.4754 |
| PAX5 -> CD22 | 134.3404 |
| PAX5 -> CD79A | 133.0548 |
| EBF1 -> MS4A1 | 133.0394 |
| TCF4 -> RUNX2 | 132.6735 |
| LEF1 -> PRKCA | 131.4753 |
| NEAT1 -> VCAN | 128.9541 |
| MYBL1 -> CCL5 | 128.6566 |
| CREB5 -> VCAN | 128.1028 |
| TCF4 -> CLEC4C | 127.5453 |
| NEAT1 -> NAMPT | 125.2810 |
| PAX5 -> FCER2 | 122.5395 |
| PAX5 -> BLK | 121.1737 |
| ZEB2 -> PLEK | 119.0654 |
| LEF1 -> FHIT | 116.9675 |
| LEF1 -> CCR7 | 114.4464 |
| NEAT1 -> VMP1 | 114.3166 |
| NEAT1 -> RPL3 | 111.9062 |
| LEF1 -> TCF7 | 110.8090 |
| MEF2C -> HLA-DRB1 | 110.3182 |

## 5. Context: same tf_only regulators, other baselines & workflows

All rows use the `tf_only.txt` regulator list (816 present). Cell count and target count differ per row, so this is context, not a controlled comparison.

| Run | Cells | Targets | Inferred edges | PBMC-TRRUST recovered | PBMC-Blood recovered |
|---|---|---|---|---|---|
| RNA-only 3k+10k multiome (this run) | 14609 | 2827 | 616,865 | 2920/8751 (33.4%) | 4210/96846 (4.3%) |
| RNA-only 3k+4k (tf_only, earlier) | 7051 | 2808 | 539629 | 2655/8751 (30.3%) | 3661/96846 (3.8%) |

*3k+4k values from `reports/rna-baseline/REPORT_tfonly.md`. Full-10k integration + reverse-imputeKNN runs (expA/expB1/expB2, 29,218 cells) are in `reports/full_10k_3k_onlyTf/comparison_report_tfonly.md`; unpaired external-ATAC run in `reports/unpaired-10k/comparison_report_tfonly_unpaired.md`.*

---
