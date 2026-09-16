# RNA-only baseline — tf_only.txt regulators — GRNBoost2 on stacked PBMC 3k + PBMC 4k RNA

**Baseline:** no integration, no imputation — the PBMC 3k multiome RNA and the 10x PBMC 4k RNA matrices are stacked and handed directly to Arboreto.
**Pipeline:** identical to the original RNA-only baseline run; the single change is the transcription-factor list.
**GRN:** Arboreto GRNBoost2; targets = `trrust_tf.txt` genes present in the matrix, **regulators = `tf_only.txt` TFs present in the matrix**.
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated).

---

## 1. Inputs (identical to the original RNA-only baseline)

- **3k RNA:** PBMC multiome 3k granulocyte-sorted, Gene Expression (36601 genes, GRCh38-2020-A), 2711 cells.
- **4k RNA:** 10x PBMC 4k from a healthy donor (v2 / Cell Ranger 2.1.0, 33694 genes, GRCh38-3.0.0), 4340 cells.
- **Stacking:** gene symbols common to both references (**21932** genes); duplicate symbols summed per reference. Each matrix log1p(CPM)-normalized independently, then cells stacked (3k rows first, then 4k).
- **All-cells matrix:** 7051 cells x 21932 genes (bit-identical to the matrix used by the original run).

## 2. Regulator / target specification — the one changed variable

| | Original run (`grn_trrust`) | This run (`grn_tfonly_reg`) |
|---|---|---|
| Regulator list | `data/trrust_tf.txt` | `data/tf_only.txt` |
| Regulator list size | 2862 | 829 |
| Regulators present in matrix | 2808 | 816 |
| Regulators appearing in the network | 2447 | 721 |
| Target genes (`trrust_tf.txt` present) | 2808 | 2808 |
| Expression subset handed to GRNBoost2 | 2808 columns | 2833 columns (union of targets + tf_only regulators) |
| Cells | 7051 | 7051 |
| Inferred edges | 751,233 | 539,629 |

- Regulators dropped relative to the original run: 1992 (non-TF genes of the mixed TRRUST list).
- 816 of the 829 `tf_only.txt` entries are present as columns in the matrix; the remaining 13 are absent from both reference gene sets.

## 3. Evaluation vs ground truth (deduplicated)

### Edges recovered

| Ground truth | Regulators | Dedup. GT edges | Recovered |
|---|---|---|---|
| PBMC-TRRUST | tf_only (816) | 8751 | 2655/8751 (30.3%) |
| PBMC-TRRUST | TRRUST (2808) | 8751 | 1433/8751 (16.4%) |
| PBMC-Blood | tf_only (816) | 96846 | 3661/96846 (3.8%) |
| PBMC-Blood | TRRUST (2808) | 96846 | 1841/96846 (1.9%) |

### Precision@K / Recall@K — PBMC-TRRUST

| Top-K | Prec@K (tf_only) | Recall@K (tf_only) | Prec@K (TRRUST) | Recall@K (TRRUST) |
|---|---|---|---|---|
| 100 | 0.0200 | 0.0002 | 0.0100 | 0.0001 |
| 500 | 0.0380 | 0.0022 | 0.0040 | 0.0002 |
| 1000 | 0.0310 | 0.0035 | 0.0050 | 0.0006 |
| 5000 | 0.0160 | 0.0091 | 0.0054 | 0.0031 |
| 10000 | 0.0130 | 0.0149 | 0.0042 | 0.0048 |
| 539629 | 0.0049 | 0.3034 | 0.0019 | 0.1638 |

### Precision@K / Recall@K — PBMC-Blood

| Top-K | Prec@K (tf_only) | Recall@K (tf_only) | Prec@K (TRRUST) | Recall@K (TRRUST) |
|---|---|---|---|---|
| 100 | 0.0100 | 0.0000 | 0.0100 | 0.0000 |
| 500 | 0.0120 | 0.0001 | 0.0040 | 0.0000 |
| 1000 | 0.0160 | 0.0002 | 0.0030 | 0.0000 |
| 5000 | 0.0148 | 0.0008 | 0.0038 | 0.0002 |
| 10000 | 0.0129 | 0.0013 | 0.0036 | 0.0004 |
| 539629 | 0.0068 | 0.0378 | 0.0025 | 0.0190 |

## 4. Top 25 edges (tf_only regulators)

| TF -> target | importance |
|---|---|
| HOPX -> GNLY | 175.5087 |
| LYAR -> CCL5 | 156.6694 |
| FOS -> DUSP1 | 155.7187 |
| HOPX -> PRF1 | 145.5202 |
| PAX5 -> EBF1 | 132.0980 |
| EBF1 -> PAX5 | 128.9517 |
| HOPX -> IL2RB | 126.3544 |
| MEF2C -> HLA-DRA | 125.0827 |
| HOPX -> GZMB | 123.9629 |
| FOS -> JUN | 123.0006 |
| NEAT1 -> LYZ | 121.7125 |
| NEAT1 -> S100A9 | 114.8055 |
| HOPX -> CCL5 | 113.6789 |
| JUNB -> JUN | 113.0770 |
| FOS -> ZFP36 | 113.0298 |
| HOPX -> RHOC | 111.4303 |
| NEAT1 -> CTSS | 108.8151 |
| NEAT1 -> VCAN | 107.4495 |
| NEAT1 -> MNDA | 106.0488 |
| FOS -> JUNB | 104.9289 |
| NEAT1 -> SAT1 | 102.3271 |
| MEF2C -> CD74 | 100.9747 |
| NEAT1 -> CD14 | 100.1105 |
| HMGB2 -> MKI67 | 99.6758 |
| MEF2C -> HLA-DPA1 | 99.5625 |

## 5. Context: same tf_only regulators under integration + imputation

For reference, the same regulator list (`tf_only.txt`, 816 present) was used on the
full-10k integration + reverse-imputeKNN runs of this repo. Those runs hold 29,218 cells
(14,609 real multiome RNA [2,711 from 3k + 11,898 from 10k] + 14,609 imputed ATAC) versus
the 7,051 real RNA cells here, and 2,827 targets (2,852 union columns) versus 2,808 (2,833):

| Run | Cells | Inferred edges | PBMC-TRRUST recovered | PBMC-Blood recovered |
|---|---|---|---|---|
| RNA-only baseline, tf_only regs (this run) | 7051 | 539,629 | 2655/8751 (30.3%) | 3661/96846 (3.8%) |
| expA (paired full-10k, SCEMENT ref) | 29218 | 829982 | 3660/8751 (41.8%) | 5269/96846 (5.4%) |
| expB1 (paired full-10k, 3k RNA ref) | 29218 | 779701 | 3454/8751 (39.5%) | 5132/96846 (5.3%) |
| expB2 (paired full-10k, 10k RNA ref) | 29218 | 832664 | 3664/8751 (41.9%) | 5362/96846 (5.5%) |
| unpaired full-10k (external 10x ATAC) | 20059 | 659314 | 3110/8751 (35.5%) | 4486/96846 (4.6%) |

*Paired/unpaired values from `reports/full_10k_3k_onlyTf/comparison_report_tfonly.md` and
`reports/unpaired-10k/comparison_report_tfonly_unpaired.md`; targets there are 2,827 vs 2,808
here because the all-cells matrices carry different gene sets. Cell count, imputation and
target count all differ, so this is context, not a controlled comparison.*

---
