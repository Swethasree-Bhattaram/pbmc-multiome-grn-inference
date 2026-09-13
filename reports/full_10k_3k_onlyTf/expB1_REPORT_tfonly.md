# Experiment B1 — 3k RNA only (2,711 cells) (full-10k, TF-only regulators)

**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, 10k-RNA, 10k-ATAC)
**Reverse-imputeKNN reference:** 3k RNA only (2,711 cells)
**GRN:** Arboreto GRNBoost2 — regulators = tf_only.txt (816 present), targets = TRRUST (2,827)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets: rna3k (2,711), atac3k (2,711), rna10k (11,898), atac10k (11,898). The 10k RNA/ATAC are the SAME 11,898 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (29,218 x 30).
3. **Reverse-imputeKNN**: reverse-imputeKNN from 3k RNA only.
4. **All-cells matrix:** 29,218 x 36,601 (rows: rna3k, rna10k, atac3k-imputed, atac10k-imputed).
5. **GRN inference:** Arboreto GRNBoost2 — regulators restricted to the 816 TF-only list present in the matrix (tf_only.txt), targets = the 2,827 TRRUST genes present.
6. **Evaluation:** against PBMC-TRRUST and PBMC-Blood ground truth.

## 2. Integration (joint embedding H)

- Joint embedding **H**: 29,218 x 30 (blocks: rna3k, atac3k, rna10k, atac10k).
- Global alignment score: 0.0984 (4 datasets, 3k anchor).
- Pairwise scores vs anchor (rna3k): atac3k 0.4027, rna10k 0.2240, atac10k 0.1236.

## 3. Imputation (reverse-imputeKNN, k=20, softmax-of-negative-distance weights)

- Imputed expression: **36601 genes x 14609 ATAC cells**.
- All-cells matrix: **29218 cells x 36601 genes** (rna3k + rna10k real + 14,609 imputed ATAC).

## 4. Arboreto GRNBoost2 results (tf_only regulators, TRRUST-only targets)

- Regulators (tf_only TFs present in expression): 816
- Target genes: 2827 (TRRUST-only); union columns (targets + regulator-only): 2852
- Edges inferred: **779,701**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| EBF1 → PAX5 | 208.1 |
| PAX5 → MS4A1 | 181.7 |
| PAX5 → EBF1 | 168.4 |
| PAX5 → CD79A | 166.4 |
| FOSB → JUN | 143.8 |
| LEF1 → PRKCA | 139.7 |
| PAX5 → CD22 | 138.6 |
| LEF1 → NELL2 | 136.7 |
| FOXP3 → IL2RA | 134.0 |
| TCF7L2 → CDKN1C | 133.9 |
| PAX5 → BLK | 133.2 |
| ZEB2 → LTB | 130.7 |
| EBF1 → MS4A1 | 129.8 |
| TCF4 → RUNX2 | 129.0 |
| MYBL1 → CCL5 | 127.3 |
| PAX5 → FCER2 | 120.5 |
| BACH2 → SELL | 116.7 |
| LEF1 → FHIT | 116.5 |
| LEF1 → BACH2 | 115.9 |
| LEF1 → CCR7 | 112.2 |
| ZEB2 → PLEK | 111.1 |
| PAX5 → BIRC3 | 110.9 |
| CREB5 → VCAN | 110.8 |
| BACH2 → FOXP1 | 108.9 |
| BACH2 → S100A4 | 108.5 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=3454 (39.5%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.05 | 0.0006 |
| 500 | 0.04 | 0.0023 |
| 1,000 | 0.026 | 0.003 |
| 5,000 | 0.018 | 0.0103 |
| 10,000 | 0.0152 | 0.0174 |
| 779,701 | 0.0044 | 0.3947 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=5132 (5.3%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.03 | 0.0 |
| 500 | 0.03 | 0.0002 |
| 1,000 | 0.023 | 0.0002 |
| 5,000 | 0.0164 | 0.0008 |
| 10,000 | 0.0161 | 0.0017 |
| 779,701 | 0.0066 | 0.053 |

---
