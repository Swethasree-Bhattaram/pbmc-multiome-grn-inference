# Experiment B2 — 10k RNA only (11,898 cells) (full-10k, TF-only regulators)

**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, 10k-RNA, 10k-ATAC)
**Reverse-imputeKNN reference:** 10k RNA only (11,898 cells)
**GRN:** Arboreto GRNBoost2 — regulators = tf_only.txt (816 present), targets = TRRUST (2,827)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets: rna3k (2,711), atac3k (2,711), rna10k (11,898), atac10k (11,898). The 10k RNA/ATAC are the SAME 11,898 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (29,218 x 30).
3. **Reverse-imputeKNN**: reverse-imputeKNN from 10k RNA only.
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
- Edges inferred: **832,664**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| SPIB → DUX4 | 233.5 |
| EBF1 → PAX5 | 200.1 |
| HDAC3 → DUX4 | 197.1 |
| PAX5 → MS4A1 | 192.4 |
| PAX5 → CD79A | 180.2 |
| PAX5 → EBF1 | 165.7 |
| PAX5 → BLK | 151.3 |
| FOSB → JUN | 148.6 |
| TSG101 → DUX4 | 147.3 |
| LEF1 → PRKCA | 146.5 |
| PAX5 → CD22 | 137.4 |
| CREB5 → VCAN | 136.8 |
| TCF7L2 → CDKN1C | 136.6 |
| ZEB2 → LTB | 135.7 |
| HEYL → DACT2 | 134.8 |
| LEF1 → NELL2 | 128.8 |
| ERMAP → DUX4 | 127.0 |
| MYBL1 → CCL5 | 126.7 |
| TCF4 → RUNX2 | 126.7 |
| EBF1 → MS4A1 | 126.6 |
| PAX5 → FCER2 | 125.8 |
| LEF1 → BACH2 | 119.9 |
| LEF1 → FHIT | 118.0 |
| BACH2 → SELL | 117.8 |
| FOXP3 → IL2RA | 117.4 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=3664 (41.9%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.04 | 0.0005 |
| 500 | 0.038 | 0.0022 |
| 1,000 | 0.027 | 0.0031 |
| 5,000 | 0.016 | 0.0091 |
| 10,000 | 0.0134 | 0.0153 |
| 832,664 | 0.0044 | 0.4187 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=5362 (5.5%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.02 | 0.0 |
| 500 | 0.026 | 0.0001 |
| 1,000 | 0.024 | 0.0002 |
| 5,000 | 0.0164 | 0.0008 |
| 10,000 | 0.0153 | 0.0016 |
| 832,664 | 0.0064 | 0.0554 |

---
