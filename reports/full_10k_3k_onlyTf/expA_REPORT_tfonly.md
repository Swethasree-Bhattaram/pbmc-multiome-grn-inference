# Experiment A — SCEMENT-integrated 3k+10k RNA (14,609 cells) (full-10k, TF-only regulators)

**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, 10k-RNA, 10k-ATAC)
**Reverse-imputeKNN reference:** SCEMENT-integrated 3k+10k RNA (14,609 cells)
**GRN:** Arboreto GRNBoost2 — regulators = tf_only.txt (816 present), targets = TRRUST (2,827)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets: rna3k (2,711), atac3k (2,711), rna10k (11,898), atac10k (11,898). The 10k RNA/ATAC are the SAME 11,898 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (29,218 x 30).
3. **Reverse-imputeKNN**: reverse-imputeKNN from SCEMENT-integrated 3k+10k RNA.
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
- Edges inferred: **829,982**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| GSC → B2M | 205.5 |
| EBF1 → PAX5 | 198.0 |
| PAX5 → EBF1 | 184.3 |
| PAX5 → MS4A1 | 181.8 |
| HDAC3 → DUX4 | 168.4 |
| LEF1 → NELL2 | 165.6 |
| TCF7L2 → CDKN1C | 159.5 |
| LEF1 → BACH2 | 159.2 |
| LEF1 → PRKCA | 154.5 |
| PAX5 → CD79A | 148.5 |
| SOX2 → B2M | 144.6 |
| EBF1 → MS4A1 | 135.2 |
| HOXC13 → B2M | 133.2 |
| PAX5 → CD22 | 132.2 |
| LEF1 → MLLT3 | 131.2 |
| OTX1 → B2M | 130.8 |
| MYBL1 → CCL5 | 130.6 |
| LEF1 → CCR7 | 129.4 |
| PAX5 → BLK | 127.5 |
| FOSB → JUN | 126.9 |
| TCF4 → RUNX2 | 126.1 |
| LEF1 → TCF7 | 123.9 |
| HEYL → DACT2 | 123.4 |
| LEF1 → FHIT | 121.2 |
| GSC → RPL3 | 120.3 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=3660 (41.8%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.05 | 0.0006 |
| 500 | 0.038 | 0.0022 |
| 1,000 | 0.024 | 0.0027 |
| 5,000 | 0.0162 | 0.0093 |
| 10,000 | 0.0136 | 0.0155 |
| 829,982 | 0.0044 | 0.4182 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=5269 (5.4%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.02 | 0.0 |
| 500 | 0.022 | 0.0001 |
| 1,000 | 0.02 | 0.0002 |
| 5,000 | 0.0144 | 0.0007 |
| 10,000 | 0.0145 | 0.0015 |
| 829,982 | 0.0063 | 0.0544 |

---
