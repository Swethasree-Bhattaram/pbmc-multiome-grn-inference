# Experiment B1 — reverse-imputeKNN from the single 3k RNA reference (full-10k)

**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, 10k-RNA, 10k-ATAC)
**Reverse-imputeKNN reference:** 3k RNA only (2,711 cells)
**GRN:** Arboreto GRNBoost2 (TRRUST regulators, TRRUST-only targets)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets: rna3k (2,711), atac3k (2,711), rna10k (11,898), atac10k (11,898). The 10k RNA/ATAC are the SAME 11,898 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (29,218 x 30).
3. **Reverse-imputeKNN**: reverse-imputeKNN from the single 3k RNA reference.
4. **All-cells matrix:** 29,218 x 36,601 (rows: rna3k, rna10k, atac3k-imputed, atac10k-imputed).
5. **GRN inference:** Arboreto GRNBoost2 with TRRUST TFs as regulators, TRRUST-only targets.
6. **Evaluation:** against PBMC-TRRUST and PBMC-Blood ground truth.

## 2. Integration (joint embedding H)

- Joint embedding **H**: 29,218 x 30 (blocks: rna3k, atac3k, rna10k, atac10k).
- Global alignment score: 0.0984 (4 datasets, 3k anchor).
- Pairwise scores vs anchor (rna3k): atac3k 0.4027, rna10k 0.2240, atac10k 0.1236.

## 3. Imputation (reverse-imputeKNN, k=20, softmax-of-negative-distance weights)

- Imputed expression: **36601 genes x 14609 ATAC cells**.
- All-cells matrix: **29218 cells x 36601 genes** (rna3k + rna10k real + 14,609 imputed ATAC).

## 4. Arboreto GRNBoost2 results (TRRUST regulators, TRRUST-only targets)

- Regulators (TFs present in expression): 2827
- Target genes: 2827 (TRRUST-only)
- Edges inferred: **1,095,408**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| MCPH1 → ANGPT2 | 232.3 |
| ANGPT2 → MCPH1 | 188.0 |
| CD8A → CD8B | 187.9 |
| CD8B → CD8A | 181.4 |
| FOSB → JUN | 145.0 |
| BCL2 → CDK6 | 141.2 |
| LEF1 → PRKCA | 110.5 |
| B2M → HLA-B | 110.3 |
| JCHAIN → SDC1 | 109.5 |
| B2M → HLA-A | 107.3 |
| GNLY → PRF1 | 107.2 |
| LEF1 → FHIT | 103.4 |
| B2M → HLA-C | 100.7 |
| MS4A1 → PAX5 | 97.8 |
| HLA-B → HLA-A | 96.1 |
| TNC → DEC1 | 94.0 |
| PRF1 → GNLY | 93.4 |
| LEF1 → CCR7 | 93.1 |
| HLA-DRA → HLA-DQA1 | 93.1 |
| CD74 → HLA-DQA1 | 92.2 |
| LEF1 → NELL2 | 90.6 |
| HLA-B → HLA-C | 89.9 |
| S100A4 → S100A6 | 88.3 |
| CD36 → VCAN | 88.0 |
| IFITM3 → CDKN1C | 87.8 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=1908 (21.8%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0100 | 0.0001 |
| 500 | 0.0120 | 0.0007 |
| 1000 | 0.0110 | 0.0013 |
| 5000 | 0.0052 | 0.0030 |
| 10000 | 0.0057 | 0.0065 |
| 1095408 | 0.0017 | 0.2180 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=2528 (2.6%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0100 | 0.0000 |
| 500 | 0.0080 | 0.0000 |
| 1000 | 0.0110 | 0.0001 |
| 5000 | 0.0050 | 0.0003 |
| 10000 | 0.0053 | 0.0005 |
| 1095408 | 0.0023 | 0.0261 |

---
