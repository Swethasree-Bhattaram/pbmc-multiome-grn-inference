# Experiment B2 — reverse-imputeKNN from the single 10k RNA reference (full-10k)

**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, 10k-RNA, 10k-ATAC)
**Reverse-imputeKNN reference:** 10k RNA only (11,898 cells)
**GRN:** Arboreto GRNBoost2 (TRRUST regulators, TRRUST-only targets)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets: rna3k (2,711), atac3k (2,711), rna10k (11,898), atac10k (11,898). The 10k RNA/ATAC are the SAME 11,898 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (29,218 x 30).
3. **Reverse-imputeKNN**: reverse-imputeKNN from the single 10k RNA reference.
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
- Edges inferred: **1,189,398**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| ANGPT2 → MCPH1 | 239.5 |
| MCPH1 → ANGPT2 | 211.8 |
| MECOM → TBX4 | 203.2 |
| CD8B → CD8A | 203.1 |
| CD8A → CD8B | 194.6 |
| FOSB → JUN | 151.6 |
| WFS1 → TBX4 | 149.0 |
| PTPN13 → TBX4 | 137.6 |
| B2M → HLA-B | 135.4 |
| LEF1 → PRKCA | 126.8 |
| BCL2 → CDK6 | 123.1 |
| PRF1 → GNLY | 112.3 |
| B2M → HLA-A | 112.1 |
| HLA-B → HLA-A | 110.2 |
| GNLY → PRF1 | 108.2 |
| LEF1 → CCR7 | 106.7 |
| S100A4 → S100A6 | 106.3 |
| CD8A → KLRK1 | 104.7 |
| HLA-DRA → HLA-DQA1 | 103.6 |
| SCG3 → TBX4 | 101.5 |
| CCL5 → A2M | 100.0 |
| LYZ → S100A9 | 99.8 |
| CD74 → HLA-DPB1 | 99.1 |
| ZEB2 → LTB | 98.0 |
| CD74 → HLA-DPA1 | 97.6 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=1939 (22.2%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0000 | 0.0000 |
| 500 | 0.0140 | 0.0008 |
| 1000 | 0.0100 | 0.0011 |
| 5000 | 0.0054 | 0.0031 |
| 10000 | 0.0050 | 0.0057 |
| 1189398 | 0.0016 | 0.2216 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=2694 (2.8%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0100 | 0.0000 |
| 500 | 0.0100 | 0.0001 |
| 1000 | 0.0080 | 0.0001 |
| 5000 | 0.0060 | 0.0003 |
| 10000 | 0.0051 | 0.0005 |
| 1189398 | 0.0023 | 0.0278 |

---
