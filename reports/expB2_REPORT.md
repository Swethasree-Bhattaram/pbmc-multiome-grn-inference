# Experiment B2 — reverse-imputeKNN from the single subsampled_3k RNA reference

**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, subsampled_3k-RNA, subsampled_3k-ATAC)
**Reverse-imputeKNN reference:** subsampled_3k RNA only (2,711 cells)
**GRN:** Arboreto GRNBoost2 (TRRUST regulators)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets (2,711 cells each): rna3k, atac3k, rna_sub3k (subsampled from 10k), atac_sub3k (subsampled from 10k). The subsampled_3k RNA/ATAC are the SAME 2,711 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (10,844 x 30).
3. **Reverse-imputeKNN**: Only the subsampled_3k RNA dataset (2,711 cells, subsampled from the 10k) is used as reference; both ATAC populations (3k + subsampled_3k) are imputed from it.
4. **All-cells matrix:** 10,844 x 36,601 (rows: rna3k, rna_sub3k, atac3k-imputed, atac_sub3k-imputed).
5. **GRN inference:** Arboreto GRNBoost2 with TRRUST TFs as regulators.
6. **Evaluation:** against PBMC-TRRUST and PBMC-Blood ground truth.

## 2. Integration (joint embedding H)

- Joint embedding **H**: 10,844 x 30 (blocks: rna3k, atac3k, rna_sub3k, atac_sub3k, each 2,711 x 30).
- Global alignment score: 0.8201 (4 datasets, 3k anchor).
- Pairwise scores vs anchor (rna3k): atac3k 0.8430, rna_sub3k 0.7876, atac_sub3k 0.8104.

## 3. Imputation (reverse-imputeKNN, k=20, softmax-of-negative-distance weights)

- Imputed expression: **36601 genes x 5422 ATAC cells**.
- Mean imputed expression: 0.0757; non-zero fraction: 0.251.
- All-cells matrix: **10844 cells x 36601 genes** (rna3k + rna_sub3k real + 5,422 imputed ATAC).

## 4. Arboreto GRNBoost2 results (TRRUST regulators)

- Regulators (TFs present in expression): 2827
- Target genes: 4389 (top-2000 HVG + TRRUST TFs)
- Edges inferred: **1,932,030**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| CD8B → CD8A | 159.3 |
| ANGPT2 → MCPH1 | 149.6 |
| CD8A → CD8B | 147.7 |
| LTC4S → PLIN1 | 138.9 |
| RS1 → BMP4 | 133.9 |
| POLD1 → BMP4 | 111.0 |
| RORA → RORA-AS1 | 101.8 |
| S100A9 → S100A8 | 99.5 |
| MCPH1 → ANGPT2 | 99.1 |
| PHLDA1 → PLIN1 | 95.8 |
| RPS3A → EEF1A1 | 94.4 |
| B2M → HLA-B | 92.8 |
| ITGB3 → SNCB | 92.1 |
| GNLY → PRF1 | 91.9 |
| PRF1 → GNLY | 89.6 |
| RPL3 → RPLP0 | 85.5 |
| FOSB → JUN | 85.1 |
| CCL5 → GZMA | 82.7 |
| RPL3 → RPS3A | 80.7 |
| MS4A1 → PAX5 | 80.1 |
| RPL3 → EEF1A1 | 78.1 |
| BCL2 → CDK6 | 78.0 |
| PAX5 → MS4A1 | 77.9 |
| RPS3A → RPS8 | 77.6 |
| PIK3R3 → PLIN1 | 77.0 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=1721 (19.7%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0000 | 0.0000 |
| 500 | 0.0060 | 0.0003 |
| 1000 | 0.0070 | 0.0008 |
| 5000 | 0.0036 | 0.0021 |
| 10000 | 0.0029 | 0.0033 |
| 1932030 | 0.0009 | 0.1967 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=4721 (4.9%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0000 | 0.0000 |
| 500 | 0.0080 | 0.0000 |
| 1000 | 0.0090 | 0.0001 |
| 5000 | 0.0080 | 0.0004 |
| 10000 | 0.0059 | 0.0006 |
| 1932030 | 0.0024 | 0.0487 |

---
