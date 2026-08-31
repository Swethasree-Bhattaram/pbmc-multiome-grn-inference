# Experiment B1 — reverse-imputeKNN from the single 3k RNA reference

**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, subsampled_3k-RNA, subsampled_3k-ATAC)
**Reverse-imputeKNN reference:** 3k RNA only (2,711 cells)
**GRN:** Arboreto GRNBoost2 (TRRUST regulators)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets (2,711 cells each): rna3k, atac3k, rna_sub3k (subsampled from 10k), atac_sub3k (subsampled from 10k). The subsampled_3k RNA/ATAC are the SAME 2,711 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (10,844 x 30).
3. **Reverse-imputeKNN**: Only the original 3k RNA dataset (2,711 cells) is used as reference; both ATAC populations (3k + subsampled_3k) are imputed from it.
4. **All-cells matrix:** 10,844 x 36,601 (rows: rna3k, rna_sub3k, atac3k-imputed, atac_sub3k-imputed).
5. **GRN inference:** Arboreto GRNBoost2 with TRRUST TFs as regulators.
6. **Evaluation:** against PBMC-TRRUST and PBMC-Blood ground truth.

## 2. Integration (joint embedding H)

- Joint embedding **H**: 10,844 x 30 (blocks: rna3k, atac3k, rna_sub3k, atac_sub3k, each 2,711 x 30).
- Global alignment score: 0.8201 (4 datasets, 3k anchor).
- Pairwise scores vs anchor (rna3k): atac3k 0.8430, rna_sub3k 0.7876, atac_sub3k 0.8104.

## 3. Imputation (reverse-imputeKNN, k=20, softmax-of-negative-distance weights)

- Imputed expression: **36601 genes x 5422 ATAC cells**.
- Mean imputed expression: 0.0756; non-zero fraction: 0.264.
- All-cells matrix: **10844 cells x 36601 genes** (rna3k + rna_sub3k real + 5,422 imputed ATAC).

## 4. Arboreto GRNBoost2 results (TRRUST regulators)

- Regulators (TFs present in expression): 2827
- Target genes: 4390 (top-2000 HVG + TRRUST TFs)
- Edges inferred: **1,880,170**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| ANGPT2 → MCPH1 | 139.8 |
| CD8A → CD8B | 134.3 |
| CD8B → CD8A | 132.5 |
| MCPH1 → ANGPT2 | 105.2 |
| S100A9 → S100A8 | 100.3 |
| RORA → RORA-AS1 | 99.8 |
| PRF1 → GNLY | 99.2 |
| ACE → FOXI1 | 96.6 |
| B2M → HLA-B | 90.5 |
| S100A4 → S100A6 | 89.2 |
| BCL2 → CDK6 | 88.2 |
| GNLY → NKG7 | 85.9 |
| PAX5 → MS4A1 | 85.7 |
| FOSB → JUN | 84.3 |
| GNLY → PRF1 | 83.2 |
| DENND4A → RAB11A | 80.0 |
| HLA-B → B2M | 79.2 |
| B2M → HLA-C | 76.3 |
| RPS3A → EEF1A1 | 76.3 |
| B2M → TMSB4X | 75.2 |
| HLA-B → HLA-A | 74.3 |
| LYZ → S100A9 | 74.1 |
| S100A4 → LGALS1 | 72.3 |
| MS4A1 → BANK1 | 71.5 |
| B2M → MYL12A | 71.5 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=1739 (19.9%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0100 | 0.0001 |
| 500 | 0.0120 | 0.0007 |
| 1000 | 0.0070 | 0.0008 |
| 5000 | 0.0050 | 0.0029 |
| 10000 | 0.0033 | 0.0038 |
| 1880170 | 0.0009 | 0.1987 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=4555 (4.7%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0100 | 0.0000 |
| 500 | 0.0100 | 0.0001 |
| 1000 | 0.0050 | 0.0001 |
| 5000 | 0.0056 | 0.0003 |
| 10000 | 0.0042 | 0.0004 |
| 1880170 | 0.0024 | 0.0470 |

---
