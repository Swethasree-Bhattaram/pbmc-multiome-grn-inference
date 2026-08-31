# Experiment A — reverse-imputeKNN from a SCEMENT-integrated combined RNA reference

**Date:** 2026-08-31
**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, subsampled_3k-RNA, subsampled_3k-ATAC)
**Reverse-imputeKNN reference:** SCEMENT-combined RNA (3k+subsampled_3k, 5,422 cells)
**GRN:** Arboreto GRNBoost2 (TRRUST regulators)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets (2,711 cells each): rna3k (real), atac3k (peaks), rna_sub3k (subsampled 10k RNA, real), atac_sub3k (subsampled 10k peaks). The subsampled_3k RNA/ATAC are the SAME 2,711 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (10,844 x 30).
3. **Reverse-imputeKNN**: SCEMENT (AluruLab) batch-integrates the two real RNA datasets (3k + subsampled_3k) into a single 5,422-cell reference; both ATAC populations are imputed from that combined reference.
4. **All-cells matrix:** 10,844 x 36,601 (rows: rna3k, rna_sub3k, atac3k-imputed, atac_sub3k-imputed).
5. **GRN inference:** Arboreto GRNBoost2 with TRRUST TFs as regulators.
6. **Evaluation:** against PBMC-TRRUST and PBMC-Blood ground truth.

## 2. Integration (joint embedding H)

- Joint embedding **H**: 10,844 x 30 (blocks: rna3k, atac3k, rna_sub3k, atac_sub3k, each 2,711 x 30).
- Global alignment score: 0.8201 (4 datasets, 3k anchor).
- Pairwise scores vs anchor (rna3k): atac3k 0.8430, rna_sub3k 0.7876, atac_sub3k 0.8104.

## 3. Imputation (reverse-imputeKNN, k=20, softmax-of-negative-distance weights)

- Imputed expression: **36601 genes x 5422 ATAC cells**.
- Mean imputed expression: 0.0911; non-zero fraction: 0.490.
- All-cells matrix: **10844 cells x 36601 genes** (rna3k + rna_sub3k real + 5,422 imputed ATAC).

## 4. Arboreto GRNBoost2 results (TRRUST regulators)

- Regulators (TFs present in expression): 2827
- Target genes: 4389 (top-2000 HVG + TRRUST TFs)
- Edges inferred: **1,982,020**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| HLA-B → MT-CYB | 215.5 |
| B2M → MALAT1 | 201.4 |
| HLA-B → MT-ATP6 | 198.6 |
| HLA-B → B2M | 172.2 |
| HLA-B → MT-CO2 | 167.7 |
| RPL3 → EEF1A1 | 155.1 |
| HLA-B → EEF1A1 | 153.9 |
| CD8A → CD8B | 149.6 |
| CD8B → CD8A | 140.7 |
| RPL3 → B2M | 137.1 |
| B2M → RPL3 | 135.4 |
| ANGPT2 → MCPH1 | 130.5 |
| HLA-B → MT-ND4 | 129.6 |
| GNLY → PRF1 | 120.9 |
| B2M → HLA-B | 120.9 |
| B2M → MT-ND3 | 119.9 |
| RORA → RORA-AS1 | 117.6 |
| TPT1 → EEF1A1 | 115.9 |
| MCPH1 → ANGPT2 | 114.6 |
| GNLY → NKG7 | 114.3 |
| HLA-B → MT-ND1 | 109.5 |
| HLA-B → RPS27 | 108.3 |
| RPL3 → RPS27 | 107.7 |
| HLA-B → RPL41 | 103.2 |
| RS1 → BMP4 | 100.7 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=1690 (19.3%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0000 | 0.0000 |
| 500 | 0.0040 | 0.0002 |
| 1000 | 0.0060 | 0.0007 |
| 5000 | 0.0028 | 0.0016 |
| 10000 | 0.0025 | 0.0029 |
| 1982020 | 0.0009 | 0.1931 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=4891 (5.1%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0000 | 0.0000 |
| 500 | 0.0040 | 0.0000 |
| 1000 | 0.0050 | 0.0001 |
| 5000 | 0.0050 | 0.0003 |
| 10000 | 0.0044 | 0.0005 |
| 1982020 | 0.0025 | 0.0505 |

---
