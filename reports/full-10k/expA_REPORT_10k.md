# Experiment A — reverse-imputeKNN from a SCEMENT-integrated combined RNA reference (full-10k)

**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, 10k-RNA, 10k-ATAC)
**Reverse-imputeKNN reference:** SCEMENT-combined RNA (3k+10k, 14,609 cells)
**GRN:** Arboreto GRNBoost2 (TRRUST regulators, TRRUST-only targets)
**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)

---

## 1. Pipeline

1. **Datasets:** 4 single-modality datasets: rna3k (2,711), atac3k (2,711), rna10k (11,898), atac10k (11,898). The 10k RNA/ATAC are the SAME 11,898 physical cells (paired multiome).
2. **Integration:** scSAGA into a shared joint embedding **H** (29,218 x 30).
3. **Reverse-imputeKNN**: reverse-imputeKNN from a SCEMENT-integrated combined RNA reference.
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
- Edges inferred: **1,177,622**

**Top 25 edges:**

| TF → target | importance |
|---|---|
| ANGPT2 → MCPH1 | 242.5 |
| MCPH1 → ANGPT2 | 233.0 |
| CD8A → CD8B | 213.9 |
| CD8B → CD8A | 207.1 |
| RPL3 → B2M | 194.6 |
| MECOM → TBX4 | 178.0 |
| HLA-A → B2M | 137.7 |
| B2M → HLA-B | 136.5 |
| LEF1 → PRKCA | 127.3 |
| GNLY → PRF1 | 126.9 |
| FOSB → JUN | 124.6 |
| PTPN13 → TBX4 | 118.9 |
| B2M → LTB | 117.5 |
| WFS1 → TBX4 | 116.2 |
| LEF1 → BACH2 | 115.1 |
| PRF1 → GNLY | 113.1 |
| B2M → HLA-E | 111.7 |
| CD8A → KLRK1 | 110.7 |
| RPL3 → TPT1 | 110.6 |
| B2M → HLA-A | 109.8 |
| B2M → HLA-C | 106.8 |
| HLA-B → HLA-C | 104.8 |
| BCL2 → CDK6 | 101.1 |
| RPL3 → RPLP0 | 101.0 |
| KLK3 → B2M | 100.1 |

## 5. Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, deduplicated=8751, recovered by GRN=1957 (22.4%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0100 | 0.0001 |
| 500 | 0.0100 | 0.0006 |
| 1000 | 0.0100 | 0.0011 |
| 5000 | 0.0054 | 0.0031 |
| 10000 | 0.0046 | 0.0053 |
| 1177622 | 0.0017 | 0.2236 |

**PBMC-Blood:** total=98338, deduplicated=96846, recovered by GRN=2642 (2.7%)

| Top-K | Precision@K | Recall@K |
|---|---|---|
| 100 | 0.0000 | 0.0000 |
| 500 | 0.0040 | 0.0000 |
| 1000 | 0.0080 | 0.0001 |
| 5000 | 0.0048 | 0.0002 |
| 10000 | 0.0041 | 0.0004 |
| 1177622 | 0.0022 | 0.0273 |

---
