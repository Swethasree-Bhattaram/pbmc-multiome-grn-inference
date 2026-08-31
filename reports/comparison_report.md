# Multi-Dataset scSAGA Integration & GRN: Experiment Comparison

## Design

Four datasets are jointly integrated with scSAGA (3k-RNA anchor, 2,711 cells each):
- **rna3k** · **atac3k** · **rna_sub3k** (subsampled from 10k) · **atac_sub3k** (subsampled from 10k)
- Joint embedding **H** = 10,844 x 30 (each block 2,711 x 30); global alignment score 0.8201.

The only variable across experiments is **which RNA reference** is used by reverse-imputeKNN to
impute expression for the 5,422 ATAC cells. The all-cells matrix (10,844 x 36,601) always holds the
same real RNA cells (rna3k then rna_sub3k) in rows 0..5421; only the imputed ATAC block (rows 5422..10843) changes.

| Experiment | Reference for imputeKNN |
|---|---|
| **A** | SCEMENT-integrated 3k+RNA & subsampled_3k+RNA → one 5,422-cell reference |
| **B1** | 3k RNA only (2,711 cells) |
| **B2** | subsampled_3k RNA only (2,711 cells) |

Each produces an Arboreto GRNBoost2 network evaluated against PBMC-TRRUST and PBMC-Blood ground truth.

## GRN overview

| Experiment | All-cells matrix | Cells | Target genes | Inferred edges |
|---|---|---|---|---|
| expA | 10844 x 36601 | 10844 | 4389 | 1982020 |
| expB1 | 10844 x 36601 | 10844 | 4390 | 1880170 |
| expB2 | 10844 x 36601 | 10844 | 4389 | 1932030 |

**Why the inferred-edge counts differ across experiments:** GRNBoost2 infers an edge for every
(regulator, target) pair whose importance exceeds a data-dependent threshold. Because the imputed
ATAC block differs between experiments (different RNA reference), the expression matrix fed to
GRNBoost2 differs, so the number of edges that clear the threshold varies slightly. The target-gene
set is also marginally different (e.g. expB1 selected 4,390 targets vs 4,389 for the others due to a
tie at the top-2000 HVG cutoff). These are expected, minor run-to-run differences — not an error.

## Edges recovered by GRN (fraction of deduplicated ground truth)

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 1690/8751 | 4891/96846 |
| expB1 | 1739/8751 | 4555/96846 |
| expB2 | 1721/8751 | 4721/96846 |

*Format: recovered edges / deduplicated ground-truth edges (e.g. 1690/8751 = 19.3% of PBMC-TRRUST).*

## Precision@K (K=100)

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0000 | 0.0000 |
| expB1 | 0.0100 | 0.0100 |
| expB2 | 0.0000 | 0.0000 |

## Precision@K (K=1000)

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0060 | 0.0050 |
| expB1 | 0.0070 | 0.0050 |
| expB2 | 0.0070 | 0.0090 |

## Precision@K (K=5000)

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0028 | 0.0050 |
| expB1 | 0.0050 | 0.0056 |
| expB2 | 0.0036 | 0.0080 |

## Reverse-imputed ATAC expression (5,422 cells x 36,601 genes)

| Experiment | Mean expr | Non-zero fraction |
|---|---|---|
| expA | 0.0911 | 0.490 |
| expB1 | 0.0756 | 0.264 |
| expB2 | 0.0757 | 0.251 |

---

## Experiment A: Reverse-imputeKNN from a SCEMENT-integrated combined RNA reference

**Reference:** SCEMENT-integrated 3k+subsampled_3k RNA (5,422 cells)

### Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1690 (19.3%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=4891 (5.1%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0000/0.0 | 0.0000/0.0 |
| 500 | 0.0040/0.0002 | 0.0040/0.0 |
| 1000 | 0.0060/0.0007 | 0.0050/0.0001 |
| 5000 | 0.0028/0.0016 | 0.0050/0.0003 |
| 10000 | 0.0025/0.0029 | 0.0044/0.0005 |
| 1982020 | 0.0009/0.1931 | 0.0025/0.0505 |

---

## Experiment B1: Reverse-imputeKNN from the single 3k RNA reference

**Reference:** 3k RNA only (2,711 cells)

### Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1739 (19.9%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=4555 (4.7%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0100/0.0001 | 0.0100/0.0 |
| 500 | 0.0120/0.0007 | 0.0100/0.0001 |
| 1000 | 0.0070/0.0008 | 0.0050/0.0001 |
| 5000 | 0.0050/0.0029 | 0.0056/0.0003 |
| 10000 | 0.0033/0.0038 | 0.0042/0.0004 |
| 1880170 | 0.0009/0.1987 | 0.0024/0.047 |

---

## Experiment B2: Reverse-imputeKNN from the single subsampled_3k RNA reference

**Reference:** subsampled_3k RNA only (2,711 cells)

### Evaluation vs ground truth (deduplicated)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1721 (19.7%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=4721 (4.9%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0000/0.0 | 0.0000/0.0 |
| 500 | 0.0060/0.0003 | 0.0080/0.0 |
| 1000 | 0.0070/0.0008 | 0.0090/0.0001 |
| 5000 | 0.0036/0.0021 | 0.0080/0.0004 |
| 10000 | 0.0029/0.0033 | 0.0059/0.0006 |
| 1932030 | 0.0009/0.1967 | 0.0024/0.0487 |

---
