# Multi-Dataset scSAGA Integration & GRN: TRRUST-only Target Comparison (full-10k)

## Design

Same four-dataset scSAGA integration and reverse-imputeKNN experiments as the main
comparison, but using the **full 10k multiome (11,898 cells)** instead of a 2,711-cell
subsample, and the Arboreto GRNBoost2 **target genes are restricted to ONLY the TRRUST
TFs present in the expression matrix** (2,827 genes) instead of top-2000 HVG + TRRUST TFs.

| Experiment | Reference for imputeKNN |
|---|---|
| **A** | SCEMENT-integrated 3k+10k RNA (14,609 cells) |
| **B1** | 3k RNA only (2,711 cells) |
| **B2** | 10k RNA only (11,898 cells) |

## GRN overview

| Experiment | All-cells matrix | Cells | Target genes | Inferred edges |
|---|---|---|---|---|
| expA | 29218 x 36601 | 29218 | 2827 | 1177622 |
| expB1 | 29218 x 36601 | 29218 | 2827 | 1095408 |
| expB2 | 29218 x 36601 | 29218 | 2827 | 1189398 |

## Edges recovered by GRN (fraction of deduplicated ground truth)

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 1957/8751 | 2642/96846 |
| expB1 | 1908/8751 | 2528/96846 |
| expB2 | 1939/8751 | 2694/96846 |

*Format: recovered edges / deduplicated ground-truth edges.*

## Precision@K

### K=100

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0100 | 0.0000 |
| expB1 | 0.0100 | 0.0100 |
| expB2 | 0.0000 | 0.0100 |

### K=1000

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0100 | 0.0080 |
| expB1 | 0.0110 | 0.0110 |
| expB2 | 0.0100 | 0.0080 |

### K=5000

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0054 | 0.0048 |
| expB1 | 0.0052 | 0.0050 |
| expB2 | 0.0054 | 0.0060 |

---

## Experiment A: Reverse-imputeKNN from a SCEMENT-integrated combined RNA reference

**Reference:** SCEMENT-integrated 3k+10k RNA (14,609 cells)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1957 (22.4%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=2642 (2.7%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0100/0.0001 | 0.0000/0.0 |
| 500 | 0.0100/0.0006 | 0.0040/0.0 |
| 1000 | 0.0100/0.0011 | 0.0080/0.0001 |
| 5000 | 0.0054/0.0031 | 0.0048/0.0002 |
| 10000 | 0.0046/0.0053 | 0.0041/0.0004 |
| 1177622 | 0.0017/0.2236 | 0.0022/0.0273 |

---

## Experiment B1: Reverse-imputeKNN from the single 3k RNA reference

**Reference:** 3k RNA only (2,711 cells)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1908 (21.8%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=2528 (2.6%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0100/0.0001 | 0.0100/0.0 |
| 500 | 0.0120/0.0007 | 0.0080/0.0 |
| 1000 | 0.0110/0.0013 | 0.0110/0.0001 |
| 5000 | 0.0052/0.003 | 0.0050/0.0003 |
| 10000 | 0.0057/0.0065 | 0.0053/0.0005 |
| 1095408 | 0.0017/0.218 | 0.0023/0.0261 |

---

## Experiment B2: Reverse-imputeKNN from the single 10k RNA reference

**Reference:** 10k RNA only (11,898 cells)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1939 (22.2%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=2694 (2.8%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0000/0.0 | 0.0100/0.0 |
| 500 | 0.0140/0.0008 | 0.0100/0.0001 |
| 1000 | 0.0100/0.0011 | 0.0080/0.0001 |
| 5000 | 0.0054/0.0031 | 0.0060/0.0003 |
| 10000 | 0.0050/0.0057 | 0.0051/0.0005 |
| 1189398 | 0.0016/0.2216 | 0.0023/0.0278 |

---
