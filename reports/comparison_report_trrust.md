# Multi-Dataset scSAGA Integration & GRN: TRRUST-only Target Comparison

## Design

Same four-dataset scSAGA integration and reverse-imputeKNN experiments as the main
comparison, but the Arboreto GRNBoost2 **target genes are restricted to ONLY the TRRUST
TFs present in the expression matrix** (2,827 genes) instead of top-2000 HVG + TRRUST TFs
(~4,389 genes). This makes the GRN inference much faster and focuses the network on
TF-to-TF regulatory edges.

| Experiment | Reference for imputeKNN |
|---|---|
| **A** | SCEMENT-integrated 3k+subsampled_3k RNA (5,422 cells) |
| **B1** | 3k RNA only (2,711 cells) |
| **B2** | subsampled_3k RNA only (2,711 cells) |

## GRN overview

| Experiment | All-cells matrix | Cells | Target genes | Inferred edges |
|---|---|---|---|---|
| expA | 10844 x 36601 | 10844 | 2827 | 905890 |
| expB1 | 10844 x 36601 | 10844 | 2827 | 898699 |
| expB2 | 10844 x 36601 | 10844 | 2827 | 920605 |

## Edges recovered by GRN (fraction of deduplicated ground truth)

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 1603/8751 | 2117/96846 |
| expB1 | 1657/8751 | 2105/96846 |
| expB2 | 1652/8751 | 2177/96846 |

*Format: recovered edges / deduplicated ground-truth edges.*

## Precision@K

### K=100

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0200 | 0.0100 |
| expB1 | 0.0200 | 0.0000 |
| expB2 | 0.0100 | 0.0000 |

### K=1000

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0100 | 0.0060 |
| expB1 | 0.0100 | 0.0090 |
| expB2 | 0.0100 | 0.0080 |

### K=5000

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0054 | 0.0040 |
| expB1 | 0.0064 | 0.0052 |
| expB2 | 0.0050 | 0.0048 |

---

## Experiment A: Reverse-imputeKNN from a SCEMENT-integrated combined RNA reference

**Reference:** SCEMENT-integrated 3k+subsampled_3k RNA (5,422 cells)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1603 (18.3%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=2117 (2.2%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0200/0.0002 | 0.0100/0.0 |
| 500 | 0.0140/0.0008 | 0.0080/0.0 |
| 1000 | 0.0100/0.0011 | 0.0060/0.0001 |
| 5000 | 0.0054/0.0031 | 0.0040/0.0002 |
| 10000 | 0.0049/0.0056 | 0.0038/0.0004 |
| 905890 | 0.0018/0.1832 | 0.0023/0.0219 |

---

## Experiment B1: Reverse-imputeKNN from the single 3k RNA reference

**Reference:** 3k RNA only (2,711 cells)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1657 (18.9%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=2105 (2.2%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0200/0.0002 | 0.0000/0.0 |
| 500 | 0.0140/0.0008 | 0.0100/0.0001 |
| 1000 | 0.0100/0.0011 | 0.0090/0.0001 |
| 5000 | 0.0064/0.0037 | 0.0052/0.0003 |
| 10000 | 0.0050/0.0057 | 0.0044/0.0005 |
| 898699 | 0.0018/0.1893 | 0.0023/0.0217 |

---

## Experiment B2: Reverse-imputeKNN from the single subsampled_3k RNA reference

**Reference:** subsampled_3k RNA only (2,711 cells)

**PBMC-TRRUST:** total=9720, dedup=8751, recovered=1652 (18.9%)
**PBMC-Blood:** total=98338, dedup=96846, recovered=2177 (2.2%)

| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |
|---|---|---|
| 100 | 0.0100/0.0001 | 0.0000/0.0 |
| 500 | 0.0140/0.0008 | 0.0100/0.0001 |
| 1000 | 0.0100/0.0011 | 0.0080/0.0001 |
| 5000 | 0.0050/0.0029 | 0.0048/0.0002 |
| 10000 | 0.0044/0.005 | 0.0040/0.0004 |
| 920605 | 0.0018/0.1888 | 0.0024/0.0225 |

---
