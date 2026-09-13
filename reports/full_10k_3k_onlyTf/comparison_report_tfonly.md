# Multi-Dataset scSAGA Integration & GRN: TF-only Regulators Comparison (full-10k)

## Design

Same four-dataset scSAGA integration and reverse-imputeKNN experiments as the TRRUST-only
full-10k runs, but the Arboreto GRNBoost2 **regulators are restricted to ONLY the tf_only.txt
TFs present in the expression matrix (816)** instead of the full TRRUST regulator set (2,827).
Target genes remain the TRRUST TFs present (2,827).
This isolates the effect of restricting regulators to true transcription factors.

| Experiment | Reference for imputeKNN |
|---|---|
| **A** | SCEMENT-integrated 3k+10k RNA (14,609 cells) |
| **B1** | 3k RNA only (2,711 cells) |
| **B2** | 10k RNA only (11,898 cells) |

## GRN overview

| Experiment | Cells | Regulators (tf_only) | Targets (TRRUST) | Union cols | Inferred edges |
|---|---|---|---|---|---|
| expA | 29218 | 816 | 2827 | 2852 | 829982 |
| expB1 | 29218 | 816 | 2827 | 2852 | 779701 |
| expB2 | 29218 | 816 | 2827 | 2852 | 832664 |

## Edges recovered by GRN (fraction of deduplicated ground truth)

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 3660/8751 (41.8%) | 5269/96846 (5.4%) |
| expB1 | 3454/8751 (39.5%) | 5132/96846 (5.3%) |
| expB2 | 3664/8751 (41.9%) | 5362/96846 (5.5%) |

*Format: recovered edges / deduplicated ground-truth edges (%).*

## Precision@K

### K=100

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.05 | 0.02 |
| expB1 | 0.05 | 0.03 |
| expB2 | 0.04 | 0.02 |

### K=1000

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.024 | 0.02 |
| expB1 | 0.026 | 0.023 |
| expB2 | 0.027 | 0.024 |

### K=5000

| Experiment | PBMC-TRRUST | PBMC-Blood |
|---|---|---|
| expA | 0.0162 | 0.0144 |
| expB1 | 0.018 | 0.0164 |
| expB2 | 0.016 | 0.0164 |

## Comparison vs TRRUST-only (full-10k): effect of TF-only regulators

Restricting regulators from 2,827 (TRRUST, mixed targets+TFs) to 816 (tf_only) roughly
doubles ground-truth edge recovery while cutting inferred edges ~30%:

| Experiment | PBMC-TRRUST recovered (trrust regs) | PBMC-TRRUST recovered (tfonly regs) |
|---|---|---|
| A  | 1,957 / 8,751 (22.4%) | 3,660 / 8,751 (41.8%) |
| B1 | 1,908 / 8,751 (21.8%) | 3,454 / 8,751 (39.5%) |
| B2 | 1,939 / 8,751 (22.2%) | 3,664 / 8,751 (41.9%) |

*TRRUST-regulator baseline values from `comparison_report_trrust_10k.md`.*

## Notes / caveat

The tf_only.txt regulator list was curated from the TF columns of these same PBMC
ground-truth files. Restricting regulators to benchmark-defined TFs removes non-TF
targets from the predictor set, matching the ground-truth schema (all edges originate
at a TF). The recovered-edge and precision gains should therefore be interpreted as
"regulator set aligned to benchmark TFs", which partly reflects feature selection
using the evaluation key itself.

---
