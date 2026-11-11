# report_extrapolated: adding unpaired PBMC datasets to multi11k

Everything in this file is a **projection** unless a row is marked *measured*. Basis for every number:
the measured per-cell cost of 1.585 ms/cell/epoch (8 CPU threads, G+P = 126,397) and the measured
scSAGA timing breakdown (`results/integration_10k/saga_runtimes.txt`). No run of this configuration
has been executed.

## 1. What is being added

| dataset | modality | cells | features | paired with anchor | enters |
|---|---|---|---|---|---|
| rna11k | RNA | 11,000 | 36,601 genes | yes (anchor) | MultiVI + scSAGA + GRN |
| atac11k | ATAC | 11,000 | 89,796 peaks | yes | MultiVI + scSAGA (imputed) |
| rna11k_unp | RNA | 11,000 | ~36,601 genes | no | MultiVI + scSAGA + GRN |
| atac11k_unp | ATAC | 11,000 | ~89,796 peaks | no | MultiVI + scSAGA (imputed) |
| rna3k | RNA | 3,000 | ~36,601 genes | no | MultiVI + scSAGA + GRN |
| atac3k | ATAC | 3,000 | ~89,796 peaks | no | MultiVI + scSAGA (imputed) |
| atac8161 | ATAC | 8,161 | 87,863 peaks | no | MultiVI + scSAGA (imputed) |

Totals, with the 11,000 co-assayed barcodes collapsed to one row each (MultiVI row space):

```
rows                    47,161     (11,000 paired + 14,000 RNA-only + 22,161 ATAC-only)
paired fraction         11,000/47,161 = 23.3%      (now: 11,000/11,000 = 100%)
ATAC cells to impute    33,161     (now: 11,000)
RNA cells for the GRN   25,000     (now: 11,000)
GRN matrix rows         58,161     (now: 22,000)
```

## 2. Assumptions

Every projected number below depends on these. Each one is a knob, not a fact.

| id | assumption | effect if wrong |
|---|---|---|
| A1 | gene axis is the union of the RNA gene sets = 36,601 (same annotation) | a different annotation adds 0-2,000 genes, ~2% of cost |
| A2 | peak axis is the **union**: 89,796 ∪ 87,863 ∪ 11k/3k peak sets = 120k-178k, base case **160,000** | halves the MultiVI cost if the sets are intersected instead (see §6) |
| A3 | h=128, L=30, 200 epochs, 2 fits (full + held-out excluded), batch 128, as in multi11k | early stopping cuts the fit time, not the per-epoch cost |
| A4 | MultiVI gets **no** batch/technology covariate, as now | a covariate would absorb the technology axis but needs batch labels |
| A5 | per-cell cost is linear in (G+P) and in epochs — the measured slope, extrapolated | sublinear only if the first layer dominates less at wider P |
| A6 | scSAGA cost ∝ Σ N_anchor·N_query over the pair solves, one pair at a time | joint refinement across 7 blocks may add a superlinear term |
| A7 | the 8,161 × 87,863 ATAC set gets a 50-d PCA computed first (scSAGA never reads peaks) | PCA on 8,161 × 87,863 sparse is minutes, not hours |

## 3. Cost projection

### MultiVI

| quantity | now (multi11k) | projected (P=160,000) | low (P=120,000) | high (P=178,000) |
|---|---|---|---|---|
| parameters | 58,441,506 *measured* | 85,540,250 | 70,100,250 | 92,488,250 |
| model + Adam | 0.94 GB | 1.37 GB | 1.12 GB | 1.48 GB |
| ms/cell/epoch | 1.585 *measured* | 2.47 | 1.96 | 2.69 |
| s/epoch at 47,161 rows | 17.4 *derived* | 116 | 93 | 127 |
| one 200-epoch fit | 58 min | 6.5 h | 5.2 h | 7.1 h |
| two fits | 1.9 h | 12.9 h | 10.3 h | 14.1 h |
| imputed matrix export | 1.6 GB | 4.9 GB | 4.9 GB | 4.9 GB |
| minibatch activations | 65 MB | 101 MB | 82 MB | 111 MB |

The 6.5x jump in fit time is 4.3x more rows times 1.56x wider feature axis. It is not a surprise
result, it is arithmetic — but it is the difference between a 2 h job and an overnight job.

### scSAGA

| quantity | now | projected |
|---|---|---|
| pair solves | 1 (atac11k → rna11k) | 6 (five ATAC + one RNA block) |
| Σ N_anchor·N_query | 1.21e8 | 5.19e8 (4.29x) |
| integration wall time | ~183 s *derived* | ~13 min |
| peak memory, largest pair | 11k x 11k = 484 MB per matrix | unchanged (pairs run one at a time) |
| cost of the 87,863-peak block | — | zero in peaks; PCA only |

scSAGA is indifferent to the peak axis: it never reads the peak matrix. The 8,161 x 87,863 dataset is
the cheapest thing in this table for scSAGA and the most expensive thing for MultiVI.

### GRN (Arboreto / GRNBoost2)

| quantity | now | projected |
|---|---|---|
| matrix rows | 22,000 | 58,161 (2.64x) |
| dense input | 3.2 GB | 8.5 GB |
| wall time | t | ~2.6x t (linear in rows at fixed TF count and trees) |
| imputed rows in the matrix | 11,000/22,000 (50%) | 33,161/58,161 (57%) |

## 4. Quality projection, step by step

### Step 1 — integration

| metric | now | projected | basis |
|---|---|---|---|
| MultiVI FOSCTTM (co-assayed block) | 0.008 *measured* | 0.02-0.06, up to 0.10 if a technology axis forms | paired fraction 100% -> 23.3%; the paper's 1-99% unpairing sweep degrades monotonically |
| scSAGA FOSCTTM (co-assayed block) | 0.241 *measured* | 0.22-0.28, i.e. unchanged | the rna11k<->atac11k pair solve is untouched by the new blocks; only the joint refinement shifts it |
| FOSCTTM on the 5 new blocks | — | **undefined** | requires a cell measured in both modalities; these cells have one |
| ARI / NMI / label transfer / isolated F1 | n/a | n/a unless `labels.txt` is supplied | unchanged blocker, not a scaling effect |
| iLISI / kBET / batch ASW | 0.269 / 0.115 / 0.015 (scSAGA) | **not comparable** | batch count goes 2 -> 7, so the scale of iLISI changes; must be recomputed per block |

The 0.008 -> 0.02-0.06 range is a degradation of the *measurement you currently report*, not
necessarily of usefulness: 0.06 is still far inside the 0.14-0.22 "good" band quoted in `report.md`.
The bigger risk is A4: the standalone 3' scRNA-seq and the multiome RNA are different technologies with
different depth profiles. Without a batch covariate, that difference is absorbed into `z` and shows up as
a technology axis, which is what pushes the projection to 0.10 rather than 0.06.

### Step 2 — imputation

| metric | now | projected | basis |
|---|---|---|---|
| MSE, MultiVI | 0.172 | 0.14-0.17 | more training cells -> better amortised posterior; still above control 0.097 |
| MSE, scSAGA | 0.124 | 0.11-0.124 | 25,000 reference cells instead of 11,000 -> finer local averages; control is 0.124 |
| densification, MultiVI | 22.5x | 15-25x | unchanged mechanism, slightly diluted by denser training input |
| densification, scSAGA | 6.4x | 5-7x | weighted average of real cells; bounded by the reference's own density |
| gene mean shift, MultiVI | 0.336 sd | 0.28-0.34 sd | |
| gene mean shift, scSAGA | 0.013 sd | ~0.01 sd | averaging preserves marginals by construction |
| MSE / rho on the new ATAC blocks | — | **unverifiable** | no held-out RNA exists for cells that were never assayed for RNA |

The step-2 gate (`MSE < control`) is projected to **fail for both methods again**. Neither mechanism
changes with cell count: MultiVI decodes a smooth latent (over-dense, mean-shifted), scSAGA averages
neighbours (marginal-preserving, cell-agnostic). Adding cells sharpens both slightly without changing
which side of the control they land on.

For the 22,161 new ATAC cells the only available checks are the ones that need no truth: pFCS against
the gene-permuted null, densification vs the reference, per-cluster marker pseudobulk, and the gene mean
shift. Say so explicitly in the run report; do not report an MSE for them.

### Step 3 — GRN

| metric | now | projected | basis |
|---|---|---|---|
| cells in the matrix | 22,000 | 58,161 | 2.64x |
| seed Jaccard, top-10k, MultiVI | 0.546 | 0.60-0.70 | more rows -> more stable trees |
| seed Jaccard, top-10k, scSAGA | 0.405 | 0.45-0.55 | |
| AUPRC vs TRRUST, MultiVI | 0.006 | 0.008-0.015 | prevalence 0.005 |
| AUPRC vs TRRUST, scSAGA | 0.010 | 0.012-0.020 | prevalence 0.006 |
| early precision, MultiVI | 0.000 | 0.000-0.005 | its top-ranked edges carry no known edge; nothing here fixes that |
| hub overlap vs truth | 0.000 / 0.000 | likely still ~0 | needs a usable network first |
| out-degree Spearman | 0.041 / 0.190 | 0.05-0.10 / 0.20-0.25 | |

The projection that matters is the negative one: **adding 22,161 imputed ATAC rows may lower MultiVI's
GRN precision.** GRNBoost2 infers co-variation, and MultiVI's imputed rows are 22x denser than measured
cells with a 0.336 sd gene-mean shift — that is exactly the structure that manufactures spurious edges.
scSAGA's imputed rows preserve marginals, so its dilution is milder. If the extrapolated run shows
MultiVI's AUPRC falling while scSAGA's rises, this is the mechanism to check first, not a bug.

## 5. What breaks first

1. **Missing labels.** 12 of the 21 step-1/step-2 correctness metrics are `n/a` for lack of
   `labels.txt`, and the 5 new blocks cannot be judged by FOSCTTM at all. At 47,161 cells this is the
   binding constraint on knowing whether anything worked. Fix before scaling.
2. **MultiVI peak axis.** 89,796 -> 160,000 peaks is +46% parameters and +56% per-cell cost, and it
   buys columns that are mostly empty (union of two different peak callings).
3. **Unmodelled technology axis.** Six datasets, at least two technologies, no batch covariate.
4. **Overnight fit.** 12.9 h for two fits on CPU; a GPU turns that into ~1 h, which is the difference
   between one experiment and an ablation grid.
5. **GRN dilution.** See step 3.
6. Not a problem at this scale: scSAGA memory (largest pair unchanged), imputed export (4.9 GB),
   GRN input (8.5 GB).

## 6. What to change before running

| change | why | effect on cost |
|---|---|---|
| supply `labels.txt` (cell types for the 11k co-assayed cells) | unlocks ARI, NMI, label transfer, isolated F1, pDES, regulon specificity | none |
| **intersect** peaks across the ATAC datasets instead of unioning | a union of two peak callings is mostly empty columns; intersection is where the accessibility signal is | P 160,000 -> ~70,000, MultiVI fit 12.9 h -> 7.0 h, params 85.5M -> 50.8M |
| add a batch/technology label per dataset (`n_batch`) | without it the decoder cannot separate technology from biology | none, but needs labels |
| run it as 4 ablated configs, not 1 (§7) | one number cannot attribute cause across 6 datasets | 4x the compute; do it on a GPU |
| cap the GRN matrix to measured RNA + only the clusters with confident imputation | halves the imputed dilution and the GRN time | GRN rows 58,161 -> ~36,000 |

## 7. Ablation grid worth running instead of one big run

| config | datasets | rows | paired frac | isolates |
|---|---|---|---|---|
| A: now | rna11k, atac11k | 11,000 | 100% | the current report |
| B: + ATAC-only | A + atac11k_unp, atac3k, atac8161 | 33,161 | 33.2% | does ATAC-only data help the latent and the imputation? |
| C: + RNA-only | A + rna11k_unp, rna3k | 25,000 | 44.0% | does a bigger reference help the GRN? |
| D: all | 7 datasets | 47,161 | 23.3% | the target configuration |

B and C are the two halves of D. If D beats B and C but neither B nor C beats A on its own, the gain is
interaction, not data volume. If D is worse than A, the cause is almost certainly the unmodelled
technology axis (§4, A4), and it is diagnosable by coloring the latent by dataset.

## 8. What would falsify this projection

- MultiVI's FOSCTTM on the co-assayed block staying below 0.02 at 23.3% paired — would mean the
  alignment term is stronger per cell than the paper's unpairing sweep suggests.
- scSAGA's MSE dropping below its control — would mean the kNN average carries cell-specific signal at
  larger reference sizes, contradicting the current MSE == control result.
- GRN AUPRC rising above ~0.03 — would mean 58k cells with imputed rows is enough for network recovery,
  which the current 22k rows do not suggest.
- Step-2 densification falling to ~1-2x for MultiVI — would mean the over-densification is a
  small-sample artefact rather than the decoder's smoothness.
