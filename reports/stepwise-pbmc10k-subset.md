# Step-wise measurement on a PBMC 10k multiome subset — first real-data run

A starting point, not a result. This is the whole experiment (integration, imputation,
GRN, every metric in `workflow.pptx`) on **450 cells and 900 genes** of real paired
PBMC multiome data, run to prove the machinery produces sane numbers on real input and
to have a table to compare later runs against. Nothing here should be quoted as a
finding about either method: at this size the network step is below its own seed noise
(see the last section).

## What was run

| item | value |
|---|---|
| source | `10x-Multiome-Pbmc10k-{RNA,ATAC}.h5ad` — 9,631 cells, identical order in both files, raw integer counts, 19 annotated `cell_type`s, ATAC given as 13,205 genomic bins |
| cells | 450 of 9,631, stratified over `cell_type` (all 19 types kept, ≥20 cells where available) |
| gene axis | 900 of 29,095: top-400 HVG, then the `tf_only` TFs, then the most frequent `trrust_tf` targets |
| ATAC axis | 1,953 of 13,205 bins, detected in ≥40 of the 450 cells |
| answer key | all 450 barcodes are shared with the anchor → 90 frozen cells (20%), their RNA removed from the reference and from MultiVI's training data |
| GRN | `tf_only.txt` ∩ axis → 544 regulators; `trrust_tf.txt` ∩ axis → 613 targets, union 613 columns |
| seeds / epochs | GRNBoost2 seeds 666, 667, 668; MultiVI 250 max epochs, early stopping |
| command | `stepwise_eval.py --config config_subset.yml --methods multivi,scsaga --epochs 250 --seeds 666,667,668 --workers 6` |
| wall clock | 80 s for integrate+impute+GRN+metrics, 100 s for the two extra seeds |

FOSCTTM variant per method, as recorded in `metrics.json`: multivi = `z_expr` vs `z_acc`
of the same co-assayed cell (paired cells share one row); scsaga = RNA row vs ATAC row
of the same cell.

## Results

| step | metric | multivi | scsaga | what good looks like |
|---|---|---|---|---|
| 1 | FOSCTTM | 0.088 | 0.503 | 0.14-0.22; 0.5 = random |
| 1 | FOSCTTM, shuffled pairing | 0.508 | 0.504 | ≈0.5 |
| 1 | ARI (labels vs kmeans of H) | 0.341 | 0.061 | high, above control |
| 1 | NMI | 0.578 | 0.186 | high, above control |
| 1 | isolated-label F1 | 0.289 | 0.103 | >0.8 (capped by k/n_label) |
| 1 | graph connectivity | 0.929 | 0.809 | >0.8 |
| 1 | control ARI, labels shuffled | -0.000 | -0.000 | ≈0 |
| 1 | iLISI | n/a | 0.205 | →1 (mixed) |
| 1 | kBET-style acceptance | n/a | 0.032 | →0.95 |
| 1 | batch \|ASW\| | n/a | 0.010 | →0 |
| 1 | label transfer accuracy | 0.489 | 0.189 | >0.85 |
| 2 | MSE | 1.866 | 1.400 | low, below control |
| 2 | MSE, 1000 HVG | 2.145 | 1.688 | low |
| 2 | pFCS (top-100 HVG) | 0.735 | 0.902 | high, above null |
| 2 | pFCS, cell-permuted control | 0.930 | 0.923 | below pFCS |
| 2 | pFCS, gene-permuted null | -0.005 | -0.015 | ≈0 |
| 2 | pDES (Welch t) | 0.556 | 0.497 | high |
| 2 | gene-wise rho, mean | 0.225 | 0.066 | 0.155-0.22; >0.4 = leak |
| 2 | gene-wise rho, median | 0.216 | 0.059 | 0.155-0.22 |
| 2 | control rho, cells permuted | 0.515 | 0.057 | ≈0 |
| 2 | densification vs measured rows | 2.049 | 2.260 | ≈1 |
| 2 | gene mean shift (median, sd) | 0.154 | 0.048 | small |
| 2 | clf accuracy on imputed | 0.256 | 0.256 | near truth row |
| 2 | clustering ARI on imputed | 0.257 | 0.075 | near truth row |
| 2 | clf accuracy on truth | 0.422 | 0.422 | reference |
| 2 | clustering ARI on truth | 0.395 | 0.395 | reference |
| 3 | AUPRC vs PBMC-TRRUST | 0.005 | 0.006 | > prevalence |
| 3 | prevalence, random baseline vs PBMC-TRRUST | 0.004 | 0.005 | reference |
| 3 | early precision vs PBMC-TRRUST | 0.010 | 0.016 | high |
| 3 | recovered vs PBMC-TRRUST (run.py) | 206 | 249 | run.py number |
| 3 | AUPRC vs PBMC-Blood | 0.005 | 0.012 | > prevalence |
| 3 | prevalence, random baseline vs PBMC-Blood | 0.006 | 0.007 | reference |
| 3 | early precision vs PBMC-Blood | 0.000 | 0.016 | high |
| 3 | recovered vs PBMC-Blood (run.py) | 288 | 374 | run.py number |
| 3 | edges written | 100,148 | 50,748 | reported |
| 3 | top-10k Jaccard, mean over 3 seeds | 0.489 | 0.397 | differences must exceed it |
| 3 | edges in every seed | 5,150 | 4,389 | reported |
| 3 | hub overlap @20 vs truth | 0.050 | 0.050 | high |
| 3 | out-degree Spearman vs truth | 0.029 | 0.151 | high |
| 3 | symmetric edge fraction | 0.338 | 0.298 | low; high = smoothing |
| 3 | regulon specificity, mean | 0.014 | 0.004 | cell-type specific |

`n_hvg` = 100 for both, `n_groups` (pDES) = 2 of 19 cell types — only two types have
≥10 of the 90 held-out cells. `n_genes_scored` for rho = 771 / 769 of 900 (the rest are
constant or undetected in the held-out cells).

Gates, as `metrics.json` records them:

| gate | multivi | scsaga |
|---|---|---|
| step 1 FOSCTTM (≤0.22 and < 0.75× control) | pass | fail |
| step 1 label transfer (≥0.85) | fail | fail |
| step 2 pFCS > gene-permuted null | pass | pass |
| step 2 MSE < cell-permuted control | fail | fail |
| step 2 rho in 0.10-0.40 | pass | fail |
| step 2 no densification (\|ratio-1\| < 0.5) | fail | fail |

## What can and cannot be read off this table

Trustworthy at this size:

- The two workflows run on identical inputs and are scored by the same code, so the
  comparison is apples to apples: same gene axis, same regulator and target lists, same
  seeds, same scorer.
- **Step 1 separates the two methods by a wide margin.** MultiVI's per-cell latent pairs
  agree (FOSCTTM 0.088, control 0.508) and the embedding keeps cell types (NMI 0.578).
  scSAGA's alignment does not survive a 450-cell subset at all: FOSCTTM 0.503, iLISI 0.205,
  kBET acceptance 0.032 — the ATAC cells do not reach the RNA cells in `H`.
- **Step 2 agrees with the published ceiling.** MultiVI's gene-wise rho is 0.225, at the
  top of the Polarbear/BABEL 0.155-0.22 range for ATAC→RNA; scSAGA's is 0.066.
- **pFCS is measuring real correlation structure, not averaging.** Both methods sit far
  above the gene-permuted null (-0.005, -0.015).

Not readable, and why:

- **MSE and pFCS are both beaten by the cell-permuted control** (pFCS 0.930 / 0.923 vs
  0.735 / 0.902; MSE 0.912 / 1.399 vs 1.866 / 1.400). Permuting reference cells makes every
  guess a near-constant average, and on 90 held-out cells a near-constant block reproduces
  the gene-gene correlation matrix better than either imputation does. So at this scale
  neither metric separates a good imputation from a blur — the design doc's own warning.
  pFCS is gated on the gene-permuted null for that reason, and MSE is not a discriminator
  here at all.
- **Densification is 2.0x / 2.3x**: the imputed rows are twice as dense as measured rows.
  That is what averaging does (a 20-cell mean is rarely zero), and it is the reason step 2's
  "one scale" gate fails for both methods.
- **Step 3 is noise.** Top-10k edge Jaccard across seeds is 0.397-0.489 (pairwise 0.480,
  0.482, 0.504 for multivi and 0.403, 0.399, 0.389 for scsaga), i.e. half the evaluated
  network changes between seeds, and AUPRC sits within 0.001-0.005 of the prevalence
  baseline. No method comparison should be read from these columns. The symmetric-edge
  fraction of 0.30-0.34 is the deck's documented signature of smoothing-inflated
  correlation, and it is visible here exactly as predicted.
- **Regulon identities are a hint, not a result.** scSAGA's five top regulons are
  lineage-coherent for PBMC (PAX5, EBF1 → B; IRF8, ZEB2 → myeloid); MultiVI's are not
  (PAX3, PGR, ONECUT2, NKX3-1, HOXC13). Both have specificity ≈0.01, so neither is
  cell-type specific in the sense slide 8 asks for. With 613 target columns this is the
  network being under-determined, not a property of the imputation.

## What changes at full scale

This subset is small in the two places that matter for the network: the gene axis (900
vs the reference build's 36,601) and the cell count (450 vs ~2,700-5,400 per modality).
On PACE, with the repo's own datasets and no `--cap`:

- the gene axis is the anchor's full `features.txt`, so `trrust_tf` targets present go
  from 613 to ~2,850, and GRNBoost2 has ~4,000-5,000 columns instead of 613;
- ~2,700 cells instead of 450, which is what makes the held-out block (540 frozen cells
  instead of 90) large enough for pDES across all 19 cell types and for a 5-fold
  downstream evaluation that is not 4-5 cells per class;
- the step-3 numbers become interpretable only once the 3-seed Jaccard is well below the
  AUPRC differences you care about — check that first, before reading any method
  difference.

## Known limitations of the numbers in this file

- pDES is a Welch t statistic per gene, not limma (R-only); `kbet_acceptance` is the
  chi-square acceptance on the kNN composition, not kBET's own estimator. Both are
  labelled as such in `metrics.json`.
- Step 1 batch metrics are `n/a` for MultiVI: its row space is collapsed (paired cells share
  one row), so there is one batch per dataset and nothing to mix. They return as soon as the
  experiment mixes a multiome dataset with a single-modality one.
- The `truth` rows are the hidden measured RNA, which carries its own dropout; a denoising
  imputer can beat them on a downstream task, so 0.422 / 0.395 are reference points, not
  ceilings.
- 450 cells, 1,953 bins and 900 genes are far below the regime either method was designed
  for. Treat every number above as "the pipeline runs and the arithmetic is right", not as
  an estimate of what the full run will give.
