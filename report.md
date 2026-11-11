# multi11k: step-wise measurement, multivi vs scsaga

datasets: rna11k, atac11k  |  anchor: rna11k  |  ATAC queries: atac11k

FOSCTTM variant: multivi = z_expr vs z_acc of the same co-assayed cell; scsaga = RNA row vs ATAC row of the same cell

| step | metric | multivi | scsaga | what good looks like |
|---|---|---|---|---|
| 1 | FOSCTTM | 0.008 | 0.241 | 0.14-0.22; 0.5 = random |
| 1 | FOSCTTM, shuffled pairing | 0.500 | 0.502 | ≈0.5 |
| 1 | ARI (labels vs kmeans of H) | n/a | n/a | high, above control |
| 1 | NMI | n/a | n/a | high, above control |
| 1 | isolated-label F1 | n/a | n/a | > 0.8 |
| 1 | graph connectivity | n/a | n/a | > 0.8 |
| 1 | control ARI (labels shuffled) | n/a | n/a | ≈0 |
| 1 | iLISI | n/a | 0.269 | → 1 (mixed) |
| 1 | kBET-style acceptance | n/a | 0.115 | → 0.95 |
| 1 | batch |ASW| | n/a | 0.015 | → 0 |
| 1 | label transfer accuracy | n/a | n/a | > 0.85 |
| 2 | MSE | 0.172 | 0.124 | low, below control |
| 2 | MSE (1000 HVG) | 1.625 | 1.283 | low |
| 2 | pFCS (top-100 HVG correlation) | 0.928 | 0.950 | high, above control |
| 2 | control pFCS (shuffled reference) | 0.957 | 0.949 | clearly below pFCS |
| 2 | control MSE | 0.097 | 0.124 | clearly above MSE |
| 2 | control pFCS (gene labels shuffled) | -0.008 | 0.007 | the null pFCS: must be below pFCS |
| 2 | pDES (Welch t, not limma) | status:unavailable, reason:needs >=2 cell types on the held-out cells (labels.txt) | status:unavailable, reason:needs >=2 cell types on the held-out cells (labels.txt) | high |
| 2 | gene-wise rho, mean | 0.333 | 0.292 | 0.155-0.22; >0.4 = leak |
| 2 | gene-wise rho, median | 0.300 | 0.245 | 0.155-0.22 |
| 2 | control gene-wise rho | 0.573 | 0.299 | ≈0 |
| 2 | densification vs measured rows | 22.506 | 6.365 | ≈1, not inflated |
| 2 | gene mean shift (median, sd units) | 0.336 | 0.013 | small |
| 2 | clf accuracy on imputed | n/a | n/a | close to truth row |
| 2 | clf macro-F1 on imputed | n/a | n/a | close to truth row |
| 2 | clustering ARI on imputed | n/a | n/a | close to truth row |
| 2 | clf accuracy on truth | n/a | n/a | ceiling |
| 2 | clustering ARI on truth | n/a | n/a | ceiling |
| 3 | AUPRC vs PBMC-TRRUST | 0.006 | 0.010 | > prevalence |
| 3 | random baseline (prevalence) vs PBMC-TRRUST | 0.005 | 0.006 | reference |
| 3 | early precision vs PBMC-TRRUST | 0.000 | 0.031 | high |
| 3 | precision@10000 vs PBMC-TRRUST | 0.007 | 0.009 | high |
| 3 | recall@10000 vs PBMC-TRRUST | 0.008 | 0.011 | high |
| 3 | recovered vs PBMC-TRRUST (run.py) | 274 | 320 | run.py number |
| 3 | AUPRC vs PBMC-Blood | 0.006 | 0.014 | > prevalence |
| 3 | random baseline (prevalence) vs PBMC-Blood | 0.007 | 0.010 | reference |
| 3 | early precision vs PBMC-Blood | 0.000 | 0.023 | high |
| 3 | precision@10000 vs PBMC-Blood | 0.005 | 0.015 | high |
| 3 | recall@10000 vs PBMC-Blood | 0.001 | 0.002 | high |
| 3 | recovered vs PBMC-Blood (run.py) | 353 | 484 | run.py number |
| 3 | edge stability, mean Jaccard (top-10k, seeds) | 0.546 | 0.405 | differences must exceed this |
| 3 | edges in every seed | 6036 | 4309 | reported |
| 3 | hub overlap vs truth | 0.000 | 0.000 | high |
| 3 | out-degree Spearman vs truth | 0.041 | 0.190 | high |
| 3 | symmetric edge fraction | 0.114 | 0.069 | low; high = smoothing |
| 3 | regulon specificity, mean | n/a | n/a | cell-type specific |
| 3 | regulons peaking per cell type | n/a | n/a | matches lineages |

The `truth` rows are the hidden measured RNA, which carries its own dropout: a denoising imputer can beat them on a downstream task, so read them as a reference point rather than a hard ceiling.

Two negative controls are reported for step 2.  `shuffled reference` permutes the reference cells (the design doc's control, and the one MSE is gated on); it makes every guess a near-constant average. `gene labels shuffled` permutes the reference's gene rows, which leaves every marginal intact and destroys only the gene-gene structure -- that is the null pFCS is gated on, because a near-constant block can otherwise score well on pFCS.

Gates (design doc sections 3-5):

| gate | multivi | scsaga |
|---|---|---||
| step 1 FOSCTTM | pass | fail |
| step 1 label transfer | n/a | n/a |
| step 2 pFCS > gene-permuted null | pass | pass |
| step 2 MSE < control | fail | fail |
| step 2 rho in range | pass | pass |
| step 2 no densification | no | no |

`n/a` means the metric is unavailable for this method, or the input it needs is absent.  metrics.json records the reason for every one of them.  Nothing is approximated silently.

## multivi

- imputed row density vs hidden answer: 22.676x
- genes whose imputed mean is off by >0.5 sd: 23.4%
- seed-to-seed top-10000 Jaccard: 0.552, 0.537, 0.550

## scsaga

- imputed row density vs hidden answer: 6.413x
- genes whose imputed mean is off by >0.5 sd: 0.0%
- seed-to-seed top-10000 Jaccard: 0.404, 0.391, 0.418

