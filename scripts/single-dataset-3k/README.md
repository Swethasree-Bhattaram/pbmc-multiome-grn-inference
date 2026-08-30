# Original single-dataset (3k) workflow

This preserves the original pbmc_3k multiome GRN pipeline (2,711 shared cells):

```
extract_data.py          # split 3k h5 -> scsaga_input (RNA + ATAC)
run_scsaga.py            # scSAGA integrate RNA + ATAC -> joint embedding H (5,422 x 30)
reverse_impute.py        # reverse-imputeKNN: impute RNA expression onto ATAC cells
run_arboreto.py          # GRNBoost2 on all-cells matrix (5,422 x 36,601)
evaluate_grn.py          # eval vs PBMC-TRRUST / PBMC-Blood
```

See the full workflow and results in the repo's `reports/` (from the original
`pbmc3k_analysis`).
