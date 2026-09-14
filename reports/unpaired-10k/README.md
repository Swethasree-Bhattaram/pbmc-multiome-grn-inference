# Unpaired 10k PBMC GRN (multiome RNA × external ATAC v1.1)

Integration + reverse-imputeKNN + Arboreto GRNBoost2 on an **unpaired** pairing:
the 10k multiome RNA (11,898 cells) and the external 10x 10k PBMC ATAC v1.1
"cells by peaks" dataset (8,161 cells). All-cells matrix 20,059 × 36,601.

- Targets: `trrust_tf.txt` genes present = 2,624
- Regulators: `tf_only.txt` TFs present = 749
- Edges: 659,314
- Recovered: PBMC-TRRUST 3,110/8,751 (35.5%), PBMC-Blood 4,486/96,846 (4.6%)

See **REPORT.md** for the full pipeline, evaluation tables and QC (including the
unpaired caveats and the block-segregation comparison against the repo's paired run).
Raw QC logs: `qc_imputation.txt`, `qc_mixing.txt`, `qc_biology.txt`.
