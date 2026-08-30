# Shared scripts

Common steps used by both the single-3k and multi-dataset workflows.

- `extract_data.py` — split a 10x multiome h5 into per-modality scSAGA-format
  files (`counts.mtx`, `barcodes.txt`, `features.txt`, `pca_50.txt`).
  RNA = `Gene Expression` features, ATAC = `Peaks`.
- `run_scsaga.py` — run the scSAGA integration from a config yaml.
- `scsaga_config_single_3k.yml` — example config for the original 2-dataset 3k run.
- `run_arboreto_single.py` — GRNBoost2 on a single all-cells matrix (the 3k run).
- `evaluate_grn_single.py` — evaluate a GRNBoost2 network vs PBMC ground truth.
