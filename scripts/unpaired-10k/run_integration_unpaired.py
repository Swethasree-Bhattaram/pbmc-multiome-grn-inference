#!/usr/bin/env python
"""UNPAIRED scSAGA integration: 10k multiome RNA  x  external 10k ATAC v1.1.

Two datasets, no shared barcodes:
  - rna10k      : 10x pbmc_granulocyte_sorted_10k multiome RNA (11,898 cells, GRCh38)
  - atac10k_ext : 10x "10k Human PBMCs ATAC v1.1" cells-by-peaks (8,161 cells, hg19)

RNA is the anchor (matching the repo's convention). scSAGA aligns the ATAC PCA onto
the RNA PCA; `s_shared_cells` is the estimated size of the shared cell population,
set to min(n_rna, n_atac) = 8,161 since this is an unpaired / partial alignment.

Writes to results/integration_unpaired:
  T_atac10k_ext_to_rna10k.npy   transport plan (11898 x 8161 coupling)
  aligned_rna10k.npy            (11898 x 30)
  aligned_atac10k_ext.npy       (8161 x 30)
  joint_embedding_H.npy         (20059 x 30), anchor block first
  saga_runtimes.txt             timing + alignment scores

Run with the scSAGA venv (py3.12, needs torch/pot/geosketch/pydantic):
  /Volumes/samsung_ssd/tmp/scSAGA/.venv/bin/python run_integration_unpaired.py
"""
import os, sys, time, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
# Repo root: $UNPAIRED_ROOT, else two levels up (scripts/<workflow>/ -> repo).
# NOTE: one os.pardir would resolve to scripts/, which is a bug when the
# env var is unset (results/ and data/ then point inside scripts/).
PROJ = os.environ.get('UNPAIRED_ROOT', os.path.abspath(os.path.join(
    HERE, os.pardir, os.pardir)))
SCAGA_REPO = os.environ.get('SCAGA_REPO') or os.path.abspath(
    os.path.join(PROJ, os.pardir, 'scSAGA'))
OUT = f'{PROJ}/results/integration_unpaired'

sys.path.insert(0, SCAGA_REPO)
from scmint.scsaga_saveH import main as saga_main


def build_config():
    return {
        'anchor': 'rna10k',
        'datasets': [
            dict(name='rna10k', modality='rna',
                 counts=f'{PROJ}/data/10k_rna/counts.mtx',
                 barcodes=f'{PROJ}/data/10k_rna/barcodes.txt',
                 features=f'{PROJ}/data/10k_rna/features.txt',
                 pca=f'{PROJ}/data/10k_rna/pca_50.txt'),
            dict(name='atac10k_ext', modality='atac',
                 counts=f'{PROJ}/data/atac10k_ext/counts.mtx',
                 barcodes=f'{PROJ}/data/atac10k_ext/barcodes.txt',
                 features=f'{PROJ}/data/atac10k_ext/features.txt',
                 pca=f'{PROJ}/data/atac10k_ext/pca_50.txt'),
        ],
        'output_dir': OUT,
        's_shared_cells': 8161,   # min(n_rna, n_atac) -> full partial alignment
        'M_samples': 2000,
        'alpha': 0.75,
        'S_iterations': 25,
        'gw_epsilon': 1e-5,
        'gw_reg': 0.001,
    }


if __name__ == '__main__':
    cfg = build_config()
    cfg_path = os.path.join(HERE, 'scsaga_unpaired_config.yml')
    with open(cfg_path, 'w') as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    print('Wrote resolved config to', cfg_path, flush=True)

    t0 = time.time()
    saga_main(cfg)
    print(f"\n[RUNNER] scSAGA UNPAIRED integration finished in {time.time()-t0:.1f}s", flush=True)
