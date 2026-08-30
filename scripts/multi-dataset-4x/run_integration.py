#!/usr/bin/env python
"""Run the 4-dataset scSAGA integration (3k-RNA anchor, 3k-ATAC, 6k-RNA, 6k-ATAC).

Reuses AluruLab scSAGA's modified scsaga_saveH.py so the full 30-dim joint
embedding H is saved directly to results/integration. Runs on CPU.

Set PBSC4K_ROOT to the repo root (defaults to two parents above this script).
Set SCAGA_REPO to the AluruLab scSAGA checkout that contains scmint/
(default: /Volumes/samsung_ssd/tmp/scSAGA on this machine).
"""
import os, sys, time, yaml

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
SCAGA_REPO = os.environ.get('SCAGA_REPO', '/Volumes/samsung_ssd/tmp/scSAGA')
OUT = f'{PROJ}/results/integration'

sys.path.insert(0, SCAGA_REPO)
from scmint.scsaga_saveH import main as saga_main

def build_config():
    return {
        'anchor': 'rna3k',
        'datasets': [
            dict(name='rna3k', modality='rna',
                 counts=f'{PROJ}/data/3k_rna/counts.mtx',
                 barcodes=f'{PROJ}/data/3k_rna/barcodes.txt',
                 features=f'{PROJ}/data/3k_rna/features.txt',
                 pca=f'{PROJ}/data/3k_rna/pca_50.txt'),
            dict(name='atac3k', modality='atac',
                 counts=f'{PROJ}/data/3k_atac/counts.mtx',
                 barcodes=f'{PROJ}/data/3k_atac/barcodes.txt',
                 features=f'{PROJ}/data/3k_atac/features.txt',
                 pca=f'{PROJ}/data/3k_atac/pca_50.txt'),
            dict(name='rna6k', modality='rna',
                 counts=f'{PROJ}/data/6k_rna/counts.mtx',
                 barcodes=f'{PROJ}/data/6k_rna/barcodes.txt',
                 features=f'{PROJ}/data/6k_rna/features.txt',
                 pca=f'{PROJ}/data/6k_rna/pca_50.txt'),
            dict(name='atac6k', modality='atac',
                 counts=f'{PROJ}/data/6k_atac/counts.mtx',
                 barcodes=f'{PROJ}/data/6k_atac/barcodes.txt',
                 features=f'{PROJ}/data/6k_atac/features.txt',
                 pca=f'{PROJ}/data/6k_atac/pca_50.txt'),
        ],
        'output_dir': OUT,
        's_shared_cells': 2711,
        'M_samples': 2000,
        'alpha': 0.75,
        'S_iterations': 25,
        'gw_epsilon': 1e-5,
        'gw_reg': 0.001,
    }

if __name__ == '__main__':
    import subprocess
    # Emit the resolved config into the multi-dataset scripts dir so it documents paths
    cfg = build_config()
    cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'scsaga_4dataset_config.yml')
    with open(cfg_path, 'w') as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    print('Wrote resolved config to', cfg_path)

    t0 = time.time()
    saga_main(cfg)
    print(f"\n[RUNNER] scSAGA 4-dataset integration finished in {time.time()-t0:.1f}s")
