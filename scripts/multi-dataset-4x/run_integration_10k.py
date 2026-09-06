#!/usr/bin/env python
"""Run the 4-dataset scSAGA integration (3k-RNA anchor, 3k-ATAC, 10k-RNA, 10k-ATAC).

The 10k datasets use ALL 11,898 cells (no subsampling). Anchor = 3k RNA.
Writes to results/integration_10k.
"""
import os, sys, time, yaml

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir)))
SCAGA_REPO = os.environ.get('SCAGA_REPO', '/Volumes/samsung_ssd/tmp/scSAGA')
OUT = f'{PROJ}/results/integration_10k'

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
            dict(name='rna10k', modality='rna',
                 counts=f'{PROJ}/data/10k_rna/counts.mtx',
                 barcodes=f'{PROJ}/data/10k_rna/barcodes.txt',
                 features=f'{PROJ}/data/10k_rna/features.txt',
                 pca=f'{PROJ}/data/10k_rna/pca_50.txt'),
            dict(name='atac10k', modality='atac',
                 counts=f'{PROJ}/data/10k_atac/counts.mtx',
                 barcodes=f'{PROJ}/data/10k_atac/barcodes.txt',
                 features=f'{PROJ}/data/10k_atac/features.txt',
                 pca=f'{PROJ}/data/10k_atac/pca_50.txt'),
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
    cfg = build_config()
    cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'scsaga_4dataset_10k_config.yml')
    with open(cfg_path, 'w') as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    print('Wrote resolved config to', cfg_path)

    t0 = time.time()
    saga_main(cfg)
    print(f"\n[RUNNER] scSAGA 4-dataset (full-10k) integration finished in {time.time()-t0:.1f}s")
