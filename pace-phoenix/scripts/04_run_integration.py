#!/usr/bin/env python
"""Run the scSAGA integration for either the PAIRED or the UNPAIRED experiment.

Both experiments are config-driven; the only differences are which datasets are
loaded, the anchor, and s_shared_cells.  Every path comes from $PACE_ROOT +
$PACE_DATA so nothing is hardcoded.

  paired    rna10k, atac10k                       (10k multiome RNA + its own ATAC)
  unpaired  rna10k, atac10k_ext                   (10k multiome RNA + external ATAC v1.1)

scSAGA integration parameters are the ones used for the published runs:
  M_samples 2000, alpha 0.75, S_iterations 25, gw_epsilon 1e-5, gw_reg 0.001.

s_shared_cells defaults match the published runs:
  paired   -> min(n_rna, n_atac); full 10k multiome has 11,898 in both modalities,
              so this is 11,898 (see results/integration_10k for the 4-dataset pair)
  unpaired -> min(11,898, 8,161) = 8,161 (a partial alignment by construction)

Requires the scSAGA checkout with the SAVE-H patch applied
(see envs/03_patch_scsaga.py); $SCAGA_REPO points at it.

Writes into $PACE_ROOT/results/integration_<experiment>/:
  joint_embedding_H.npy            (n_cells x 30), anchor block first
  joint_embedding_barcodes.txt
  aligned_<name>.npy / aligned_<name>_barcodes.txt
  T_<name>_to_<anchor>.npy         transport plan
  saga_runtimes.txt, run.log

Usage:
  python 04_run_integration.py --experiment unpaired
  python 04_run_integration.py --experiment paired
"""
import argparse
import os
import subprocess
import sys
import time

import yaml


def build_config(exp, proj, data):
    if exp == 'unpaired':
        datasets = [
            dict(name='rna10k', modality='rna', dir=f'{data}/10k_rna'),
            dict(name='atac10k_ext', modality='atac', dir=f'{data}/atac10k_ext'),
        ]
        anchor, shared = 'rna10k', 8161
    elif exp == 'paired':
        datasets = [
            dict(name='rna10k', modality='rna', dir=f'{data}/10k_rna'),
            dict(name='atac10k', modality='atac', dir=f'{data}/10k_atac'),
        ]
        anchor, shared = 'rna10k', 11898
    else:
        raise SystemExit(f'unknown experiment {exp!r} (use paired|unpaired)')

    specs = []
    for d in datasets:
        counts = f"{d['dir']}/counts.mtx"
        if not os.path.exists(counts):
            raise SystemExit(f'{d["name"]}: missing {counts} -- run 00_split_10x_h5.py first')
        specs.append(dict(name=d['name'], modality=d['modality'], dir=d['dir'],
                          counts=counts,
                          barcodes=f"{d['dir']}/barcodes.txt",
                          features=f"{d['dir']}/features.txt",
                          pca=f"{d['dir']}/pca_50.txt"))

    return {
        'anchor': anchor,
        'datasets': specs,
        'output_dir': f'{proj}/results/integration_{exp}',
        's_shared_cells': shared,
        'M_samples': 2000,
        'alpha': 0.75,
        'S_iterations': 25,
        'gw_epsilon': 1e-5,
        'gw_reg': 0.001,
    }, anchor, shared


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--experiment', choices=['paired', 'unpaired'], required=True)
    ap.add_argument('--scsaga-repo', default=os.environ.get('SCAGA_REPO'))
    args = ap.parse_args()

    proj = os.environ.get('PACE_ROOT', os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    data = os.environ.get('PACE_DATA', f'{proj}/data')
    repo = args.scsaga_repo
    if not repo or not os.path.isdir(os.path.join(repo, 'scmint')):
        raise SystemExit(f'scSAGA checkout not found (SCAGA_REPO={repo!r}); '
                         f'run envs/02_setup_envs.sh first')

    cfg, anchor, shared = build_config(args.experiment, proj, data)
    cfg_path = f"{proj}/results/integration_{args.experiment}/scsaga_config.yml"
    os.makedirs(cfg['output_dir'], exist_ok=True)
    with open(cfg_path, 'w') as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    print(f'wrote resolved config -> {cfg_path}', flush=True)
    print(f'experiment={args.experiment} anchor={anchor} s_shared_cells={shared}',
          flush=True)

    sys.path.insert(0, repo)
    try:
        from scmint.scsaga_saveH import main as saga_main   # patched copy
        patched = True
    except ImportError:
        from scmint.scsaga import main as saga_main
        patched = False
    if not patched:
        print('WARNING: using unpatched scmint.scsaga -- it will NOT write '
              'joint_embedding_H.npy. Apply envs/03_patch_scsaga.py first.',
              flush=True)

    t0 = time.time()
    saga_main(cfg)
    dt = time.time() - t0

    H = f"{cfg['output_dir']}/joint_embedding_H.npy"
    if not os.path.exists(H):
        raise SystemExit(f'scSAGA finished but {H} is missing; '
                         f'the SAVE-H patch is not applied or failed')
    import numpy as np
    h = np.load(H, mmap_mode='r')
    print(f'[OK] joint_embedding_H.npy {h.shape} written in {dt:.1f}s', flush=True)


if __name__ == '__main__':
    sys.exit(main())
