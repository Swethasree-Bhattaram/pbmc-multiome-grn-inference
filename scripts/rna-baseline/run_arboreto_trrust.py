#!/usr/bin/env python
"""Run Arboreto GRNBoost2 on the stacked two-RNA baseline with TRRUST-only targets.

RNA-only baseline (no integration, no imputation): the stacked PBMC 3k + 4k RNA
matrix is handed directly to GRNBoost2. TARGET genes = ONLY the TRRUST TFs present
in the expression matrix (no top-2000 HVG).

Reads  results/rna_baseline/all_cells_gene_expression.npy
Writes results/rna_baseline/grn_trrust/

Run in the .venv39 (python 3.9, dask 2021.10, arboreto 0.1.6).
Must be wrapped in `if __name__ == '__main__'` for macOS multiprocessing.

Worker count via N_WORKERS (default 5). The TRRUST-only input is small
(~7k cells x ~2.8k genes), so the RAM-thrash constraint that forces 2 workers on
the full-10k runs does not apply here.
"""
import os, sys, numpy as np, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
TF_F = f'{PROJ}/data/trrust_tf.txt'


def main():
    DOWN = f'{PROJ}/results/rna_baseline'
    OUT = f'{DOWN}/grn_trrust'
    os.makedirs(OUT, exist_ok=True)

    X = np.load(f'{DOWN}/all_cells_gene_expression.npy').astype(np.float32)
    genes = np.load(f'{DOWN}/genes.npy', allow_pickle=True).astype(str)
    print('Expression matrix:', X.shape, 'dtype:', X.dtype, 'genes:', len(genes), flush=True)

    tfs = [l.strip() for l in open(TF_F)]
    tf_set = set(tfs)
    present_tfs = [g for g in genes if g in tf_set]
    print('TRRUST TFs present in expression:', len(present_tfs), flush=True)

    target_genes = present_tfs
    print('Target genes (TRRUST-only):', len(target_genes), flush=True)

    n_cells = X.shape[0]
    gene_pos = {g: i for i, g in enumerate(genes)}
    col_idx = [gene_pos[g] for g in target_genes]
    expr_sub = pd.DataFrame(X[:, col_idx], columns=target_genes)
    del X
    print('Subset expression shape:', expr_sub.shape, flush=True)

    from dask.distributed import Client, LocalCluster
    from arboreto.algo import grnboost2

    n_workers = int(os.environ.get('N_WORKERS', '5'))
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=2, processes=True)
    client = Client(cluster)
    print(f'Dask client ready (n_workers={n_workers}). Workers:',
          len(client.scheduler_info()['workers']), flush=True)

    network = grnboost2(expression_data=expr_sub, tf_names=present_tfs,
                        client_or_address=client, verbose=False, seed=666)
    print('GRNBoost2 done. Edges:', len(network), flush=True)

    network.to_csv(f'{OUT}/grnboost2_network.tsv', sep='\t', index=False)
    network.to_csv(f'{OUT}/grnboost2_network.csv', index=False)
    with open(f'{OUT}/tf_regulators.txt', 'w') as f:
        f.write('\n'.join(present_tfs))
    with open(f'{OUT}/target_genes.txt', 'w') as f:
        f.write('\n'.join(target_genes))
    with open(f'{OUT}/n_cells.txt', 'w') as f:
        f.write(str(n_cells))
    print('Saved grnboost2_network to', OUT, flush=True)
    client.close(); cluster.close()


if __name__ == '__main__':
    main()
