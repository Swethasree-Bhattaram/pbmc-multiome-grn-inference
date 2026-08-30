#!/usr/bin/env python
"""Run Arboreto GRNBoost2 on an experiment's all-cells expression matrix.

Consumes results/<exp>/all_cells_gene_expression.npy -> results/<exp>/grn.
Run in the .venv39 (python 3.9, dask 2021.10, arboreto 0.1.6).
Must be wrapped in `if __name__ == '__main__'` for macOS multiprocessing.
"""
import os, sys, numpy as np, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
TF_F = f'{PROJ}/trrust_tf.txt'

def main():
    exp = sys.argv[1]
    DOWN = f'{PROJ}/results/{exp}'
    OUT = f'{PROJ}/results/{exp}/grn'
    os.makedirs(OUT, exist_ok=True)

    X = np.load(f'{DOWN}/all_cells_gene_expression.npy')
    genes = np.load(f'{DOWN}/genes.npy', allow_pickle=True).astype(str)
    print('Expression matrix:', X.shape, 'genes:', len(genes), flush=True)

    expr = pd.DataFrame(X, columns=genes)

    tfs = [l.strip() for l in open(TF_F)]
    tf_set = set(tfs)
    present_tfs = [g for g in genes if g in tf_set]
    print('TRRUST TFs present in expression:', len(present_tfs), flush=True)

    var = expr.var(axis=0).sort_values(ascending=False)
    hvg = list(var.index[:2000])
    target_genes = list(dict.fromkeys(hvg + present_tfs))
    print('Target genes:', len(target_genes), flush=True)
    expr_sub = expr[target_genes]
    print('Subset expression shape:', expr_sub.shape, flush=True)

    from dask.distributed import Client, LocalCluster
    from arboreto.algo import grnboost2

    cluster = LocalCluster(n_workers=2, threads_per_worker=2, processes=True)
    client = Client(cluster)
    print('Dask client ready. Workers:', len(client.scheduler_info()['workers']), flush=True)

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
        f.write(str(X.shape[0]))
    print('Saved grnboost2_network to', OUT, flush=True)
    client.close(); cluster.close()

if __name__ == '__main__':
    main()
