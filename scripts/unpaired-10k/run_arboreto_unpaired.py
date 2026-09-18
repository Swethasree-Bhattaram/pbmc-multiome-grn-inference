#!/usr/bin/env python
"""Arboreto GRNBoost2 for the UNPAIRED 10k run (RNA10k multiome + external 10k ATAC v1.1).

Regulator / target specification (as requested, same shape as the repo's tf_only
regulator variant `run_arboreto_tfonly_reg_10k.py`):

  * TARGET genes  = every gene in the all-cells matrix that appears in
                    data/trrust_tf.txt (the mixed TRRUST list).  ~2,827 genes.
  * REGULATORS    = only the transcription factors in data/tf_only.txt that are
                    present as columns (tf_names fed to GRNBoost2).      ~816 TFs.

Every regulator must be a matrix column, so the expression subset is the UNION of
the TRRUST target columns and the present tf_only regulator columns (~2,852).

Expression source: results/unpaired_grn/all_cells_gene_expression.npy
  rows  = [ rna10k real (11,898) , atac10k_ext imputed (8,161) ]
  cols  = 36,601 genes

Writes to results/unpaired_grn/grn_tfonly_reg/ .

Run in the old-stack Arboreto env (python 3.9, numpy 1.21, dask 2021.10,
arboreto 0.1.6).  Must be wrapped in `if __name__ == '__main__'` on macOS
(multiprocessing spawn).
"""
import os, sys, numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.environ.get('UNPAIRED_ROOT', os.path.abspath(os.path.join(HERE, os.pardir)))
TR_F = os.environ.get('TRRUST_FILE', f'{PROJ}/data/trrust_tf.txt')  # mixed TRRUST list -> TARGETS
TF_F = os.environ.get('TF_ONLY_FILE', f'{PROJ}/data/tf_only.txt')   # purified TF list -> REGULATORS
DOWN = f'{PROJ}/results/unpaired_grn'
OUT = f'{DOWN}/grn_tfonly_reg'
N_WORKERS = int(os.environ.get('GRN_WORKERS', '4'))


def main():
    os.makedirs(OUT, exist_ok=True)

    genes = np.load(f'{DOWN}/genes.npy', allow_pickle=True).astype(str)
    # memory-safe: map the 2.9 GB matrix read-only, subset columns, then materialise
    Xmm = np.load(f'{DOWN}/all_cells_gene_expression.npy', mmap_mode='r')
    print('Expression matrix:', Xmm.shape, Xmm.dtype, 'genes:', len(genes), flush=True)

    trrust = [l.strip() for l in open(TR_F) if l.strip()]
    tfonly = [l.strip() for l in open(TF_F) if l.strip()]

    gidx = {g: i for i, g in enumerate(genes)}
    target_genes = [g for g in genes if g in set(trrust)]
    present_tfs = [g for g in genes if g in set(tfonly)]
    cols = list(dict.fromkeys(target_genes + present_tfs))   # union, order-preserving
    print('TRRUST list:', len(trrust), '-> targets present:', len(target_genes), flush=True)
    print('tf_only list:', len(tfonly), '-> regulators present:', len(present_tfs), flush=True)
    print('Union columns:', len(cols), flush=True)

    col_idx = [gidx[g] for g in cols]
    expr_sub = pd.DataFrame(np.asarray(Xmm[:, col_idx], dtype=np.float32), columns=cols)
    del Xmm
    n_cells = expr_sub.shape[0]
    print('Subset expression shape:', expr_sub.shape, flush=True)

    from dask.distributed import Client, LocalCluster
    from arboreto.algo import grnboost2

    cluster = LocalCluster(n_workers=N_WORKERS, threads_per_worker=2, processes=True)
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
    with open(f'{OUT}/all_columns.txt', 'w') as f:
        f.write('\n'.join(cols))
    with open(f'{OUT}/n_cells.txt', 'w') as f:
        f.write(str(n_cells))
    print('Saved grnboost2_network to', OUT, flush=True)
    client.close(); cluster.close()


if __name__ == '__main__':
    main()
