#!/usr/bin/env python
"""RNA-only baseline with the *purified* TF list as regulators (tf_only.txt).

Same pipeline and same input matrix as scripts/rna-baseline/run_arboreto_trrust.py
(stacked PBMC 3k multiome RNA + 10x PBMC 4k RNA, 7051 cells x 21932 genes, no
integration, no imputation). The ONE change is the transcription-factor list fed
to GRNBoost2:

  * run_arboreto_trrust.py (original) : REGULATORS = every gene of
        data/trrust_tf.txt present in the matrix -> 2,808 regulators, and the
        expression subset handed to GRNBoost2 is those 2,808 columns.
  * this script (tf_only)            : REGULATORS = only the transcription
        factors of data/tf_only.txt present in the matrix -> 816 regulators.

TARGET genes are the genes of data/trrust_tf.txt present in the matrix (2,808),
exactly as in the original run, so the two networks are directly comparable:
the only variable is the regulator list.

Every regulator must be a matrix column, so the expression subset is the UNION of
the TRRUST target columns and the present tf_only regulator columns (2,833),
whereas the original run handed GRNBoost2 only the TRRUST columns (2,808).

Inputs:
  results/rna_baseline/all_cells_gene_expression.npy   (7051 x 21932 float32)
  results/rna_baseline/genes.npy
Output:
  results/rna_baseline/grn_tfonly_reg/

Run in the .venv39 (python 3.9, numpy 1.21, dask 2021.10, arboreto 0.1.6); must
be wrapped in `if __name__ == '__main__'` for macOS multiprocessing (spawn).
Worker count via N_WORKERS (default 5): the tf_only input is small
(7051 x 2833 float32 = 80 MB), so the RAM-thrash constraint that forces 2 workers
on the full-10k runs does not apply here.
"""
import os, numpy as np, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
TR_F = f'{PROJ}/data/trrust_tf.txt'   # mixed TRRUST list  -> TARGET genes
TF_F = f'{PROJ}/data/tf_only.txt'     # purified TF list   -> REGULATORS


def main():
    DOWN = f'{PROJ}/results/rna_baseline'
    OUT = f'{DOWN}/grn_tfonly_reg'
    os.makedirs(OUT, exist_ok=True)

    # memory-safe: map the matrix read-only, subset columns, then materialise
    Xmm = np.load(f'{DOWN}/all_cells_gene_expression.npy', mmap_mode='r')
    genes = np.load(f'{DOWN}/genes.npy', allow_pickle=True).astype(str)
    print('Expression matrix:', Xmm.shape, Xmm.dtype, 'genes:', len(genes), flush=True)

    trrust = [l.strip() for l in open(TR_F) if l.strip()]
    tfonly = [l.strip() for l in open(TF_F) if l.strip()]

    target_genes = [g for g in genes if g in set(trrust)]
    present_tfs = [g for g in genes if g in set(tfonly)]
    cols = list(dict.fromkeys(target_genes + present_tfs))   # union, order-preserving
    print(f'TRRUST list: {len(trrust)} -> targets present: {len(target_genes)}', flush=True)
    print(f'tf_only list: {len(tfonly)} -> regulators present: {len(present_tfs)}', flush=True)
    print('Union columns:', len(cols), flush=True)

    gidx = {g: i for i, g in enumerate(genes)}
    col_idx = [gidx[g] for g in cols]
    expr_sub = pd.DataFrame(np.asarray(Xmm[:, col_idx], dtype=np.float32), columns=cols)
    del Xmm
    n_cells = expr_sub.shape[0]
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
    with open(f'{OUT}/all_columns.txt', 'w') as f:
        f.write('\n'.join(cols))
    with open(f'{OUT}/n_cells.txt', 'w') as f:
        f.write(str(n_cells))
    print('Saved grnboost2_network to', OUT, flush=True)
    client.close(); cluster.close()


if __name__ == '__main__':
    main()
