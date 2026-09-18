#!/usr/bin/env python
"""Prove that Arboreto GRNBoost2 accepts a pre-built Dask Client with >1 worker.

Builds a small but REAL expression matrix from the prepared 10k multiome RNA
counts (log1p CPM), then runs grnboost2 with a custom LocalCluster and prints the
scheduler's worker count and the resulting edge count.

This is the pattern the PACE GRN step uses; see pipeline/scripts/04_run_grn.py.
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

ROOT = os.environ.get('PACE_ROOT', '/Volumes/samsung_ssd/tmp/pp/test')
DATA = f'{ROOT}/data'
N_WORKERS = int(os.environ.get('GRN_WORKERS', '2'))
THREADS = int(os.environ.get('GRN_THREADS_PER_WORKER', '1'))


def main():
    X = scipy.io.mmread(f'{DATA}/10k_rna/counts.mtx').tocsr()      # cells x genes
    genes = np.array([l.strip() for l in open(f'{DATA}/10k_rna/features.txt')])
    # deterministic 400-cell, 300-gene slice
    rng = np.random.RandomState(0)
    cells = np.sort(rng.choice(X.shape[0], 400, replace=False))
    gcols = np.sort(rng.choice(X.shape[1], 300, replace=False))
    Xs = X[cells][:, gcols].astype(np.float64)
    tot = np.asarray(Xs.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    Xs = Xs.multiply(1e4 / tot[:, None]).tocsr()
    Xs.data = np.log1p(Xs.data)
    expr = pd.DataFrame(np.asarray(Xs.todense(), dtype=np.float32), columns=genes[gcols])
    print(f'expression for the test: {expr.shape}', flush=True)

    from dask.distributed import Client, LocalCluster
    from arboreto.algo import grnboost2

    cluster = LocalCluster(n_workers=N_WORKERS, threads_per_worker=THREADS,
                           processes=True)
    client = Client(cluster)
    n_actual = len(client.scheduler_info()['workers'])
    print(f'dask client ready: requested {N_WORKERS} workers, '
          f'scheduler reports {n_actual}', flush=True)
    print('worker ncores:', client.ncores(), flush=True)
    assert n_actual == N_WORKERS, f'expected {N_WORKERS} workers, got {n_actual}'

    t0 = time.time()
    net = grnboost2(expression_data=expr, tf_names=list(expr.columns),
                    client_or_address=client, verbose=False, seed=666)
    print(f'GRNBoost2 done in {time.time()-t0:.1f}s: {len(net)} edges', flush=True)
    print(net.head(5).to_string(index=False), flush=True)
    assert len(net) > 0, 'no edges returned'
    print('CUSTOM CLIENT MULTI-WORKER PATH OK', flush=True)
    client.close()
    cluster.close()


if __name__ == '__main__':
    main()
