#!/usr/bin/env python
"""Benchmark GRNBoost2 wall time vs number of dask workers on one node.

Runs the same real expression slice with n_workers = 1, 2, 4, 8 and reports wall
time.  Use the numbers to size --cpus-per-task in the GRN .slurm script.

Usage:
  GRN_THREADS_PER_WORKER=1 python _bench_arboreto_workers.py
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.io

ROOT = os.environ.get('PACE_ROOT', '/Volumes/samsung_ssd/tmp/pp/test')
DATA = f'{ROOT}/data'
N_CELLS = int(os.environ.get('BENCH_CELLS', '2000'))
N_GENES = int(os.environ.get('BENCH_GENES', '400'))
THREADS = int(os.environ.get('GRN_THREADS_PER_WORKER', '1'))
WORKER_LIST = [int(x) for x in os.environ.get('BENCH_WORKERS', '1,2,4,8').split(',')]


def build():
    X = scipy.io.mmread(f'{DATA}/10k_rna/counts.mtx').tocsr()
    genes = np.array([l.strip() for l in open(f'{DATA}/10k_rna/features.txt')])
    rng = np.random.RandomState(0)
    cells = np.sort(rng.choice(X.shape[0], N_CELLS, replace=False))
    gcols = np.sort(rng.choice(X.shape[1], N_GENES, replace=False))
    Xs = X[cells][:, gcols].astype(np.float64)
    tot = np.asarray(Xs.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    Xs = Xs.multiply(1e4 / tot[:, None]).tocsr()
    Xs.data = np.log1p(Xs.data)
    return pd.DataFrame(np.asarray(Xs.todense(), dtype=np.float32), columns=genes[gcols])


def main():
    expr = build()
    print(f'expression: {expr.shape}, threads_per_worker={THREADS}', flush=True)
    from dask.distributed import Client, LocalCluster
    from arboreto.algo import grnboost2

    print(f'{"workers":>8} {"ncores":>7} {"wall_s":>9} {"edges":>8}', flush=True)
    for nw in WORKER_LIST:
        cluster = LocalCluster(n_workers=nw, threads_per_worker=THREADS, processes=True)
        client = Client(cluster)
        ncores = sum(client.ncores().values())
        t0 = time.time()
        net = grnboost2(expression_data=expr, tf_names=list(expr.columns),
                        client_or_address=client, verbose=False, seed=666)
        dt = time.time() - t0
        print(f'{nw:>8} {ncores:>7} {dt:>9.1f} {len(net):>8}', flush=True)
        client.close()
        cluster.close()


if __name__ == '__main__':
    main()
