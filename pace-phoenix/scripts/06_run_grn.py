#!/usr/bin/env python
"""Arboreto GRNBoost2 on the all-cells matrix, with a configurable multi-worker
Dask client.

Regulator / target definition (matches the repo's tf_only-regulator runs):
  * REGULATORS (tf_names) = TFs of data/tf_only.txt present as columns   (~816)
  * TARGETS              = genes of data/trrust_tf.txt present as columns (~2,827)
  * Expression subset    = UNION of target columns and present regulator columns
                           (~2,852), because every regulator must also be a column.

Dask parallelisation  -- the answer to "can the Arboreto step run in parallel with
more workers?"  Yes.  Arboreto accepts a pre-built `distributed.Client` via
`client_or_address`, so the worker count is a runtime choice rather than a code
change (https://arboreto.readthedocs.io/en/latest/userguide.html#running-with-a-custom-dask-client).
Three modes are supported here:

  local (default)   LocalCluster(n_workers=$GRN_WORKERS, threads_per_worker=$GRN_THREADS_PER_WORKER)
                    -> single node, one process per worker. On a SLURM node set
                    GRN_WORKERS so that GRN_WORKERS * GRN_THREADS_PER_WORKER
                    equals --cpus-per-task.
  external          connect to an already-running scheduler: GRN_SCHEDULER=tcp://host:8786
                    -> lets you start `dask-scheduler` + N x `dask-worker` yourself
                    (or use dask-jobqueue, see --jobqueue below).
  jobqueue          SLURMCluster via dask_jobqueue: the client submits its own
                    worker jobs to SLURM, so GRN_NODE_WORKERS workers each get a
                    node (multi-node GRN).

Requires the py3.9 arboreto env (numpy 1.21, pandas 1.4, dask/distributed 2021.10,
arboreto 0.1.6) -- see envs/02_setup_envs.sh.  Must run under
`if __name__ == '__main__'` for spawn-based multiprocessing.

Usage:
  GRN_WORKERS=8 python 06_run_grn.py --experiment unpaired
  GRN_SCHEDULER=tcp://127.0.0.1:8786 python 06_run_grn.py --experiment unpaired
  GRN_NODE_WORKERS=2 python 06_run_grn.py --experiment unpaired --jobqueue \
       --slurm-account GT-xxxx --slurm-qos embers --slurm-partition cpu-small

Writes $PACE_ROOT/results/<exp>_grn/grn_tfonly_reg/.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

DEMON_SEED = 666


def resolve_root():
    return os.environ.get('PACE_ROOT', os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))


def make_client(args, n_workers, threads_per_worker):
    """Return (client, cluster_or_None, mode_label)."""
    from dask.distributed import Client, LocalCluster

    if args.scheduler:
        client = Client(args.scheduler)
        n = len(client.scheduler_info()['workers'])
        print(f'[dask] connected to external scheduler {args.scheduler} '
              f'({n} workers)', flush=True)
        return client, None, 'external'

    if args.jobqueue:
        try:
            from dask_jobqueue import SLURMCluster
        except ImportError:
            raise SystemExit(
                'dask-jobqueue is not installed in this env. Install it with\n'
                '  pip install "dask-jobqueue==0.7.2"\n'
                '(compatible with dask/distributed 2021.10) or run the pipeline\n'
                'in envs/02_setup_envs.sh, which installs it.'
            )
        job_extra = []
        for flag, val in [('--qos', args.slurm_qos),
                          ('--account', args.slurm_account),
                          ('--partition', args.slurm_partition)]:
            if val:
                job_extra.append(f'{flag}={val}')
        cluster = SLURMCluster(
            cores=args.slurm_cores,
            memory=args.slurm_memory,
            processes=1,                       # 1 dask worker process per job
            walltime=args.slurm_time,
            job_name='grn-worker',
            job_extra=job_extra or None,
            local_directory=args.local_directory,
        )
        cluster.scale(jobs=args.jobqueue_workers)
        client = Client(cluster)
        print(f'[dask] SLURMCluster scaled to {args.jobqueue_workers} worker jobs x '
              f'{args.slurm_cores} cores ({args.slurm_memory}); waiting for workers...',
              flush=True)
        client.wait_for_workers(args.jobqueue_workers, timeout=args.worker_timeout)
        print(f'[dask] {len(client.scheduler_info()["workers"])} workers online',
              flush=True)
        return client, cluster, 'jobqueue'

    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker,
                           processes=True, local_directory=args.local_directory)
    client = Client(cluster)
    n = len(client.scheduler_info()['workers'])
    ncores = sum(client.ncores().values())
    print(f'[dask] LocalCluster: {n} workers x {threads_per_worker} threads '
          f'= {ncores} cores', flush=True)
    if ncores != args.cpus:
        print(f'[dask] WARNING: {ncores} dask cores but --cpus {args.cpus}; '
              f'set GRN_WORKERS*GRN_THREADS_PER_WORKER == --cpus-per-task',
              flush=True)
    return client, cluster, 'local'


def main():
    proj = resolve_root()
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--experiment', choices=['paired', 'unpaired'], required=True)
    ap.add_argument('--cpus', type=int, default=int(os.environ.get('GRN_CPUS',
                                                                   os.environ.get('SLURM_CPUS_PER_TASK', '8'))))
    ap.add_argument('--workers', type=int, default=int(os.environ.get('GRN_WORKERS', '0')),
                    help='local mode worker processes (0 = derive from --cpus)')
    ap.add_argument('--threads-per-worker', type=int,
                    default=int(os.environ.get('GRN_THREADS_PER_WORKER', '1')))
    ap.add_argument('--scheduler', default=os.environ.get('GRN_SCHEDULER'))
    ap.add_argument('--jobqueue', action='store_true',
                    default=bool(int(os.environ.get('GRN_JOBQUEUE', '0'))))
    ap.add_argument('--jobqueue-workers', type=int,
                    default=int(os.environ.get('GRN_NODE_WORKERS', '2')))
    ap.add_argument('--slurm-cores', type=int, default=int(os.environ.get('GRN_SLURM_CORES', '8')))
    ap.add_argument('--slurm-memory', default=os.environ.get('GRN_SLURM_MEMORY', '32GB'))
    ap.add_argument('--slurm-time', default=os.environ.get('GRN_SLURM_TIME', '02:00:00'))
    ap.add_argument('--slurm-account', default=os.environ.get('GRN_SLURM_ACCOUNT'))
    ap.add_argument('--slurm-qos', default=os.environ.get('GRN_SLURM_QOS'))
    ap.add_argument('--slurm-partition', default=os.environ.get('GRN_SLURM_PARTITION'))
    ap.add_argument('--worker-timeout', type=int,
                    default=int(os.environ.get('GRN_WORKER_TIMEOUT', '600')))
    ap.add_argument('--local-directory', default=os.environ.get('GRN_LOCAL_DIR'))
    ap.add_argument('--seed', type=int, default=DEMON_SEED)
    # Smoke-test knobs: cap the target-gene count and/or the cell count so the
    # whole job can be validated end to end in minutes before spending hours on
    # the real run.  Leave unset for a real run.
    ap.add_argument('--max-targets', type=int,
                    default=int(os.environ.get('GRN_MAX_TARGETS', '0')) or None)
    ap.add_argument('--max-cells', type=int,
                    default=int(os.environ.get('GRN_MAX_CELLS', '0')) or None)
    args = ap.parse_args()

    n_workers = args.workers or max(1, args.cpus // args.threads_per_worker)

    exp = args.experiment
    down = f'{proj}/results/{exp}_grn'
    out = f'{down}/grn_tfonly_reg'
    os.makedirs(out, exist_ok=True)

    # Regulator / target lists live in the repo (data/), not at the repo root.
    tr_f = os.environ.get('TRRUST_FILE', f'{proj}/data/trrust_tf.txt')
    tf_f = os.environ.get('TF_ONLY_FILE', f'{proj}/data/tf_only.txt')
    for p in (tr_f, tf_f):
        if not os.path.exists(p):
            raise SystemExit(f'missing {p} (clone the repo including data/)')

    genes = np.load(f'{down}/genes.npy', allow_pickle=True).astype(str)
    Xmm = np.load(f'{down}/all_cells_gene_expression.npy', mmap_mode='r')
    print(f'expression matrix: {Xmm.shape} {Xmm.dtype}  genes: {len(genes)}', flush=True)

    trrust = [l.strip() for l in open(tr_f) if l.strip()]
    tfonly = [l.strip() for l in open(tf_f) if l.strip()]
    gset = set(genes)
    target_genes = [g for g in genes if g in set(trrust)]
    present_tfs = [g for g in genes if g in set(tfonly)]
    cols = list(dict.fromkeys(target_genes + present_tfs))     # union, order-preserving
    print(f'TRRUST list {len(trrust)} -> targets present {len(target_genes)}', flush=True)
    print(f'tf_only list {len(tfonly)} -> regulators present {len(present_tfs)}', flush=True)
    print(f'union columns -> {len(cols)}', flush=True)
    if not present_tfs:
        raise SystemExit('no regulators present -- wrong gene list or matrix')

    gidx = {g: i for i, g in enumerate(genes)}
    expr_sub = pd.DataFrame(np.asarray(Xmm[:, [gidx[g] for g in cols]], dtype=np.float32),
                            columns=cols)
    del Xmm
    if args.max_cells:
        expr_sub = expr_sub.iloc[:args.max_cells].reset_index(drop=True)
        print(f'[smoke] capped to {args.max_cells} cells', flush=True)
    if args.max_targets:
        keep_genes = set(target_genes[:args.max_targets]) | set(present_tfs)
        expr_sub = expr_sub[[c for c in cols if c in keep_genes]]
        print(f'[smoke] capped to {args.max_targets} target genes '
              f'(columns now {expr_sub.shape[1]})', flush=True)
    n_cells = expr_sub.shape[0]
    print(f'subset expression: {expr_sub.shape} '
          f'({expr_sub.memory_usage(deep=True).sum()/1e9:.2f} GB)', flush=True)

    from arboreto.algo import grnboost2

    client, cluster, mode = make_client(args, n_workers, args.threads_per_worker)
    n_online = len(client.scheduler_info()['workers'])
    try:
        t0 = time.time()
        network = grnboost2(expression_data=expr_sub, tf_names=present_tfs,
                            client_or_address=client, verbose=False, seed=args.seed)
        dt = time.time() - t0
    finally:
        client.close()
        if cluster is not None:
            cluster.close()

    print(f'GRNBoost2 done: {len(network)} edges in {dt:.1f}s '
          f'({dt/60:.1f} min) with {n_online} workers [{mode}]', flush=True)

    network.to_csv(f'{out}/grnboost2_network.tsv', sep='\t', index=False)
    network.to_csv(f'{out}/grnboost2_network.csv', index=False)
    with open(f'{out}/tf_regulators.txt', 'w') as f:
        f.write('\n'.join(present_tfs))
    with open(f'{out}/target_genes.txt', 'w') as f:
        f.write('\n'.join(target_genes))
    with open(f'{out}/all_columns.txt', 'w') as f:
        f.write('\n'.join(cols))
    with open(f'{out}/n_cells.txt', 'w') as f:
        f.write(str(n_cells))
    with open(f'{out}/run_info.txt', 'w') as f:
        f.write(f'experiment: {exp}\n')
        f.write(f'cells: {n_cells}\n')
        f.write(f'regulators (tf_only present): {len(present_tfs)}\n')
        f.write(f'targets (trrust present): {len(target_genes)}\n')
        f.write(f'union columns: {len(cols)}\n')
        f.write(f'edges: {len(network)}\n')
        f.write(f'dask mode: {mode}\n')
        f.write(f'dask workers online: {n_online}\n')
        f.write(f'threads_per_worker: {args.threads_per_worker}\n')
        f.write(f'seed: {args.seed}\n')
        f.write(f'wall_seconds: {dt:.1f}\n')
    print('DONE ->', out, flush=True)


if __name__ == '__main__':
    sys.exit(main())
