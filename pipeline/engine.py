#!/usr/bin/env python
"""Core engine for the multimodal GRN pipeline.  Imported by run.py.

Pipeline
    1. integrate()          scSAGA joint embedding H over a set of datasets
    2. build_reference()    "combined" (SCEMENT) or single-dataset RNA reference
    3. impute()             reverse-imputeKNN: propagate RNA expression to ATAC
    4. run_grn()            Arboreto GRNBoost2 on the all-cells matrix
    5. evaluate()           score the network against ground-truth edge lists

Nothing here is specific to a particular dataset, modality or organism: every
step reads its inputs from a config (see config.yml).
"""
from __future__ import annotations

import csv
import importlib.util
import json
import os
import sys
import time

import numpy as np
import scipy.io
import scipy.sparse as sp

K_NEIGHBOURS = 20          # reverse-imputeKNN k
GENE_CHUNK = 4000          # dense chunk size when building the RNA reference


# --------------------------------------------------------------------------- #
# config helpers
# --------------------------------------------------------------------------- #
def proj_root(explicit=None):
    return os.path.abspath(explicit or os.environ.get('PIPELINE_ROOT') or
                           os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))


def resolve(root, path):
    """Config paths are relative to PIPELINE_ROOT unless already absolute."""
    return path if os.path.isabs(path) else os.path.join(root, path)


def read_lines(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


class Dataset:
    """One dataset: a PCA matrix, plus (optionally) its raw expression."""

    def __init__(self, name, spec, root):
        self.name = name
        self.modality = spec.get('modality', 'rna')
        d = spec.get('dir')
        if d:
            d = resolve(root, d)
        self.pca = resolve(root, spec.get('pca') or f'{d}/pca_50.txt')
        self.counts = resolve(root, spec.get('counts') or f'{d}/counts.mtx') if d or spec.get('counts') else None
        self.barcodes = resolve(root, spec.get('barcodes') or f'{d}/barcodes.txt') if d or spec.get('barcodes') else None
        self.features = resolve(root, spec.get('features') or f'{d}/features.txt') if d or spec.get('features') else None

    # ---- introspection -----------------------------------------------------
    def has_pca(self):
        return os.path.exists(self.pca)

    def has_counts(self):
        return bool(self.counts) and os.path.exists(self.counts)

    def has_barcodes(self):
        return bool(self.barcodes) and os.path.exists(self.barcodes)

    def has_features(self):
        return bool(self.features) and os.path.exists(self.features)

    def has_expression(self):
        """True when counts+barcodes+features are all present."""
        return self.has_counts() and self.has_barcodes() and self.has_features()

    def n_cells(self):
        if self.barcodes and os.path.exists(self.barcodes):
            return sum(1 for _ in open(self.barcodes))
        if self.has_pca():
            return np.loadtxt(self.pca, comments='#', dtype=np.float32).shape[0]
        return None

    def check(self, need_expression=False):
        problems = []
        if not self.has_pca():
            problems.append(f'no PCA file ({self.pca})')
        if need_expression and not self.has_expression():
            miss = [p for p in (self.counts, self.barcodes, self.features)
                    if not p or not os.path.exists(p)]
            problems.append('missing expression files: ' + ', '.join(miss))
        return problems

    def load_pca(self):
        return np.loadtxt(self.pca, comments='#', dtype=np.float32)

    def load_barcodes(self):
        """Barcodes for this dataset.

        An ATAC dataset may legitimately have no barcodes file -- only its PCA
        is required.  In that case fall back to the PCA row order, so config
        needs nothing but `pca:`.
        """
        if self.barcodes and os.path.exists(self.barcodes):
            return read_lines(self.barcodes)
        if self.has_pca():
            n = np.loadtxt(self.pca, comments='#', dtype=np.float32).shape[0]
            prefix = self.name or 'cell'
            return [f'{prefix}_cell{i}' for i in range(n)]
        raise SystemExit(
            f'{self.name}: no barcodes file and no PCA to derive them from '
            f'(barcodes={self.barcodes}, pca={self.pca})')

    def load_features(self):
        return read_lines(self.features)

    def load_expression_genes_by_cells(self, n_features=None):
        """counts.mtx (cells x features) -> genes x cells float32 log1p(CPM/1e4).

        This is the transform applied to every reverse-imputeKNN reference (and to
        the real RNA rows of the all-cells matrix).  Chunked over genes to bound
        peak memory.
        """
        X = scipy.io.mmread(self.counts).tocsr()          # cells x features
        if n_features is None:
            n_features = X.shape[1]
        tot = np.asarray(X.sum(axis=1)).ravel()
        tot[tot == 0] = 1
        Xcsc = X.multiply((1e4 / tot)[:, None]).tocsc()
        out = np.empty((n_features, X.shape[0]), dtype=np.float32)
        for lo in range(0, n_features, GENE_CHUNK):
            hi = min(lo + GENE_CHUNK, n_features)
            sub = Xcsc[:, lo:hi].toarray().astype(np.float32)
            np.log1p(sub, out=sub)
            out[lo:hi, :] = sub.T
        return out


# --------------------------------------------------------------------------- #
# 1. scSAGA integration
# --------------------------------------------------------------------------- #
def integrate(datasets, anchor, outdir, scsaga_repo, params=None, scsaga_python=None):
    """Run scSAGA over `datasets` and write joint_embedding_H.npy.

    The anchor block is first in H, then the other datasets in `datasets` order.

    scSAGA's own CLI discards the joint embedding after computing it, so we call
    the Saga class directly and persist H ourselves.  That means NO patched
    scSAGA checkout is required -- any clone of the upstream repo works.
    """
    os.makedirs(outdir, exist_ok=True)
    params = dict(params or {})
    p = dict(M_samples=2000, alpha=0.75, S_iterations=25, gw_epsilon=1e-5, gw_reg=0.001)
    p.update({k: v for k, v in params.items() if v is not None})

    # s_shared_cells defaults to the smallest dataset (a partial alignment when
    # the modalities do not share cells).
    sizes = {d.name: d.n_cells() for d in datasets}
    n_shared = params.get('s_shared_cells') or min(s for s in sizes.values() if s)
    p['s_shared_cells'] = int(n_shared)

    import torch
    import sys
    # Prefer an installed scsaga package (pip install "git+https://github.com/
    # AluruLab/scSAGA.git"); fall back to a source checkout if config/env points
    # at one.  Either works -- the pipeline only needs `from scmint.scsaga import
    # Saga`, so there is no required manual clone.
    _saga = None
    if scsaga_repo and os.path.isdir(scsaga_repo):
        sys.path.insert(0, scsaga_repo)
    try:
        from scmint.scsaga import Saga as _saga
    except ImportError:
        _saga = None
    if _saga is None:
        raise SystemExit(
            'cannot import scmint.scsaga.\n'
            'Install it (recommended):\n'
            '  pip install "scsaga @ git+https://github.com/AluruLab/scSAGA.git"\n'
            'or point config.yml at a source checkout:\n'
            '  git clone https://github.com/AluruLab/scSAGA.git tools/scSAGA\n'
            '  scsaga_repo: tools/scSAGA'
        )
    Saga = _saga

    data_by_name = {}
    barcodes_by_name = {}
    for d in datasets:
        X = d.load_pca()
        data_by_name[d.name] = X
        if d.barcodes and os.path.exists(d.barcodes):
            bc = d.load_barcodes()
            if len(bc) == X.shape[0]:
                barcodes_by_name[d.name] = bc
        print(f'  {d.name}: PCA {X.shape}', flush=True)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'  device: {device}', flush=True)
    saga = Saga(device=device)

    t0 = time.time()
    T_dict, aligned = saga.run_multi(anchor_name=anchor, data_by_name=data_by_name,
                                    outdir=outdir, params=p)
    dt = time.time() - t0

    order = [anchor] + [d.name for d in datasets if d.name != anchor]
    H = np.vstack([aligned[n] for n in order]).astype(np.float32)
    np.save(f'{outdir}/joint_embedding_H.npy', H)
    all_bc = []
    for n in order:
        if n in barcodes_by_name:
            all_bc.extend(barcodes_by_name[n])
        np.save(f'{outdir}/aligned_{n}.npy', aligned[n].astype(np.float32))
    if len(all_bc) == H.shape[0]:
        with open(f'{outdir}/joint_embedding_barcodes.txt', 'w') as fh:
            fh.write('\n'.join(all_bc) + '\n')

    with open(f'{outdir}/integration_info.json', 'w') as fh:
        json.dump({'anchor': anchor, 'order': order, 'sizes': sizes,
                   'H_shape': list(H.shape), 's_shared_cells': p['s_shared_cells'],
                   'params': {k: v for k, v in p.items()}, 'seconds': round(dt, 1)},
                  fh, indent=2)
    print(f'  [OK] joint_embedding_H.npy {H.shape} in {dt:.1f}s', flush=True)
    return H, order, sizes


# --------------------------------------------------------------------------- #
# 2. RNA reference (single dataset, or SCEMENT-combined)
# --------------------------------------------------------------------------- #
def scsaga_paths(root):
    return os.path.join(root, 'pipeline')


def _load_scement(root):
    """Load the vendored pure-Python sct_sparse (numpy2-safe)."""
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scement_sparse.py')
    if not os.path.exists(src):
        raise SystemExit(f'SCEMENT helper not found: {src}')
    spec = importlib.util.spec_from_file_location('scement_sparse', src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_reference(spec, rna_datasets, root, outdir):
    """Return (expr genes x cells, H cells x 30, barcodes) for one strategy.

    spec is either {'dataset': name} or {'combine': 'all_rna'}.
    """
    if 'dataset' in spec:
        name = spec['dataset']
        d = rna_datasets[name]
        genes = d.load_features()
        expr = d.load_expression_genes_by_cells(len(genes))
        H = np.load(f'{outdir}/aligned_{name}.npy')
        bc = d.load_barcodes()
        print(f'  [ref] {name}: {expr.shape[0]} genes x {expr.shape[1]} cells', flush=True)
        return expr, H, bc, list(genes)

    # ---- SCEMENT-combined reference ---------------------------------------
    names = list(rna_datasets)
    if len(names) < 2:
        raise SystemExit(f"'combine' needs at least 2 RNA datasets, got {names}")
    print(f'  [ref] SCEMENT-combining {names}', flush=True)
    genes = rna_datasets[names[0]].load_features()
    Hs, bcs = [], []
    in_args = []
    for n in names:
        d = rna_datasets[n]
        g = d.load_features()
        if g != genes:
            raise SystemExit(
                f"SCEMENT combine needs identical gene ordering; {n} differs from "
                f"{names[0]}.  Reference datasets must share a feature list."
            )
        Hs.append(np.load(f'{outdir}/aligned_{n}.npy'))
        bcs.append(d.load_barcodes())
        in_args += ['--in', f'{n}={d.counts}:{d.features}:{d.barcodes}']

    # SCEMENT needs anndata/scanpy, which the main pipeline env does not carry,
    # so it runs as a subprocess in its own environment.
    npz = f'{outdir}/_scement_combined.npz'
    import subprocess
    here = os.path.dirname(os.path.abspath(__file__))
    py = os.environ.get('SCEMENT_PYTHON')
    cmd = ([py] if py else [sys.executable]) + [
        os.path.join(here, 'scement_combine.py'), '--out', npz] + in_args
    print(f'  [ref] $ {" ".join(cmd)}', flush=True)
    env = dict(os.environ)
    r = subprocess.run(cmd, env=env)
    if r.returncode != 0:
        raise SystemExit(
            'SCEMENT combining failed.  It needs anndata + scanpy; run that step '
            'with an interpreter that has them by setting SCEMENT_PYTHON, e.g.\n'
            '  export SCEMENT_PYTHON=/path/to/scement/env/bin/python'
        )
    z = np.load(npz, allow_pickle=True)
    expr = z['expression'].astype(np.float32)
    bc = [str(x) for x in z['barcodes']]
    genes = [str(x) for x in z['genes']]
    H = np.vstack(Hs).astype(np.float32)
    print(f'  [ref] combined: {expr.shape[0]} genes x {expr.shape[1]} cells',
          flush=True)
    return expr, H, bc, genes


# --------------------------------------------------------------------------- #
# 3. reverse-imputeKNN
# --------------------------------------------------------------------------- #
def impute(H, order, sizes, ref_expr, ref_H, ref_bc, query_names, outdir,
           query_barcodes, rna_names):
    """Propagate reference RNA expression onto the query (ATAC) cells.

    For each query cell: k=20 nearest reference cells in H, weights exp(-dist)
    normalised per row, then imputed = ref_expr @ W.T.

    all-cells rows = [ every real RNA cell in H order, then every query cell ].
    The real RNA expression is taken from whichever datasets are RNA in the
    integration; the `ref_*` arguments only decide what gets propagated to ATAC.

    Returns (imputed genes x n_query, all_expr n_all x genes, all_barcodes, genes)
    """
    os.makedirs(outdir, exist_ok=True)
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse import coo_matrix

    # H block slice for a dataset name
    def block(name):
        start = sum(sizes[n] for n in order[:order.index(name)])
        return H[start:start + sizes[name]]

    Hq = np.vstack([block(q) for q in query_names])
    q_bc = [b for q in query_names for b in query_barcodes[q]]
    assert Hq.shape[0] == len(q_bc), (Hq.shape, len(q_bc))

    nbrs = NearestNeighbors(n_neighbors=K_NEIGHBOURS, metric='euclidean').fit(ref_H)
    dist, idx = nbrs.kneighbors(Hq)
    e = np.exp(-dist)
    Wk = e / e.sum(axis=1, keepdims=True)
    rows = np.repeat(np.arange(Hq.shape[0]), K_NEIGHBOURS)
    W = coo_matrix((Wk.ravel(), (rows, idx.ravel())),
                   shape=(Hq.shape[0], ref_H.shape[0])).tocsr()
    imputed = (ref_expr @ W.T).astype(np.float32)           # genes x queries
    print(f'  imputed {imputed.shape[0]} genes x {imputed.shape[1]} ATAC cells '
          f'(mean nn dist {float(dist.mean()):.4f}, '
          f'ref cells used {len(np.unique(idx))}/{ref_H.shape[0]})', flush=True)

    np.save(f'{outdir}/imputed_expression_genes_x_atac.npy', imputed)
    with open(f'{outdir}/imputed_atac_barcodes.txt', 'w') as fh:
        fh.write('\n'.join(q_bc) + '\n')
    return imputed, q_bc


# --------------------------------------------------------------------------- #
# 4. Arboreto GRNBoost2
# --------------------------------------------------------------------------- #
def run_grn(expr_all, genes, reg_file, tgt_file, outdir, n_workers=8,
            threads_per_worker=1, seed=666, scheduler=None,
            max_targets=None, max_regulators=None):
    """GRNBoost2 on the all-cells matrix.

    Regulators (tf_names) = TFs from `reg_file` present in the gene list.
    Targets               = genes from `tgt_file` present in the gene list.
    The matrix handed to GRNBoost2 is the UNION of both, because every regulator
    must also be a column.

    IMPORTANT -- how much work this actually is:
        grnboost2() does NOT accept a target-gene list.  It forwards to
        create_graph(), whose `target_genes` defaults to 'all', so a regression
        is fitted for EVERY COLUMN of the matrix -- including the regulator
        columns.  Runtime therefore scales with the number of COLUMNS
        (targets + regulators), not with the number of targets.

        In this pipeline that is ~2,852 regressions (2,827 targets + 816
        regulators, union 2,852), not 2,827.

    max_targets / max_regulators exist only to make smoke tests cheap: they prune
    the column set, which is the only way to reduce the regression count.  Leave
    both unset for a real run.

    Parallelism: Arboreto takes a pre-built dask Client via client_or_address,
    so the worker count is a runtime choice.
    """
    os.makedirs(outdir, exist_ok=True)
    import pandas as pd
    from arboreto.algo import grnboost2
    from dask.distributed import Client, LocalCluster

    regs = [l.strip() for l in open(reg_file) if l.strip()]
    tgts = [l.strip() for l in open(tgt_file) if l.strip()]
    genes = list(genes)
    gset = set(genes)
    target_genes = [g for g in genes if g in set(tgts)]
    present_tfs = [g for g in genes if g in set(regs)]
    if max_targets:
        target_genes = target_genes[:max_targets]
    if max_regulators:
        present_tfs = present_tfs[:max_regulators]
    cols = list(dict.fromkeys(target_genes + present_tfs))
    print(f'  regulators present {len(present_tfs)}, targets present '
          f'{len(target_genes)}, union columns {len(cols)}', flush=True)
    print(f'  NOTE: GRNBoost2 fits one regression per COLUMN -> '
          f'{len(cols)} regressions (regulator columns are fitted too)', flush=True)
    if not present_tfs:
        raise SystemExit('no regulators present -- check the regulators file')

    gene_order = {g: i for i, g in enumerate(genes)}
    X = np.load(expr_all, mmap_mode='r') if isinstance(expr_all, str) else expr_all
    sub = pd.DataFrame(np.asarray(X[:, [gene_order[g] for g in cols]], dtype=np.float32),
                       columns=cols)
    del X
    n_cells = sub.shape[0]
    print(f'  GRN input: {sub.shape} ({sub.memory_usage(deep=True).sum()/1e9:.2f} GB)',
          flush=True)

    close_cluster = None
    if scheduler:
        client = Client(scheduler)
        print(f'  dask: external scheduler {scheduler}', flush=True)
    else:
        close_cluster = LocalCluster(n_workers=n_workers,
                                     threads_per_worker=threads_per_worker,
                                     processes=True)
        client = Client(close_cluster)
        print(f'  dask: LocalCluster {n_workers} workers x {threads_per_worker} '
              f'threads', flush=True)

    try:
        t0 = time.time()
        net = grnboost2(expression_data=sub, tf_names=present_tfs,
                        client_or_address=client, verbose=False, seed=seed)
        dt = time.time() - t0
    finally:
        client.close()
        if close_cluster is not None:
            close_cluster.close()

    n_online = len(sub.columns)
    print(f'  GRNBoost2: {len(net)} edges in {dt:.1f}s ({dt/60:.1f} min)', flush=True)
    net.to_csv(f'{outdir}/network.tsv', sep='\t', index=False)
    for name, vals in [('tf_regulators.txt', present_tfs),
                       ('target_genes.txt', target_genes),
                       ('all_columns.txt', cols)]:
        with open(f'{outdir}/{name}', 'w') as fh:
            fh.write('\n'.join(vals))
    with open(f'{outdir}/n_cells.txt', 'w') as fh:
        fh.write(str(n_cells))
    with open(f'{outdir}/run_info.txt', 'w') as fh:
        fh.write(f'cells: {n_cells}\ncolumns: {len(cols)}\n'
                 f'regulators: {len(present_tfs)}\ntargets: {len(target_genes)}\n'
                 f'edges: {len(net)}\nseed: {seed}\nworkers: {n_workers}\n'
                 f'wall_seconds: {dt:.1f}\n')
    return net


# --------------------------------------------------------------------------- #
# 5. evaluation
# --------------------------------------------------------------------------- #
def evaluate(network_tsv, gt_files, outdir):
    os.makedirs(outdir, exist_ok=True)
    import pandas as pd
    net = pd.read_csv(network_tsv, sep='\t')
    tfcol = 'TF' if 'TF' in net.columns else net.columns[0]
    tgcol = 'target' if 'target' in net.columns else net.columns[1]
    net = net.sort_values('importance', ascending=False).reset_index(drop=True)
    edge_set = set(zip(net[tfcol], net[tgcol]))
    results = {}
    with open(f'{outdir}/evaluation_summary.txt', 'w') as fh:
        fh.write(f'inferred edges: {len(net)}\n')
        fh.write('=' * 70 + '\n')
        for name, path in gt_files.items():
            if not path or not os.path.exists(path):
                print(f'  skipping {name}: not found', flush=True)
                continue
            gt = [(r['TF'], r['TARGET']) for r in csv.DictReader(open(path))]
            dedup = set(gt)
            recovered = len(edge_set & dedup)
            pct = 100.0 * recovered / len(dedup) if dedup else 0.0
            fh.write(f'\n===== vs {name} =====\n')
            fh.write(f'total edges: {len(gt)}\ndeduplicated: {len(dedup)}\n')
            fh.write(f'recovered: {recovered} ({pct:.1f}%)\n')
            for K in sorted({100, 1000, 10000, len(net)}):
                top = set(zip(net.head(K)[tfcol], net.head(K)[tgcol]))
                fh.write(f'  top-{K:<7} prec@K {len(top & dedup)/K:.4f}  '
                         f'recall@K {len(top & dedup)/len(dedup):.4f}\n')
            results[name] = {'dedup': len(dedup), 'recovered': recovered,
                             'pct': round(pct, 1)}
            print(f'  {name}: {recovered}/{len(dedup)} ({pct:.1f}%)', flush=True)
    net.head(25).to_csv(f'{outdir}/top_edges.csv', index=False)
    with open(f'{outdir}/evaluation.json', 'w') as fh:
        json.dump({'edges': len(net), 'ground_truths': results}, fh, indent=2)
    return results


# --------------------------------------------------------------------------- #
# 5. Compare reference strategies
# --------------------------------------------------------------------------- #
def compare_strategies(exp_dir, references, root, ground_truth):
    """Compare the networks from every reference strategy of one experiment.

    An experiment with N RNA datasets produces N+1 strategies (one SCEMENT
    combination, one per individual RNA dataset).  This writes

        results/<exp>/comparison.md     human-readable table
        results/<exp>/comparison.json   machine-readable, for plotting

    covering, per strategy:
      * network size
      * edge overlap with every other strategy (top-100 shared + Jaccard)
      * recovery against each ground truth, when ground truth is configured

    Ground truth is optional: without it the comparison still reports network
    size and pairwise overlap, which is what shows whether the reference choice
    changes the inferred network at all.
    """
    import pandas as pd
    out = {}
    for strat in references:
        net_path = f'{exp_dir}/{strat}/grn_tfonly_reg/network.tsv'
        if not os.path.exists(net_path):
            print(f'  [compare] skipping {strat}: no network.tsv', flush=True)
            continue
        net = pd.read_csv(net_path, sep='\t')
        tfcol = 'TF' if 'TF' in net.columns else net.columns[0]
        tgcol = 'target' if 'target' in net.columns else net.columns[1]
        net = net.sort_values('importance', ascending=False).reset_index(drop=True)
        edges = list(zip(net[tfcol], net[tgcol]))
        gj = f'{exp_dir}/{strat}/grn_tfonly_reg/evaluation/evaluation.json'
        ev = json.load(open(gj)) if os.path.exists(gj) else {}
        out[strat] = {'edges': len(edges), 'edge_set': set(edges),
                      'top100': set(edges[:100]),
                      'eval': ev.get('ground_truths', {})}

    if not out:
        print('  [compare] no networks found to compare', flush=True)
        return {}

    names = list(out)

    # ---- markdown ---------------------------------------------------------
    lines = ['# Reference-strategy comparison', '',
             f'Strategies: {", ".join(names)}', '']
    lines += ['## Network size', '', '| strategy | edges |', '|---|---|']
    for s in names:
        lines.append(f'| {s} | {out[s]["edges"]:,} |')

    lines += ['', '## Overlap between strategies (top-100 edges)', '',
              '| A | B | shared | Jaccard |', '|---|---|---|---|']
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            inter = len(out[a]['top100'] & out[b]['top100'])
            union = len(out[a]['top100'] | out[b]['top100'])
            jac = inter / union if union else 0.0
            lines.append(f'| {a} | {b} | {inter} | {jac:.3f} |')

    gt_names = set()
    for s in names:
        gt_names.update(out[s]['eval'])
    for g in sorted(gt_names):
        lines += ['', f'## Recovery vs {g}', '',
                  '| strategy | recovered | % |', '|---|---|---|']
        for s in names:
            r = out[s]['eval'].get(g)
            if r:
                lines.append(f'| {s} | {r["recovered"]}/{r["dedup"]} | {r["pct"]}% |')
        best = max((out[s]['eval'].get(g, {}).get('pct', -1) for s in names),
                   default=None)
        if best is not None and best >= 0:
            winners = [s for s in names
                       if out[s]['eval'].get(g, {}).get('pct', -1) == best]
            lines += ['', f'-> best: {", ".join(winners)} at {best}%']

    with open(f'{exp_dir}/comparison.md', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    # ---- json -------------------------------------------------------------
    js = {'strategies': {s: {'edges': out[s]['edges'],
                             'recovered': {g: out[s]['eval'].get(g, {}).get('recovered')
                                           for g in sorted(gt_names)}}
                         for s in names},
          'overlap_top100': {f'{a}|{b}': len(out[a]['top100'] & out[b]['top100'])
                             for i, a in enumerate(names)
                             for b in names[i + 1:]}}
    with open(f'{exp_dir}/comparison.json', 'w') as fh:
        json.dump(js, fh, indent=2)

    print(f'  [compare] wrote {exp_dir}/comparison.md and comparison.json',
          flush=True)
    return js
