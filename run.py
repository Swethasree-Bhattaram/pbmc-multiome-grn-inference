#!/usr/bin/env python
"""Multimodal GRN inference: integration -> imputation -> GRN -> evaluation.

One entry point, one config file.

    python run.py                 # every experiment in config.yml
    python run.py --check         # validate paths only, no compute
    python run.py --list          # show what config declares
    python run.py -e my_exp       # one experiment
    python run.py -e my_exp --only grn,evaluate    # reuse earlier steps

The five steps, once per experiment:

    integrate   scSAGA joint embedding over all datasets
    reference   build the RNA reference (one dataset, or all RNA combined)
    impute      reverse-imputeKNN: propagate reference RNA onto the ATAC cells
    grn         Arboreto GRNBoost2 on the all-cells matrix
    evaluate    score the network against ground-truth edge lists
    compare     (automatic when an experiment has >1 reference)

An experiment with N RNA datasets gets N+1 reference strategies for free: one
combining every RNA dataset, plus one per individual dataset.  Nothing here is
specific to a dataset, an organism or a modality -- all paths come from config.

Everything runs in ONE environment (see requirements.txt / setup.sh).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

import numpy as np
import scipy.io
import scipy.sparse as sp
import yaml

K_NEIGHBOURS = 20          # reverse-imputeKNN k
GENE_CHUNK = 4000          # dense chunk when building an RNA reference
STEPS = ['integrate', 'impute', 'grn', 'evaluate']

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def read_lines(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


def feature_symbols(path):
    """Feature file -> gene symbols (one per feature, in file order).

    Cell Ranger writes feature lists in three shapes and all three show up in
    real data:

        SYMBOL                                  old cellranger mtx dirs
        ENSG00000000003<TAB>SYMBOL              mtx gene_id/gene_name
        ENSG00000000003<TAB>SYMBOL<TAB>Type     features.tsv.gz, cellranger >= 3
                                                (the Type column is
                                                 'Gene Expression' for RNA)

    Gene lists (tf_only.txt, trrust_tf.txt) hold BARE SYMBOLS, so the symbol
    field is the only one that can ever match.  Using the raw line instead
    matches nothing for the 2- and 3-column forms, which makes every regulator
    and target look absent.
    """
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split('\t') if '\t' in line else line.split()
            out.append(parts[1] if len(parts) >= 2 else parts[0])
    return out


def abspath(root, path):
    return path if os.path.isabs(path) else os.path.join(root, path)


class Dataset:
    """One dataset: a PCA matrix, plus (optionally) its raw expression."""

    def __init__(self, name, spec, root):
        self.name = name
        self.modality = (spec.get('modality') or 'rna').lower()
        if self.modality not in ('rna', 'atac'):
            raise SystemExit(f'{name}: modality must be "rna" or "atac", '
                             f'got {self.modality!r}')
        # `dir` is a shorthand: the four files live together and are named
        # pca_50.txt / counts.mtx / barcodes.txt / features.txt.  Any of them
        # can still be overridden individually (e.g. a PCA stored elsewhere).
        d = abspath(root, spec['dir']) if spec.get('dir') else None

        def path(key, filename):
            if spec.get(key):
                return abspath(root, spec[key])
            return os.path.join(d, filename) if d else None

        self.pca = path('pca', 'pca_50.txt')
        self.counts = path('counts', 'counts.mtx')
        self.barcodes = path('barcodes', 'barcodes.txt')
        self.features = path('features', 'features.txt')
        if not self.pca:
            raise SystemExit(f'{name}: give either "dir", or "pca" (path to a '
                             f'cells x N PCA matrix)')

    # ---- presence ---------------------------------------------------------
    def has_pca(self):
        return os.path.exists(self.pca)

    def has_expression(self):
        return all(p and os.path.exists(p)
                   for p in (self.counts, self.barcodes, self.features))

    def needs_expression(self):
        """Only RNA datasets are read as expression.  ATAC expression is
        imputed, so an ATAC dataset needs nothing but its PCA."""
        return self.modality == 'rna'

    def missing(self):
        """Missing files that this dataset actually needs."""
        out = []
        if not self.has_pca():
            out.append(f'pca ({self.pca})')
        if self.needs_expression() and not self.has_expression():
            for label, p in (('counts', self.counts), ('barcodes', self.barcodes),
                             ('features', self.features)):
                if not p:
                    out.append(f'{label} (not set)')
                elif not os.path.exists(p):
                    out.append(f'{label} ({p})')
        return out

    # ---- load -------------------------------------------------------------
    def load_pca(self):
        return np.loadtxt(self.pca, comments='#', dtype=np.float32)

    def n_cells(self):
        if self.barcodes and os.path.exists(self.barcodes):
            return sum(1 for _ in open(self.barcodes))
        return self.load_pca().shape[0]

    def load_barcodes(self):
        if self.barcodes and os.path.exists(self.barcodes):
            bc = read_lines(self.barcodes)
            n = self.n_cells_pca()
            if len(bc) != n:
                raise SystemExit(f'{self.name}: {len(bc)} barcodes but '
                                 f'{n} rows in {self.pca}')
            return bc
        # ATAC without barcodes is fine -- derive stable ids from PCA order.
        return [f'{self.name}_cell{i}' for i in range(self.n_cells_pca())]

    def n_cells_pca(self):
        return self.load_pca().shape[0]

    def mtx_shape(self):
        """Read just the MatrixMarket header -> (rows, cols, nnz).

        Our layout is cells x features, so rows == cells and cols == features.
        """
        with open(self.counts) as fh:
            for line in fh:
                if line.startswith('%'):
                    continue
                r, c, nnz = line.split()[:3]
                return int(r), int(c), int(nnz)
        raise SystemExit(f'{self.name}: {self.counts} has no dimension line')

    def alignment_problems(self):
        """Cheap consistency check: are the four files describing the same cells?

        A dataset is only meaningful if the PCA rows, the counts rows and the
        barcodes all refer to the SAME cells in the SAME order.  Getting this
        wrong silently misaligns the integration, so it is worth checking.

        ATAC datasets are exempt: their expression is imputed, so counts are
        never read and their orientation is irrelevant.  Flagging them here
        would report a problem the pipeline cannot actually be affected by.
        """
        problems = []
        if not self.needs_expression():
            return problems
        if not (self.has_pca() and self.has_expression()):
            return problems
        rows, cols, _ = self.mtx_shape()
        n_pca = self.n_cells_pca()
        n_bc = sum(1 for _ in open(self.barcodes))
        n_feat = sum(1 for _ in open(self.features))
        if rows != n_pca:
            problems.append(f'counts has {rows} rows but PCA has {n_pca} cells')
        if rows != n_bc:
            problems.append(f'counts has {rows} rows but {n_bc} barcodes')
        if cols != n_feat:
            problems.append(f'counts has {cols} columns but {n_feat} features')
        return problems

    def load_features(self):
        """Gene SYMBOLS, in file order (see feature_symbols).

        These names index the matrix columns and are matched against
        tf_only.txt / trrust_tf.txt, so they must be bare symbols -- a feature
        file with gene_id/type columns would otherwise match nothing.
        """
        if not self.features or not os.path.exists(self.features):
            raise SystemExit(f'{self.name}: needs "features" to be used as an '
                             f'RNA reference')
        return feature_symbols(self.features)

    def load_expression(self, n_features=None):
        """counts.mtx (cells x features) -> genes x cells, log1p(CPM/1e4).

        Identical transform for every reverse-imputeKNN reference and for the
        real RNA rows of the all-cells matrix.  Chunked over genes to bound
        peak memory.
        """
        X = scipy.io.mmread(self.counts).tocsr()          # cells x features
        if n_features is None:
            n_features = X.shape[1]
        if X.shape[1] != n_features:
            raise SystemExit(f'{self.name}: counts has {X.shape[1]} columns but '
                             f'{n_features} features were listed')
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
# integrate: scSAGA joint embedding
# --------------------------------------------------------------------------- #
def integrate(datasets, anchor, outdir, params):
    """Run scSAGA over `datasets`; write joint_embedding_H.npy.

    scSAGA's CLI computes the joint embedding and discards it, so we call the
    Saga class directly and save H ourselves -- any upstream clone works, no
    patched fork.
    """
    import torch
    from scmint.scsaga import Saga

    try:
        os.makedirs(outdir, exist_ok=True)
        p = dict(M_samples=2000, alpha=0.75, S_iterations=25,
                 gw_epsilon=1e-5, gw_reg=0.001)
        p.update({k: v for k, v in (params or {}).items() if v is not None})

        # s_shared_cells defaults to the smallest dataset (a partial alignment
        # when modalities do not share cells).  Set it explicitly to align on
        # more cells, e.g. the full cell count for paired multiome data.
        data_by_name, barcodes_by_name = {}, {}
        for d in datasets:
            X = d.load_pca()
            data_by_name[d.name] = X
            if d.barcodes and os.path.exists(d.barcodes):
                bc = d.load_barcodes()
                if len(bc) == X.shape[0]:
                    barcodes_by_name[d.name] = bc
            print(f'  {d.name}: PCA {X.shape}', flush=True)

        # Block sizes come from the matrices actually handed to scSAGA, so H
        # slicing downstream can never disagree with the data.
        sizes = {n: X.shape[0] for n, X in data_by_name.items()}
        p['s_shared_cells'] = int(params.get('s_shared_cells')
                                  or min(sizes.values()))

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(f'  device: {device}', flush=True)
        saga = Saga(device=device)

        t0 = time.time()
        _T, aligned = saga.run_multi(anchor_name=anchor, data_by_name=data_by_name,
                                     outdir=outdir, params=p)
        dt = time.time() - t0

        # H: anchor block first, then the other datasets in declaration order.
        order = [anchor] + [d.name for d in datasets if d.name != anchor]
        H = np.vstack([aligned[n] for n in order]).astype(np.float32)
        np.save(f'{outdir}/joint_embedding_H.npy', H)
        for n in order:
            np.save(f'{outdir}/aligned_{n}.npy', aligned[n].astype(np.float32))

        all_bc = [b for n in order for b in barcodes_by_name.get(n, [])]
        if len(all_bc) == H.shape[0]:
            with open(f'{outdir}/joint_embedding_barcodes.txt', 'w') as fh:
                fh.write('\n'.join(all_bc) + '\n')

        with open(f'{outdir}/integration_info.json', 'w') as fh:
            json.dump({'anchor': anchor, 'order': order, 'sizes': sizes,
                       'H_shape': list(H.shape), 'params': p,
                       'seconds': round(dt, 1)}, fh, indent=2)
        print(f'  [OK] joint_embedding_H.npy {H.shape} in {dt:.1f}s', flush=True)
        return H, order, sizes
    except Exception as exc:                      # noqa: BLE001
        raise SystemExit(f'scSAGA integration failed: {exc}') from exc


# --------------------------------------------------------------------------- #
# gene-axis alignment
# --------------------------------------------------------------------------- #
def shared_gene_axis(gene_lists, mode='union'):
    """One gene ordering covering several datasets.

    SCEMENT (like ComBat generally) needs every dataset on the SAME gene axis
    before the counts can be stacked.  Its own evaluation pipelines do this with
    AnnData's `concat(..., merge='same')`, i.e. a union, and ship both variants:

        combat_union.json / combat_intersect.json

    union     every gene seen in any dataset, in first-seen order.  A dataset
              missing a gene contributes zeros there.  Keeps the widest gene
              coverage; matches merge='same'.
    intersect only genes present in ALL datasets.  No invented zeros, but a
              smaller gene set.
    """
    if mode not in ('union', 'intersect'):
        raise SystemExit(f'gene axis mode must be union|intersect, got {mode!r}')
    if mode == 'union':
        seen, out = set(), []
        for genes in gene_lists:
            for g in genes:
                if g not in seen:
                    seen.add(g)
                    out.append(g)
        return out
    common = set(gene_lists[0])
    for genes in gene_lists[1:]:
        common &= set(genes)
    return [g for g in gene_lists[0] if g in common]


def reindex_expression(expr, src_genes, dst_genes):
    """genes x cells -> dst_genes x cells.

    Genes present in dst_genes but absent from src_genes come out as zeros, so
    datasets measured against different references can share one axis.
    """
    out = np.zeros((len(dst_genes), expr.shape[1]), dtype=np.float32)
    pos = {g: i for i, g in enumerate(src_genes)}
    dst_at, src_at = [], []
    for j, g in enumerate(dst_genes):
        i = pos.get(g)
        if i is not None:
            dst_at.append(j)
            src_at.append(i)
    if dst_at:
        out[dst_at, :] = expr[src_at, :]
    return out


def reindex_counts_columns(X, src_genes, dst_genes):
    """cells x genes sparse -> cells x dst_genes, gene-aligned (zeros added)."""
    pos = {g: i for i, g in enumerate(src_genes)}
    dst_at, src_at = [], []
    for j, g in enumerate(dst_genes):
        i = pos.get(g)
        if i is not None:
            dst_at.append(j)
            src_at.append(i)
    if not dst_at:
        return sp.csr_matrix((X.shape[0], len(dst_genes)), dtype=X.dtype)
    sub = X[:, src_at].tocoo()
    col = np.asarray([dst_at[c] for c in sub.col], dtype=np.int32)
    return sp.csr_matrix((sub.data, (sub.row, col)),
                         shape=(X.shape[0], len(dst_genes)), dtype=X.dtype)



def _combined_reference(rna_datasets, counts_by_name=None, gene_axis='union'):
    """ComBat-combine several RNA datasets into one batch-corrected reference.

    genes x cells, clipped at 0.  Uses the vendored SCEMENT sct_sparse (see
    scement.py), which reproduces the published SCEMENT numbers exactly.
    scanpy's pp.combat is a DIFFERENT implementation: on real PBMC counts it
    differs by mean |delta| 9.5e-4, with 10.6% of nonzeros off by more than
    1e-3 -- so switching implementation would change published results.

    Datasets measured against different references do NOT need identical feature
    lists: they are first put on one shared gene axis (`gene_axis`), exactly as
    SCEMENT's own evaluation pipelines do with AnnData's
    `concat(..., merge='same')`.  Genes a dataset lacks contribute zeros.
    """
    import anndata as ad
    import pandas as pd

    import scement

    names = list(rna_datasets)
    if len(names) < 2:
        raise SystemExit(f"'combine' needs at least 2 RNA datasets, got {names}")

    # ---- put every dataset on one gene axis ------------------------------
    gene_lists = [rna_datasets[n].load_features() for n in names]
    axis = shared_gene_axis(gene_lists, gene_axis)
    if not axis:
        raise SystemExit("'combine' found no genes in common across "
                         f"{names} (gene_axis={gene_axis})")
    dropped = [len(set(g) - set(axis)) for g in gene_lists]
    print(f'  [ref] gene axis: {gene_axis}, {len(axis)} genes '
          f'(dataset sizes {[len(g) for g in gene_lists]}, '
          f'dropped per dataset {dropped})', flush=True)

    blocks, bcs, batches = [], [], []
    for n, genes in zip(names, gene_lists):
        d = rna_datasets[n]
        X = scipy.io.mmread(d.counts).tocsr()             # cells x genes
        if X.shape[1] != len(genes):
            raise SystemExit(f'{n}: counts has {X.shape[1]} columns but '
                             f'{len(genes)} features were listed')
        bc = d.load_barcodes()
        if len(bc) != X.shape[0]:
            raise SystemExit(f'{n}: {X.shape[0]} count rows but {len(bc)} barcodes')
        Xa = reindex_counts_columns(X.astype(np.float32), genes, axis)
        blocks.append(Xa.tocsr())
        bcs.extend(bc)
        batches.extend([n] * X.shape[0])
        print(f'  [ref] {n}: {X.shape[0]} cells x {len(genes)} genes '
              f'-> {Xa.shape[1]} on the shared axis', flush=True)

    X = sp.vstack(blocks).tocsr()
    adata = ad.AnnData(X=X, obs=pd.DataFrame({'batch': pd.Categorical(batches)}))
    print(f'  [ref] ComBat on {adata.shape} '
          f'{dict(zip(*np.unique(batches, return_counts=True)))}', flush=True)
    corrected = np.clip(np.asarray(scement.sct_sparse(adata, key='batch',
                                                      inplace=False)), 0, None)
    expr = np.ascontiguousarray(corrected.T, dtype=np.float32)   # genes x cells
    return expr, list(axis), bcs


def build_reference(ref_spec, rna_datasets, integ_dir):
    """-> (expr genes x cells, H cells x d, barcodes, genes) for one strategy.

    ref_spec is {'dataset': name} or {'combine': ...}.
    """
    if 'dataset' in ref_spec:
        name = ref_spec['dataset']
        if name not in rna_datasets:
            raise SystemExit(f'strategy references {name!r}, which is not an RNA '
                             f'dataset of this experiment')
        d = rna_datasets[name]
        genes = d.load_features()
        expr = d.load_expression(len(genes))
        H = np.load(f'{integ_dir}/aligned_{name}.npy')
        print(f'  [ref] {name}: {expr.shape[1]} cells x {expr.shape[0]} genes',
              flush=True)
        return expr, H, d.load_barcodes(), genes

    expr, genes, bc = _combined_reference(rna_datasets)
    H = np.vstack([np.load(f'{integ_dir}/aligned_{n}.npy')
                   for n in rna_datasets]).astype(np.float32)
    return expr, H, bc, genes


# --------------------------------------------------------------------------- #
# impute: reverse-imputeKNN
# --------------------------------------------------------------------------- #
def impute(H, order, sizes, ref_expr, ref_H, query_names, outdir, query_barcodes):
    """Propagate reference RNA expression onto the query (ATAC) cells.

    Per query cell: k=20 nearest reference cells in H, weights exp(-distance)
    normalised per row, imputed = ref_expr @ W.T.
    """
    os.makedirs(outdir, exist_ok=True)
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse import coo_matrix

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
    imputed = (ref_expr @ W.T).astype(np.float32)             # genes x queries
    print(f'  imputed {imputed.shape[0]} genes x {imputed.shape[1]} ATAC cells '
          f'(mean nn dist {float(dist.mean()):.4f}, ref cells used '
          f'{len(np.unique(idx))}/{ref_H.shape[0]})', flush=True)

    np.save(f'{outdir}/imputed_expression_genes_x_atac.npy', imputed)
    with open(f'{outdir}/imputed_atac_barcodes.txt', 'w') as fh:
        fh.write('\n'.join(q_bc) + '\n')
    return imputed, q_bc


# --------------------------------------------------------------------------- #
# grn: Arboreto GRNBoost2
# --------------------------------------------------------------------------- #
def run_grn(expr_all, genes, reg_file, tgt_file, outdir, n_workers=8,
            threads_per_worker=1, seed=666, scheduler=None,
            max_columns=None):
    """GRNBoost2 on the all-cells matrix.

    Regulators (tf_names) = TFs from `reg_file` present in the gene list.
    Targets               = genes from `tgt_file` present in the gene list.
    The matrix is the UNION of both, because every regulator must also be a
    column.

    grnboost2() takes no target-gene list: it forwards to create_graph(), whose
    target_genes defaults to 'all', so ONE REGRESSION IS FITTED PER COLUMN --
    regulator columns included.  Runtime therefore scales with the number of
    COLUMNS, not targets.  `max_columns` exists only to make smoke tests cheap.
    """
    import grn_compat
    grn_compat.apply()
    import pandas as pd
    from arboreto.algo import grnboost2
    from dask.distributed import Client, LocalCluster

    os.makedirs(outdir, exist_ok=True)
    genes = list(genes)
    regs = set(read_lines(reg_file))
    tgts = set(read_lines(tgt_file))
    target_genes = [g for g in genes if g in tgts]
    present_tfs = [g for g in genes if g in regs]
    cols = list(dict.fromkeys(target_genes + present_tfs))
    if max_columns:
        cols = cols[:max_columns]
    print(f'  regulators present {len(present_tfs)}/{len(regs)}, targets present '
          f'{len(target_genes)}/{len(tgts)}, union columns {len(cols)}',
          flush=True)
    # Report what the gene axis cost us.  The matrix axis comes from the
    # reference dataset, so a target or regulator that exists in a dataset but
    # not on that axis is silently unusable -- say so rather than let it vanish.
    absent_t = sorted(tgts - set(genes))
    absent_r = sorted(regs - set(genes))
    if absent_t:
        print(f'  [!] {len(absent_t)} targets absent from the gene axis: '
              f'{absent_t[:8]}{" ..." if len(absent_t) > 8 else ""}', flush=True)
    if absent_r:
        print(f'  [!] {len(absent_r)} regulators absent from the gene axis: '
              f'{absent_r[:8]}{" ..." if len(absent_r) > 8 else ""}', flush=True)
    print(f'  NOTE: one regression per COLUMN -> {len(cols)} regressions '
          f'(regulator columns are fitted too)', flush=True)
    if not present_tfs:
        raise SystemExit(f'no regulators from {reg_file} are present in the '
                         f'matrix -- check the gene names')

    gene_order = {g: i for i, g in enumerate(genes)}
    X = np.load(expr_all, mmap_mode='r') if isinstance(expr_all, str) else expr_all
    sub = pd.DataFrame(np.asarray(X[:, [gene_order[g] for g in cols]],
                                  dtype=np.float32), columns=cols)
    del X
    n_cells = sub.shape[0]
    print(f'  GRN input: {sub.shape} '
          f'({sub.memory_usage(deep=True).sum()/1e9:.2f} GB)', flush=True)

    close = None
    if scheduler:
        client = Client(scheduler)
        print(f'  dask: external scheduler {scheduler}', flush=True)
    else:
        close = LocalCluster(n_workers=n_workers,
                             threads_per_worker=threads_per_worker,
                             processes=True, dashboard_address=None)
        client = Client(close)
        print(f'  dask: {n_workers} workers x {threads_per_worker} threads',
              flush=True)
    try:
        t0 = time.time()
        net = grnboost2(expression_data=sub, tf_names=present_tfs,
                        client_or_address=client, verbose=False, seed=seed)
        dt = time.time() - t0
    finally:
        client.close()
        if close is not None:
            close.close()

    print(f'  GRNBoost2: {len(net)} edges in {dt:.1f}s ({dt/60:.1f} min)',
          flush=True)
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
# evaluate + compare
# --------------------------------------------------------------------------- #
def evaluate(network_tsv, gt_files, outdir):
    """Score a network against ground-truth edge lists (missing files skipped)."""
    import pandas as pd

    os.makedirs(outdir, exist_ok=True)
    net = pd.read_csv(network_tsv, sep='\t')
    tfcol = 'TF' if 'TF' in net.columns else net.columns[0]
    tgcol = 'target' if 'target' in net.columns else net.columns[1]
    net = net.sort_values('importance', ascending=False).reset_index(drop=True)
    edge_set = set(zip(net[tfcol], net[tgcol]))
    results = {}
    with open(f'{outdir}/evaluation_summary.txt', 'w') as fh:
        fh.write(f'inferred edges: {len(net)}\n')
        fh.write('=' * 70 + '\n')
        for name, path in (gt_files or {}).items():
            if not path or not os.path.exists(path):
                print(f'  skipping {name}: not found ({path})', flush=True)
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


def _network(exp_dir, strat):
    path = f'{exp_dir}/{strat}/grn/network.tsv'
    if not os.path.exists(path):
        return None
    import pandas as pd
    net = pd.read_csv(path, sep='\t')
    tfcol = 'TF' if 'TF' in net.columns else net.columns[0]
    tgcol = 'target' if 'target' in net.columns else net.columns[1]
    net = net.sort_values('importance', ascending=False).reset_index(drop=True)
    edges = list(zip(net[tfcol], net[tgcol]))
    ev = f'{exp_dir}/{strat}/grn/evaluation/evaluation.json'
    return {'edges': edges,
            'eval': json.load(open(ev)).get('ground_truths', {})
            if os.path.exists(ev) else {}}


def compare_strategies(exp_dir, references):
    """Compare the networks from every reference strategy of one experiment.

    Answers the question the experiment exists to ask: does the choice of RNA
    reference change the inferred network?  Writes comparison.md / .json.
    """
    got = {}
    for strat in references:
        n = _network(exp_dir, strat)
        if n is None:
            print(f'  [compare] skipping {strat}: no network.tsv', flush=True)
            continue
        got[strat] = n
    if not got:
        return {}

    names = list(got)
    top = {s: set(got[s]['edges'][:100]) for s in names}
    lines = ['# Reference-strategy comparison', '',
             f'Strategies: {", ".join(names)}', '',
             '## Network size', '', '| strategy | edges |', '|---|---|']
    for s in names:
        lines.append(f'| {s} | {len(got[s]["edges"]):,} |')

    lines += ['', '## Overlap between strategies (top-100 edges)', '',
              '| A | B | shared | Jaccard |', '|---|---|---|---|']
    overlap = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            inter = len(top[a] & top[b])
            union = len(top[a] | top[b])
            jac = inter / union if union else 0.0
            lines.append(f'| {a} | {b} | {inter} | {jac:.3f} |')
            overlap[f'{a}|{b}'] = inter

    gt_names = sorted({g for s in names for g in got[s]['eval']})
    for g in gt_names:
        lines += ['', f'## Recovery vs {g}', '',
                  '| strategy | recovered | % |', '|---|---|---|']
        for s in names:
            r = got[s]['eval'].get(g)
            if r:
                lines.append(f'| {s} | {r["recovered"]}/{r["dedup"]} | {r["pct"]}% |')
        best = max((got[s]['eval'].get(g, {}).get('pct', -1) for s in names),
                   default=None)
        if best is not None and best >= 0:
            winners = [s for s in names
                       if got[s]['eval'].get(g, {}).get('pct', -1) == best]
            lines += ['', f'-> best: {", ".join(winners)} at {best}%']

    with open(f'{exp_dir}/comparison.md', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    js = {'strategies': {s: {'edges': len(got[s]['edges']),
                             'recovered': {g: got[s]['eval'].get(g, {})
                                           .get('recovered') for g in gt_names}}
                         for s in names},
          'overlap_top100': overlap}
    with open(f'{exp_dir}/comparison.json', 'w') as fh:
        json.dump(js, fh, indent=2)
    print(f'  [compare] wrote {exp_dir}/comparison.md', flush=True)
    return js


# --------------------------------------------------------------------------- #
# one strategy: reference -> impute -> grn   (isolated per strategy)
# --------------------------------------------------------------------------- #
def _run_strategy(cfg, datasets, rna, rna_expr, rna_genes, H, order, sizes,
                  queries, integ_dir, exp_dir, sdir, grn_dir,
                  all_path, strat, ref_spec, steps, args, root):
    """Build one reference strategy's all-cells matrix, then its network.

    Extracted from run_experiment so the caller can isolate failures: raising a
    SystemExit here only loses THIS strategy, not the whole experiment.
    """
    genes = None
    if 'impute' in steps:
        print('--- reference + reverse-imputeKNN')
        ref_expr, ref_H, ref_bc, genes = build_reference(
            ref_spec, rna, integ_dir)
        q_bc = {q: datasets[q].load_barcodes() for q in queries}
        imputed, q_bc_flat = impute(H, order, sizes, ref_expr, ref_H,
                                    queries, sdir, q_bc)
        np.save(f'{sdir}/genes.npy', np.array(genes))
        with open(f'{sdir}/genes.txt', 'w') as fh:
            fh.write('\n'.join(genes) + '\n')

        # all-cells rows: every real RNA cell in H order, then the queries.
        # The reference choice only decides what is propagated to ATAC.
        #
        # Every RNA dataset is put on the REFERENCE's gene axis before stacking.
        # A dataset from a different 10x reference has a different gene list,
        # and stacking it unaligned would silently label its rows with the
        # wrong genes (row j would hold a different gene in each block).
        rna_order = [o for o in order if o in rna]
        blocks = []
        for n in rna_order:
            src_genes = rna_genes[n]
            e = rna_expr[n]
            if src_genes == genes:
                blocks.append(e.T)
            else:
                n_absent = len(set(genes) - set(src_genes))
                if n_absent:
                    print(f'    [ref] {n}: {n_absent} genes of the reference '
                          f'axis absent -> zero-filled', flush=True)
                blocks.append(reindex_expression(e, src_genes, genes).T)
        real = np.vstack(blocks).astype(np.float32)
        all_expr = np.vstack([real, imputed.T]).astype(np.float32)
        all_bc = [b for n in rna_order for b in datasets[n].load_barcodes()]
        all_bc += q_bc_flat
        np.save(all_path, all_expr)
        with open(f'{sdir}/all_cells_barcodes.txt', 'w') as fh:
            fh.write('\n'.join(all_bc) + '\n')
        with open(f'{sdir}/n_cells.txt', 'w') as fh:
            fh.write(str(all_expr.shape[0]))
        print(f'    all-cells matrix {all_expr.shape} (rows: '
              f'{len(all_bc) - len(q_bc_flat)} real RNA + '
              f'{len(q_bc_flat)} imputed ATAC)')
    elif {'grn', 'evaluate'} & set(steps):
        if not os.path.exists(all_path):
            raise SystemExit(f'no {all_path}; run --only impute first')
        genes = read_lines(f'{sdir}/genes.txt')
        print(f'--- reuse all-cells matrix '
              f'{np.load(all_path, mmap_mode="r").shape}')

    if 'grn' in steps:
        if genes is None:
            raise SystemExit(f'no genes list for {strat}; '
                             f'run --only impute first')
        print(f'--- Arboreto GRNBoost2 ({args.workers} workers)')
        g = cfg.get('grn') or {}
        run_grn(all_path, genes,
                abspath(root, g.get('regulators') or 'data/tf_only.txt'),
                abspath(root, g.get('targets') or 'data/trrust_tf.txt'),
                grn_dir, n_workers=args.workers,
                threads_per_worker=args.threads_per_worker,
                seed=args.seed, scheduler=args.scheduler,
                max_columns=args.max_columns)

    net_path = f'{grn_dir}/network.tsv'
    if 'evaluate' in steps:
        if os.path.exists(net_path):
            print('--- evaluation')
            gt = {k: abspath(root, v)
                  for k, v in (cfg.get('ground_truth') or {}).items()}
            evaluate(net_path, gt, f'{grn_dir}/evaluation')
        else:
            print(f'    [!] no {net_path} -- evaluation SKIPPED for '
                  f'strategy {strat}', flush=True)


def run_experiment(cfg, datasets, exp_name, steps, args):
    spec = cfg['experiments'][exp_name]
    names = spec['datasets']
    anchor = spec['anchor']
    queries = spec.get('queries') or []
    references = spec.get('references') or {}
    root = cfg['_root']
    exp_dir = f'{root}/results/{exp_name}'
    integ_dir = f'{exp_dir}/integration'

    ds = [datasets[n] for n in names]
    rna = {d.name: d for d in ds if d.modality == 'rna'}
    if anchor not in rna:
        raise SystemExit(f'{exp_name}: anchor {anchor!r} must be an RNA dataset '
                         f'(got {rna.keys()})')
    os.makedirs(exp_dir, exist_ok=True)

    print(f'\n{"="*72}\nEXPERIMENT {exp_name}\n{"="*72}')
    print(f'  datasets   : {names}')
    print(f'  anchor     : {anchor}')
    print(f'  queries    : {queries}')
    print(f'  strategies : {list(references)}')

    # ---------------------------------------------------------- integrate
    H_path = f'{integ_dir}/joint_embedding_H.npy'
    if 'integrate' in steps:
        print(f'\n--- integrate ({len(ds)} datasets, anchor {anchor})')
        H, order, sizes = integrate(ds, anchor, integ_dir,
                                    cfg.get('integration') or {})
    else:
        if not os.path.exists(H_path):
            raise SystemExit(f'no {H_path}; run --only integrate first')
        info = json.load(open(f'{integ_dir}/integration_info.json'))
        order = info['order']
        sizes = {k: int(v) for k, v in info['sizes'].items()}
        H = np.load(H_path)
        print(f'\n--- reuse integration H {H.shape}')

    if not ({'impute', 'grn', 'evaluate'} & set(steps)):
        print(f'\nintegration written to {integ_dir}')
        return

    # real RNA expression, loaded once per experiment (these rows are the top
    # of the all-cells matrix, whichever reference strategy is used)
    rna_expr = {}
    rna_genes = {}
    for n, d in rna.items():
        genes_n = d.load_features()
        rna_genes[n] = genes_n
        rna_expr[n] = d.load_expression(len(genes_n))

    failed = {}
    missing_networks = []
    for strat, ref_spec in references.items():
        sdir = f'{exp_dir}/{strat}'
        grn_dir = f'{sdir}/grn'
        os.makedirs(sdir, exist_ok=True)
        all_path = f'{sdir}/all_cells_gene_expression.npy'
        genes = None
        print(f'\n===== strategy {strat}: {ref_spec}')

        # Each strategy is independent: the reference choice is the whole point
        # of the experiment, so one reference that cannot be built (e.g. a
        # 'combine' whose datasets do not share a feature list) must not take
        # the other strategies down with it.  build_reference() reports its
        # problems as SystemExit, which is NOT an Exception subclass, so both
        # are caught here.
        try:
            _run_strategy(cfg, datasets, rna, rna_expr, rna_genes, H, order,
                          sizes, queries, integ_dir, exp_dir, sdir, grn_dir,
                          all_path, strat, ref_spec, steps, args, root)
        except (Exception, SystemExit) as exc:              # noqa: BLE001
            print(f'    [!] strategy {strat} FAILED: {exc}', flush=True)
            failed[strat] = str(exc)
            continue

        if 'evaluate' in steps and not os.path.exists(f'{grn_dir}/network.tsv'):
            missing_networks.append(strat)

    if failed:
        print(f'\n{exp_name}: {len(failed)} strategy(ies) produced no results:')
        for s, why in failed.items():
            print(f'  - {s}: {why}')
    if missing_networks:
        raise SystemExit(
            f'{exp_name}: no network.tsv for {missing_networks}; the run is '
            f'incomplete (impute or grn did not produce a network)')

    if len(references) > 1:
        compare_strategies(exp_dir, [s for s in references if s not in failed])


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #
def check(cfg, datasets, only=None):
    exps = cfg.get('experiments') or {}
    todo = [only] if only else list(exps)
    ok = True
    print(f'root: {cfg["_root"]}\n')
    print(f'datasets declared: {len(datasets)}')
    for n, d in datasets.items():
        miss = d.missing()
        align = d.alignment_problems() if not miss else []
        if miss:
            state = 'MISSING: ' + ', '.join(miss)
        elif align:
            state = 'INCONSISTENT: ' + ', '.join(align)
        else:
            state = 'OK'
        n_cells = d.n_cells() if d.has_pca() or (d.barcodes and
                                                 os.path.exists(d.barcodes)) else None
        print(f'  {n:16s} {d.modality:5s} {state}'
              f'{"  " + str(n_cells) + " cells" if n_cells else ""}')
        if miss or align:
            ok = False

    print('\nexperiments:')
    for e in todo:
        spec = exps.get(e)
        if not spec:
            print(f'  {e}: NOT IN CONFIG')
            ok = False
            continue
        names = spec.get('datasets') or []
        unknown = [n for n in names if n not in datasets]
        if unknown:
            print(f'  {e}: unknown datasets {unknown}')
            ok = False
            continue
        flags = []
        anchor = spec.get('anchor')
        if anchor not in rna_of(spec, datasets):
            flags.append(f'anchor {anchor!r} is not an RNA dataset')
        for q in spec.get('queries') or []:
            if datasets[q].modality != 'atac':
                flags.append(f'query {q!r} is not an atac dataset')
        n_rna = len(rna_of(spec, datasets))
        refs = spec.get('references') or {}
        if not refs:
            flags.append('no references declared')
        for label, r in refs.items():
            if 'combine' in r and n_rna < 2:
                flags.append(f'strategy {label}: "combine" needs >=2 RNA datasets')
            if 'dataset' in r and r['dataset'] not in rna_of(spec, datasets):
                flags.append(f'strategy {label}: {r["dataset"]!r} is not an RNA '
                             f'dataset here')
        bad = [n for n in names if datasets[n].missing()]
        if bad:
            flags.append(f'missing files: {bad}')
        print(f'  {e}: {" + ".join(names)}  anchor={anchor}')
        print(f'      strategies: {", ".join(refs)}'
              f'  ({n_rna + (1 if n_rna > 1 else 0)} expected)'
              f'   -> {"OK" if not flags else " | ".join(flags)}')
        if flags:
            ok = False

    g = cfg.get('grn') or {}
    for label, key, default in (('regulators', 'regulators', 'data/tf_only.txt'),
                                ('targets', 'targets', 'data/trrust_tf.txt')):
        p = abspath(cfg['_root'], g.get(key) or default)
        exists = os.path.exists(p)
        print(f'  {label:11s} {p}  {"OK" if exists else "MISSING"}')
        if not exists:
            ok = False
    for name, p in (cfg.get('ground_truth') or {}).items():
        full = abspath(cfg['_root'], p)
        exists = os.path.exists(full)
        print(f'  truth {name:9s} {full}  {"OK" if exists else "MISSING (skipped)"}')
    return ok


def rna_of(spec, datasets):
    return [n for n in (spec.get('datasets') or [])
            if n in datasets and datasets[n].modality == 'rna']


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default='config.yml')
    ap.add_argument('-e', '--experiment', action='append', default=[],
                    help='experiment name; repeatable (default: all)')
    ap.add_argument('--list', action='store_true', help='list config and exit')
    ap.add_argument('--check', action='store_true', help='validate paths and exit')
    ap.add_argument('--only', default=','.join(STEPS),
                    help=f'comma-separated subset of {",".join(STEPS)}')
    ap.add_argument('--workers', type=int,
                    default=int(os.environ.get('GRN_WORKERS', '0')) or None,
                    help='dask workers for GRN (default: SLURM_CPUS_PER_TASK, '
                         'else min(cores, 8))')
    ap.add_argument('--threads-per-worker', type=int,
                    default=int(os.environ.get('GRN_THREADS_PER_WORKER', '1')))
    ap.add_argument('--scheduler', default=os.environ.get('GRN_SCHEDULER'),
                    help='attach to an existing dask scheduler')
    ap.add_argument('--seed', type=int, default=666)
    ap.add_argument('--max-columns', type=int,
                    default=int(os.environ.get('GRN_MAX_COLUMNS', '0')) or None,
                    help='smoke test: keep only N columns (default: all)')
    args = ap.parse_args()

    config_path = os.path.abspath(args.config)
    if not os.path.exists(config_path):
        raise SystemExit(f'no config at {config_path}')
    cfg = yaml.safe_load(open(config_path))
    root = os.path.abspath(os.environ.get('PIPELINE_ROOT')
                           or os.path.dirname(config_path))
    cfg['_root'] = root
    os.environ['PIPELINE_ROOT'] = root

    datasets = {n: Dataset(n, s or {}, root)
                for n, s in (cfg.get('datasets') or {}).items()}
    exps = cfg.get('experiments') or {}

    if args.list:
        print(f'config : {config_path}\nroot   : {root}\n\nexperiments:')
        for name, spec in exps.items():
            print(f'  {name:14s} {" + ".join(spec.get("datasets") or [])}  '
                  f'(anchor {spec.get("anchor")}, strategies: '
                  f'{", ".join(spec.get("references") or {})})')
        return 0

    if args.check:
        return 0 if check(cfg, datasets, args.experiment[0] if args.experiment
                          else None) else 1

    steps = [s.strip() for s in args.only.split(',') if s.strip()]
    bad = [s for s in steps if s not in STEPS]
    if bad:
        raise SystemExit(f'unknown step(s) {bad}; choose from {STEPS}')

    if args.workers is None:
        import multiprocessing
        args.workers = int(os.environ.get('SLURM_CPUS_PER_TASK')
                           or min(multiprocessing.cpu_count(), 8))

    targets = args.experiment or list(exps)
    unknown = [t for t in targets if t not in exps]
    if unknown:
        raise SystemExit(f'unknown experiment(s) {unknown}; '
                         f'available: {list(exps)}')

    print(f'config      : {config_path}')
    print(f'root        : {root}')
    print(f'experiments : {targets}')
    print(f'steps       : {steps}')
    print(f'workers     : {args.workers} '
          f'(threads/worker {args.threads_per_worker})')

    t0 = time.time()
    for e in targets:
        run_experiment(cfg, datasets, e, steps, args)
    print(f'\n{"="*72}\nDONE in {time.time()-t0:.1f}s\n'
          f'results in {root}/results/\n{"="*72}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
