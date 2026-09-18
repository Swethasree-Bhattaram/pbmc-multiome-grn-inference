#!/usr/bin/env python
"""Reverse-imputeKNN: propagate real RNA expression onto the ATAC cells using the
scSAGA joint embedding H.

For each ATAC (query) cell, take its k=20 nearest RNA (reference) cells in H,
weight them by exp(-distance) (row-normalised), and propagate:

    imputed(genes x query) = ref_expr(genes x ref) @ W.T(ref x query)

Identical recipe to the repo's scripts/unpaired-10k/build_imputation_unpaired.py
(K=20, exp(-dist) weights).  Cell counts are read from the barcode files instead
of being hardcoded, so the same script serves both experiments.

  unpaired : reference = 10k multiome RNA (11,898), query = atac10k_ext (8,161)
  paired   : reference = 10k multiome RNA (11,898), query = 10k multiome ATAC (11,898)

Reads  $PACE_ROOT/results/integration_<exp>/joint_embedding_H.npy
Writes $PACE_ROOT/results/<exp>_grn/{imputed_expression_genes_x_atac.npy,
                                     all_cells_gene_expression.npy,
                                     genes.txt, genes.npy,
                                     all_cells_barcodes.txt, n_cells.txt}

all_cells rows = [ reference RNA (real), query ATAC (imputed) ].

Usage: python 05_build_imputation.py --experiment unpaired
"""
import argparse
import os
import sys
import time

import numpy as np
import scipy.io
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

K = 20
GENE_CHUNK = 4000


def read_lines(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]


def log1p_cpm_genes_by_cells(path, n_features):
    """cells x features MTX -> genes x cells float32 log1p(CPM/1e4), chunked over genes."""
    X = scipy.io.mmread(path).tocsr()                 # cells x features
    n_cells = X.shape[0]
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    scale = (1e4 / tot).astype(np.float64)
    Xcsc = X.multiply(scale[:, None]).tocsc()
    out = np.empty((n_features, n_cells), dtype=np.float32)
    for lo in range(0, n_features, GENE_CHUNK):
        hi = min(lo + GENE_CHUNK, n_features)
        sub = Xcsc[:, lo:hi].toarray().astype(np.float32)
        np.log1p(sub, out=sub)
        out[lo:hi, :] = sub.T
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--experiment', choices=['paired', 'unpaired'], required=True)
    ap.add_argument('--k', type=int, default=K)
    args = ap.parse_args()
    exp = args.experiment

    proj = os.environ.get('PACE_ROOT', os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    data = os.environ.get('PACE_DATA', f'{proj}/data')
    res = f'{proj}/results/integration_{exp}'
    out = f'{proj}/results/{exp}_grn'
    os.makedirs(out, exist_ok=True)

    ref_dir = f'{data}/10k_rna'
    qry_dir = f'{data}/atac10k_ext' if exp == 'unpaired' else f'{data}/10k_atac'

    genes = read_lines(f'{ref_dir}/features.txt')
    bc_ref = read_lines(f'{ref_dir}/barcodes.txt')
    bc_qry = read_lines(f'{qry_dir}/barcodes.txt')
    n_ref, n_qry = len(bc_ref), len(bc_qry)
    print(f'{exp}: reference {n_ref} RNA cells, query {n_qry} ATAC cells, '
          f'{len(genes)} genes', flush=True)

    H_path = f'{res}/joint_embedding_H.npy'
    if not os.path.exists(H_path):
        raise SystemExit(f'missing {H_path}; run 04_run_integration.py first')
    H = np.load(H_path)
    print('joint embedding H:', H.shape, flush=True)
    if H.shape[0] != n_ref + n_qry:
        raise SystemExit(f'H has {H.shape[0]} rows but expected {n_ref}+{n_qry}='
                         f'{n_ref+n_qry}; is the anchor/order right?')
    H_ref, H_qry = H[:n_ref], H[n_ref:]

    t0 = time.time()
    rna = log1p_cpm_genes_by_cells(f'{ref_dir}/counts.mtx', len(genes))
    assert rna.shape == (len(genes), n_ref), rna.shape
    print(f'RNA reference expression (genes x cells): {rna.shape} in '
          f'{time.time()-t0:.1f}s', flush=True)

    nbrs = NearestNeighbors(n_neighbors=args.k, metric='euclidean').fit(H_ref)
    dist, idx = nbrs.kneighbors(H_qry)
    e = np.exp(-dist)
    Wk = e / e.sum(axis=1, keepdims=True)
    rows = np.repeat(np.arange(H_qry.shape[0]), args.k)
    W = sp.coo_matrix((Wk.ravel(), (rows, idx.ravel())),
                      shape=(H_qry.shape[0], H_ref.shape[0])).tocsr()
    print(f'mean nn distance: {float(dist.mean()):.4f}  '
          f'distinct reference cells used: {len(np.unique(idx))}/{n_ref}', flush=True)

    imputed = (rna @ W.T).astype(np.float32)          # genes x n_qry
    print('imputed expression (genes x ATAC):', imputed.shape, flush=True)
    del W
    np.save(f'{out}/imputed_expression_genes_x_atac.npy', imputed)
    with open(f'{out}/imputed_atac_barcodes.txt', 'w') as f:
        f.write('\n'.join(bc_qry) + '\n')
    with open(f'{out}/genes.txt', 'w') as f:
        f.write('\n'.join(genes) + '\n')
    np.save(f'{out}/genes.npy', np.array(genes))

    all_expr = np.vstack([rna.T, imputed.T]).astype(np.float32)
    all_bc = bc_ref + bc_qry
    print('all-cells matrix:', all_expr.shape, all_expr.dtype, flush=True)
    np.save(f'{out}/all_cells_gene_expression.npy', all_expr)
    with open(f'{out}/all_cells_barcodes.txt', 'w') as f:
        f.write('\n'.join(all_bc) + '\n')
    with open(f'{out}/n_cells.txt', 'w') as f:
        f.write(str(all_expr.shape[0]))
    print('DONE ->', out, flush=True)


if __name__ == '__main__':
    sys.exit(main())
