#!/usr/bin/env python
"""Split a 10x Genomics multiome / ATAC / RNA h5 into scSAGA-format per-modality files.

One splitter for every input in this repo, so the PACE pipeline does not have to
carry three near-identical one-off scripts.

Produces, per modality (cells x features MTX is the on-disk orientation):
  counts.mtx    cells x features
  barcodes.txt  one per line
  features.txt  gene symbols / peak intervals, one per line
  pca_50.txt    50 PCs on log1p(CPM/1e4) of the top-2000 variable features

The PCA recipe is numerically identical to both recipes used in the repo
(`scripts/common/extract_data.py` densifies everything then takes the top-2000
by variance; `preprocess_10k_full.py` computes the exact log1p variance from the
sparse structure then densifies only the top-2000 columns).  Verified on the 3k
RNA modality: identical feature selection and max |PC difference| = 4.3e-12.
This script uses the memory-safe sparse-variance form.

Modes
-----
multiome   (default)  a Cell Ranger ARC feature-bc matrix with both
                      'Gene Expression' and 'Peaks' features. Writes both
                      modalities into --out-subdir/{tag}.
atac                  an ATAC-only "cells by peaks" matrix (single feature type).
rna                   an RNA-only matrix.

Examples
--------
# 3k and 10k multiome -> data/3k_{rna,atac}, data/10k_{rna,atac}
python split_10x_h5.py --h5 raw/3k_multiome.h5 --out-subdir 3k
python split_10x_h5.py --h5 raw/10k_multiome.h5 --out-subdir 10k

# external 10x ATAC v1.1 -> data/atac10k_ext
python split_10x_h5.py --h5 raw/atac10k_v1.1.h5 --mode atac --out-subdir atac10k_ext

# 6k = deterministic 2,711-cell subsample of the 10k multiome (both modalities,
# the SAME cells, seed 0 -- reproduces scripts/multi-dataset-4x/preprocess_10k.py)
python split_10x_h5.py --h5 raw/10k_multiome.h5 --out-subdir 6k --subsample 2711 --seed 0

The '6k' files are only needed for the 3k+6k workflow; the unpaired-10k and
3k+10k workflows do not use them.
"""
import argparse
import os
import sys

import h5py
import numpy as np
import scipy.io
import scipy.sparse as sp
from sklearn.decomposition import PCA

GEX_LABELS = {b'Gene Expression', b'Gene Expression ', b'RNA'}
PEAK_LABELS = {b'Peaks', b'Peaks ', b'ATAC'}


def read_h5(path):
    """Return (features x cells CSR, barcodes, feature_names, feature_types)."""
    with h5py.File(path, 'r') as f:
        if 'matrix' in f:
            grp = f['matrix']
        else:
            raise SystemExit(f'{path}: no /matrix group (not a 10x h5)')
        data = grp['data'][:]
        indices = grp['indices'][:]
        indptr = grp['indptr'][:]
        shape = tuple(int(x) for x in grp['shape'][:])       # (n_features, n_cells)
        barcodes = [b.decode() for b in grp['barcodes'][:]]
        feats = grp['features']
        names = [x.decode() for x in feats['name'][:]]
        if 'feature_type' in feats:
            types = [x for x in feats['feature_type'][:]]
        else:
            types = None
    # 10x h5 is cell-major in the sense that indptr has n_cells+1 entries
    M = sp.csr_matrix((data, indices, indptr), shape=(shape[1], shape[0]))  # cells x features
    return M.transpose().tocsr(), barcodes, names, types


def pca_50(X):
    """50 PCs on log1p(CPM/1e4) of the top-2000 variable features.

    X is features x cells (sparse). Exact log1p variance from the sparse
    structure; densify only the top-2000 columns.
    """
    Xc = X.transpose().tocsr()                      # cells x features
    tot = np.asarray(Xc.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    scale = (1e4 / tot).astype(np.float64)
    Xs = Xc.multiply(scale[:, None]).tocsc()        # log1p(CPM/1e4), still sparse
    n = Xs.shape[0]
    Xs = Xs.tocoo()
    lx = np.log1p(Xs.data)
    sum_lx = np.bincount(Xs.col, weights=lx, minlength=Xs.shape[1])
    sum_lx2 = np.bincount(Xs.col, weights=lx * lx, minlength=Xs.shape[1])
    mean = sum_lx / n
    mean_sq = sum_lx2 / n
    var = mean_sq - mean ** 2
    n_keep = min(2000, Xs.shape[1])
    keepf = np.argsort(var)[-n_keep:]
    Xt = np.log1p(Xs.tocsr()[:, keepf].toarray())   # cells x <=2000
    return PCA(n_components=50, random_state=0).fit_transform(Xt)


def emit(outdir, X, barcodes, names, tag):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'counts.mtx'), 'wb') as fo:
        scipy.io.mmwrite(fo, X.transpose().tocsr())   # cells x features
    with open(os.path.join(outdir, 'barcodes.txt'), 'w') as fo:
        fo.write('\n'.join(barcodes) + '\n')
    with open(os.path.join(outdir, 'features.txt'), 'w') as fo:
        fo.write('\n'.join(names) + '\n')
    np.savetxt(os.path.join(outdir, 'pca_50.txt'), pca_50(X), fmt='%.6f')
    print(f'  {tag:>12s}: {X.shape[0]:>6d} features x {X.shape[1]:>6d} cells -> {outdir}',
          flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5', required=True, help='input 10x .h5')
    ap.add_argument('--data-dir', default=None,
                    help='output data dir (default <workspace>/data)')
    ap.add_argument('--out-subdir', required=True,
                    help='dataset name; files go to <data-dir>/<out-subdir>[_rna|_atac]')
    ap.add_argument('--mode', choices=['multiome', 'atac', 'rna'], default='multiome')
    ap.add_argument('--subsample', type=int, default=0,
                    help='if >0, deterministically subsample this many cells '
                         '(same cells in both modalities)')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = args.data_dir or os.environ.get('PACE_DATA') or \
        os.path.join(os.environ.get('PIPELINE_ROOT') or
                     os.path.abspath(os.path.join(here, os.pardir)), 'data')

    M, barcodes, names, types = read_h5(args.h5)
    print(f'{args.h5}: {M.shape[0]} features x {M.shape[1]} cells, nnz={M.nnz}', flush=True)
    if types is not None:
        print('feature types present:', sorted({t.decode() for t in types}), flush=True)

    keep = None
    if args.subsample:
        if args.subsample > M.shape[1]:
            raise SystemExit(f'--subsample {args.subsample} > {M.shape[1]} cells')
        rng = np.random.RandomState(args.seed)
        keep = np.sort(rng.choice(M.shape[1], size=args.subsample, replace=False))
        print(f'deterministic subsample ({args.seed}): {args.subsample} of {M.shape[1]} cells',
              flush=True)

    subsets = []          # (tag, row_mask into features)
    if args.mode == 'multiome':
        if types is None:
            raise SystemExit('--mode multiome but the h5 has no feature_type field')
        rna_mask = np.array([t in GEX_LABELS for t in types])
        atac_mask = np.array([t in PEAK_LABELS for t in types])
        if not rna_mask.any() or not atac_mask.any():
            raise SystemExit('--mode multiome: could not find both Gene Expression and Peaks')
        subsets = [('rna', rna_mask), ('atac', atac_mask)]
    else:
        subsets = [('rna' if args.mode == 'rna' else 'atac', np.ones(M.shape[0], dtype=bool))]

    for tag, mask in subsets:
        X = M[mask, :].tocsr()
        bc = barcodes
        nm = [names[i] for i in np.nonzero(mask)[0]]
        if keep is not None:
            X = X[:, keep].tocsr()
            bc = [barcodes[i] for i in keep]
        sub = args.out_subdir if len(subsets) == 1 and args.mode != 'multiome' \
            else f'{args.out_subdir}_{tag}'
        if args.mode == 'multiome':
            sub = f'{args.out_subdir}_{tag}'
        emit(os.path.join(data_dir, sub), X, bc, nm, sub)

    print('Done.', flush=True)


if __name__ == '__main__':
    sys.exit(main())
