#!/usr/bin/env python
"""Preprocess the 10x PBMC granulocyte-sorted 10k multiome h5 into scSAGA-format
per-modality files, subsampling the SAME 2711 cells across both RNA and ATAC.

We pick a deterministic subset of 2711 cells (same index set for both modalities,
since 10k multiome barcodes are shared/paired) so the '6k' RNA and '6k' ATAC
datasets are the same physical cells — matching how the original pbmc_3k dataset
is one multiome object with 2711 shared cells.

Outputs per modality (data/6k_rna, data/6k_atac):
  counts.mtx    (cells x features, MTX)
  barcodes.txt  (one per line)
  features.txt  (one per line)
  pca_50.txt    (50 PCs on log1p(CPM) of top-2000 variable features)
"""
import os, h5py, numpy as np, scipy.sparse as sp, scipy.io
from sklearn.decomposition import PCA

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)   # project root
RAW = os.path.join(PROJ, 'raw', 'pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5')
OUT = os.path.join(PROJ, 'data')
N_SUB = 2711          # same number of cells per modality as the original 3k
SEED = 0              # deterministic subsample

f = h5py.File(RAW, 'r')
grp = f['matrix']
data = grp['data'][:]
indices = grp['indices'][:]
indptr = grp['indptr'][:]
shape = tuple(grp['shape'][:])
barcodes = grp['barcodes'][:]
feat_type = grp['features']['feature_type'][:]
feat_names = grp['features']['name'][:]

M = sp.csr_matrix((data, indices, indptr), shape=(shape[1], shape[0]))  # cells x features
M = M.transpose().tocsr()  # features x cells

rna_mask = feat_type == b'Gene Expression'
atac_mask = feat_type == b'Peaks'
print('RNA features:', rna_mask.sum(), 'ATAC features:', atac_mask.sum(),
      'cells:', M.shape[1])

# Choose a fixed set of 2711 cells. Same indices for both modalities (paired).
rng = np.random.RandomState(SEED)
cell_idx = np.sort(rng.choice(M.shape[1], size=N_SUB, replace=False))
keep = np.zeros(M.shape[1], dtype=bool)
keep[cell_idx] = True
print('Subsampled cells:', int(keep.sum()))


def emit(tag, mask):
    X = M[mask, :][:, keep].tocsr()          # features x n_subsampled_cells
    feat_col = [x.decode() for x in feat_names[mask]]
    outdir = os.path.join(OUT, f'6k_{tag}')
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'counts.mtx'), 'wb') as fo:
        scipy.io.mmwrite(fo, X.transpose().tocsr())  # cells x features
    with open(os.path.join(outdir, 'barcodes.txt'), 'w') as fo:
        fo.write('\n'.join(barcodes[i].decode() for i in cell_idx) + '\n')
    with open(os.path.join(outdir, 'features.txt'), 'w') as fo:
        fo.write('\n'.join(feat_col) + '\n')
    # PCA on cells-as-rows, log1p(CPM), top-2000 variable features
    Xt = X.transpose().tocsr().toarray()     # cells x features
    Xt = np.log1p(Xt / Xt.sum(axis=1, keepdims=True) * 1e4)
    var = Xt.var(axis=0)
    keepf = np.argsort(var)[-2000:]
    pca = PCA(n_components=50, random_state=0)
    P = pca.fit_transform(Xt[:, keepf])
    np.savetxt(os.path.join(outdir, 'pca_50.txt'), P, fmt='%.6f')
    return X.shape


for tag, mask in [('rna', rna_mask), ('atac', atac_mask)]:
    d = emit(tag, mask)
    print(f'6k_{tag}: features={d[0]} cells={d[1]}')

print('Done.')
