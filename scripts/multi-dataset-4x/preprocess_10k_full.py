#!/usr/bin/env python
"""Preprocess the 10x PBMC granulocyte-sorted 10k multiome h5 into scSAGA-format
per-modality files using ALL cells (no subsampling).

The full 10k multiome has 11,898 cells. We split RNA/ATAC and keep every cell
(paired multiome, so RNA and ATAC share the same barcodes).

Outputs per modality (data/10k_rna, data/10k_atac):
  counts.mtx    (cells x features, MTX)
  barcodes.txt  (one per line)
  features.txt  (one per line)
  pca_50.txt    (50 PCs on log1p(CPM) of top-2000 variable features)
"""
import os, h5py, numpy as np, scipy.sparse as sp, scipy.io
from sklearn.decomposition import PCA

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)   # project root
RAW = os.path.join(PROJ, 'raw', '10k.h5')
OUT = os.path.join(PROJ, 'data')

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

def emit(tag, mask):
    X = M[mask, :].tocsr()          # features x n_cells
    feat_col = [x.decode() for x in feat_names[mask]]
    outdir = os.path.join(OUT, f'10k_{tag}')
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'counts.mtx'), 'wb') as fo:
        scipy.io.mmwrite(fo, X.transpose().tocsr())  # cells x features
    with open(os.path.join(outdir, 'barcodes.txt'), 'w') as fo:
        fo.write('\n'.join(barcodes[i].decode() for i in range(M.shape[1])) + '\n')
    with open(os.path.join(outdir, 'features.txt'), 'w') as fo:
        fo.write('\n'.join(feat_col) + '\n')
    # PCA on cells-as-rows, log1p(CPM), top-2000 variable features.
    # Memory-safe: compute exact log1p variance from the sparse structure, then
    # densify only the top-2000 columns.
    Xc = X.transpose().tocsr()          # cells x features (sparse)
    tot = np.asarray(Xc.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    scale = (1e4 / tot).astype(np.float64)
    Xs = Xc.multiply(scale[:, None]).tocsc()   # scaled counts (CPM/1e4)
    n = Xs.shape[0]
    # exact E[log1p(x)] and E[log1p(x)^2] per column from sparse data
    Xs = Xs.tocoo()
    lx = np.log1p(Xs.data)
    sum_lx = np.bincount(Xs.col, weights=lx, minlength=Xs.shape[1])
    sum_lx2 = np.bincount(Xs.col, weights=lx * lx, minlength=Xs.shape[1])
    mean = sum_lx / n
    mean_sq = sum_lx2 / n
    var = mean_sq - mean ** 2
    keepf = np.argsort(var)[-2000:]
    Xt = np.log1p(Xs.tocsr()[:, keepf].toarray())   # cells x 2000
    pca = PCA(n_components=50, random_state=0)
    P = pca.fit_transform(Xt)
    np.savetxt(os.path.join(outdir, 'pca_50.txt'), P, fmt='%.6f')
    return X.shape

for tag, mask in [('rna', rna_mask), ('atac', atac_mask)]:
    d = emit(tag, mask)
    print(f'10k_{tag}: features={d[0]} cells={d[1]}')

print('Done.')
