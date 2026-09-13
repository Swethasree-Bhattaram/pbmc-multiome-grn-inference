#!/usr/bin/env python
"""Preprocess the external 10x PBMC 10k ATAC v1.1 (hg19, Cell Ranger ATAC 1.1.0)
"cells by peaks" filtered matrix into scSAGA-format files.

Dataset page:
https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-atac-v-1-1-chromium-controller-1-1-standard-2-0-0
Downloaded from:
https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_v1/atac_pbmc_10k_v1_filtered_peak_bc_matrix.h5

This is an UNPAIRED companion to the 10k multiome RNA: different cells, different
donor/protocol, and hg19 peaks (the multiome is GRCh38). Only the per-dataset PCA
is used by scSAGA for integration, so peak-coordinate build differences do not
enter the alignment directly.

Outputs (data/atac10k_ext):
  counts.mtx    (cells x features, MTX)
  barcodes.txt  (one per line)
  features.txt  (one per line)
  pca_50.txt    (50 PCs on log1p(CPM/1e4) of top-2000 variable features)

Same recipe as scripts/preprocess_10k_full.py so the two modalities are treated
identically. Memory-safe: exact log1p variance from the sparse structure, then
densify only the top-2000 columns.
"""
import os, h5py, numpy as np, scipy.sparse as sp, scipy.io
from sklearn.decomposition import PCA

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.environ.get('UNPAIRED_ROOT',
                      os.path.abspath(os.path.join(HERE, os.pardir)))
RAW = os.path.join(PROJ, 'raw', 'atac10k_v1.1_peaks.h5')
OUT = os.path.join(PROJ, 'data', 'atac10k_ext')

f = h5py.File(RAW, 'r')
grp = f['matrix']
data = grp['data'][:]
indices = grp['indices'][:]
indptr = grp['indptr'][:]
shape = tuple(int(x) for x in grp['shape'][:])   # (n_features, n_cells)
barcodes = grp['barcodes'][:]
feat_type = grp['features']['feature_type'][:]
feat_names = grp['features']['name'][:]

# 10x h5 is cell-major: indptr length == n_cells + 1
M = sp.csr_matrix((data, indices, indptr), shape=(shape[1], shape[0]))  # cells x features
M = M.transpose().tocsr()                                              # features x cells

print('shape (features, cells):', shape)
print('feature types present:', sorted({x.decode() for x in feat_type}))
print('cells:', M.shape[1], 'features:', M.shape[0], 'nnz:', M.nnz)

os.makedirs(OUT, exist_ok=True)
X = M.tocsr()                                     # features x cells
feat_col = [x.decode() for x in feat_names]

with open(os.path.join(OUT, 'counts.mtx'), 'wb') as fo:
    scipy.io.mmwrite(fo, X.transpose().tocsr())   # cells x features
with open(os.path.join(OUT, 'barcodes.txt'), 'w') as fo:
    fo.write('\n'.join(barcodes[i].decode() for i in range(M.shape[1])) + '\n')
with open(os.path.join(OUT, 'features.txt'), 'w') as fo:
    fo.write('\n'.join(feat_col) + '\n')

# --- PCA (identical recipe to preprocess_10k_full.py) ---
Xc = X.transpose().tocsr()                    # cells x features (sparse)
tot = np.asarray(Xc.sum(axis=1)).ravel()
tot[tot == 0] = 1
scale = (1e4 / tot).astype(np.float64)
Xs = Xc.multiply(scale[:, None]).tocsc()      # CPM/1e4
n = Xs.shape[0]
Xs = Xs.tocoo()
lx = np.log1p(Xs.data)
sum_lx = np.bincount(Xs.col, weights=lx, minlength=Xs.shape[1])
sum_lx2 = np.bincount(Xs.col, weights=lx * lx, minlength=Xs.shape[1])
mean = sum_lx / n
mean_sq = sum_lx2 / n
var = mean_sq - mean ** 2
keepf = np.argsort(var)[-2000:]
Xt = np.log1p(Xs.tocsr()[:, keepf].toarray())  # cells x 2000
print('top-2000 variable-feature matrix for PCA:', Xt.shape)

P = PCA(n_components=50, random_state=0).fit_transform(Xt)
np.savetxt(os.path.join(OUT, 'pca_50.txt'), P, fmt='%.6f')
print('wrote pca_50.txt:', P.shape)
print('explained variance (50 PCs):', float(P.var(axis=0).sum()))
print('ATAC-ext cells/features:', M.shape[1], M.shape[0])
print('Done ->', OUT)
