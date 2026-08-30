#!/usr/bin/env python
"""Split the 10x PBMC 3k granulocyte-sorted multiome h5 into per-modality files
matching scSAGA's config format:
  counts.mtx, barcodes.txt, features.txt  (the 3 files per modality)
  plus a PCA file (scSAGA config requires it; code reads np.loadtxt pca).
RNA = 'Gene Expression' features, ATAC = 'Peaks' features.
"""
import os, h5py, numpy as np, scipy.sparse as sp, scipy.io
from sklearn.decomposition import PCA

RAW = '/Users/sbhattaram/pbmc3k_analysis/raw/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5'
OUT = '/Users/sbhattaram/pbmc3k_analysis/scsaga_input'
os.makedirs(OUT, exist_ok=True)

f = h5py.File(RAW, 'r')
grp = f['matrix']
data = grp['data'][:]
indices = grp['indices'][:]
indptr = grp['indptr'][:]
shape = tuple(grp['shape'][:])
barcodes = grp['barcodes'][:]
feat_type = grp['features']['feature_type'][:]
feat_names = grp['features']['name'][:]
feat_ids = grp['features']['id'][:]

M = sp.csr_matrix((data, indices, indptr), shape=(shape[1], shape[0]))  # cells x features
M = M.transpose().tocsr()  # features x cells

rna_mask = feat_type == b'Gene Expression'
atac_mask = feat_type == b'Peaks'

def emit(tag, mask, feat_names, feat_ids, M, barcodes):
    X = M[mask, :].tocsr()          # features x cells
    # features file: use name (gene symbol / peak interval). scSAGA features col used
    # for reference only; we write id for uniqueness and name for human readability.
    if tag == 'rna':
        feat_col = feat_names[mask]
    else:
        feat_col = feat_names[mask]  # peak names like chr:start-end are in 'name'
    with open(os.path.join(OUT, f'{tag}_counts.mtx'), 'wb') as fo:
        scipy.io.mmwrite(fo, X.transpose().tocsr())  # cells x features (mtx standard)
    with open(os.path.join(OUT, f'{tag}_barcodes.txt'), 'w') as fo:
        fo.write('\n'.join(b.decode() for b in barcodes) + '\n')
    with open(os.path.join(OUT, f'{tag}_features.txt'), 'w') as fo:
        fo.write('\n'.join(x.decode() for x in feat_col) + '\n')
    # PCA (50 comps) on the cells-as-rows log1p-normalized counts
    Xt = X.transpose().tocsr().toarray()   # cells x features
    Xt = np.log1p(Xt / Xt.sum(axis=1, keepdims=True) * 1e4)
    # restrict to top variable features to be quick & robust
    var = Xt.var(axis=0)
    keep = np.argsort(var)[-2000:]
    pca = PCA(n_components=50, random_state=0)
    P = pca.fit_transform(Xt[:, keep])
    np.savetxt(os.path.join(OUT, f'{tag}_pca_50.txt'), P, fmt='%.6f')
    return X.shape[0], X.shape[1]

for tag, mask in [('rna', rna_mask), ('atac', atac_mask)]:
    nfeat, ncell = emit(tag, mask, feat_names, feat_ids, M, barcodes)
    print(f'{tag}: features={nfeat} cells={ncell}')

print('Done. Files in', OUT)
