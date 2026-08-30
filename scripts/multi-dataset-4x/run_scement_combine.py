#!/usr/bin/env python
"""Build the SCEMENT-integrated combined RNA reference for Experiment A.

SCEMENT (ComBat-style sparse batch integration) merges the two real RNA datasets
into a single batch-corrected reference: genes x 5422 cells (3k RNA then 6k RNA).

We use SCEMENT's pure-Python `sct_sparse` (patched for numpy2/scipy2, self-contained
copy at scripts/scement_py/scement_sparse.py). The reference embedding in joint
scSAGA space is the
vstack of the two RNA H-blocks (rna3k then rna6k), matched to the expression order.

Outputs (results/expA):
  combined_ref_expression.npy   genes x 5422  (batch-corrected counts)
  combined_ref_embedding.npy    5422 x 30     (vstack of H_rna3k, H_rna6k)
  combined_ref_barcodes.txt     5422 (3k then 6k)

Run with the .venv-scement python.
"""
import os, sys, importlib.util
import numpy as np, scipy.io, scipy.sparse as sp, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
DATA = f'{PROJ}/data'
RES = f'{PROJ}/results/integration'
OUT = f'{PROJ}/results/expA'
os.makedirs(OUT, exist_ok=True)

# --- Load patched pure-Python SCEMENT (self-contained copy next to this script) ---
_SCEMENT_PY = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'scement_py', 'scement_sparse.py')
spec = importlib.util.spec_from_file_location(
    'scement_sparse', _SCEMENT_PY)
scm = importlib.util.module_from_spec(spec); spec.loader.exec_module(scm)

# --- Load raw (un-normalized) RNA counts, cells x genes ---
def load(d):
    X = scipy.io.mmread(f'{DATA}/{d}/counts.mtx').tocsr()
    bc = [l.strip() for l in open(f'{DATA}/{d}/barcodes.txt')]
    feat = [l.strip() for l in open(f'{DATA}/{d}/features.txt')]
    return X, bc, feat

r3, bc3, feat = load('3k_rna')
r6, bc6, feat2 = load('6k_rna')
assert feat == feat2 and r3.shape[1] == r6.shape[1]

import anndata as ad
X = sp.vstack([r3, r6]).tocsr()           # 5422 cells x genes (raw counts)
obs = pd.DataFrame({'batch': ['3k'] * r3.shape[0] + ['6k'] * r6.shape[0]})
obs['batch'] = obs['batch'].astype('category')
adata = ad.AnnData(X=X, obs=obs)
print('SCEMENT input AnnData:', adata.shape, 'batches:', obs['batch'].value_counts().to_dict())

# --- Run SCEMENT batch integration (out-of-place) ---
corrected = scm.sct_sparse(adata, key='batch', inplace=False)   # cells x genes ndarray
print('SCEMENT output:', corrected.shape, corrected.dtype)
# clip negatives to 0 (ComBat can produce small negatives)
corrected = np.clip(corrected, 0, None).astype(np.float32)

# --- Combined reference embedding: vstack of the two RNA H-blocks ---
H = np.load(f'{RES}/joint_embedding_H.npy')
n = 2711
H_rna3k = H[0:n]
H_rna6k = H[2*n:3*n]
ref_H = np.vstack([H_rna3k, H_rna6k])
assert ref_H.shape[0] == corrected.shape[0] == 5422

np.save(f'{OUT}/combined_ref_expression.npy', corrected.T)   # genes x 5422
np.save(f'{OUT}/combined_ref_embedding.npy', ref_H)
with open(f'{OUT}/combined_ref_barcodes.txt', 'w') as f:
    f.write('\n'.join(bc3 + bc6) + '\n')
with open(f'{OUT}/combined_ref_genes.txt', 'w') as f:
    f.write('\n'.join(feat) + '\n')

print('combined_ref_expression (genes x cells):', corrected.T.shape)
print('combined_ref_embedding:', ref_H.shape)
print('DONE')
