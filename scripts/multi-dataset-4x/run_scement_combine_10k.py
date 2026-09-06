#!/usr/bin/env python
"""Build the SCEMENT-integrated combined RNA reference for Experiment A (full-10k).

SCEMENT (ComBat-style sparse batch integration) merges the two real RNA datasets
(3k RNA, 2,711 cells + 10k RNA, 11,898 cells) into a single batch-corrected
reference: genes x 14,609 cells (3k RNA then 10k RNA).

The reference embedding in joint scSAGA space is the vstack of the two RNA
H-blocks (rna3k then rna10k), matched to the expression order.

Outputs (results/expA_10k):
  combined_ref_expression.npy   genes x 14609  (batch-corrected counts, float32)
  combined_ref_embedding.npy    14609 x 30     (vstack of H_rna3k, H_rna10k)
  combined_ref_barcodes.txt     14609 (3k then 10k)

Run with the .venv-scement python.
"""
import os, sys, importlib.util
import numpy as np, scipy.io, scipy.sparse as sp, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir)))
DATA = f'{PROJ}/data'
RES = f'{PROJ}/results/integration_10k'
OUT = f'{PROJ}/results/expA_10k'
os.makedirs(OUT, exist_ok=True)

# --- Load patched pure-Python SCEMENT (self-contained copy) ---
_SCEMENT_PY = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'scement_py', 'scement_sparse.py')
spec = importlib.util.spec_from_file_location('scement_sparse', _SCEMENT_PY)
scm = importlib.util.module_from_spec(spec); spec.loader.exec_module(scm)

# --- Load raw (un-normalized) RNA counts, cells x genes ---
def load(d):
    X = scipy.io.mmread(f'{DATA}/{d}/counts.mtx').tocsr()
    bc = [l.strip() for l in open(f'{DATA}/{d}/barcodes.txt')]
    feat = [l.strip() for l in open(f'{DATA}/{d}/features.txt')]
    return X, bc, feat

r3, bc3, feat = load('3k_rna')
r10, bc10, feat2 = load('10k_rna')
assert feat == feat2 and r3.shape[1] == r10.shape[1]
print('3k RNA:', r3.shape, '10k RNA:', r10.shape)

import anndata as ad
X = sp.vstack([r3, r10]).tocsr()           # 14609 cells x genes (raw counts)
obs = pd.DataFrame({'batch': ['3k'] * r3.shape[0] + ['10k'] * r10.shape[0]})
obs['batch'] = obs['batch'].astype('category')
adata = ad.AnnData(X=X, obs=obs)
print('SCEMENT input AnnData:', adata.shape, 'batches:', obs['batch'].value_counts().to_dict())

# --- Run SCEMENT batch integration (out-of-place) ---
corrected = scm.sct_sparse(adata, key='batch', inplace=False)   # cells x genes ndarray
print('SCEMENT output:', corrected.shape, corrected.dtype)
corrected = np.clip(corrected, 0, None).astype(np.float32)

# --- Combined reference embedding: vstack of the two RNA H-blocks ---
H = np.load(f'{RES}/joint_embedding_H.npy')
n3 = 2711
n10 = 11898
H_rna3k = H[0:n3]
H_rna10k = H[n3 + n3: n3 + n3 + n10]   # blocks: rna3k, atac3k, rna10k, atac10k
ref_H = np.vstack([H_rna3k, H_rna10k])
assert ref_H.shape[0] == corrected.shape[0] == 14609, (ref_H.shape, corrected.shape)

np.save(f'{OUT}/combined_ref_expression.npy', corrected.T)   # genes x 14609
np.save(f'{OUT}/combined_ref_embedding.npy', ref_H)
with open(f'{OUT}/combined_ref_barcodes.txt', 'w') as f:
    f.write('\n'.join(bc3 + bc10) + '\n')
with open(f'{OUT}/combined_ref_genes.txt', 'w') as f:
    f.write('\n'.join(feat) + '\n')

print('combined_ref_expression (genes x cells):', corrected.T.shape)
print('combined_ref_embedding:', ref_H.shape)
print('DONE')
