#!/usr/bin/env python
"""Reverse imputeKNN: propagate RNA gene expression onto ATAC cells.

Identical to reverse_impute.py but consumes the NEW integration output that saved
H directly (results/scsaga_output), writing to results/downstream.
No separate extract_embedding.py step was used to produce H.
"""
import numpy as np, scipy.io, scipy.sparse as sp

RES  = '/Users/sbhattaram/pbmc3k_analysis/results/scsaga_output'
IN   = '/Users/sbhattaram/pbmc3k_analysis/scsaga_input'
DOWN = '/Users/sbhattaram/pbmc3k_analysis/results/downstream'
import os; os.makedirs(DOWN, exist_ok=True)

H = np.load(f'{RES}/joint_embedding_H.npy')           # (5422, 30)  rna then atac
all_bc = [l.strip() for l in open(f'{RES}/joint_embedding_barcodes.txt')]
n_rna = 2711
H_rna, H_atac = H[:n_rna], H[n_rna:]
bc_rna, bc_atac = all_bc[:n_rna], all_bc[n_rna:]

rna_mtx = scipy.io.mmread(f'{IN}/rna_counts.mtx').tocsr()  # cells x genes
rna_genes = [l.strip() for l in open(f'{IN}/rna_features.txt')]
r = rna_mtx.toarray()
tot = r.sum(axis=1, keepdims=True); tot[tot==0]=1
r = np.log1p(r / tot * 1e4).T  # genes x cells (reference RNA cells)
print('RNA expr genes x cells', r.shape)

K = 20
from sklearn.neighbors import NearestNeighbors
nbrs = NearestNeighbors(n_neighbors=K, metric='euclidean').fit(H_rna)
dist, idx = nbrs.kneighbors(H_atac)
e = np.exp(-dist)
Wk = e / e.sum(axis=1, keepdims=True)
from scipy.sparse import coo_matrix
rows = np.repeat(np.arange(n_atac := H_atac.shape[0]), K)
cols = idx.ravel()
W = coo_matrix((Wk.ravel(), (rows, cols)), shape=(n_atac, H_rna.shape[0])).tocsr()
imputed = r @ W.T
print('imputed expression genes x ATAC', imputed.shape)

np.save(f'{DOWN}/imputed_expression_genes_x_atac.npy', imputed)
with open(f'{DOWN}/imputed_atac_barcodes.txt','w') as f:
    f.write('\n'.join(bc_atac)+'\n')
with open(f'{DOWN}/imputed_atac_genes.txt','w') as f:
    f.write('\n'.join(rna_genes)+'\n')

all_expr = np.vstack([r.T, imputed.T])
np.save(f'{DOWN}/all_cells_gene_expression.npy', all_expr)
with open(f'{DOWN}/all_cells_barcodes.txt','w') as f:
    f.write('\n'.join(all_bc)+'\n')
np.save(f'{DOWN}/genes.npy', np.array(rna_genes))

import csv
with open(f'{DOWN}/all_cells_gene_expression.tsv','w',newline='') as f:
    w=csv.writer(f, delimiter='\t')
    w.writerow(['barcode']+rna_genes)
    for i,b in enumerate(all_bc):
        w.writerow([b]+[f'{v:.4f}' for v in all_expr[i]])
print('Wrote all_cells_gene_expression.npy and .tsv (results/downstream)')
print('all_cells shape', all_expr.shape)
