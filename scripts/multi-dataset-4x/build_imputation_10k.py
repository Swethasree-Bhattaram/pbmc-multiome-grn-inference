#!/usr/bin/env python
"""Reverse-imputeKNN for the full-10k 4-dataset scenario.

Given the 4-dataset joint embedding H (blocks: rna3k, atac3k, rna10k, atac10k)
and a chosen RNA reference, propagate RNA expression onto both ATAC populations
(atac3k + atac10k).

Reference expression can be:
  - a single RNA dataset (3k only, or 10k only)  -> Exp B1 / B2
  - a SCEMENT-integrated combination of 3k+10k RNA -> Exp A  (produced by
    scripts/run_scement_combine_10k.py; reference has 14,609 cells)

Query = atac3k cells + atac10k cells (both imputed).

Output (to results/<exp>_10k):
  imputed_expression_genes_x_atac.npy  (genes x 14609; ATAC3k rows then ATAC10k)
  imputed_atac_barcodes.txt            (14609 barcodes: ATAC3k then ATAC10k)
  imputed_atac_genes.txt
  all_cells_gene_expression.npy        (29218 x genes, float32)
  all_cells_barcodes.txt
  genes.npy
  ref_block_barcodes.txt

All-cells matrix rows: [RNA3k real, RNA10k real, ATAC3k imputed, ATAC10k imputed].

Usage: python build_imputation_10k.py <exp> <ref_mode>
  exp       : expA | expB1 | expB2
  ref_mode  : combined | rna3k | rna10k
"""
import os, sys, numpy as np, scipy.io, scipy.sparse as sp

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir)))
RES = f'{PROJ}/results/integration_10k'
DATA = f'{PROJ}/data'
K = 20
N3 = 2711
N10 = 11898

def load_mtx_counts(d):
    X = scipy.io.mmread(f'{DATA}/{d}/counts.mtx').tocsr()  # cells x features
    bc = [l.strip() for l in open(f'{DATA}/{d}/barcodes.txt')]
    feat = [l.strip() for l in open(f'{DATA}/{d}/features.txt')]
    return X, bc, feat

def log1p_cpm(X):
    # X is cells x features (sparse). Compute per-cell totals, then log1p(CPM).
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    Xc = X.tocsr().multiply((1e4 / tot)[:, None]).tocsc()
    return np.log1p(Xc.toarray()).T.astype(np.float32)  # genes x cells

def main():
    exp, ref_mode = sys.argv[1], sys.argv[2]
    outdir = f'{PROJ}/results/{exp}_10k'
    os.makedirs(outdir, exist_ok=True)

    # --- RNA real expression (genes x cells) ---
    rna3k, bc_rna3k, genes = load_mtx_counts('3k_rna')
    rna10k, bc_rna10k, _ = load_mtx_counts('10k_rna')
    assert rna3k.shape[1] == rna10k.shape[1] == len(genes)
    r3 = log1p_cpm(rna3k).astype(np.float32)   # genes x 2711
    r10 = log1p_cpm(rna10k).astype(np.float32) # genes x 11898
    print('RNA expr genes x cells: 3k', r3.shape, '10k', r10.shape)

    # --- Load joint embedding blocks ---
    H = np.load(f'{RES}/joint_embedding_H.npy')
    Hrna3k, Hatac3k, Hrna10k, Hatac10k = H[0:N3], H[N3:2*N3], H[2*N3:2*N3+N10], H[2*N3+N10:]

    # --- Build reference expression + embedding ---
    if ref_mode == 'rna3k':
        ref_expr, ref_H, ref_bc = r3, Hrna3k, bc_rna3k
    elif ref_mode == 'rna10k':
        ref_expr, ref_H, ref_bc = r10, Hrna10k, bc_rna10k
    elif ref_mode == 'combined':
        ref_expr = np.load(f'{PROJ}/results/expA_10k/combined_ref_expression.npy').astype(np.float32)
        ref_H = np.load(f'{PROJ}/results/expA_10k/combined_ref_embedding.npy')
        ref_bc = [l.strip() for l in open(f'{PROJ}/results/expA_10k/combined_ref_barcodes.txt')]
    else:
        raise SystemExit(f'unknown ref_mode {ref_mode}')
    print(f'[ref] {ref_mode}: expression {ref_expr.shape}, embedding {ref_H.shape}, cells {len(ref_bc)}')
    assert ref_expr.shape[1] == ref_H.shape[0] == len(ref_bc)

    # --- Query: both ATAC populations (ATAC3k then ATAC10k) ---
    atac_bc = bc_rna3k + bc_rna10k
    H_query = np.vstack([Hatac3k, Hatac10k])
    print(f'[query] ATAC cells: {H_query.shape[0]}')

    # --- Reverse imputeKNN ---
    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=K, metric='euclidean').fit(ref_H)
    dist, idx = nbrs.kneighbors(H_query)
    e = np.exp(-dist)
    Wk = e / e.sum(axis=1, keepdims=True)
    from scipy.sparse import coo_matrix
    rows = np.repeat(np.arange(H_query.shape[0]), K)
    cols = idx.ravel()
    W = coo_matrix((Wk.ravel(), (rows, cols)), shape=(H_query.shape[0], ref_H.shape[0])).tocsr()
    imputed = (ref_expr @ W.T).astype(np.float32)   # genes x n_atac_query
    print('imputed expression genes x ATAC:', imputed.shape)

    np.save(f'{outdir}/imputed_expression_genes_x_atac.npy', imputed)
    with open(f'{outdir}/imputed_atac_barcodes.txt', 'w') as f:
        f.write('\n'.join(atac_bc) + '\n')
    with open(f'{outdir}/genes.txt', 'w') as f:
        f.write('\n'.join(genes) + '\n')
    with open(f'{outdir}/ref_block_barcodes.txt', 'w') as f:
        f.write('\n'.join(ref_bc) + '\n')

    # --- All-cells matrix: rows = [RNA3k real, RNA10k real, ATAC3k imputed, ATAC10k imputed] ---
    all_expr = np.vstack([r3.T, r10.T, imputed.T]).astype(np.float32)
    all_bc = bc_rna3k + bc_rna10k + atac_bc
    np.save(f'{outdir}/all_cells_gene_expression.npy', all_expr)
    np.save(f'{outdir}/genes.npy', np.array(genes))
    with open(f'{outdir}/all_cells_barcodes.txt', 'w') as f:
        f.write('\n'.join(all_bc) + '\n')
    print('all_cells_gene_expression shape:', all_expr.shape, 'dtype', all_expr.dtype)
    print(f'[ref] {ref_mode}  [exp] {exp}  DONE')

if __name__ == '__main__':
    main()
