#!/usr/bin/env python
"""Reverse-imputeKNN for the 4-dataset scenario.

Given the 4-dataset joint embedding H (blocks: rna3k, atac3k, rna6k, atac6k, each
2711x30) and a chosen RNA reference, propagate RNA expression onto both ATAC
populations (atac3k + atac6k).

The reference expression can be:
  - a single RNA dataset (3k only, or 6k only)  -> Exp B1 / B2
  - a SCEMENT-integrated combination of 3k+6k RNA -> Exp A  (produced by
    scripts/run_scement_combine.py; reference has 5422 cells, one block per orig)

Query = atac3k cells + atac6k cells (both imputed).

Output (to results/<exp>):
  imputed_expression_genes_x_atac.npy  (genes x 5422; ATAC3k rows then ATAC6k)
  imputed_atac_barcodes.txt            (5422 barcodes: ATAC3k then ATAC6k)
  imputed_atac_genes.txt
  all_cells_gene_expression.npy        (n_cells x genes)
  all_cells_barcodes.txt
  genes.npy
  ref_block_barcodes.txt               (reference cell order)

The all-cells matrix always contains the same real-RNA cells (RNA3k then RNA6k)
as rows 0..5421, then the imputed ATAC cells (ATAC3k then ATAC6k) as rows
5422..10843. Only the imputed ATAC values vary by experiment, which isolates the
reverse-imputeKNN reference strategy as the sole variable.

Usage: python build_imputation.py <exp> <ref_mode>
  exp       : expA | expB1 | expB2
  ref_mode  : combined | rna3k | rna6k
    combined -> load SCEMENT-combined reference (results/expA/combined_ref*)
    rna3k    -> use only the 3k RNA dataset as reference
    rna6k    -> use only the 6k RNA dataset as reference
"""
import os, sys, numpy as np, scipy.io, scipy.sparse as sp

# Project root: set via env PBSC4K_ROOT, else 2 parents above this script
PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
RES = f'{PROJ}/results/integration'
DATA = f'{PROJ}/data'
K = 20

def load_mtx_counts(d):
    X = scipy.io.mmread(f'{DATA}/{d}/counts.mtx').tocsr().toarray()  # cells x features
    bc = [l.strip() for l in open(f'{DATA}/{d}/barcodes.txt')]
    feat = [l.strip() for l in open(f'{DATA}/{d}/features.txt')]
    return X, bc, feat

def log1p_cpm(X):
    tot = X.sum(axis=1, keepdims=True); tot[tot == 0] = 1
    return np.log1p(X / tot * 1e4).T  # genes x cells

def main():
    exp, ref_mode = sys.argv[1], sys.argv[2]
    outdir = f'{PROJ}/results/{exp}'
    os.makedirs(outdir, exist_ok=True)

    # --- RNA real expression (genes x cells), shared for all experiments ---
    rna3k, bc_rna3k, genes = load_mtx_counts('3k_rna')
    rna6k, bc_rna6k, _ = load_mtx_counts('6k_rna')
    assert rna3k.shape[1] == rna6k.shape[1] == len(genes), f"gene mismatch {rna3k.shape[1]} {len(genes)}"
    # keep only common gene set (both use identical feature list)
    r3 = log1p_cpm(rna3k)   # genes x 2711
    r6 = log1p_cpm(rna6k)   # genes x 2711

    # --- Load joint embedding blocks ---
    H = np.load(f'{RES}/joint_embedding_H.npy')
    n = 2711
    Hrna3k, Hatac3k, Hrna6k, Hatac6k = H[0:n], H[n:2*n], H[2*n:3*n], H[3*n:4*n]

    # --- Build reference expression + reference embedding ---
    if ref_mode == 'rna3k':
        ref_expr = r3
        ref_H = Hrna3k
        ref_bc = bc_rna3k
        ref_tag = 'rna3k'
    elif ref_mode == 'rna6k':
        ref_expr = r6
        ref_H = Hrna6k
        ref_bc = bc_rna6k
        ref_tag = 'rna6k'
    elif ref_mode == 'combined':
        # SCEMENT-combined reference (genes x 5422: 3k then 6k)
        ref_expr = np.load(f'{PROJ}/results/expA/combined_ref_expression.npy')
        ref_H = np.load(f'{PROJ}/results/expA/combined_ref_embedding.npy')
        ref_bc = [l.strip() for l in open(f'{PROJ}/results/expA/combined_ref_barcodes.txt')]
        ref_tag = 'combined'
    else:
        raise SystemExit(f'unknown ref_mode {ref_mode}')
    print(f'[ref] {ref_mode}: reference expression {ref_expr.shape}, embedding {ref_H.shape}, cells {len(ref_bc)}')
    assert ref_expr.shape[1] == ref_H.shape[0] == len(ref_bc)

    # --- Query: both ATAC populations (ATAC3k then ATAC6k) ---
    atac_bc = bc_rna3k + bc_rna6k
    H_query = np.vstack([Hatac3k, Hatac6k])
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
    imputed = ref_expr @ W.T   # genes x n_atac_query
    print('imputed expression genes x ATAC:', imputed.shape)

    np.save(f'{outdir}/imputed_expression_genes_x_atac.npy', imputed)
    with open(f'{outdir}/imputed_atac_barcodes.txt', 'w') as f:
        f.write('\n'.join(atac_bc) + '\n')
    with open(f'{outdir}/genes.txt', 'w') as f:
        f.write('\n'.join(genes) + '\n')
    with open(f'{outdir}/ref_block_barcodes.txt', 'w') as f:
        f.write('\n'.join(ref_bc) + '\n')

    # --- All-cells matrix: rows = [RNA3k real, RNA6k real, ATAC3k imputed, ATAC6k imputed] ---
    all_expr = np.vstack([r3.T, r6.T, imputed.T])
    all_bc = bc_rna3k + bc_rna6k + atac_bc
    np.save(f'{outdir}/all_cells_gene_expression.npy', all_expr)
    np.save(f'{outdir}/genes.npy', np.array(genes))
    with open(f'{outdir}/all_cells_barcodes.txt', 'w') as f:
        f.write('\n'.join(all_bc) + '\n')
    print('all_cells_gene_expression shape:', all_expr.shape)
    print(f'[ref] {ref_tag}  [exp] {exp}  DONE')

if __name__ == '__main__':
    main()
