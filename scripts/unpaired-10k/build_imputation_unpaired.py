#!/usr/bin/env python
"""Reverse-imputeKNN for the UNPAIRED 10k scenario.

Joint embedding H (results/integration_unpaired/joint_embedding_H.npy) has blocks
  [ rna10k (11,898 real RNA cells), atac10k_ext (8,161 external ATAC cells) ].
H is the anchor-first ordering written by scmint.scsaga_saveH.main.

RNA is the reference (its real expression is known); the external ATAC cells are
the query. For each query cell we take its k=20 nearest reference cells in H,
weight them by softmax(-distance) (here exp(-dist) normalised per row, matching
the repo's build_imputation_10k.py), and propagate RNA expression:
    imputed(genes x query) = ref_expr(genes x ref) @ W.T(ref x query)

Output (results/unpaired_grn/):
  imputed_expression_genes_x_atac.npy  (36601 x 8161)
  imputed_atac_barcodes.txt            (8161)
  genes.txt                            (36601)
  all_cells_gene_expression.npy        (20059 x 36601, float32)
  all_cells_barcodes.txt               (20059)
  genes.npy
  n_cells.txt

All-cells rows: [RNA10k real (11898), ATAC10k-ext imputed (8161)].
"""
import os, numpy as np, scipy.io, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
# Repo root: $UNPAIRED_ROOT, else two levels up (scripts/<workflow>/ -> repo).
# NOTE: one os.pardir would resolve to scripts/, which is a bug when the
# env var is unset (results/ and data/ then point inside scripts/).
PROJ = os.environ.get('UNPAIRED_ROOT', os.path.abspath(os.path.join(
    HERE, os.pardir, os.pardir)))
RES = f'{PROJ}/results/integration_unpaired'
DATA = f'{PROJ}/data'
OUT = f'{PROJ}/results/unpaired_grn'
K = 20
N_RNA = 11898
N_ATAC = 8161


def log1p_cpm_genes_by_cells(path, gene_chunk=4000):
    """cells x features MTX -> genes x cells float32 log1p(CPM/1e4), chunked over genes."""
    X = scipy.io.mmread(path).tocsr()           # cells x features
    n_cells, n_feat = X.shape
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    scale = (1e4 / tot).astype(np.float64)
    Xcsc = X.multiply(scale[:, None]).tocsc()
    out = np.empty((n_feat, n_cells), dtype=np.float32)
    for lo in range(0, n_feat, gene_chunk):
        hi = min(lo + gene_chunk, n_feat)
        sub = Xcsc[:, lo:hi].toarray().astype(np.float32)   # cells x chunk
        np.log1p(sub, out=sub)
        out[lo:hi, :] = sub.T
    return out


def main():
    os.makedirs(OUT, exist_ok=True)

    genes = [l.strip() for l in open(f'{DATA}/10k_rna/features.txt') if l.strip()]
    bc_rna = [l.strip() for l in open(f'{DATA}/10k_rna/barcodes.txt') if l.strip()]
    bc_atac = [l.strip() for l in open(f'{DATA}/atac10k_ext/barcodes.txt') if l.strip()]
    assert len(bc_rna) == N_RNA and len(bc_atac) == N_ATAC

    rna = log1p_cpm_genes_by_cells(f'{DATA}/10k_rna/counts.mtx')
    assert rna.shape == (len(genes), N_RNA), rna.shape
    print('RNA reference expression (genes x cells):', rna.shape, rna.dtype, flush=True)

    H = np.load(f'{RES}/joint_embedding_H.npy')
    print('joint embedding H:', H.shape, flush=True)
    assert H.shape[0] == N_RNA + N_ATAC
    H_ref, H_query = H[:N_RNA], H[N_RNA:]

    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse import coo_matrix
    nbrs = NearestNeighbors(n_neighbors=K, metric='euclidean').fit(H_ref)
    dist, idx = nbrs.kneighbors(H_query)
    e = np.exp(-dist)
    Wk = e / e.sum(axis=1, keepdims=True)
    rows = np.repeat(np.arange(H_query.shape[0]), K)
    W = coo_matrix((Wk.ravel(), (rows, idx.ravel())),
                   shape=(H_query.shape[0], H_ref.shape[0])).tocsr()
    print('query/reference cells:', H_query.shape[0], H_ref.shape[0],
          'mean nn dist:', float(dist.mean()), flush=True)

    imputed = (rna @ W.T).astype(np.float32)      # genes x n_atac
    print('imputed expression (genes x ATAC):', imputed.shape, flush=True)
    del W
    np.save(f'{OUT}/imputed_expression_genes_x_atac.npy', imputed)
    with open(f'{OUT}/imputed_atac_barcodes.txt', 'w') as f:
        f.write('\n'.join(bc_atac) + '\n')
    with open(f'{OUT}/genes.txt', 'w') as f:
        f.write('\n'.join(genes) + '\n')

    all_expr = np.vstack([rna.T, imputed.T]).astype(np.float32)
    all_bc = bc_rna + bc_atac
    print('all-cells matrix:', all_expr.shape, all_expr.dtype, flush=True)
    np.save(f'{OUT}/all_cells_gene_expression.npy', all_expr)
    np.save(f'{OUT}/genes.npy', np.array(genes))
    with open(f'{OUT}/all_cells_barcodes.txt', 'w') as f:
        f.write('\n'.join(all_bc) + '\n')
    with open(f'{OUT}/n_cells.txt', 'w') as f:
        f.write(str(all_expr.shape[0]))
    print('DONE ->', OUT, flush=True)


if __name__ == '__main__':
    main()
