#!/usr/bin/env python
"""RNA-only baseline: stack PBMC 3k (multiome RNA) + PBMC 4k (10x v2) gene
expression matrices for Arboreto.

No integration / no imputation: this is a pure two-RNA-matrix baseline. The two
matrices come from different Cell Ranger references (3k = GRCh38-2020-A,
36,601 genes; 4k = GRCh38-3.0.0, 33,694 genes), so we restrict to the gene
symbols common to both and stack cells (3k rows first, then 4k).

Gene symbols that occur more than once in a reference (multiple Ensembl IDs
sharing a symbol) are aggregated by summing their counts, so each symbol has
exactly one column in each matrix before intersecting.

Each matrix is log1p(CPM) normalized independently (per-cell library-size
scaling), the same transform the reverse-imputeKNN pipeline applies to its
reference RNA.

Inputs  (produced as for the other workflows; see PIPELINE / ENVIRONMENT):
  data/3k_rna/counts.mtx barcodes.txt features.txt   # from the 3k multiome h5
  data/4k_rna/matrix.mtx barcodes.txt genes.tsv      # 10x v2 4k filtered matrices

Output (results/rna_baseline):
  all_cells_gene_expression.npy  (7051 x 21932) float32
  genes.npy
  all_cells_barcodes.txt
  n_cells.txt
  dataset_sizes.txt

Run with the scSAGA .venv (numpy/scipy).
"""
import os, numpy as np, scipy.io, scipy.sparse as sp

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
DATA = f'{PROJ}/data'
OUT = f'{PROJ}/results/rna_baseline'
os.makedirs(OUT, exist_ok=True)


def load_3k():
    """3k multiome RNA: counts.mtx is cells x features (gene symbols)."""
    X = scipy.io.mmread(f'{DATA}/3k_rna/counts.mtx').tocsr()          # cells x features
    bc = [l.strip() for l in open(f'{DATA}/3k_rna/barcodes.txt')]
    feat = [l.strip() for l in open(f'{DATA}/3k_rna/features.txt')]
    assert X.shape == (len(bc), len(feat)), (X.shape, len(bc), len(feat))
    return X, bc, feat


def load_4k():
    """10x v2 4k: matrix.mtx is genes x cells -> transpose to cells x genes."""
    M = scipy.io.mmread(f'{DATA}/4k_rna/matrix.mtx').tocsr()          # genes x cells
    bc = [l.strip() for l in open(f'{DATA}/4k_rna/barcodes.txt')]
    genes = [l.rstrip('\n').split('\t')[1] for l in open(f'{DATA}/4k_rna/genes.tsv')]
    assert M.shape == (len(genes), len(bc)), (M.shape, len(genes), len(bc))
    return M.T.tocsr(), bc, genes                                     # cells x genes


def aggregate_by_symbol(X, symbols):
    """Sum columns that share a gene symbol. Returns (X_agg, unique_symbols, n_merged)."""
    uniq, first_idx = [], {}
    for i, s in enumerate(symbols):
        if s not in first_idx:
            first_idx[s] = len(uniq)
            uniq.append(s)
    if len(uniq) == len(symbols):
        return X, uniq, 0
    rows = np.array([first_idx[s] for s in symbols])
    cols = np.arange(len(symbols))
    A = sp.coo_matrix((np.ones(len(symbols)), (rows, cols)),
                      shape=(len(uniq), len(symbols))).tocsr()
    return (X @ A.T).tocsr(), uniq, len(symbols) - len(uniq)


def log1p_cpm_sparse(X):
    """Per-cell CPM scaling + log1p, preserving sparsity."""
    X = X.tocsr().astype(np.float64)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0
    X = X.multiply(1e4 / tot[:, None]).tocsr()
    X.data = np.log1p(X.data)
    return X


def main():
    X3, bc3, g3 = load_3k()
    X4, bc4, g4 = load_4k()
    print(f'3k: {X3.shape} cells x genes (symbols)')
    print(f'4k: {X4.shape} cells x genes (symbols)')

    X3, u3, d3 = aggregate_by_symbol(X3, g3)
    X4, u4, d4 = aggregate_by_symbol(X4, g4)
    print(f'3k: {len(g3)} -> {len(u3)} unique symbols ({d3} merged)')
    print(f'4k: {len(g4)} -> {len(u4)} unique symbols ({d4} merged)')

    set4 = set(u4)
    gene_order = [g for g in u3 if g in set4]
    print(f'common symbols: {len(gene_order)}')

    idx3 = {g: i for i, g in enumerate(u3)}
    idx4 = {g: i for i, g in enumerate(u4)}
    c3 = [idx3[g] for g in gene_order]
    c4 = [idx4[g] for g in gene_order]

    N3 = log1p_cpm_sparse(X3)[:, c3].toarray().astype(np.float32)
    N4 = log1p_cpm_sparse(X4)[:, c4].toarray().astype(np.float32)
    print(f'normalized+subset  3k: {N3.shape}  4k: {N4.shape}')

    all_expr = np.vstack([N3, N4]).astype(np.float32)
    all_bc = bc3 + bc4
    print(f'all_cells_gene_expression: {all_expr.shape}')

    np.save(f'{OUT}/all_cells_gene_expression.npy', all_expr)
    np.save(f'{OUT}/genes.npy', np.array(gene_order))
    with open(f'{OUT}/all_cells_barcodes.txt', 'w') as f:
        f.write('\n'.join(all_bc) + '\n')
    with open(f'{OUT}/n_cells.txt', 'w') as f:
        f.write(str(all_expr.shape[0]))
    with open(f'{OUT}/dataset_sizes.txt', 'w') as f:
        f.write(f'3k cells: {len(bc3)}\n4k cells: {len(bc4)}\n'
                f'3k genes (GRCh38-2020-A): {len(g3)}\n'
                f'4k genes (GRCh38-3.0.0): {len(g4)}\n'
                f'common gene symbols: {len(gene_order)}\n')
    print('Saved to', OUT)


if __name__ == '__main__':
    main()
