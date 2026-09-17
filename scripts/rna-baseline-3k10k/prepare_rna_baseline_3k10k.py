#!/usr/bin/env python
"""RNA-only baseline #2: stack PBMC multiome 3k RNA + PBMC multiome 10k RNA gene
expression matrices for Arboreto.

Same philosophy as scripts/rna-baseline/prepare_rna_baseline.py (the 3k + 4k
baseline): no integration, no imputation -- a pure two-RNA-matrix baseline.

Here BOTH matrices come from the 10x PBMC granulocyte-sorted multiome libraries
and use the SAME Cell Ranger reference (GRCh38-2020-A, 36,601 genes), so they
share the same feature list and gene order. We keep the gene symbols common to
both (identical order), collapse any duplicate symbols by summing counts, and
stack cells (3k rows first, then 10k).

Each matrix is log1p(CPM) normalized independently (per-cell library-size
scaling), the same transform the reverse-imputeKNN pipeline applies to its
reference RNA.

Inputs (multiome Gene Expression; produced by scripts/common/extract_data.py
and the 10k full-multiome preprocessor):
  data/3k_rna/counts.mtx barcodes.txt features.txt   # 2,711 cells x 36,601 genes
  data/10k_rna/counts.mtx barcodes.txt features.txt  # 11,898 cells x 36,601 genes
  (both under the workspace root; override with RNA_3K_DIR / RNA_10K_DIR)

Output (results/rna_baseline_3k10k):
  all_cells_gene_expression.npy  (14609 x 36591) float32
  genes.npy
  all_cells_barcodes.txt
  n_cells.txt
  dataset_sizes.txt

Run with the scSAGA .venv (numpy/scipy).
"""
import os, numpy as np, scipy.io, scipy.sparse as sp

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
D3 = os.environ.get('RNA_3K_DIR', f'{PROJ}/data/3k_rna')
D10 = os.environ.get('RNA_10K_DIR', f'{PROJ}/data/10k_rna')
OUT = f'{PROJ}/results/rna_baseline_3k10k'
os.makedirs(OUT, exist_ok=True)


def load(dir_):
    """Load a multiome RNA matrix: cells x features (gene symbols)."""
    X = scipy.io.mmread(f'{dir_}/counts.mtx').tocsr()          # cells x features
    bc = [l.strip() for l in open(f'{dir_}/barcodes.txt')]
    feat = [l.strip() for l in open(f'{dir_}/features.txt')]
    assert X.shape == (len(bc), len(feat)), (X.shape, len(bc), len(feat))
    return X, bc, feat


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
    X3, bc3, g3 = load(D3)
    X10, bc10, g10 = load(D10)
    print(f'3k : {X3.shape} cells x genes (symbols)')
    print(f'10k: {X10.shape} cells x genes (symbols)')

    X3, u3, d3 = aggregate_by_symbol(X3, g3)
    X10, u10, d10 = aggregate_by_symbol(X10, g10)
    print(f'3k : {len(g3)} -> {len(u3)} unique symbols ({d3} merged)')
    print(f'10k: {len(g10)} -> {len(u10)} unique symbols ({d10} merged)')

    set10 = set(u10)
    gene_order = [g for g in u3 if g in set10]
    print(f'common symbols: {len(gene_order)}  (shared order: {u3 == u10})')

    idx3 = {g: i for i, g in enumerate(u3)}
    idx10 = {g: i for i, g in enumerate(u10)}
    c3 = [idx3[g] for g in gene_order]
    c10 = [idx10[g] for g in gene_order]

    N3 = log1p_cpm_sparse(X3)[:, c3].toarray().astype(np.float32)
    N10 = log1p_cpm_sparse(X10)[:, c10].toarray().astype(np.float32)
    print(f'normalized+subset  3k: {N3.shape}  10k: {N10.shape}')

    all_expr = np.vstack([N3, N10]).astype(np.float32)
    all_bc = bc3 + bc10
    print(f'all_cells_gene_expression: {all_expr.shape}')

    np.save(f'{OUT}/all_cells_gene_expression.npy', all_expr)
    np.save(f'{OUT}/genes.npy', np.array(gene_order))
    with open(f'{OUT}/all_cells_barcodes.txt', 'w') as f:
        f.write('\n'.join(all_bc) + '\n')
    with open(f'{OUT}/n_cells.txt', 'w') as f:
        f.write(str(all_expr.shape[0]))
    with open(f'{OUT}/dataset_sizes.txt', 'w') as f:
        f.write(f'3k cells: {len(bc3)}\n10k cells: {len(bc10)}\n'
                f'3k genes (GRCh38-2020-A): {len(g3)}\n'
                f'10k genes (GRCh38-2020-A): {len(g10)}\n'
                f'common gene symbols: {len(gene_order)}\n')
    print('Saved to', OUT)


if __name__ == '__main__':
    main()
