#!/usr/bin/env python
"""QC for the UNPAIRED 10k imputation + integration.

Checks that the reverse-imputeKNN output is not degenerate:
  1. Real RNA block vs imputed ATAC block: mean/std/nonzero per cell.
  2. Correlation of each imputed ATAC cell with its assigned RNA neighbours vs a
     random set of RNA cells (neighbour assignment must be informative).
  3. Lineage-marker sanity: imputed ATAC profiles should recover a plausible
     PBMC marker structure (CD3D/MS4A1/CD79A/LYZ/NKG7/CD14).
  4. Integration embedding: is the ATAC block distributed across RNA (not collapsed
     into one corner)?  kNN label mixing statistic.
Writes results/unpaired_grn/qc_imputation.txt and qc_report.json.
"""
import os, json, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.environ.get('UNPAIRED_ROOT', os.path.abspath(os.path.join(HERE, os.pardir)))
RES = f'{PROJ}/results/integration_unpaired'
OUTD = f'{PROJ}/results/unpaired_grn'
K = 20
N_RNA = 11898
N_ATAC = 8161
MARKERS = ['CD3D', 'CD3E', 'MS4A1', 'CD79A', 'LYZ', 'NKG7', 'GNLY', 'CD14', 'FCGR3A', 'PPBP']

lines, payload = [], {}


def p(s):
    print(s, flush=True)
    lines.append(s)


def main():
    genes = [l.strip() for l in open(f'{OUTD}/genes.txt') if l.strip()]
    gi = {g: i for i, g in enumerate(genes)}

    all_expr = np.load(f'{OUTD}/all_cells_gene_expression.npy', mmap_mode='r')
    rna = np.asarray(all_expr[:N_RNA], dtype=np.float32)
    atac = np.asarray(all_expr[N_RNA:], dtype=np.float32)
    p(f'all-cells matrix: {all_expr.shape} (RNA real {rna.shape}, ATAC imputed {atac.shape})')

    p('\n--- 1. Per-cell distributional summary (log1p CPM/1e4 units) ---')
    for name, M in [('RNA real', rna), ('ATAC imputed', atac)]:
        mean = M.mean(axis=1); std = M.std(axis=1)
        nz = (M > 0).mean(axis=1)
        p(f'{name:14s} mean-of-cell-means={mean.mean():.4f}  mean-std={std.mean():.4f}  '
          f'mean-frac-nonzero={nz.mean():.4f}  library(mean sum)={M.sum(axis=1).mean():.1f}')
    p(f'global: RNA real mean={rna.mean():.4f}  ATAC imputed mean={atac.mean():.4f}')

    p('\n--- 2. Is the kNN assignment informative? ---')
    H = np.load(f'{RES}/joint_embedding_H.npy')
    H_ref, H_q = H[:N_RNA], H[N_RNA:]
    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=K, metric='euclidean').fit(H_ref)
    dist, idx = nbrs.kneighbors(H_q)

    # correlation of each imputed cell with (a) its true neighbours, (b) random RNA cells
    rng = np.random.default_rng(0)
    samp = rng.choice(N_ATAC, size=200, replace=False)
    r_true, r_rand = [], []
    for s in samp:
        q = atac[s]
        ref = rna[idx[s]].mean(axis=0)
        rnd = rna[rng.choice(N_RNA, size=K, replace=False)].mean(axis=0)
        r_true.append(np.corrcoef(q, ref)[0, 1])
        r_rand.append(np.corrcoef(q, rnd)[0, 1])
    p(f'mean corr(imputed cell, its {K} RNA neighbours) = {np.nanmean(r_true):.4f}')
    p(f'mean corr(imputed cell, {K} random RNA cells)     = {np.nanmean(r_rand):.4f}')
    payload['corr_true_nbr'] = float(np.nanmean(r_true))
    payload['corr_rand_nbr'] = float(np.nanmean(r_rand))
    p(f'mean NN distance of ATAC queries to RNA: {dist.mean():.5f} (median {np.median(dist):.5f})')

    p('\n--- 3. Lineage-marker recovery in imputed ATAC vs real RNA ---')
    rows = []
    for g in MARKERS:
        if g not in gi:
            continue
        a = float(atac[:, gi[g]].mean()); r = float(rna[:, gi[g]].mean())
        rows.append((g, r, a, a / r if r > 1e-9 else float('nan')))
        p(f'  {g:8s} RNA-real mean={r:.4f}   ATAC-imputed mean={a:.4f}   ratio={a/r if r>1e-9 else float("nan"):.2f}')
    payload['markers'] = [{'gene': g, 'rna': r, 'atac': a, 'ratio': rt} for g, r, a, rt in rows]

    p('\n--- 4. Integration: is the ATAC block spread across the RNA manifold? ---')
    comb = np.vstack([H_ref, H_q])
    labels = np.array([0] * N_RNA + [1] * N_ATAC)
    nb = NearestNeighbors(n_neighbors=31, metric='euclidean').fit(comb)
    ii = nb.kneighbors(comb, return_distance=False)[:, 1:]
    frac_same = (labels[ii] == labels[:, None]).mean(axis=1)
    p(f'mean fraction of same-modality neighbours: ATAC={frac_same[N_RNA:].mean():.4f} '
      f'RNA={frac_same[:N_RNA].mean():.4f}  (0.5 = perfectly mixed)')
    p(f'ATAC block: {N_ATAC} cells = {N_ATAC/(N_RNA+N_ATAC)*100:.1f}% of H')
    payload['mix_atac'] = float(frac_same[N_RNA:].mean())
    payload['mix_rna'] = float(frac_same[:N_RNA].mean())

    with open(f'{OUTD}/qc_imputation.txt', 'w') as f:
        f.write('\n'.join(lines) + '\n')
    with open(f'{OUTD}/qc_report.json', 'w') as f:
        json.dump(payload, f, indent=2)
    p(f'\nwrote {OUTD}/qc_imputation.txt')


if __name__ == '__main__':
    main()
