#!/usr/bin/env python
"""Decisive QC: does the UNPAIRED alignment carry real biological signal?

The mixing statistic is block-segregated in BOTH this run and the repo's own paired
4-dataset run (see qc_mixing.txt), so mixing alone cannot tell us whether the
alignment is meaningful.  This script tests biology instead:

  1. Cluster the 8,161 external ATAC cells independently, on their OWN peak-accessible
     PCA (data/atac10k_ext/pca_50.txt) -- this uses no RNA and no alignment.
  2. Cluster the 11,898 real RNA cells on their own PCA.
  3. For each ATAC cluster, report the mean IMPUTED marker profile.
     For each RNA cluster, report the mean REAL marker profile.
  4. Test whether the two cluster-marker structures agree (Hungarian-matched
     correlation between cluster-centroid marker profiles).

If the imputation assigned RNA neighbours essentially at random, imputed marker
levels would be nearly identical across ATAC clusters.  Coherent, PBMC-plausible
marker structure across independently-derived ATAC clusters is evidence the
alignment/imputation carries real signal.
"""
import os, numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.environ.get('UNPAIRED_ROOT', os.path.abspath(os.path.join(HERE, os.pardir)))
OUTD = f'{PROJ}/results/unpaired_grn'
DATA = f'{PROJ}/data'
N_RNA, N_ATAC, KC = 11898, 8161, 6
MARKERS = ['CD3D', 'CD3E', 'IL7R', 'CCR7', 'MS4A1', 'CD79A', 'CD79B', 'LYZ',
           'S100A8', 'S100A9', 'CD14', 'FCGR3A', 'NKG7', 'GNLY', 'PPBP', 'FCER1A']
lines = []


def p(s):
    print(s, flush=True)
    lines.append(s)


def main():
    genes = [l.strip() for l in open(f'{OUTD}/genes.txt') if l.strip()]
    gi = {g: i for i, g in enumerate(genes)}
    mk = [g for g in MARKERS if g in gi]

    atac_pca = np.loadtxt(f'{DATA}/atac10k_ext/pca_50.txt', dtype=np.float32)
    rna_pca = np.loadtxt(f'{DATA}/10k_rna/pca_50.txt', dtype=np.float32)
    p(f'independent PCAs: ATAC {atac_pca.shape}, RNA {rna_pca.shape}')

    ca = KMeans(n_clusters=KC, n_init=10, random_state=0).fit_predict(atac_pca)
    cr = KMeans(n_clusters=KC, n_init=10, random_state=0).fit_predict(rna_pca)
    p(f'kmeans k={KC}: ATAC cluster sizes {np.bincount(ca).tolist()}')
    p(f'kmeans k={KC}: RNA  cluster sizes {np.bincount(cr).tolist()}')

    expr = np.load(f'{OUTD}/all_cells_gene_expression.npy', mmap_mode='r')
    midx = [gi[g] for g in mk]
    rna_m = np.asarray(expr[:N_RNA][:, midx], dtype=np.float32)   # RNA x markers
    atac_m = np.asarray(expr[N_RNA:][:, midx], dtype=np.float32)  # ATAC x markers

    A = np.vstack([atac_m[ca == c].mean(axis=0) for c in range(KC)])   # KC x markers (imputed)
    R = np.vstack([rna_m[cr == c].mean(axis=0) for c in range(KC)])    # KC x markers (real)

    def table(M, title):
        p(f'\n{title}')
        p('cluster  n     ' + '  '.join(f'{g:>7s}' for g in mk))
        for c in range(KC):
            n = (ca == c).sum() if title.startswith('ATAC') else (cr == c).sum()
            p(f'  C{c}    {n:5d}  ' + '  '.join(f'{v:7.3f}' for v in M[c]))

    table(R, f'REAL RNA marker profile per RNA cluster (n={N_RNA})')
    table(A, f'IMPUTED marker profile per ATAC cluster (n={N_ATAC})')
    p('\n(ATAC clusters come from the ATAC peak PCA alone; the values are the imputed'
      '\n RNA expression of those cells -- coherence = alignment carries signal)')

    # variability across ATAC clusters relative to within-cluster spread
    p('\n--- signal-to-noise of imputed markers across ATAC clusters ---')
    for j, g in enumerate(mk):
        between = A[:, j].std()
        within = np.mean([atac_m[ca == c][:, j].std() for c in range(KC)])
        p(f'  {g:8s} between-cluster sd={between:.4f}  mean within-cluster sd={within:.4f}  '
          f'ratio={between/within if within > 1e-9 else float("nan"):.3f}')

    # Do the two marker structures agree after optimal cluster matching?
    Az = StandardScaler().fit_transform(A)
    Rz = StandardScaler().fit_transform(R)
    C = 1 - np.corrcoef(Az, Rz)[:KC, KC:]
    ri, ci = linear_sum_assignment(C)
    matched = [np.corrcoef(Az[i], Rz[j])[0, 1] for i, j in zip(ri, ci)]
    p('\n--- agreement of ATAC-cluster imputed profile with best-matched RNA cluster ---')
    for i, j in zip(ri, ci):
        p(f'  ATAC C{i} <-> RNA C{j}: marker-profile r = {np.corrcoef(Az[i], Rz[j])[0,1]:+.3f}')
    p(f'mean matched r = {np.mean(matched):+.3f}')

    # Random-neighbour control: rebuild imputation with shuffled reference assignment
    rng = np.random.default_rng(0)
    perm = rng.permutation(N_ATAC)
    A_ctrl = np.vstack([atac_m[perm][ca == c].mean(axis=0) for c in range(KC)])
    Cz = 1 - np.corrcoef(StandardScaler().fit_transform(A_ctrl),
                         StandardScaler().fit_transform(R))[:KC, KC:]
    ri2, ci2 = linear_sum_assignment(Cz)
    ctrl = [np.corrcoef(StandardScaler().fit_transform(A_ctrl)[i],
                        StandardScaler().fit_transform(R)[j])[0, 1] for i, j in zip(ri2, ci2)]
    p(f'\nCONTROL (imputed profiles shuffled across ATAC clusters), mean matched r = {np.mean(ctrl):+.3f}')
    p('  -> if the real value clearly exceeds this control, the cluster structure of the')
    p('     imputation is driven by the ATAC cells themselves, not by chance.')

    with open(f'{OUTD}/qc_biology.txt', 'w') as f:
        f.write('\n'.join(lines) + '\n')
    p(f'\nwrote {OUTD}/qc_biology.txt')


if __name__ == '__main__':
    main()
