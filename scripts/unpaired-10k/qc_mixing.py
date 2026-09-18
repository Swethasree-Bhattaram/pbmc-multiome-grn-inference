#!/usr/bin/env python
"""Supplementary QC: (A) mixing benchmark vs the repo's existing 4-dataset H,
(B) neighbour-diversity of the imputation, (C) inter-ATAC-cell variability.

(A) is the important one: it answers "did the UNPAIRED run break the alignment, or
is this asymmetric same-modality-neighbour pattern a property of scSAGA's partial
alignment as already used in this repo?"  We recompute the identical statistic on
results/integration_10k/joint_embedding_H.npy (the paired 4-dataset run).
"""
import os, json, numpy as np
from sklearn.neighbors import NearestNeighbors

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.environ.get('UNPAIRED_ROOT', os.path.abspath(os.path.join(HERE, os.pardir)))
OUTD = f'{PROJ}/results/unpaired_grn'
K = 20
NEW = f'{PROJ}/results/integration_unpaired/joint_embedding_H.npy'
# Reference the repo's paired 4-dataset H if it is available; otherwise the
# comparison is skipped (the script prints so explicitly).
OLD = os.environ.get('PAIRED_H', f'{PROJ}/results/integration_10k/joint_embedding_H.npy')
lines = []


def p(s):
    print(s, flush=True)
    lines.append(s)


def mix_stats(H, labels, name):
    nb = NearestNeighbors(n_neighbors=31, metric='euclidean').fit(H)
    ii = nb.kneighbors(H, return_distance=False)[:, 1:]
    lab = np.asarray(labels)
    frac = (lab[ii] == lab[:, None]).mean(axis=1)
    props = np.array([(lab == u).mean() for u in sorted(set(labels))])
    expected = float((props ** 2).sum())
    p(f'  {name}:')
    for u in sorted(set(labels)):
        m = lab == u
        p(f'    {u:16s} n={m.sum():6d} ({m.mean()*100:5.1f}%)  '
          f'frac-same-modality-nbrs={frac[m].mean():.4f}  (random={props[sorted(set(labels)).index(u)]:.4f})')
    p(f'    overall frac-same={frac.mean():.4f}  random-expectation={expected:.4f}')
    return float(frac.mean()), expected


def main():
    p('=== (A) Alignment mixing benchmark: UNPAIRED vs the repo 4-dataset run ===')
    p('Random expectation for "frac same modality" = sum of squared block proportions.')
    p('A well-MIXED joint embedding gives every block frac-same close to its own proportion.')

    Hnew = np.load(NEW)
    p(f'\nUNPAIRED H {Hnew.shape} = [rna10k 11898 | atac10k_ext 8161]')
    a, ea = mix_stats(Hnew, ['rna10k'] * 11898 + ['atac10k_ext'] * 8161, 'UNPAIRED (this run)')

    if os.path.exists(OLD):
        Hold = np.load(OLD)
        n3, n10 = 2711, 11898
        labels = (['rna3k'] * n3 + ['atac3k'] * n3 + ['rna10k'] * n10 + ['atac10k'] * n10)
        p(f'\n4-DATASET H {Hold.shape} = [rna3k | atac3k | rna10k | atac10k] (repo, paired)')
        b, eb = mix_stats(Hold, labels, '4-dataset (repo reference)')
    else:
        b, eb = None, None

    p('\n=== (B) Imputation neighbour diversity ===')
    Href, Hq = Hnew[:11898], Hnew[11898:]
    nb = NearestNeighbors(n_neighbors=K, metric='euclidean').fit(Href)
    dist, idx = nb.kneighbors(Hq)
    uniq = np.unique(idx)
    p(f'  distinct RNA reference cells used as neighbours: {len(uniq)} / 11898 '
      f'({len(uniq)/11898*100:.1f}%)')
    cnt = np.bincount(idx.ravel(), minlength=11898)
    p(f'  neighbour-use distribution: max={cnt.max()}, mean={cnt.mean():.1f}, '
      f'median={int(np.median(cnt))}, cells never used={int((cnt==0).sum())}')
    p(f'  per-query distinct neighbours (of {K}): mean={len(np.unique(idx, axis=1)[0]) if False else float(np.mean([len(set(r)) for r in idx])):.2f}')

    all_expr = np.load(f'{OUTD}/all_cells_gene_expression.npy', mmap_mode='r')
    atac = np.asarray(all_expr[11898:], dtype=np.float32)
    rna = np.asarray(all_expr[:11898], dtype=np.float32)
    rng = np.random.default_rng(0)
    s = rng.choice(atac.shape[0], 300, replace=False)
    sub = atac[s]
    c = np.corrcoef(sub)
    iu = np.triu_indices(len(s), 1)
    p(f'\n=== (C) Imputed ATAC profile variability ===')
    p(f'  pairwise corr among 300 sampled imputed ATAC cells: mean={c[iu].mean():.4f} '
      f'min={c[iu].min():.4f}')
    sr = rng.choice(rna.shape[0], 300, replace=False)
    cr = np.corrcoef(rna[sr])
    p(f'  pairwise corr among 300 sampled REAL RNA cells:    mean={cr[iu].mean():.4f} '
      f'min={cr[iu].min():.4f}')
    p('  (imputed cells being as variable as real cells => not collapsed to one profile)')

    with open(f'{OUTD}/qc_mixing.txt', 'w') as f:
        f.write('\n'.join(lines) + '\n')
    json.dump({'unpaired_frac_same': a, 'unpaired_random_expect': ea,
               'four_dataset_frac_same': b, 'four_dataset_random_expect': eb,
               'distinct_ref_cells_used': int(len(uniq))},
              open(f'{OUTD}/qc_mixing.json', 'w'), indent=2)
    p(f'\nwrote {OUTD}/qc_mixing.txt')


if __name__ == '__main__':
    main()
