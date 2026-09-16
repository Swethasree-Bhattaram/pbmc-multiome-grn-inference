#!/usr/bin/env python
"""Write reports/rna-baseline/REPORT_tfonly.md — RNA-only baseline with the
tf_only.txt regulator list, side by side with the original TRRUST-regulator run.

Reads:
  results/rna_baseline/grn_tfonly_reg/                 (this run)
  /Volumes/samsung_ssd/tmp/pbmc-4k-baseline/results/rna_baseline/grn_trrust/  (original)
Writes:
  <repo>/reports/rna-baseline/REPORT_tfonly.md

Run with the scSAGA .venv (numpy/pandas).
"""
import os, re
import numpy as np, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', '/Users/sbhattaram/pbmc-multiome-grn-inference')
ORIG = '/Volumes/samsung_ssd/tmp/pbmc-4k-baseline/results/rna_baseline/grn_trrust'
# the network .tsv is gitignored, so this run's GRN dir may live outside the repo
GRN = os.environ.get('GRN_DIR', f'{PROJ}/results/rna_baseline/grn_tfonly_reg')
REP = f'{PROJ}/reports/rna-baseline'
DOWN = f'{PROJ}/results/rna_baseline'


def parse_eval(path):
    ev, cur = {}, None
    for line in open(path):
        m = re.match(r'========== EVALUATION vs (.*) ==========', line.strip())
        if m:
            cur = m.group(1); ev[cur] = {'precrec': []}; continue
        if cur is None:
            continue
        m = re.match(r'Total edges: (\d+)', line)
        if m: ev[cur]['total'] = int(m.group(1)); continue
        m = re.match(r'Deduplicated edges: (\d+)', line)
        if m: ev[cur]['dedup'] = int(m.group(1)); continue
        m = re.match(r'Edges recovered by GRN: (\d+) \(([\d.]+)%\)', line)
        if m: ev[cur]['rec'] = int(m.group(1)); ev[cur]['recpct'] = float(m.group(2)); continue
        m = re.match(r'\s*(\d+)\s+([\d.]+)\s+([\d.]+)', line)
        if m: ev[cur]['precrec'].append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    return ev


def load_run(base):
    net = pd.read_csv(f'{base}/grnboost2_network.tsv', sep='\t')
    net = net.sort_values('importance', ascending=False).reset_index(drop=True)
    n_cells = int(open(f'{base}/n_cells.txt').read().strip())
    regs = len(open(f'{base}/tf_regulators.txt').read().splitlines())
    tgts = len(open(f'{base}/target_genes.txt').read().splitlines())
    ev = parse_eval(f'{base}/evaluation/evaluation_summary.txt')
    return dict(net=net, n_cells=n_cells, regs=regs, tgts=tgts, ev=ev)


def main():
    sizes = dict(l.strip().split(': ') for l in open(f'{DOWN}/dataset_sizes.txt'))
    # the all-cells matrix is gitignored (*.npy); allow an out-of-repo matrix dir
    mroot = os.environ.get('MATRIX_DIR', DOWN)
    X = np.load(f'{mroot}/all_cells_gene_expression.npy', mmap_mode='r')
    genes = np.load(f'{mroot}/genes.npy', allow_pickle=True).astype(str)

    T = load_run(GRN)     # tf_only regulators (this run)
    R = load_run(ORIG)    # TRRUST regulators (original run)
    union_cols = len(open(f'{GRN}/all_columns.txt').read().splitlines())

    L = []
    A = L.append
    A('# RNA-only baseline — tf_only.txt regulators — GRNBoost2 on stacked PBMC 3k + PBMC 4k RNA')
    A('')
    A('**Baseline:** no integration, no imputation — the PBMC 3k multiome RNA and the 10x PBMC 4k RNA '
      'matrices are stacked and handed directly to Arboreto.')
    A('**Pipeline:** identical to the original RNA-only baseline run; the single change is the '
      'transcription-factor list.')
    A('**GRN:** Arboreto GRNBoost2; targets = `trrust_tf.txt` genes present in the matrix, '
      '**regulators = `tf_only.txt` TFs present in the matrix**.')
    A('**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated).')
    A('')
    A('---')
    A('')
    A('## 1. Inputs (identical to the original RNA-only baseline)')
    A('')
    A(f'- **3k RNA:** PBMC multiome 3k granulocyte-sorted, Gene Expression '
      f'({sizes["3k genes (GRCh38-2020-A)"]} genes, GRCh38-2020-A), {sizes["3k cells"]} cells.')
    A(f'- **4k RNA:** 10x PBMC 4k from a healthy donor (v2 / Cell Ranger 2.1.0, '
      f'{sizes["4k genes (GRCh38-3.0.0)"]} genes, GRCh38-3.0.0), {sizes["4k cells"]} cells.')
    A(f'- **Stacking:** gene symbols common to both references (**{sizes["common gene symbols"]}** genes); '
      f'duplicate symbols summed per reference. Each matrix log1p(CPM)-normalized independently, then '
      f'cells stacked (3k rows first, then 4k).')
    A(f'- **All-cells matrix:** {X.shape[0]} cells x {X.shape[1]} genes (bit-identical to the matrix used '
      f'by the original run).')
    A('')
    A('## 2. Regulator / target specification — the one changed variable')
    A('')
    A('| | Original run (`grn_trrust`) | This run (`grn_tfonly_reg`) |')
    A('|---|---|---|')
    A(f'| Regulator list | `data/trrust_tf.txt` | `data/tf_only.txt` |')
    A(f'| Regulator list size | 2862 | 829 |')
    A(f'| Regulators present in matrix | {R["regs"]} | {T["regs"]} |')
    A(f'| Regulators appearing in the network | {R["net"]["TF"].nunique()} | {T["net"]["TF"].nunique()} |')
    A(f'| Target genes (`trrust_tf.txt` present) | {R["tgts"]} | {T["tgts"]} |')
    A(f'| Expression subset handed to GRNBoost2 | {R["regs"]} columns | '
      f'{union_cols} columns (union of targets + tf_only regulators) |')
    A(f'| Cells | {T["n_cells"]} | {T["n_cells"]} |')
    A(f'| Inferred edges | {len(R["net"]):,} | {len(T["net"]):,} |')
    A('')
    A(f'- Regulators dropped relative to the original run: {R["regs"] - T["regs"]} '
      f'(non-TF genes of the mixed TRRUST list).')
    A(f'- {T["regs"]} of the 829 `tf_only.txt` entries are present as columns in the matrix; the remaining '
      f'{829 - T["regs"]} are absent from both reference gene sets.')
    A('')
    A('## 3. Evaluation vs ground truth (deduplicated)')
    A('')
    A('### Edges recovered')
    A('')
    A('| Ground truth | Regulators | Dedup. GT edges | Recovered |')
    A('|---|---|---|---|')
    for gt in ['PBMC-TRRUST', 'PBMC-Blood']:
        e, r = T['ev'][gt], R['ev'][gt]
        A(f'| {gt} | tf_only ({T["regs"]}) | {e["dedup"]} | {e["rec"]}/{e["dedup"]} ({e["recpct"]}%) |')
        A(f'| {gt} | TRRUST ({R["regs"]}) | {r["dedup"]} | {r["rec"]}/{r["dedup"]} ({r["recpct"]}%) |')
    A('')
    A('### Precision@K / Recall@K — PBMC-TRRUST')
    A('')
    A('| Top-K | Prec@K (tf_only) | Recall@K (tf_only) | Prec@K (TRRUST) | Recall@K (TRRUST) |')
    A('|---|---|---|---|---|')
    for (K, p, rc), (K2, p2, rc2) in zip(T['ev']['PBMC-TRRUST']['precrec'],
                                         R['ev']['PBMC-TRRUST']['precrec']):
        A(f'| {K} | {p:.4f} | {rc:.4f} | {p2:.4f} | {rc2:.4f} |')
    A('')
    A('### Precision@K / Recall@K — PBMC-Blood')
    A('')
    A('| Top-K | Prec@K (tf_only) | Recall@K (tf_only) | Prec@K (TRRUST) | Recall@K (TRRUST) |')
    A('|---|---|---|---|---|')
    for (K, p, rc), (K2, p2, rc2) in zip(T['ev']['PBMC-Blood']['precrec'],
                                         R['ev']['PBMC-Blood']['precrec']):
        A(f'| {K} | {p:.4f} | {rc:.4f} | {p2:.4f} | {rc2:.4f} |')
    A('')
    A('## 4. Top 25 edges (tf_only regulators)')
    A('')
    A('| TF -> target | importance |')
    A('|---|---|')
    for _, r in T['net'].head(25).iterrows():
        A(f'| {r["TF"]} -> {r["target"]} | {r["importance"]:.4f} |')
    A('')
    A('## 5. Context: same tf_only regulators under integration + imputation')
    A('')
    A('For reference, the same regulator list (`tf_only.txt`, 816 present) was used on the')
    A('full-10k integration + reverse-imputeKNN runs of this repo. Those runs hold 29,218 cells')
    A('(14,609 real multiome RNA [2,711 from 3k + 11,898 from 10k] + 14,609 imputed ATAC) versus')
    A('the 7,051 real RNA cells here, and 2,827 targets (2,852 union columns) versus 2,808 (2,833):')
    A('')
    A('| Run | Cells | Inferred edges | PBMC-TRRUST recovered | PBMC-Blood recovered |')
    A('|---|---|---|---|---|')
    A(f'| RNA-only baseline, tf_only regs (this run) | {T["n_cells"]} | {len(T["net"]):,} | '
      f'{T["ev"]["PBMC-TRRUST"]["rec"]}/{T["ev"]["PBMC-TRRUST"]["dedup"]} '
      f'({T["ev"]["PBMC-TRRUST"]["recpct"]}%) | '
      f'{T["ev"]["PBMC-Blood"]["rec"]}/{T["ev"]["PBMC-Blood"]["dedup"]} '
      f'({T["ev"]["PBMC-Blood"]["recpct"]}%) |')
    A('| expA (paired full-10k, SCEMENT ref) | 29218 | 829982 | 3660/8751 (41.8%) | 5269/96846 (5.4%) |')
    A('| expB1 (paired full-10k, 3k RNA ref) | 29218 | 779701 | 3454/8751 (39.5%) | 5132/96846 (5.3%) |')
    A('| expB2 (paired full-10k, 10k RNA ref) | 29218 | 832664 | 3664/8751 (41.9%) | 5362/96846 (5.5%) |')
    A('| unpaired full-10k (external 10x ATAC) | 20059 | 659314 | 3110/8751 (35.5%) | 4486/96846 (4.6%) |')
    A('')
    A('*Paired/unpaired values from `reports/full_10k_3k_onlyTf/comparison_report_tfonly.md` and')
    A('`reports/unpaired-10k/comparison_report_tfonly_unpaired.md`; targets there are 2,827 vs 2,808')
    A('here because the all-cells matrices carry different gene sets. Cell count, imputation and')
    A('target count all differ, so this is context, not a controlled comparison.*')
    A('')
    A('---')

    os.makedirs(REP, exist_ok=True)
    out = f'{REP}/REPORT_tfonly.md'
    with open(out, 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('Wrote', out)


if __name__ == '__main__':
    main()
