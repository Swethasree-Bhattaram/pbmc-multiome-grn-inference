#!/usr/bin/env python
"""Write reports/rna-baseline-3k10k/REPORT_tfonly.md — the 3k+10k multiome
RNA-only baseline (tf_only.txt regulators), with the 3k+4k tf_only baseline and
the full-10k integration runs as context.

Reads:
  results/rna_baseline_3k10k/grn_tfonly_reg/      (this run)
Writes:
  <repo>/reports/rna-baseline-3k10k/REPORT_tfonly.md

Run with the scSAGA .venv (numpy/pandas).
"""
import os, re
import numpy as np, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', '/Users/sbhattaram/pbmc-multiome-grn-inference')
GRN = os.environ.get('GRN_DIR', f'{PROJ}/results/rna_baseline_3k10k/grn_tfonly_reg')
DOWN = f'{PROJ}/results/rna_baseline_3k10k'
REP = f'{PROJ}/reports/rna-baseline-3k10k'


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
    mroot = os.environ.get('MATRIX_DIR', DOWN)
    X = np.load(f'{mroot}/all_cells_gene_expression.npy', mmap_mode='r')
    genes = np.load(f'{mroot}/genes.npy', allow_pickle=True).astype(str)

    T = load_run(GRN)            # this run (3k+10k tf_only)
    union_cols = len(open(f'{GRN}/all_columns.txt').read().splitlines())

    L = []
    A = L.append
    A('# RNA-only baseline (3k + 10k multiome) — tf_only.txt regulators — GRNBoost2')
    A('')
    A('**Baseline:** no integration, no imputation — the PBMC multiome 3k Gene Expression '
      'and the PBMC multiome 10k Gene Expression matrices are stacked and handed directly '
      'to Arboreto.')
    A('**GRN:** Arboreto GRNBoost2; targets = `trrust_tf.txt` genes present in the matrix, '
      '**regulators = `tf_only.txt` TFs present in the matrix**.')
    A('**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated).')
    A('')
    A('---')
    A('')
    A('## 1. Inputs')
    A('')
    A(f'- **3k RNA:** PBMC multiome 3k granulocyte-sorted, Gene Expression '
      f'({sizes["3k genes (GRCh38-2020-A)"]} genes, GRCh38-2020-A), {sizes["3k cells"]} cells.')
    A(f'- **10k RNA:** PBMC multiome 10k granulocyte-sorted, Gene Expression '
      f'({sizes["10k genes (GRCh38-2020-A)"]} genes, GRCh38-2020-A), {sizes["10k cells"]} cells.')
    A(f'- **Stacking:** both matrices use the same Cell Ranger reference (GRCh38-2020-A) and '
      f'the identical feature list, so no reordering is needed. Duplicate symbols (10 per '
      f'reference) summed per reference; keep the **{sizes["common gene symbols"]}** common symbols. '
      f'Each matrix log1p(CPM)-normalized independently, then cells stacked '
      f'(3k rows first, then 10k).')
    A(f'- **All-cells matrix:** {X.shape[0]} cells x {X.shape[1]} genes.')
    A('')
    A('## 2. Regulator / target specification')
    A('')
    A('| Regulator list | `data/tf_only.txt` |')
    A('|---|---|')
    A(f'| Regulator list size | 829 |')
    A(f'| Regulators present in matrix | {T["regs"]} |')
    A(f'| Regulators appearing in the network | {T["net"]["TF"].nunique()} |')
    A(f'| Target genes (`trrust_tf.txt` present) | {T["tgts"]} |')
    A(f'| Expression subset handed to GRNBoost2 | {union_cols} columns (union of targets + tf_only regulators) |')
    A(f'| Cells | {T["n_cells"]} |')
    A(f'| Inferred edges | {len(T["net"]):,} |')
    A('')
    A('## 3. Evaluation vs ground truth (deduplicated)')
    A('')
    A('| Ground truth | Dedup. GT edges | Recovered |')
    A('|---|---|---|')
    for gt in ['PBMC-TRRUST', 'PBMC-Blood']:
        e = T['ev'][gt]
        A(f'| {gt} | {e["dedup"]} | {e["rec"]}/{e["dedup"]} ({e["recpct"]}%) |')
    A('')
    A('### Precision@K / Recall@K — PBMC-TRRUST')
    A('')
    A('| Top-K | Prec@K | Recall@K |')
    A('|---|---|---|')
    for (K, p, rc) in T['ev']['PBMC-TRRUST']['precrec']:
        A(f'| {K} | {p:.4f} | {rc:.4f} |')
    A('')
    A('### Precision@K / Recall@K — PBMC-Blood')
    A('')
    A('| Top-K | Prec@K | Recall@K |')
    A('|---|---|---|')
    for (K, p, rc) in T['ev']['PBMC-Blood']['precrec']:
        A(f'| {K} | {p:.4f} | {rc:.4f} |')
    A('')
    A('## 4. Top 25 edges')
    A('')
    A('| TF -> target | importance |')
    A('|---|---|')
    for _, r in T['net'].head(25).iterrows():
        A(f'| {r["TF"]} -> {r["target"]} | {r["importance"]:.4f} |')
    A('')
    A('## 5. Context: same tf_only regulators, other baselines & workflows')
    A('')
    A('All rows use the `tf_only.txt` regulator list (816 present). Cell count and target '
      'count differ per row, so this is context, not a controlled comparison.')
    A('')
    A('| Run | Cells | Targets | Inferred edges | PBMC-TRRUST recovered | PBMC-Blood recovered |')
    A('|---|---|---|---|---|---|')
    A(f'| RNA-only 3k+10k multiome (this run) | {T["n_cells"]} | {T["tgts"]} | {len(T["net"]):,} | '
      f'{T["ev"]["PBMC-TRRUST"]["rec"]}/{T["ev"]["PBMC-TRRUST"]["dedup"]} '
      f'({T["ev"]["PBMC-TRRUST"]["recpct"]}%) | '
      f'{T["ev"]["PBMC-Blood"]["rec"]}/{T["ev"]["PBMC-Blood"]["dedup"]} '
      f'({T["ev"]["PBMC-Blood"]["recpct"]}%) |')
    A('| RNA-only 3k+4k (tf_only, earlier) | 7051 | 2808 | 539629 | 2655/8751 (30.3%) | 3661/96846 (3.8%) |')
    A('')
    A('*3k+4k values from `reports/rna-baseline/REPORT_tfonly.md`. Full-10k integration + '
      'reverse-imputeKNN runs (expA/expB1/expB2, 29,218 cells) are in '
      '`reports/full_10k_3k_onlyTf/comparison_report_tfonly.md`; unpaired external-ATAC run in '
      '`reports/unpaired-10k/comparison_report_tfonly_unpaired.md`.*')
    A('')
    A('---')

    os.makedirs(REP, exist_ok=True)
    out = f'{REP}/REPORT_tfonly.md'
    with open(out, 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('Wrote', out)


if __name__ == '__main__':
    main()
