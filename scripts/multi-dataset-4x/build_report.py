#!/usr/bin/env python
"""Generate per-experiment and consolidated comparison reports + figures.

Reads results/exp{A,B1,B2}: all-cells matrices, integration metrics, GRN eval,
and renders README-style markdown reports into reports/.

Run with the scSAGA .venv after all GRN+eval completes.
"""
import os, sys, csv, json, datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
RES = f'{PROJ}/results'
REP = f'{PROJ}/reports'
os.makedirs(REP, exist_ok=True)
os.makedirs(f'{REP}/figures', exist_ok=True)

EXP_META = {
    'expA':  dict(label='Experiment A', title='Reverse-imputeKNN from a SCEMENT-integrated combined RNA reference',
                  ref='SCEMENT-integrated 3k+6k RNA (5,422 cells)', short='combined reference (SCEMENT)'),
    'expB1': dict(label='Experiment B1', title='Reverse-imputeKNN from the single 3k RNA reference',
                  ref='3k RNA only (2,711 cells)', short='3k RNA reference'),
    'expB2': dict(label='Experiment B2', title='Reverse-imputeKNN from the single 6k RNA reference',
                  ref='6k RNA only (2,711 cells)', short='6k RNA reference'),
}

def read_eval(exp):
    p = f'{RES}/{exp}/grn/evaluation/evaluation_summary.txt'
    if not os.path.exists(p): return None
    txt = open(p).read()
    return txt

def read_net(exp):
    p = f'{RES}/{exp}/grn/grnboost2_network.tsv'
    if not os.path.exists(p): return None
    import pandas as pd
    net = pd.read_csv(p, sep='\t').sort_values('importance', ascending=False)
    return net

def read_imputed(exp):
    imp = np.load(f'{RES}/{exp}/imputed_expression_genes_x_atac.npy', mmap_mode='r')
    return imp

def parse_eval(txt):
    """Return dict: gt -> {total, dedup, recovered, prec/rec by K}"""
    import re
    out = {}
    cur = None
    for line in txt.splitlines():
        m = re.match(r'EVALUATION vs (.*)', line)
        if m:
            cur = {}; out[m.group(1)] = cur; continue
        if cur is None: continue
        m = re.match(r'Total edges: (\d+)', line)
        if m: cur['total'] = int(m.group(1)); continue
        m = re.match(r'Deduplicated edges: (\d+)', line)
        if m: cur['dedup'] = int(m.group(1)); continue
        m = re.match(r'Edges recovered by GRN: (\d+) \(([\d.]+)%\)', line)
        if m: cur['recovered'] = int(m.group(1)); cur['recovered_pct'] = float(m.group(2)); continue
        m = re.match(r'\s*(\d+)\s+([\d.]+)\s+([\d.]+)', line)
        if m: cur.setdefault('precrec', {})[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    return out

def make_imputation_summary(exp):
    imp = read_imputed(exp)
    X = np.asarray(imp)
    return {
        'shape': list(X.shape),
        'mean': float(X.mean()),
        'nonzero_frac': float((X > 0).mean()),
        'median_nonzero': float(np.median(X[X > 0])) if (X > 0).any() else 0.0,
    }

def main():
    summary_rows = []
    for exp, meta in EXP_META.items():
        ev = parse_eval(read_eval(exp) or '')
        imp = make_imputation_summary(exp)
        net = read_net(exp)
        n_cells = None
        nc = f'{RES}/{exp}/grn/n_cells.txt'
        if os.path.exists(nc): n_cells = int(open(nc).read().strip())
        n_edges = len(net) if net is not None else None
        summary_rows.append({
            'exp': exp, 'label': meta['label'], 'title': meta['title'], 'ref': meta['ref'],
            'n_cells': n_cells, 'n_edges': n_edges, 'eval': ev, 'imp': imp,
        })

    # ---- Consolidated comparison table ----
    gt_names = ['PBMC-TRRUST', 'PBMC-Blood']
    lines = []
    lines.append('# Multi-Dataset scSAGA Integration & GRN: Experiment Comparison')
    lines.append('')
    lines.append(f'**Generated:** {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}')
    lines.append('')
    lines.append('## Design')
    lines.append('')
    lines.append('Four datasets are jointly integrated with scSAGA (3k-RNA anchor, 2,711 cells each):')
    lines.append('- **rna3k** (real expression) · **atac3k** (peaks) · **rna6k** (subsampled 10k RNA) · **atac6k** (subsampled 10k peaks)')
    lines.append('- Joint embedding **H** = 10,844 x 30 (each block 2,711 x 30); global alignment score 0.221.')
    lines.append('')
    lines.append('The only variable across experiments is **which RNA reference** is used by reverse-imputeKNN to')
    lines.append('impute expression for the 5,422 ATAC cells. The all-cells matrix (10,844 x 36,601) always holds the')
    lines.append('same real RNA cells (rna3k then rna6k) in rows 0..5421; only the imputed ATAC block (rows 5422..10843) changes.')
    lines.append('')
    lines.append('| Experiment | Reference for imputeKNN |')
    lines.append('|---|---|')
    lines.append('| **A** | SCEMENT-integrated 3k+RNA & 6k+RNA → one 5,422-cell reference |')
    lines.append('| **B1** | 3k RNA only (2,711 cells) |')
    lines.append('| **B2** | 6k RNA only (2,711 cells) |')
    lines.append('')
    lines.append('Each produces an Arboreto GRNBoost2 network evaluated against PBMC-TRRUST and PBMC-Blood ground truth.')
    lines.append('')
    lines.append('## GRN overview')
    lines.append('')
    lines.append('| Experiment | All-cells matrix | Cells | Target genes | Inferred edges |')
    lines.append('|---|---|---|---|---|')
    for r in summary_rows:
        lines.append(f'| {r["exp"]} | {r["imp"]["shape"][0]} x {r["imp"]["shape"][1]} | {r["n_cells"]} | — | {r["n_edges"]} |')
    lines.append('')

    # ---- Edges recovered ----
    lines.append('## Edges recovered by GRN (fraction of deduplicated ground truth)')
    lines.append('')
    lines.append('| Experiment | PBMC-TRRUST | PBMC-Blood |')
    lines.append('|---|---|---|')
    for r in summary_rows:
        t = r['eval'].get('PBMC-TRRUST', {}).get('recovered', '-')
        b = r['eval'].get('PBMC-Blood', {}).get('recovered', '-')
        lines.append(f'| {r["exp"]} | {t} | {b} |')

    # ---- Precision@K ----
    for K in [100, 1000, 5000]:
        lines.append(f'')
        lines.append(f'## Precision@K (K={K})')
        lines.append('')
        lines.append('| Experiment | PBMC-TRRUST | PBMC-Blood |')
        lines.append('|---|---|---|')
        for r in summary_rows:
            t = r['eval'].get('PBMC-TRRUST', {}).get('precrec', {}).get(K, ('-','-'))[0]
            b = r['eval'].get('PBMC-Blood', {}).get('precrec', {}).get(K, ('-','-'))[0]
            t = f'{t:.4f}' if isinstance(t, float) else t
            b = f'{b:.4f}' if isinstance(b, float) else b
            lines.append(f'| {r["exp"]} | {t} | {b} |')

    # ---- Imputation summary ----
    lines.append('')
    lines.append('## Reverse-imputed ATAC expression (5,422 cells x 36,601 genes)')
    lines.append('')
    lines.append('| Experiment | Mean expr | Non-zero fraction |')
    lines.append('|---|---|---|')
    for r in summary_rows:
        lines.append(f'| {r["exp"]} | {r["imp"]["mean"]:.4f} | {r["imp"]["nonzero_frac"]:.3f} |')
    lines.append('')

    # ---- Per-experiment detail ----
    lines.append('---')
    lines.append('')
    for r in summary_rows:
        lines.append(f'## {r["label"]}: {r["title"]}')
        lines.append('')
        lines.append(f'**Reference:** {r["ref"]}')
        lines.append('')
        if r['eval']:
            lines.append('### Evaluation vs ground truth (deduplicated)')
            lines.append('')
            for gt in gt_names:
                e = r['eval'].get(gt, {})
                lines.append(f'**{gt}:** total={e.get("total","-")}, dedup={e.get("dedup","-")}, recovered={e.get("recovered","-")} ({e.get("recovered_pct","-")}%)')
            lines.append('')
            lines.append('| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |')
            lines.append('|---|---|---|')
            Ks = sorted(set(list(r['eval'].get('PBMC-TRRUST', {}).get('precrec', {}).keys()) +
                             list(r['eval'].get('PBMC-Blood', {}).get('precrec', {}).keys())))
            for K in Ks:
                t = r['eval'].get('PBMC-TRRUST', {}).get('precrec', {}).get(K, ('-','-'))
                b = r['eval'].get('PBMC-Blood', {}).get('precrec', {}).get(K, ('-','-'))
                lines.append(f'| {K} | {t[0] if isinstance(t[0],float) else t[0]:.4f}/{t[1]} | {b[0] if isinstance(b[0],float) else b[0]:.4f}/{b[1]} |')
            lines.append('')
        lines.append('---')
        lines.append('')

    report = '\n'.join(lines)
    with open(f'{REP}/comparison_report.md', 'w') as f:
        f.write(report)
    print('Wrote', f'{REP}/comparison_report.md')

    # figures
    make_figures(summary_rows, gt_names)

def make_figures(rows, gt_names):
    # Bar: edges recovered fraction per exp per GT
    exps = [r['exp'] for r in rows]
    for gt in gt_names:
        vals = [r['eval'].get(gt, {}).get('recovered', 0) for r in rows]
        plt.figure(figsize=(7, 4))
        plt.bar(exps, vals, color=['#c0392b', '#2980b9', '#27ae60'])
        for i, v in enumerate(vals):
            plt.text(i, v, str(v), ha='center', va='bottom')
        plt.title(f'Edges recovered by GRN — {gt}')
        plt.ylabel('deduplicated GT edges found')
        plt.tight_layout()
        plt.savefig(f'{REP}/figures/recovered_{gt.replace("-","_")}.png', dpi=150)
        plt.close()
    # Imputation density
    plt.figure(figsize=(7, 4))
    means = [r['imp']['mean'] for r in rows]
    plt.bar(exps, means, color=['#c0392b', '#2980b9', '#27ae60'])
    plt.title('Mean imputed ATAC expression')
    plt.ylabel('mean log1p(CPM)')
    plt.tight_layout()
    plt.savefig(f'{REP}/figures/imputed_mean.png', dpi=150)
    plt.close()
    print('Wrote figures to', f'{REP}/figures/')

if __name__ == '__main__':
    main()
