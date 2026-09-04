#!/usr/bin/env python
"""Generate a consolidated comparison report for the TRRUST-only GRN runs.

Reads results/exp{A,B1,B2}/grn_trrust/ (networks + evaluations) and writes
reports/comparison_report_trrust.md. Reuses the same metrics/format as the
original comparison report but for the TRRUST-only target-gene runs.

Run with the scSAGA .venv.
"""
import os, sys, csv, datetime
import numpy as np, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
RES = f'{PROJ}/results'
REP = f'{PROJ}/reports'
os.makedirs(REP, exist_ok=True)

EXP_META = {
    'expA':  dict(label='Experiment A', title='Reverse-imputeKNN from a SCEMENT-integrated combined RNA reference',
                  ref='SCEMENT-integrated 3k+subsampled_3k RNA (5,422 cells)'),
    'expB1': dict(label='Experiment B1', title='Reverse-imputeKNN from the single 3k RNA reference',
                  ref='3k RNA only (2,711 cells)'),
    'expB2': dict(label='Experiment B2', title='Reverse-imputeKNN from the single subsampled_3k RNA reference',
                  ref='subsampled_3k RNA only (2,711 cells)'),
}

def parse_eval(txt):
    import re
    out = {}
    cur = None
    for line in txt.splitlines():
        m = re.search(r'EVALUATION vs (.*?) =+', line)
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

def main():
    gt_names = ['PBMC-TRRUST', 'PBMC-Blood']
    rows = []
    for exp, meta in EXP_META.items():
        grn = f'{RES}/{exp}/grn_trrust'
        net = pd.read_csv(f'{grn}/grnboost2_network.tsv', sep='\t').sort_values('importance', ascending=False)
        ev = parse_eval(open(f'{grn}/evaluation/evaluation_summary.txt').read())
        ncells = int(open(f'{grn}/n_cells.txt').read().strip())
        ntargets = len([l for l in open(f'{grn}/target_genes.txt') if l.strip()])
        rows.append({'exp': exp, 'label': meta['label'], 'title': meta['title'], 'ref': meta['ref'],
                     'n_edges': len(net), 'n_cells': ncells, 'n_targets': ntargets, 'eval': ev})

    L = []
    L.append('# Multi-Dataset scSAGA Integration & GRN: TRRUST-only Target Comparison')
    L.append('')
    L.append('## Design')
    L.append('')
    L.append('Same four-dataset scSAGA integration and reverse-imputeKNN experiments as the main')
    L.append('comparison, but the Arboreto GRNBoost2 **target genes are restricted to ONLY the TRRUST')
    L.append('TFs present in the expression matrix** (2,827 genes) instead of top-2000 HVG + TRRUST TFs')
    L.append('(~4,389 genes). This makes the GRN inference much faster and focuses the network on')
    L.append('TF-to-TF regulatory edges.')
    L.append('')
    L.append('| Experiment | Reference for imputeKNN |')
    L.append('|---|---|')
    for r in rows:
        L.append(f'| **{r["exp"].replace("exp","")}** | {r["ref"]} |')
    L.append('')
    L.append('## GRN overview')
    L.append('')
    L.append('| Experiment | All-cells matrix | Cells | Target genes | Inferred edges |')
    L.append('|---|---|---|---|---|')
    for r in rows:
        L.append(f'| {r["exp"]} | {r["n_cells"]} x 36601 | {r["n_cells"]} | {r["n_targets"]} | {r["n_edges"]} |')
    L.append('')
    L.append('## Edges recovered by GRN (fraction of deduplicated ground truth)')
    L.append('')
    L.append('| Experiment | PBMC-TRRUST | PBMC-Blood |')
    L.append('|---|---|---|')
    for r in rows:
        t = r['eval'].get('PBMC-TRRUST', {})
        b = r['eval'].get('PBMC-Blood', {})
        t_str = f"{t.get('recovered','-')}/{t.get('dedup','-')}" if t else '-'
        b_str = f"{b.get('recovered','-')}/{b.get('dedup','-')}" if b else '-'
        L.append(f'| {r["exp"]} | {t_str} | {b_str} |')
    L.append('')
    L.append('*Format: recovered edges / deduplicated ground-truth edges.*')
    L.append('')
    L.append('## Precision@K')
    L.append('')
    for K in [100, 1000, 5000]:
        L.append(f'### K={K}')
        L.append('')
        L.append('| Experiment | PBMC-TRRUST | PBMC-Blood |')
        L.append('|---|---|---|')
        for r in rows:
            t = r['eval'].get('PBMC-TRRUST', {}).get('precrec', {}).get(K, ('-','-'))[0]
            b = r['eval'].get('PBMC-Blood', {}).get('precrec', {}).get(K, ('-','-'))[0]
            t = f'{t:.4f}' if isinstance(t, float) else t
            b = f'{b:.4f}' if isinstance(b, float) else b
            L.append(f'| {r["exp"]} | {t} | {b} |')
        L.append('')
    L.append('---')
    L.append('')
    for r in rows:
        L.append(f'## {r["label"]}: {r["title"]}')
        L.append('')
        L.append(f'**Reference:** {r["ref"]}')
        L.append('')
        for gt in gt_names:
            e = r['eval'].get(gt, {})
            L.append(f'**{gt}:** total={e.get("total","-")}, dedup={e.get("dedup","-")}, recovered={e.get("recovered","-")} ({e.get("recovered_pct","-")}%)')
        L.append('')
        L.append('| Top-K | TRRUST Prec/Rec | Blood Prec/Rec |')
        L.append('|---|---|---|')
        Ks = sorted(set(list(r['eval'].get('PBMC-TRRUST', {}).get('precrec', {}).keys()) +
                         list(r['eval'].get('PBMC-Blood', {}).get('precrec', {}).keys())))
        for K in Ks:
            t = r['eval'].get('PBMC-TRRUST', {}).get('precrec', {}).get(K, ('-','-'))
            b = r['eval'].get('PBMC-Blood', {}).get('precrec', {}).get(K, ('-','-'))
            L.append(f'| {K} | {t[0] if isinstance(t[0],float) else t[0]:.4f}/{t[1]} | {b[0] if isinstance(b[0],float) else b[0]:.4f}/{b[1]} |')
        L.append('')
        L.append('---')
        L.append('')

    report = '\n'.join(L)
    with open(f'{REP}/comparison_report_trrust.md', 'w') as f:
        f.write(report)
    print('Wrote', f'{REP}/comparison_report_trrust.md')

if __name__ == '__main__':
    main()
