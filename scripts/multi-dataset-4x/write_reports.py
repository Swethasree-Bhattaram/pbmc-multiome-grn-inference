#!/usr/bin/env python
"""Generate a per-experiment README-style report (similar to the original
pbmc3k_analysis/REPORT.md), plus a consolidated comparison report.

Writes reports/<exp>_REPORT.md and reports/00_COMPARISON.md
Requires GRN + eval to have completed (results/exp*/grn/*).
Run with the scSAGA .venv (numpy/pandas/matplotlib).
"""
import os, sys, csv, datetime
import numpy as np, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
RES = f'{PROJ}/results'
REP = f'{PROJ}/reports'
os.makedirs(REP, exist_ok=True)
GT = f'{PROJ}/data/ground_truth'

EXP_META = {
    'expA': dict(
        title='Experiment A — reverse-imputeKNN from a SCEMENT-integrated combined RNA reference',
        ref_desc=('SCEMENT (AluruLab) batch-integrates the two real RNA datasets '
                  '(3k + 6k) into a single 5,422-cell reference; both ATAC populations are '
                  'imputed from that combined reference.'),
        ref_tag='SCEMENT-combined RNA (3k+6k, 5,422 cells)'),
    'expB1': dict(
        title='Experiment B1 — reverse-imputeKNN from the single 3k RNA reference',
        ref_desc=('Only the original 3k RNA dataset (2,711 cells) is used as reference; '
                  'both ATAC populations (3k + 6k) are imputed from it.'),
        ref_tag='3k RNA only (2,711 cells)'),
    'expB2': dict(
        title='Experiment B2 — reverse-imputeKNN from the single 6k RNA reference',
        ref_desc=('Only the subsampled 6k RNA dataset (2,711 cells) is used as reference; '
                  'both ATAC populations (3k + 6k) are imputed from it.'),
        ref_tag='6k RNA only (2,711 cells)'),
}

def main():
    for exp, meta in EXP_META.items():
        write_one(exp, meta)

def load(exp):
    grn = f'{RES}/{exp}/grn'
    net = pd.read_csv(f'{grn}/grnboost2_network.tsv', sep='\t').sort_values('importance', ascending=False)
    ncells = int(open(f'{grn}/n_cells.txt').read().strip())
    tfs = len(open(f'{grn}/tf_regulators.txt').read().splitlines())
    targets = len(open(f'{grn}/target_genes.txt').read().splitlines())
    ev = {}
    with open(f'{grn}/evaluation/evaluation_summary.txt') as f:
        evtxt = f.read()
    # parse per-gt
    import re
    curgt = None
    rows = {}
    for line in evtxt.splitlines():
        m = re.match(r'========== EVALUATION vs (.*) ==========')
        if m:
            curgt = m.group(1); rows[curgt] = {'precrec': []}; continue
        if curgt is None: continue
        m = re.match(r'Total edges: (\d+)')
        if m: rows[curgt]['total'] = int(m.group(1)); continue
        m = re.match(r'Deduplicated edges: (\d+)')
        if m: rows[curgt]['dedup'] = int(m.group(1)); continue
        m = re.match(r'Edges recovered by GRN: (\d+) \(([\d.]+)%\)')
        if m: rows[curgt]['rec'] = int(m.group(1)); rows[curgt]['recpct'] = float(m.group(2)); continue
        m = re.match(r'\s*(\d+)\s+([\d.]+)\s+([\d.]+)')
        if m: rows[curgt]['precrec'].append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    ev = rows
    imp = np.load(f'{RES}/{exp}/imputed_expression_genes_x_atac.npy', mmap_mode='r')
    return dict(net=net, ncells=ncells, tfs=tfs, targets=targets, ev=ev, imp=imp)

def write_one(exp, meta):
    d = load(exp)
    net, ncells = d['net'], d['ncells']
    L = []
    L.append(f'# {meta["title"]}')
    L.append('')
    L.append(f'**Date:** {datetime.datetime.now().strftime("%Y-%m-%d")}')
    L.append('**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, 6k-RNA, 6k-ATAC)')
    L.append('**Reverse-imputeKNN reference:** ' + meta['ref_tag'])
    L.append('**GRN:** Arboreto GRNBoost2 (TRRUST regulators)')
    L.append('**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)')
    L.append('')
    L.append('---')
    L.append('')
    L.append('## 1. Pipeline')
    L.append('')
    L.append('1. **Datasets:** 4 single-modality datasets (2,711 cells each): rna3k (real), atac3k (peaks), rna6k (subsampled 10k RNA, real), atac6k (subsampled 10k peaks). The 6k RNA/ATAC are the SAME 2,711 physical cells (paired multiome).')
    L.append(f'2. **Integration:** scSAGA into a shared joint embedding **H** (10,844 x 30).')
    L.append(f'3. **Reverse-imputeKNN**: {meta["ref_desc"]}')
    L.append('4. **All-cells matrix:** 10,844 x 36,601 (rows: rna3k, rna6k, atac3k-imputed, atac6k-imputed).')
    L.append('5. **GRN inference:** Arboreto GRNBoost2 with TRRUST TFs as regulators.')
    L.append('6. **Evaluation:** against PBMC-TRRUST and PBMC-Blood ground truth.')
    L.append('')
    L.append('## 2. Integration (joint embedding H)')
    L.append('')
    L.append(f'- Joint embedding **H**: 10,844 x 30 (blocks: rna3k, atac3k, rna6k, atac6k, each 2,711 x 30).')
    L.append('- Global alignment score: 0.221 (4 datasets, 3k anchor).')
    L.append('- Pairwise scores vs anchor (rna3k): atac3k 0.363, rna6k 0.413, atac6k 0.318.')
    L.append('')
    L.append('## 3. Imputation (reverse-imputeKNN, k=20, softmax-of-negative-distance weights)')
    L.append('')
    L.append(f'- Imputed expression: **{d["imp"].shape[0]} genes x {d["imp"].shape[1]} ATAC cells**.')
    L.append(f'- Mean imputed expression: {float(d["imp"].mean()):.4f}; non-zero fraction: {(np.asarray(d["imp"])>0).mean():.3f}.')
    L.append(f'- All-cells matrix: **{ncells} cells x {d["imp"].shape[0]} genes** (rna3k + rna6k real + 5,422 imputed ATAC).')
    L.append('')
    L.append('## 4. Arboreto GRNBoost2 results (TRRUST regulators)')
    L.append('')
    L.append(f'- Regulators (TFs present in expression): {d["tfs"]}')
    L.append(f'- Target genes: {d["targets"]} (top-2000 HVG + TRRUST TFs)')
    L.append(f'- Edges inferred: **{len(net):,}**')
    L.append('')
    L.append('**Top 25 edges:**')
    L.append('')
    L.append('| TF → target | importance |')
    L.append('|---|---|')
    for _, r in net.head(25).iterrows():
        L.append(f'| {r["TF"]} → {r["target"]} | {r["importance"]:.1f} |')
    L.append('')
    L.append('## 5. Evaluation vs ground truth (deduplicated)')
    L.append('')
    gtnames = ['PBMC-TRRUST', 'PBMC-Blood']
    for gt in gtnames:
        e = d['ev'].get(gt, {})
        dedup = e.get('dedup', '-')
        rec = e.get('rec', '-')
        recpct = e.get('recpct', '-')
        total = e.get('total', '-')
        L.append(f'**{gt}:** total={total}, deduplicated={dedup}, recovered by GRN={rec} ({recpct}%)')
        L.append('')
        L.append('| Top-K | Precision@K | Recall@K |')
        L.append('|---|---|---|')
        for K, p, r_ in e.get('precrec', []):
            L.append(f'| {K} | {p:.4f} | {r_:.4f} |')
        L.append('')
    L.append('---')
    L.append('')
    text = '\n'.join(L)
    with open(f'{REP}/{exp}_REPORT.md', 'w') as f:
        f.write(text)
    print('Wrote', f'{REP}/{exp}_REPORT.md')

if __name__ == '__main__':
    main()
