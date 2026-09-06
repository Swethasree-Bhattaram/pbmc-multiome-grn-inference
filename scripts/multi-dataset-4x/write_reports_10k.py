#!/usr/bin/env python
"""Generate per-experiment reports for the full-10k TRRUST-only GRN runs.

Reads results/exp{A,B1,B2}_10k/grn_trrust/ and writes reports/exp{EXP}_REPORT_10k.md.
Run with the scSAGA .venv.
"""
import os, sys, csv, re
import numpy as np, pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir)))
RES = f'{PROJ}/results'
REP = f'{PROJ}/reports'
os.makedirs(REP, exist_ok=True)

EXP_META = {
    'expA':  dict(label='Experiment A', title='reverse-imputeKNN from a SCEMENT-integrated combined RNA reference',
                  ref='SCEMENT-combined RNA (3k+10k, 14,609 cells)',
                  datasets='rna3k, atac3k, rna10k, atac10k'),
    'expB1': dict(label='Experiment B1', title='reverse-imputeKNN from the single 3k RNA reference',
                  ref='3k RNA only (2,711 cells)',
                  datasets='rna3k, atac3k, rna10k, atac10k'),
    'expB2': dict(label='Experiment B2', title='reverse-imputeKNN from the single 10k RNA reference',
                  ref='10k RNA only (11,898 cells)',
                  datasets='rna3k, atac3k, rna10k, atac10k'),
}

def parse_eval(txt):
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
    for exp, meta in EXP_META.items():
        grn = f'{RES}/{exp}_10k/grn_trrust'
        net = pd.read_csv(f'{grn}/grnboost2_network.tsv', sep='\t').sort_values('importance', ascending=False)
        ev = parse_eval(open(f'{grn}/evaluation/evaluation_summary.txt').read())
        ncells = int(open(f'{grn}/n_cells.txt').read().strip())
        ntargets = len([l for l in open(f'{grn}/target_genes.txt') if l.strip()])
        nregs = len([l for l in open(f'{grn}/tf_regulators.txt') if l.strip()])

        L = []
        L.append(f'# {meta["label"]} — {meta["title"]} (full-10k)')
        L.append('')
        L.append('**Integration:** scSAGA — 4 datasets jointly embedded (3k-RNA anchor, 3k-ATAC, 10k-RNA, 10k-ATAC)')
        L.append(f'**Reverse-imputeKNN reference:** {meta["ref"]}')
        L.append('**GRN:** Arboreto GRNBoost2 (TRRUST regulators, TRRUST-only targets)')
        L.append('**Evaluation:** PBMC-TRRUST, PBMC-Blood (deduplicated)')
        L.append('')
        L.append('---')
        L.append('')
        L.append('## 1. Pipeline')
        L.append('')
        L.append(f'1. **Datasets:** 4 single-modality datasets: rna3k (2,711), atac3k (2,711), rna10k (11,898), atac10k (11,898). The 10k RNA/ATAC are the SAME 11,898 physical cells (paired multiome).')
        L.append('2. **Integration:** scSAGA into a shared joint embedding **H** (29,218 x 30).')
        L.append(f'3. **Reverse-imputeKNN**: {meta["title"]}.')
        L.append('4. **All-cells matrix:** 29,218 x 36,601 (rows: rna3k, rna10k, atac3k-imputed, atac10k-imputed).')
        L.append('5. **GRN inference:** Arboreto GRNBoost2 with TRRUST TFs as regulators, TRRUST-only targets.')
        L.append('6. **Evaluation:** against PBMC-TRRUST and PBMC-Blood ground truth.')
        L.append('')
        L.append('## 2. Integration (joint embedding H)')
        L.append('')
        L.append('- Joint embedding **H**: 29,218 x 30 (blocks: rna3k, atac3k, rna10k, atac10k).')
        L.append('- Global alignment score: 0.0984 (4 datasets, 3k anchor).')
        L.append('- Pairwise scores vs anchor (rna3k): atac3k 0.4027, rna10k 0.2240, atac10k 0.1236.')
        L.append('')
        L.append('## 3. Imputation (reverse-imputeKNN, k=20, softmax-of-negative-distance weights)')
        L.append('')
        L.append('- Imputed expression: **36601 genes x 14609 ATAC cells**.')
        L.append('- All-cells matrix: **29218 cells x 36601 genes** (rna3k + rna10k real + 14,609 imputed ATAC).')
        L.append('')
        L.append('## 4. Arboreto GRNBoost2 results (TRRUST regulators, TRRUST-only targets)')
        L.append('')
        L.append(f'- Regulators (TFs present in expression): {nregs}')
        L.append(f'- Target genes: {ntargets} (TRRUST-only)')
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
        for gt in ['PBMC-TRRUST', 'PBMC-Blood']:
            e = ev.get(gt, {})
            L.append(f'**{gt}:** total={e.get("total","-")}, deduplicated={e.get("dedup","-")}, recovered by GRN={e.get("recovered","-")} ({e.get("recovered_pct","-")}%)')
            L.append('')
            L.append('| Top-K | Precision@K | Recall@K |')
            L.append('|---|---|---|')
            for K in sorted(e.get('precrec', {}).keys()):
                p, r = e['precrec'][K]
                L.append(f'| {K} | {p:.4f} | {r:.4f} |')
            L.append('')
        L.append('---')
        L.append('')

        report = '\n'.join(L)
        with open(f'{REP}/{exp}_REPORT_10k.md', 'w') as f:
            f.write(report)
        print('Wrote', f'{REP}/{exp}_REPORT_10k.md')

if __name__ == '__main__':
    main()
