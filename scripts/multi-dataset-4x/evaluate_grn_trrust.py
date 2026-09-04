#!/usr/bin/env python
"""Evaluate a TRRUST-only Arboreto GRNBoost2 network against PBMC ground truth.

Same metrics as evaluate_grn.py but reads from results/<exp>/grn_trrust/ and
writes evaluation to results/<exp>/grn_trrust/evaluation/.

Run with the scSAGA .venv (numpy/pandas).
Usage: python evaluate_grn_trrust.py <exp>
"""
import sys, os, csv
import pandas as pd

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
GT = f'{PROJ}/data/ground_truth'

def main():
    exp = sys.argv[1]
    GRN = f'{PROJ}/results/{exp}/grn_trrust'
    OUT = f'{PROJ}/results/{exp}/grn_trrust/evaluation'
    os.makedirs(OUT, exist_ok=True)

    net = pd.read_csv(f'{GRN}/grnboost2_network.tsv', sep='\t')
    net = net.sort_values('importance', ascending=False).reset_index(drop=True)
    print(f'[exp {exp}] Inferred edges: {len(net)}')

    edge_set = set(zip(net['TF'], net['target']))
    tfs = set(net['TF'].unique()); targets = set(net['target'].unique())
    print(f'GRN regulators (TFs): {len(tfs)}, targets: {len(targets)}')

    def evaluate(gt_name, gt_edges, out):
        total = len(gt_edges)
        dedup = set(gt_edges)
        topN = sorted(set([100, 500, 1000, 5000, 10000, len(net)]))
        out.write(f'========== EVALUATION vs {gt_name} ==========\n')
        out.write(f'Total edges: {total}\n')
        out.write(f'Deduplicated edges: {len(dedup)}\n\n')
        out.write(f'{"Top-K":>12} {"Prec@K":>10} {"Recall@K":>10}\n')
        rows = {}
        for K in topN:
            top_set = set(zip(net.head(K)['TF'], net.head(K)['target']))
            prec = len(top_set & dedup) / K
            rec = len(top_set & dedup) / len(dedup) if dedup else 0
            rows[K] = (prec, rec)
            out.write(f'{K:>12} {prec:>10.4f} {rec:>10.4f}\n')
        recovered = len(edge_set & dedup)
        frac = recovered / len(dedup) if dedup else 0
        out.write(f'\nEdges recovered by GRN: {recovered} ({frac*100:.1f}%)\n')
        out.write('-' * 70 + '\n')
        return {'total': total, 'dedup': len(dedup), 'recovered': recovered, 'edges': len(net), 'rows': rows}

    with open(f'{OUT}/evaluation_summary.txt', 'w') as out:
        out.write(f'ARBORETO GRNBOOST2 NETWORK EVALUATION SUMMARY (exp {exp}, TRRUST-only targets)\n')
        out.write('=' * 70 + '\n')
        out.write(f'Inferred edges: {len(net)}, regulators: {len(tfs)}, targets: {len(targets)}\n\n')
        summary = {}
        gt1 = []
        with open(f'{GT}/PBMC-TRRUST.csv') as f:
            for x in csv.DictReader(f): gt1.append((x['TF'], x['TARGET']))
        summary['PBMC-TRRUST'] = evaluate('PBMC-TRRUST', gt1, out)
        gt2 = []
        with open(f'{GT}/PBMC-Blood.csv') as f:
            for x in csv.DictReader(f): gt2.append((x['TF'], x['TARGET']))
        summary['PBMC-Blood'] = evaluate('PBMC-Blood', gt2, out)

    top25 = net.head(25)
    with open(f'{OUT}/top_edges.csv', 'w') as f:
        top25.to_csv(f, index=False)
    print('\nTOP 25 EDGES:')
    print(top25.head(25).to_string(index=False))
    print('\n=== DONE ===')
    print(summary)

if __name__ == '__main__':
    main()
