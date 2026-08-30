#!/usr/bin/env python
"""Evaluate the Arboreto GRNBoost2 network against PBMC ground-truth regulatory
edges (PBMC-TRRUST.csv and PBMC-Blood.csv).

For each ground-truth set we report retrieval metrics over the ranked edge list
(sorted by importance descending):
  - Total edges      : number of edge rows in the ground-truth file
  - Deduplicated edges : unique TF->target pairs (repeated pairs removed)
  - Precision@K      : fraction of the top-K inferred edges that are in the
                       deduplicated ground-truth set
  - Recall@K         : fraction of all deduplicated ground-truth edges found in
                       the top-K
  - Edges recovered by GRN : number (and fraction) of deduplicated ground-truth
                       edges that appear anywhere in the inferred network

No filtering beyond deduplication is applied to the ground truth.

Run with the scSAGA .venv (numpy/pandas).
"""
import pandas as pd, os, csv

BASE = '/Users/sbhattaram/pbmc3k_analysis'
GRN  = f'{BASE}/results/grn'
GT   = f'{BASE}/data/ground_truth'
OUT  = f'{BASE}/results/grn/evaluation'
os.makedirs(OUT, exist_ok=True)

# --- Load inferred network, ranked by importance descending ---
net = pd.read_csv(f'{GRN}/grnboost2_network.tsv', sep='\t')
net = net.sort_values('importance', ascending=False).reset_index(drop=True)
print(f'Inferred edges: {len(net)}')

edge_set = set(zip(net['TF'], net['target']))
tfs = set(net['TF'].unique())
targets = set(net['target'].unique())
print(f'GRN regulators (TFs): {len(tfs)}, targets: {len(targets)}')


def evaluate(gt_name, gt_edges, out):
    print(f'\n========== EVALUATION vs {gt_name} ==========')
    total = len(gt_edges)               # raw rows in ground-truth file
    dedup = set(gt_edges)               # unique TF->target pairs
    print(f'Total edges: {total}')
    print(f'Deduplicated edges: {len(dedup)}')

    topN = sorted(set([100, 500, 1000, 5000, 10000, len(net)]))

    out.write(f'========== EVALUATION vs {gt_name} ==========\n')
    out.write(f'Total edges: {total}\n')
    out.write(f'Deduplicated edges: {len(dedup)}\n\n')
    out.write(f'{"Top-K":>12} {"Prec@K":>10} {"Recall@K":>10}\n')
    print(f'{"Top-K":>12} {"Prec@K":>10} {"Recall@K":>10}')

    for K in topN:
        top_set = set(zip(net.head(K)['TF'], net.head(K)['target']))
        prec = len(top_set & dedup) / K
        rec = len(top_set & dedup) / len(dedup) if dedup else 0
        out.write(f'{K:>12} {prec:>10.4f} {rec:>10.4f}\n')
        print(f'{K:>12} {prec:>10.4f} {rec:>10.4f}')

    # Number of deduplicated ground-truth edges that appear anywhere in the network
    recovered = len(edge_set & dedup)
    frac_recovered = recovered / len(dedup) if dedup else 0
    out.write(f'\nEdges recovered by GRN: {recovered} ({frac_recovered*100:.1f}%)\n')
    out.write('-'*70 + '\n')
    print(f'Edges recovered by GRN: {recovered} ({frac_recovered*100:.1f}%)')
    return {'name': gt_name, 'total': total, 'dedup': len(dedup),
            'recovered': recovered, 'edges': len(net)}


summary = {}
with open(f'{OUT}/evaluation_summary.txt', 'w') as out:
    out.write('ARBORETO GRNBOOST2 NETWORK EVALUATION SUMMARY\n')
    out.write('='*70 + '\n')
    out.write(f'Inferred edges: {len(net)}, regulators (TFs): {len(tfs)}, targets: {len(targets)}\n\n')

    gt1 = []
    with open(f'{GT}/PBMC-TRRUST.csv') as f:
        for x in csv.DictReader(f):
            gt1.append((x['TF'], x['TARGET']))
    summary['PBMC-TRRUST'] = evaluate('PBMC-TRRUST', gt1, out)

    gt2 = []
    with open(f'{GT}/PBMC-Blood.csv') as f:
        for x in csv.DictReader(f):
            gt2.append((x['TF'], x['TARGET']))
    summary['PBMC-Blood'] = evaluate('PBMC-Blood', gt2, out)

print('\n=== DONE ===')
print(summary)
