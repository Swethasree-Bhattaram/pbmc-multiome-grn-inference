#!/usr/bin/env python
"""Evaluate a GRNBoost2 network against the PBMC ground-truth edge lists.

Ground truths (TF,TARGET edge lists, deduplicated before scoring):
  data/ground_truth/PBMC-TRRUST.csv   (8,751 dedup edges)
  data/ground_truth/PBMC-Blood.csv    (96,846 dedup edges)

Writes $PACE_ROOT/results/<exp>_grn/grn_tfonly_reg/evaluation/:
  evaluation_summary.txt   edges, Precision@K / Recall@K table, recovered fractions
  top_edges.csv            top 25 inferred edges
  evaluation.json          machine-readable summary

Usage: python 07_evaluate_grn.py --experiment unpaired
"""
import argparse
import csv
import json
import os
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--experiment', choices=['paired', 'unpaired'], required=True)
    args = ap.parse_args()
    exp = args.experiment

    proj = os.environ.get('PACE_ROOT', os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    grn = f'{proj}/results/{exp}_grn/grn_tfonly_reg'
    out = f'{grn}/evaluation'
    gt = f'{proj}/data/ground_truth'
    os.makedirs(out, exist_ok=True)

    net_path = f'{grn}/grnboost2_network.tsv'
    if not os.path.exists(net_path):
        raise SystemExit(f'missing {net_path}; run 06_run_grn.py first')
    net = pd.read_csv(net_path, sep='\t')
    net = net.sort_values('importance', ascending=False).reset_index(drop=True)
    print('inferred edges:', len(net), flush=True)

    edge_set = set(zip(net['TF'], net['target']))
    tfs, targets = set(net['TF'].unique()), set(net['target'].unique())

    K_LIST = sorted(set([100, 500, 1000, 5000, 10000, len(net)]))
    results = {}

    def evaluate(name, gt_edges, fh):
        dedup = set(gt_edges)
        fh.write(f'========== EVALUATION vs {name} ==========\n')
        fh.write(f'Total edges: {len(gt_edges)}\n')
        fh.write(f'Deduplicated edges: {len(dedup)}\n\n')
        fh.write(f'{"Top-K":>12} {"Prec@K":>10} {"Recall@K":>10}\n')
        rows = {}
        for K in K_LIST:
            top = set(zip(net.head(K)['TF'], net.head(K)['target']))
            prec = len(top & dedup) / K
            rec = len(top & dedup) / len(dedup) if dedup else 0.0
            rows[str(K)] = [round(prec, 4), round(rec, 4)]
            fh.write(f'{K:>12} {prec:>10.4f} {rec:>10.4f}\n')
        recovered = len(edge_set & dedup)
        frac = recovered / len(dedup) if dedup else 0.0
        fh.write(f'\nEdges recovered by GRN: {recovered} ({frac*100:.1f}%)\n')
        fh.write('-' * 70 + '\n')
        return {'total': len(gt_edges), 'dedup': len(dedup), 'recovered': recovered,
                'edges': len(net), 'pct': round(frac * 100, 1), 'rows': rows}

    regs = int(open(f'{grn}/tf_regulators.txt').read().count('\n')) if \
        os.path.exists(f'{grn}/tf_regulators.txt') else len(tfs)
    with open(f'{out}/evaluation_summary.txt', 'w') as fh:
        fh.write(f'ARBORETO GRNBOOST2 EVALUATION -- {exp}\n')
        fh.write(f'regulators = tf_only.txt TFs present ({regs}), '
                 f'targets = trrust_tf.txt genes present\n')
        fh.write('=' * 70 + '\n')
        fh.write(f'Inferred edges: {len(net)}, regulators: {len(tfs)}, '
                 f'targets: {len(targets)}\n\n')
        for name, fn in [('PBMC-TRRUST', 'PBMC-TRRUST.csv'),
                         ('PBMC-Blood', 'PBMC-Blood.csv')]:
            p = f'{gt}/{fn}'
            if not os.path.exists(p):
                print(f'skipping {name}: {p} not found', flush=True)
                continue
            edges = []
            with open(p) as f2:
                for x in csv.DictReader(f2):
                    edges.append((x['TF'], x['TARGET']))
            results[name] = evaluate(name, edges, fh)

    top25 = net.head(25)
    top25.to_csv(f'{out}/top_edges.csv', index=False)

    payload = {
        'experiment': exp,
        'n_edges': len(net),
        'n_regulators': len(tfs),
        'n_targets': len(targets),
        'n_cells': int(open(f'{grn}/n_cells.txt').read().strip()),
        'ground_truths': results,
        'top25': [[r.TF, r.target, float(r.importance)] for r in top25.itertuples()],
    }
    with open(f'{out}/evaluation.json', 'w') as f3:
        json.dump(payload, f3, indent=2)

    print('\nTOP 25 EDGES:', flush=True)
    print(top25.to_string(index=False), flush=True)
    print('\n=== DONE ===', flush=True)


if __name__ == '__main__':
    sys.exit(main())
