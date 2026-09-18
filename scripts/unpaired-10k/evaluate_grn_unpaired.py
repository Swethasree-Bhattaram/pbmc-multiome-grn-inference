#!/usr/bin/env python
"""Evaluate the UNPAIRED 10k GRNBoost2 network against PBMC ground truth.

Reads  results/unpaired_grn/grn_tfonly_reg/grnboost2_network.tsv
Writes results/unpaired_grn/grn_tfonly_reg/evaluation/{evaluation_summary.txt,top_edges.csv}

Ground truths (TF,TARGET edge lists), deduplicated before scoring:
  data/ground_truth/PBMC-TRRUST.csv
  data/ground_truth/PBMC-Blood.csv

Also emits a machine-readable evaluation.json so the report generator does not
have to re-parse text.
"""
import os, sys, csv, json
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
# Repo root: $UNPAIRED_ROOT, else two levels up (scripts/<workflow>/ -> repo).
# NOTE: one os.pardir would resolve to scripts/, which is a bug when the
# env var is unset (results/ and data/ then point inside scripts/).
PROJ = os.environ.get('UNPAIRED_ROOT', os.path.abspath(os.path.join(
    HERE, os.pardir, os.pardir)))
GRN = f'{PROJ}/results/unpaired_grn/grn_tfonly_reg'
OUT = f'{GRN}/evaluation'
GT = f'{PROJ}/data/ground_truth'


def main():
    os.makedirs(OUT, exist_ok=True)
    net = pd.read_csv(f'{GRN}/grnboost2_network.tsv', sep='\t')
    net = net.sort_values('importance', ascending=False).reset_index(drop=True)
    print('Inferred edges:', len(net))
    print('Columns:', list(net.columns))

    edge_set = set(zip(net['TF'], net['target']))
    tfs = set(net['TF'].unique())
    targets = set(net['target'].unique())
    print(f'GRN regulators: {len(tfs)}, targets: {len(targets)}')

    K_LIST = sorted(set([100, 500, 1000, 5000, 10000, len(net)]))
    results = {}

    def evaluate(gt_name, gt_edges, out):
        dedup = set(gt_edges)
        top_sets = {}
        out.write(f'========== EVALUATION vs {gt_name} ==========\n')
        out.write(f'Total edges: {len(gt_edges)}\n')
        out.write(f'Deduplicated edges: {len(dedup)}\n\n')
        out.write(f'{"Top-K":>12} {"Prec@K":>10} {"Recall@K":>10}\n')
        rows = {}
        for K in K_LIST:
            top_set = set(zip(net.head(K)['TF'], net.head(K)['target']))
            top_sets[K] = top_set
            prec = len(top_set & dedup) / K
            rec = len(top_set & dedup) / len(dedup) if dedup else 0.0
            rows[str(K)] = [round(prec, 4), round(rec, 4)]
            out.write(f'{K:>12} {prec:>10.4f} {rec:>10.4f}\n')
        recovered = len(edge_set & dedup)
        frac = recovered / len(dedup) if dedup else 0.0
        out.write(f'\nEdges recovered by GRN: {recovered} ({frac*100:.1f}%)\n')
        out.write('-' * 70 + '\n')
        return {'total': len(gt_edges), 'dedup': len(dedup), 'recovered': recovered,
                'edges': len(net), 'rows': rows}

    with open(f'{OUT}/evaluation_summary.txt', 'w') as out:
        out.write('ARBORETO GRNBOOST2 EVALUATION — UNPAIRED 10k (rna10k x atac10k_ext v1.1)\n')
        out.write('targets = trrust_tf.txt genes present (2827), regulators = tf_only.txt TFs present (816)\n')
        out.write('=' * 70 + '\n')
        out.write(f'Inferred edges: {len(net)}, regulators: {len(tfs)}, targets: {len(targets)}\n\n')
        for name, fn in [('PBMC-TRRUST', 'PBMC-TRRUST.csv'), ('PBMC-Blood', 'PBMC-Blood.csv')]:
            edges = []
            with open(f'{GT}/{fn}') as f:
                for x in csv.DictReader(f):
                    edges.append((x['TF'], x['TARGET']))
            results[name] = evaluate(name, edges, out)

    top25 = net.head(25)
    top25.to_csv(f'{OUT}/top_edges.csv', index=False)

    payload = {
        'n_edges': len(net), 'n_regulators': len(tfs), 'n_targets': len(targets),
        'n_cells': int(open(f'{GRN}/n_cells.txt').read().strip()),
        'ground_truths': results,
        'top25': [[r.TF, r.target, float(r.importance)] for r in top25.itertuples()],
    }
    with open(f'{OUT}/evaluation.json', 'w') as f:
        json.dump(payload, f, indent=2)

    print('\nTOP 25 EDGES:')
    print(top25.to_string(index=False))
    print('\n=== DONE ===')


if __name__ == '__main__':
    main()
