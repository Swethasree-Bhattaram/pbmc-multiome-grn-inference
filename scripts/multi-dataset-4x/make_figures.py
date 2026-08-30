#!/usr/bin/env python
"""Generate diagnostic figures for the multi-dataset experiments.

Reads results/exp*/imputed and grn evaluations, writes figures to reports/figures.
Run with the scSAGA .venv.
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJ = os.environ.get('PBSC4K_ROOT', os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)))
RES = f'{PROJ}/results'
REP = f'{PROJ}/reports/figures'
os.makedirs(REP, exist_ok=True)

EXPS = ['expA', 'expB1', 'expB2']
COLORS = {'expA': '#c0392b', 'expB1': '#2980b9', 'expB2': '#27ae60'}

def load_eval(exp):
    p = f'{RES}/{exp}/grn/evaluation/evaluation_summary.txt'
    d = {}
    if not os.path.exists(p): return d
    import re
    cur = None
    for line in open(p):
        m = re.match(r'EVALUATION vs (.*)', line)
        if m: cur = {}; d[m.group(1)] = cur; continue
        if cur is None: continue
        m = re.match(r'Edges recovered by GRN: (\d+) \(([\d.]+)%\)', line)
        if m: cur['recovered'] = int(m.group(1)); cur['pct'] = float(m.group(2))
    return d

def main():
    # Fig 1: Edges recovered per GT per experiment
    gts = ['PBMC-TRRUST', 'PBMC-Blood']
    for gt in gts:
        vals = [load_eval(e).get(gt, {}).get('recovered', 0) for e in EXPS]
        plt.figure(figsize=(6.5, 4))
        bars = plt.bar(EXPS, vals, color=[COLORS[e] for e in EXPS])
        for b, v in zip(bars, vals):
            plt.text(b.get_x()+b.get_width()/2, v, str(v), ha='center', va='bottom')
        plt.title(f'GRN edges recovered — {gt}')
        plt.ylabel('deduplicated ground-truth edges found')
        plt.tight_layout(); plt.savefig(f'{REP}/recovered_{gt}.png', dpi=150); plt.close()

    # Fig 2: mean imputed expression density
    means = []
    for e in EXPS:
        imp = np.load(f'{RES}/{e}/imputed_expression_genes_x_atac.npy', mmap_mode='r')
        means.append(float(np.asarray(imp).mean()))
    plt.figure(figsize=(6.5, 4))
    plt.bar(EXPS, means, color=[COLORS[e] for e in EXPS])
    plt.title('Mean reverse-imputed ATAC expression')
    plt.ylabel('mean log1p(CPM)')
    plt.tight_layout(); plt.savefig(f'{REP}/imputed_mean.png', dpi=150); plt.close()

    # Fig 3: joint embedding 2D (from integration)
    import scipy.io
    j2 = f'{RES}/integration/joint_embedding_2d.csv'
    if os.path.exists(j2):
        df = pd.read_csv(j2)
        plt.figure(figsize=(7, 6))
        for ds, clr in [('rna3k','#2c3e50'),('atac3k','#c0392b'),('rna6k','#27ae60'),('atac6k','#8e44ad')]:
            sub = df[df.dataset == ds]
            plt.scatter(sub.x, sub.y, s=3, label=ds, color=clr, alpha=0.6)
        plt.legend(markerscale=3)
        plt.title('scSAGA 4-dataset joint embedding (2D)')
        plt.tight_layout(); plt.savefig(f'{REP}/joint_embedding_4datasets.png', dpi=150); plt.close()

    print('Wrote figures to', REP)
    for f in sorted(os.listdir(REP)):
        print(' -', f)

if __name__ == '__main__':
    main()
