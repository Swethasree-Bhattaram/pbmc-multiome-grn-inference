#!/usr/bin/env python
"""Run the multimodal GRN pipeline from a config file.

    python pipeline/run.py --list                     # show experiments
    python pipeline/run.py --check                    # validate inputs only
    python pipeline/run.py --experiment ab            # integration+impute+GRN+eval
    python pipeline/run.py --experiment ab --stage integrate
    python pipeline/run.py --experiment ab --stage grn
    python pipeline/run.py --all                      # every experiment in config

Stages: integrate | impute | grn | evaluate.  Running a later stage reuses what
earlier stages already wrote, so you can iterate on GRN without redoing scSAGA.

Every path comes from config.yml (see that file).  Output layout:

    results/<experiment>/integration/            joint_embedding_H.npy + aligned_*.npy
    results/<experiment>/<strategy>/             all_cells_gene_expression.npy
                                                 grn_tfonly_reg/network.tsv
                                                 grn_tfonly_reg/evaluation/
"""
from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import yaml

import engine
from engine import Dataset, resolve

STAGES = ['integrate', 'impute', 'grn', 'evaluate']


def load_config(path, root_override=None):
    path = os.path.abspath(path)
    cfg = yaml.safe_load(open(path))
    root = os.path.abspath(root_override or os.environ.get('PIPELINE_ROOT')
                           or os.path.dirname(path))
    os.environ['PIPELINE_ROOT'] = root
    cfg['_root'] = root
    cfg['_config_path'] = path
    return cfg


def build_datasets(cfg):
    out = {}
    for name, spec in (cfg.get('datasets') or {}).items():
        out[name] = Dataset(name, spec or {}, cfg['_root'])
    return out


def check(cfg, datasets, experiment=None):
    """Print what is present/missing. Returns True if the requested work can run."""
    exps = cfg.get('experiments') or {}
    todo = [experiment] if experiment else list(exps)
    ok = True
    print(f'root: {cfg["_root"]}')
    print(f'\ndatasets declared: {len(datasets)}')
    for n, d in datasets.items():
        state = []
        state.append('pca' if d.has_pca() else 'NO-PCA')
        state.append('expr' if d.has_expression() else 'no-expr')
        n_cells = d.n_cells()
        print(f'  {n:16s} {d.modality:5s} {", ".join(state):18s}'
              f'{"  " + str(n_cells) + " cells" if n_cells else ""}')
    print('\nexperiments:')
    for e in todo:
        spec = exps.get(e)
        if not spec:
            print(f'  {e}: NOT IN CONFIG'); ok = False; continue
        names = spec['datasets']
        unknown = [n for n in names if n not in datasets]
        if unknown:
            print(f'  {e}: unknown datasets {unknown}'); ok = False; continue
        refs = spec.get('references') or {}
        missing_pca = [n for n in names if not datasets[n].has_pca()]
        need_expr = set(spec.get('queries') or [])
        for r in refs.values():
            if 'dataset' in r:
                need_expr.add(r['dataset'])
            else:
                need_expr.update(n for n in names if datasets[n].modality == 'rna')
        missing_expr = [n for n in sorted(need_expr)
                        if not datasets[n].has_expression()]
        flags = []
        if missing_pca:
            flags.append(f'MISSING PCA: {missing_pca}'); ok = False
        if missing_expr:
            flags.append(f'MISSING EXPRESSION: {missing_expr}'); ok = False
        status = 'OK' if not flags else ' | '.join(flags)
        print(f'  {e}: {" + ".join(names)}  anchor={spec["anchor"]}')
        print(f'      strategies: {", ".join(refs)}   -> {status}')
    grn = cfg.get('grn') or {}
    for label, key, default in [('regulators', 'regulators', 'data/tf_only.txt'),
                                ('targets', 'targets', 'data/trrust_tf.txt')]:
        rel = grn.get(key) or default
        p = resolve(cfg['_root'], rel)
        print(f'  {label:11s} {p}  {"OK" if os.path.exists(p) else "MISSING"}')
        if not os.path.exists(p):
            ok = False
    return ok


def run_experiment(cfg, datasets, exp_name, stages, args):
    spec = cfg['experiments'][exp_name]
    names = spec['datasets']
    anchor = spec['anchor']
    queries = spec.get('queries') or []
    references = spec.get('references') or {}
    root = cfg['_root']
    exp_dir = f'{root}/results/{exp_name}'
    integ_dir = f'{exp_dir}/integration'
    os.makedirs(exp_dir, exist_ok=True)

    ds = [datasets[n] for n in names]
    rna_names = [d.name for d in ds if d.modality == 'rna']
    rna_datasets = {d.name: d for d in ds if d.modality == 'rna'}
    scsaga_repo = os.environ.get('SCAGA_REPO') or \
        resolve(root, cfg.get('scsaga_repo', 'tools/scSAGA'))

    print(f'\n{"="*72}\nEXPERIMENT {exp_name}\n{"="*72}')
    print(f'  datasets   : {names}')
    print(f'  anchor     : {anchor}')
    print(f'  queries    : {queries}')
    print(f'  RNA for ref: {rna_names}')
    print(f'  strategies : {list(references)}')
    print(f'  scSAGA     : {scsaga_repo}')

    # -------------------------------------------------- integrate
    H_path = f'{integ_dir}/joint_embedding_H.npy'
    if 'integrate' in stages:
        print(f'\n--- integrate ({len(ds)} datasets, anchor {anchor})')
        t0 = time.time()
        H, order, sizes = engine.integrate(
            ds, anchor, integ_dir, scsaga_repo,
            params=cfg.get('integration') or {},
        )
        print(f'    ({time.time()-t0:.1f}s)')
    else:
        if not os.path.exists(H_path):
            raise SystemExit(f'no {H_path}; run --stage integrate first')
        info = yaml.safe_load(open(f'{integ_dir}/integration_info.json'))
        order, sizes = info['order'], {k: int(v) for k, v in info['sizes'].items()}
        H = np.load(H_path)
        print(f'\n--- reuse integration H {H.shape}')

    # cache real RNA expression once per experiment
    rna_expr = {}
    for n in rna_names:
        d = rna_datasets[n]
        genes = d.load_features()
        rna_expr[n] = d.load_expression_genes_by_cells(len(genes))

    # -------------------------------------------------- per-strategy
    for strat, ref_spec in references.items():
        sdir = f'{exp_dir}/{strat}'
        grn_dir = f'{sdir}/grn_tfonly_reg'
        os.makedirs(sdir, exist_ok=True)
        all_expr_path = f'{sdir}/all_cells_gene_expression.npy'
        print(f'\n===== strategy {strat}: {ref_spec}')

        genes = rna_datasets[rna_names[0]].load_features()
        if 'impute' in stages:
            print('--- build reference + reverse-imputeKNN')
            ref_expr, ref_H, ref_bc, genes = engine.build_reference(
                ref_spec, rna_datasets, root, integ_dir)
            query_bc = {q: datasets[q].load_barcodes() for q in queries}
            imputed, q_bc = engine.impute(H, order, sizes, ref_expr, ref_H,
                                          ref_bc, queries, sdir, query_bc,
                                          rna_names)
            np.save(f'{sdir}/genes.npy', np.array(genes))
            with open(f'{sdir}/genes.txt', 'w') as fh:
                fh.write('\n'.join(genes) + '\n')
            # all-cells rows: every real RNA cell in H order, then queries
            real = np.vstack([rna_expr[n].T for n in
                              [o for o in order if o in rna_names]]).astype(np.float32)
            all_expr = np.vstack([real, imputed.T]).astype(np.float32)
            all_bc = [b for n in [o for o in order if o in rna_names]
                      for b in datasets[n].load_barcodes()] + q_bc
            np.save(all_expr_path, all_expr)
            with open(f'{sdir}/all_cells_barcodes.txt', 'w') as fh:
                fh.write('\n'.join(all_bc) + '\n')
            with open(f'{sdir}/n_cells.txt', 'w') as fh:
                fh.write(str(all_expr.shape[0]))
            print(f'    all-cells matrix {all_expr.shape} '
                  f'(rows: {len(all_bc) - len(q_bc)} real RNA + {len(q_bc)} imputed ATAC)')
        else:
            if not os.path.exists(all_expr_path):
                raise SystemExit(f'no {all_expr_path}; run --stage impute first')
            genes = [l.strip() for l in open(f'{sdir}/genes.txt')]
            print(f'--- reuse all-cells matrix {np.load(all_expr_path, mmap_mode="r").shape}')

        if 'grn' in stages:
            print(f'--- Arboreto GRNBoost2 ({args.workers} workers)')
            engine.run_grn(all_expr_path, genes,
                           resolve(root, cfg['grn']['regulators']),
                           resolve(root, cfg['grn']['targets']),
                           grn_dir, n_workers=args.workers,
                           threads_per_worker=args.threads_per_worker,
                           seed=args.seed, scheduler=args.scheduler,
                           max_targets=args.max_targets,
                           max_regulators=args.max_regulators)

        if 'evaluate' in stages:
            print('--- evaluation')
            gt = {k: resolve(root, v) for k, v in (cfg.get('ground_truth') or {}).items()}
            engine.evaluate(f'{grn_dir}/network.tsv', gt, f'{grn_dir}/evaluation')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default=os.environ.get('PIPELINE_CONFIG', 'config.yml'))
    ap.add_argument('--root', default=None, help='override PIPELINE_ROOT')
    ap.add_argument('--experiment', '-e', action='append', default=[],
                    help='experiment name from config.yml (repeatable)')
    ap.add_argument('--all', action='store_true', help='every experiment in config')
    ap.add_argument('--list', action='store_true', help='list experiments and exit')
    ap.add_argument('--check', action='store_true', help='validate inputs and exit')
    ap.add_argument('--stage', action='append', default=[],
                    choices=STAGES, help='run only these stages')
    ap.add_argument('--workers', type=int,
                    default=int(os.environ.get('GRN_WORKERS', '0')) or None)
    ap.add_argument('--threads-per-worker', type=int,
                    default=int(os.environ.get('GRN_THREADS_PER_WORKER', '1')))
    ap.add_argument('--scheduler', default=os.environ.get('GRN_SCHEDULER'))
    ap.add_argument('--seed', type=int, default=666)
    # Smoke-test only.  grnboost2() fits one regression per COLUMN, so the only
    # way to make a run cheap is to prune the column set -- these truncate the
    # target and regulator lists before the union is built.  Leave unset for a
    # real run; both default off.
    ap.add_argument('--max-targets', type=int,
                    default=int(os.environ.get('GRN_MAX_TARGETS', '0')) or None)
    ap.add_argument('--max-regulators', type=int,
                    default=int(os.environ.get('GRN_MAX_REGULATORS', '0')) or None)
    args = ap.parse_args()

    if args.workers is None:
        import multiprocessing
        args.workers = int(os.environ.get('SLURM_CPUS_PER_TASK')
                           or multiprocessing.cpu_count())

    cfg = load_config(args.config, args.root)
    datasets = build_datasets(cfg)
    exps = cfg.get('experiments') or {}

    if args.list:
        print(f'config : {cfg["_config_path"]}')
        print(f'root   : {cfg["_root"]}')
        print('\nexperiments:')
        for name, spec in exps.items():
            print(f'  {name:12s} {" + ".join(spec["datasets"])}  '
                  f'(anchor {spec["anchor"]}, strategies: '
                  f'{", ".join(spec.get("references") or {})})')
        return 0

    if args.check:
        return 0 if check(cfg, datasets, args.experiment[0] if args.experiment else None) else 1

    targets = list(exps) if args.all else args.experiment
    if not targets:
        raise SystemExit('choose --experiment NAME, --all, --list or --check '
                         '(see config.yml)')
    missing = [t for t in targets if t not in exps]
    if missing:
        raise SystemExit(f'unknown experiment(s) {missing}; available: {list(exps)}')

    stages = args.stage or STAGES
    print(f'config   : {cfg["_config_path"]}')
    print(f'root     : {cfg["_root"]}')
    print(f'experiments: {targets}')
    print(f'stages   : {stages}')
    print(f'workers  : {args.workers} (threads/worker {args.threads_per_worker})')

    t0 = time.time()
    for e in targets:
        run_experiment(cfg, datasets, e, stages, args)
    print(f'\n{"="*72}\nALL DONE in {time.time()-t0:.1f}s\n{"="*72}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
