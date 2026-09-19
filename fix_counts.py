#!/usr/bin/env python
"""Repair counts.mtx orientation for every RNA dataset in config.yml.

counts.mtx must be CELLS x FEATURES (rows == barcodes, cols == features).
A dataset whose MatrixMarket header reads

    <features> <cells> <nnz>       e.g.  38606 5048 13996959

is stored the other way round.  Every RNA dataset in config.yml is checked
against its own barcodes.txt / features.txt:

    already cells x features  -> left byte-identical
    features x cells          -> rewritten transposed, atomically
    matches neither           -> refused, nothing written

ATAC datasets are skipped: their expression is imputed, never read.

    python fix_counts.py --dry-run        # report only, changes nothing
    python fix_counts.py                  # repair
    CONFIG=/path/to/config.yml python fix_counts.py

Run this before a long pipeline run: run.py only discovers the problem after
the integration step, whereas this reads file headers and finishes in seconds.
Idempotent -- safe to run twice, and safe to run on data that is already fine.
"""
from __future__ import annotations

import os
import sys

import scipy.io as io
import scipy.sparse as sp
import yaml

CONFIG = os.environ.get('CONFIG', 'config.yml')


def paths(spec, root):
    """-> (counts, barcodes, features) absolute paths for one dataset spec."""
    if spec.get('dir'):
        d = spec['dir']
        if not os.path.isabs(d):
            d = os.path.join(root, d)
        return (os.path.join(d, 'counts.mtx'),
                os.path.join(d, 'barcodes.txt'),
                os.path.join(d, 'features.txt'))

    def p(key):
        v = spec.get(key)
        if not v:
            return None
        return v if os.path.isabs(v) else os.path.join(root, v)

    return p('counts'), p('barcodes'), p('features')


def header(path):
    """MatrixMarket dimension line -> (rows, cols)."""
    with open(path) as fh:
        for line in fh:
            if not line.startswith('%'):
                r, c = line.split()[:2]
                return int(r), int(c)
    raise SystemExit(f'{path}: no dimension line')


def nlines(path):
    with open(path) as fh:
        return sum(1 for _ in fh)


def main():
    dry = '--dry-run' in sys.argv
    root = os.path.abspath(os.path.dirname(os.path.abspath(CONFIG)) or '.')
    cfg = yaml.safe_load(open(CONFIG))
    todo = {}

    print(f'config: {os.path.abspath(CONFIG)}\n')
    for name, spec in (cfg.get('datasets') or {}).items():
        spec = spec or {}
        if (spec.get('modality') or 'rna').lower() != 'rna':
            continue                       # ATAC is imputed, never read
        c, b, f = paths(spec, root)
        if not all(p and os.path.exists(p) for p in (c, b, f)):
            print(f'  {name:12s} SKIP  missing counts/barcodes/features')
            continue
        rows, cols = header(c)
        n_bc, n_ft = nlines(b), nlines(f)
        if (rows, cols) == (n_bc, n_ft):
            print(f'  {name:12s} OK    {rows} x {cols} '
                  f'(already cells x features)')
        elif (rows, cols) == (n_ft, n_bc):
            print(f'  {name:12s} FIX   {rows} x {cols} -> {n_bc} x {n_ft}')
            todo[name] = (c, n_bc, n_ft)
        else:
            print(f'  {name:12s} SKIP  {rows} x {cols} matches neither '
                  f'{n_bc} x {n_ft} nor {n_ft} x {n_bc} -- NOT touching it')

    print(f'\n{len(todo)} dataset(s) to transpose'
          f'{"" if not dry else " (dry run, nothing written)"}')
    if dry or not todo:
        return 0

    for name, (c, n_bc, n_ft) in todo.items():
        # scipy's mmwrite appends '.mtx' to any path not already ending in it,
        # so the temp file must carry the suffix for os.replace to find it.
        tmp = c + '.new.mtx'
        io.mmwrite(tmp, sp.csr_matrix(io.mmread(c).transpose()),
                   field='integer')
        os.replace(tmp, c)
        rows, cols = header(c)
        assert (rows, cols) == (n_bc, n_ft), f'{name}: wrote {rows} x {cols}'
        print(f'  {name:12s} written {rows} x {cols}')
    print('\ndone -- re-run with --dry-run to confirm, then: python run.py --check')
    return 0


if __name__ == '__main__':
    sys.exit(main())
