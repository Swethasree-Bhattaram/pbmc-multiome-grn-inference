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

ATAC datasets are included when they carry counts: MultiVI reads raw peaks, so an
ATAC counts.mtx stored features x cells breaks it (a PCA-only ATAC dir is skipped).

    python fix_counts.py --dry-run        # report only, changes nothing
    python fix_counts.py                  # repair
    CONFIG=/path/to/config.yml python fix_counts.py

Covers ATAC datasets as well as RNA ones: MultiVI reads raw peaks, so an ATAC
counts.mtx stored features x cells breaks that arm (a PCA-only ATAC dir is skipped).

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


COUNTS_NAMES = ('counts.mtx', 'counts.txt', 'counts.tsv')


def _in_dir(d, names):
    """First of `names` that exists in d, else the first name (for the message)."""
    for n in names:
        if os.path.exists(os.path.join(d, n)):
            return os.path.join(d, n)
    return os.path.join(d, names[0])


def is_matrixmarket(path):
    try:
        with open(path, 'rb') as fh:
            return fh.read(64).lstrip().startswith(b'%%MatrixMarket')
    except OSError:
        return False


def paths(spec, root):
    """-> (counts, barcodes, features) absolute paths for one dataset spec.

    Mirrors run.py's Dataset: explicit keys win, `dir` is a shorthand.  A matrix
    named counts.txt / counts.tsv is found too -- scipy.io.mmread does not care
    about the suffix, only about the content.
    """
    def p(key, default=None):
        v = spec.get(key)
        if v:
            return v if os.path.isabs(v) else os.path.join(root, v)
        return default

    d = p('dir')
    if d and os.path.isdir(d):
        return (_in_dir(d, COUNTS_NAMES) if not spec.get('counts') else p('counts'),
                p('barcodes', _in_dir(d, ('barcodes.txt',))) if not spec.get('barcodes')
                else p('barcodes'),
                p('features', _in_dir(d, ('features.txt',))) if not spec.get('features')
                else p('features'))
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

    rename = {}
    print(f'config: {os.path.abspath(CONFIG)}\n')
    for name, spec in (cfg.get('datasets') or {}).items():
        spec = spec or {}
        c, b, f = paths(spec, root)
        if not all(p and os.path.exists(p) for p in (c, b, f)):
            print(f'  {name:12s} SKIP  missing counts/barcodes/features '
                  f'(looked for {", ".join(COUNTS_NAMES)} in the dataset dir)')
            continue
        if os.path.basename(c) != 'counts.mtx':
            if is_matrixmarket(c):
                print(f'  {name:12s} NAME  {os.path.basename(c)} is MatrixMarket; '
                      f'will become counts.mtx')
                rename[name] = c
            else:
                print(f'  {name:12s} SKIP  {os.path.basename(c)} is not MatrixMarket '
                      f'(no %%MatrixMarket header) -- convert it first, see '
                      f'README_STEPWISE.md')
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

    print(f'\n{len(todo)} dataset(s) to transpose, {len(rename)} to rename'
          f'{" (dry run, nothing written)" if dry else ""}')
    if dry or not (todo or rename):
        return 0

    # rename first: a later transpose must write to the new path
    for name, src in rename.items():
        dst = os.path.join(os.path.dirname(src), 'counts.mtx')
        if os.path.exists(dst):
            print(f'  {name:12s} SKIP  {dst} already exists, not overwriting')
            continue
        os.replace(src, dst)
        print(f'  {name:12s} renamed {os.path.basename(src)} -> counts.mtx')
        if name in todo:
            todo[name] = (dst,) + todo[name][1:]

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
