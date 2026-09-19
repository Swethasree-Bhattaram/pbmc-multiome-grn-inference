#!/usr/bin/env python
"""Compare the feature lists of the RNA datasets used by a `combine` strategy.

The `combine` reference stacks several RNA counts matrices in gene space, so all
of its datasets must share an IDENTICAL feature list -- same genes, same order.
When they differ, build_reference() aborts with:

    combine needs identical feature ordering; X differs from Y

This script says HOW they differ, which decides what to do about it:

    same genes, different order   -> reordering is a safe, exact fix
    different lengths             -> different reference builds; intersect or drop
    same length, a few mismatches -> different reference builds; inspect them

Usage:
    python compare_features.py A/features.txt B/features.txt [C/features.txt ...]
    python compare_features.py --config config.yml --experiment multi
"""
from __future__ import annotations

import os
import sys


def read_features(path):
    """Feature file -> gene SYMBOLS.

    Cell Ranger feature files come in three shapes; only the symbol field can
    ever match a gene list, so always take column 2 when present:

        SYMBOL
        ENSG...<TAB>SYMBOL
        ENSG...<TAB>SYMBOL<TAB>Type      ('Gene Expression' etc.)

    Comparing raw lines instead reports 'shared genes: 0' for files that in
    fact share every gene but use different formats.
    """
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split('\t') if '\t' in line else line.split()
            out.append(parts[1] if len(parts) >= 2 else parts[0])
    return out


def report(ref_path, ref, other_path, other):
    print(f'\n===== {os.path.basename(os.path.dirname(other_path))} vs '
          f'{os.path.basename(os.path.dirname(ref_path))} =====')
    print(f'  lengths      : {len(ref)} vs {len(other)}'
          f'{"  (DIFFERENT)" if len(ref) != len(other) else ""}')

    if ref == other:
        print('  VERDICT      : identical (same genes, same order) -- fine')
        return 'same'

    set_ref, set_other = set(ref), set(other)
    common, only_ref, only_other = (set_ref & set_other,
                                    set_ref - set_other,
                                    set_other - set_ref)
    print(f'  shared genes : {len(common)}')
    print(f'  only in {os.path.basename(os.path.dirname(ref_path)):18s}: '
          f'{len(only_ref)}')
    print(f'  only in {os.path.basename(os.path.dirname(other_path)):18s}: '
          f'{len(only_other)}')

    if len(ref) == len(other) and not only_ref and not only_other:
        first = next(i for i, (a, b) in enumerate(zip(ref, other)) if a != b)
        print(f'  SAME GENE SET, DIFFERENT ORDER (first mismatch at index {first}: '
              f'{ref[first]!r} vs {other[first]!r})')
        print('  VERDICT      : reorder this dataset\'s genes/counts to match '
              'the first -> exact fix')
        return 'order'

    if only_ref:
        print(f'  e.g. only in ref   : {sorted(only_ref)[:8]}')
    if only_other:
        print(f'  e.g. only in other : {sorted(only_other)[:8]}')
    print('  VERDICT      : different gene sets -- NOT a reordering. These are '
          'different reference builds.')
    return 'geneset'


def main():
    args = sys.argv[1:]
    if not args:
        raise SystemExit(__doc__)

    if args[0] == '--config':
        import yaml
        cfg_path = args[1]
        exp = None
        if '--experiment' in args:
            exp = args[args.index('--experiment') + 1]
        root = os.path.abspath(os.path.dirname(os.path.abspath(cfg_path)) or '.')
        cfg = yaml.safe_load(open(cfg_path))
        spec = (cfg.get('experiments') or {})[exp] if exp else None
        if spec is None:
            raise SystemExit(f'no experiment {exp!r} in {cfg_path}')
        names = [n for n in (spec.get('datasets') or [])]
        paths = []
        for n in names:
            d = (cfg.get('datasets') or {}).get(n) or {}
            if (d.get('modality') or 'rna').lower() != 'rna':
                continue
            f = d.get('features')
            if not f:
                base = d.get('dir')
                if not base:
                    continue
                f = os.path.join(base if os.path.isabs(base)
                                 else os.path.join(root, base), 'features.txt')
            elif not os.path.isabs(f):
                f = os.path.join(root, f)
            if os.path.exists(f):
                paths.append((n, f))
        if not paths:
            raise SystemExit('no RNA datasets with a features.txt found')
        print(f'RNA datasets in experiment {exp!r}: '
              f'{", ".join(n for n, _ in paths)}')
    else:
        paths = [(os.path.basename(os.path.dirname(p)) or p, p) for p in args]

    loaded = [(n, p, read_features(p)) for n, p in paths]
    for n, p, feats in loaded:
        print(f'  {n:12s} {len(feats):7d} features   {p}')

    ref_name, ref_path, ref = loaded[0]
    verdicts = set()
    for n, p, feats in loaded[1:]:
        verdicts.add(report(ref_path, ref, p, feats))

    print('\n' + '=' * 68)
    if verdicts == {'same'}:
        print('All feature lists identical -> `combine` should work.')
    elif 'geneset' in verdicts:
        print('Some datasets use different gene sets. `combine` CANNOT stack '
              'them as-is.\nOptions:\n'
              '  1. drop the combined strategy from that experiment (keep the\n'
              '     per-dataset strategies), or\n'
              '  2. re-derive the datasets from one shared 10x reference, or\n'
              '  3. subset every dataset to the shared genes, in a fixed order\n'
              '     (changes the gene set the reference covers).')
    else:
        print('All lists share the same genes but in different orders.\n'
              'Reordering each dataset to the first list is an exact fix.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
