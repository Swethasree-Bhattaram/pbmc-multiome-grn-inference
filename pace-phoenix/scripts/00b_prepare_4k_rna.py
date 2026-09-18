#!/usr/bin/env python
"""Convert the 10x PBMC 4k (v2, Cell Ranger 2.1.0) filtered matrix tar.gz into the
scSAGA/4k-baseline file layout.

Input:  raw/pbmc4k_filtered_gene_bc_matrices.tar.gz
        (https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/
         pbmc4k_filtered_gene_bc_matrices.tar.gz)
        The archive holds filtered_gene_bc_matrices/GRCh38/{matrix.mtx,genes.tsv,barcodes.tsv}

Output: data/4k_rna/{matrix.mtx,genes.tsv,barcodes.txt}

Notes
-----
* matrix.mtx is genes x cells and is kept as-is (the consumer
  scripts/rna-baseline/prepare_rna_baseline.py mmread()s it and transposes).
* barcodes.tsv is copied to barcodes.txt -- only the name changes; the original
  run's 4k_rna/barcodes.txt is byte-identical to barcodes.tsv.
* genes.tsv (Ensembl id <TAB> symbol) is kept as-is and parsed for the symbol
  column by the consumer.
* No PCA file: the 4k dataset is never handed to scSAGA, only stacked with the 3k
  RNA by the RNA-only baseline.

Usage:
  python 00b_prepare_4k_rna.py [--raw-dir raw] [--data-dir data]
"""
import argparse
import os
import shutil
import subprocess
import sys
import tarfile

EXPECTED = {'matrix.mtx': 71464390, 'genes.tsv': 840938, 'barcodes.tsv': 82460}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.environ.get('PACE_ROOT', os.path.dirname(here))
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw-dir', default=f'{root}/raw')
    ap.add_argument('--data-dir', default=f'{root}/data')
    args = ap.parse_args()

    tar_path = os.path.join(args.raw_dir, 'pbmc4k_filtered_gene_bc_matrices.tar.gz')
    if not os.path.exists(tar_path):
        raise SystemExit(f'missing {tar_path}; run 01_download_data.sh first')

    out = os.path.join(args.data_dir, '4k_rna')
    os.makedirs(out, exist_ok=True)

    with tarfile.open(tar_path, 'r:gz') as tf:
        members = {m.name: m for m in tf.getmembers() if m.isfile()}
        want = {}
        for name in members:
            base = os.path.basename(name)
            if base in EXPECTED:
                want[base] = name
        missing = set(EXPECTED) - set(want)
        if missing:
            raise SystemExit(f'archive is missing {sorted(missing)}')
        for base, member in want.items():
            dest_name = 'barcodes.txt' if base == 'barcodes.tsv' else base
            dest = os.path.join(out, dest_name)
            with tf.extractfile(member) as src, open(dest, 'wb') as dst:
                shutil.copyfileobj(src, dst)
            print(f'  {base:>14s} -> {dest} ({os.path.getsize(dest)} bytes)', flush=True)

    bad = {b: (os.path.getsize(os.path.join(out, 'barcodes.txt' if b == 'barcodes.tsv' else b)), n)
           for b, n in EXPECTED.items()
           if os.path.getsize(os.path.join(out, 'barcodes.txt' if b == 'barcodes.tsv' else b)) != n}
    if bad:
        print('WARNING: unexpected sizes (expected vs got):', bad, flush=True)
    n_cells = sum(1 for _ in open(os.path.join(out, 'barcodes.txt')))
    n_genes = sum(1 for _ in open(os.path.join(out, 'genes.tsv')))
    print(f'4k_rna: {n_genes} genes x {n_cells} cells -> {out}', flush=True)
    print('Done.', flush=True)


if __name__ == '__main__':
    sys.exit(main())
