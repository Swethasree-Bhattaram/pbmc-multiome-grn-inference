#!/usr/bin/env python
"""SCEMENT-combine several RNA datasets into one batch-corrected reference.

Runs in the SCEMENT environment (needs anndata/scanpy).  Called as a subprocess
by engine.build_reference so the anndata dependency never has to live in the main
pipeline env.

    python pipeline/scement_combine.py --out ref.npz \
        --in  name=/path/counts.mtx:/path/features.txt \
        --in  name2=/path/counts.mtx:/path/features.txt

Writes an .npz with:
    expression   genes x cells  float32  (batch-corrected, clipped at 0)
    genes        (n_genes,)     unicode
    barcodes     (n_cells,)     unicode
    batches      (n_cells,)     unicode  (which input each cell came from)
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys

import numpy as np
import scipy.io
import scipy.sparse as sp


def load_sct_sparse():
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scement_sparse.py')
    spec = importlib.util.spec_from_file_location('scement_sparse', src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_lines(p):
    with open(p) as fh:
        return [l.strip() for l in fh if l.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--in', dest='inputs', action='append', required=True,
                    metavar='NAME=COUNTS:FEATURES[:BARCODES]',
                    help='one RNA dataset; barcodes optional, defaults to '
                         'a sibling barcodes.txt')
    args = ap.parse_args()

    import pandas as pd
    import anndata as ad

    blocks, genes_ref, bcs, batches = [], None, [], []
    for item in args.inputs:
        name, _, rest = item.partition('=')
        parts = rest.split(':')
        counts, features = parts[0], parts[1]
        barcodes = parts[2] if len(parts) > 2 else os.path.join(
            os.path.dirname(counts), 'barcodes.txt')

        genes = read_lines(features)
        if genes_ref is None:
            genes_ref = genes
        elif genes != genes_ref:
            raise SystemExit(
                f'{name}: feature list differs from the first dataset. '
                f'SCEMENT combining requires identical gene ordering.'
            )
        X = scipy.io.mmread(counts).tocsr()          # cells x genes
        if X.shape[1] != len(genes):
            raise SystemExit(f'{name}: counts has {X.shape[1]} columns but '
                             f'{len(genes)} features were listed')
        bc = read_lines(barcodes)
        if len(bc) != X.shape[0]:
            raise SystemExit(f'{name}: {X.shape[0]} count rows but {len(bc)} barcodes')
        blocks.append(X)
        bcs.extend(bc)
        batches.extend([name] * X.shape[0])
        print(f'  {name}: {X.shape} cells x genes', flush=True)

    X = sp.vstack(blocks).tocsr()
    obs = pd.DataFrame({'batch': pd.Categorical(batches)})
    adata = ad.AnnData(X=X, obs=obs)
    print(f'  SCEMENT input: {adata.shape}, '
          f'batches={obs["batch"].value_counts().to_dict()}', flush=True)

    scm = load_sct_sparse()
    corrected = scm.sct_sparse(adata, key='batch', inplace=False)
    corrected = np.clip(np.asarray(corrected), 0, None).astype(np.float32)
    expr = corrected.T                                # genes x cells
    print(f'  SCEMENT output: {expr.shape} (genes x cells)', flush=True)

    np.savez_compressed(args.out,
                        expression=expr,
                        genes=np.array(genes_ref),
                        barcodes=np.array(bcs),
                        batches=np.array(batches))
    print(f'  wrote {args.out}', flush=True)


if __name__ == '__main__':
    sys.exit(main())
