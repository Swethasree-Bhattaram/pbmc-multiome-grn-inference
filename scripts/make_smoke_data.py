"""Build a tiny paired multiome dataset so the whole pipeline can be smoke-run.

Two datasets (rna + atac) over the SAME cells (shared barcodes) generated from one
shared per-cell latent state, so the cell's cross-modal partner is genuinely
recoverable and every metric in stepwise_eval.py has something to measure:

  * a_c in R^K   : per-cell program activity, shared by both modalities
  * gene loading W, peak loading V : program -> feature weights
  * counts ~ Poisson(depth * a @ W) and Poisson(depth * a @ V)
  * labels.txt   : cell type = argmax of the cell's first four program activities
  * truth.csv    : TF_k -> the genes most loaded on program k (a planted network)
  * pca_50.txt   : 20-dim PCA of each modality's log1p counts

Writes the four-file convention config.yml already expects (pca_50.txt /
counts.mtx / barcodes.txt / features.txt) plus a config_smoke.yml next to it.

    python scripts/make_smoke_data.py --out results/smoke/data
    python scripts/stepwise_eval.py --config results/smoke/config_smoke.yml --check
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import scipy.io
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/smoke/data")
    ap.add_argument("--cells", type=int, default=320)
    ap.add_argument("--genes", type=int, default=240)
    ap.add_argument("--peaks", type=int, default=400)
    ap.add_argument("--programs", type=int, default=6)
    ap.add_argument("--types", type=int, default=4)
    ap.add_argument("--per-program", type=int, default=25,
                    help="genes driven by each program (= planted targets per TF)")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    os.makedirs(a.out, exist_ok=True)

    n, g, p, K, T = a.cells, a.genes, a.peaks, a.programs, a.types
    genes = [f"G{i:03d}" for i in range(g)]
    peaks = [f"P{i:03d}" for i in range(p)]
    barcodes = [f"CELL{i:04d}-1" for i in range(n)]

    # ---- shared per-cell latent state ------------------------------------
    act = rng.gamma(2.0, 1.0, size=(n, K))
    act[rng.random((n, K)) < 0.35] = 0.0             # sparse programs per cell
    # a cell-specific component shared by both modalities: without it every cell
    # of a type looks alike and no integrator can recover the true partner
    # (measured: a supervised CCA oracle scores FOSCTTM 0.30 on program-only data)
    act = np.clip(act + 0.45 * rng.normal(0, 1.0, size=(n, K)), 0, None)
    depth = rng.uniform(0.6, 1.8, size=n)

    # ---- program -> feature loadings -------------------------------------
    W = rng.gamma(1.5, 0.6, size=(K, g)) * (rng.random((K, g)) < 0.25)
    V = rng.gamma(1.5, 0.6, size=(K, p)) * (rng.random((K, p)) < 0.20)

    mu_r = (act @ W) * depth[:, None] * 1.2
    mu_p = (act @ V) * depth[:, None] * 1.2
    Xr = sp.csr_matrix(rng.poisson(mu_r).astype(np.float32))
    Xa = sp.csr_matrix(rng.poisson(mu_p).astype(np.float32))

    # ---- labels: cell type = argmax over the first T programs ------------
    y = np.argmax(act[:, :T], axis=1)

    def pca(X, k=20):
        from sklearn.decomposition import PCA
        Z = np.log1p(np.asarray(X.todense()))
        Z = (Z - Z.mean(0)) / (Z.std(0) + 1e-9)
        return PCA(n_components=k, random_state=0).fit_transform(Z).astype(np.float32)

    for name, X, feats in (("rna", Xr, genes), ("atac", Xa, peaks)):
        d = os.path.join(a.out, name)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "counts.mtx"), "wb") as fh:
            scipy.io.mmwrite(fh, X, field="integer")
        with open(os.path.join(d, "barcodes.txt"), "w") as fh:
            fh.write("\n".join(barcodes) + "\n")
        with open(os.path.join(d, "features.txt"), "w") as fh:
            fh.write("\n".join(feats) + "\n")
        np.savetxt(os.path.join(d, "pca_50.txt"), pca(X), fmt="%.6f")
        with open(os.path.join(d, "labels.txt"), "w") as fh:
            fh.write("\n".join(f"type{int(t)}" for t in y) + "\n")

    # ---- planted network: one TF per program drives its top genes -------
    truth = []
    tfs = []
    for k in range(T):
        top = np.argsort(-W[k])[:a.per_program]
        tf = genes[int(top[0])]                      # the TF is itself a gene
        tfs.append(tf)
        truth += [(tf, genes[int(j)]) for j in top[1:]]
    with open(os.path.join(a.out, "tf.txt"), "w") as fh:
        fh.write("\n".join(tfs) + "\n")
    with open(os.path.join(a.out, "targets.txt"), "w") as fh:
        fh.write("\n".join(genes) + "\n")
    with open(os.path.join(a.out, "truth.csv"), "w") as fh:
        fh.write("TF,TARGET\n")
        for tf, tg in truth:
            fh.write(f"{tf},{tg}\n")

    root = os.path.abspath(a.out)
    cfg = f"""datasets:
  rna_s:  {{modality: rna,  dir: {os.path.join(root, 'rna')}}}
  atac_s: {{modality: atac, dir: {os.path.join(root, 'atac')}}}
experiments:
  smoke:
    datasets: [rna_s, atac_s]
    anchor:   rna_s
    queries:  [atac_s]
    references:
      rna_s: {{dataset: rna_s}}
grn:
  regulators: {os.path.join(root, 'tf.txt')}
  targets:    {os.path.join(root, 'targets.txt')}
ground_truth:
  planted: {os.path.join(root, 'truth.csv')}
integration:
  s_shared_cells: {n}
"""
    path = os.path.abspath(os.path.join(a.out, "..", "config_smoke.yml"))
    with open(path, "w") as fh:
        fh.write(cfg)
    print(f"wrote {a.out}: {n} cells, {g} genes, {p} peaks, {K} programs, "
          f"{T} types; planted {len(truth)} edges from {tfs}")
    print(f"config: {path}")


if __name__ == "__main__":
    main()
