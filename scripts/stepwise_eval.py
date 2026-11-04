#!/usr/bin/env python
"""Step-wise multiome -> GRN experiment, measured at every step, two ways.

WHAT THIS DOES
--------------
workflow.pptx describes several routes from per-dataset PCA matrices to a gene
network.  This script runs two of them over the SAME datasets declared in
config.yml and scores every step with the metrics the deck lists for that step,
before the next step runs:

  step 1 integration   slides 4 / 11  FOSCTTM, label conservation (ARI, NMI,
                                      silhouette, isolated-label F1, graph
                                      connectivity), batch correction (kBET-style
                                      acceptance, iLISI, batch ASW), label transfer
  step 2 imputation    slides 6 / 11  MSE, pFCS, pDES, gene-wise rho, densification,
                                      downstream (5 classification + 5 clustering)
  step 3 GRN           slide 8        AUPRC, early precision, stability, topology,
                                      perturbation (recorded as not computable),
                                      regulon activity -- plus run.py's existing
                                      vs-ground-truth recovery (recovered,
                                      precision@K, recall@K)

The two workflows differ in exactly one place, on purpose:

  multivi   one model does integration + imputation (workflow 2)
            raw RNA counts + raw ATAC peak counts + a per-cell modality mask
            -> joint latent -> decoded RNA for the ATAC cells
  scsaga    integration to a joint embedding H (workflow 1), then this repo's
            reverse-imputeKNN (k=20, exp(-distance)) as a separate gate-able step

From imputation onwards both use the same code, so a difference in the network is
attributable to the embedding/imputation and not to the GRN step.

HOW EACH METHOD SEES THE CELLS (this matters)
---------------------------------------------
scSAGA gets one block of cells per dataset, as run.py has always done: RNA cells
and ATAC cells are separate rows of H.

MultiVI gets a COLLAPSED row space: a barcode present in both an RNA and an ATAC
dataset becomes ONE cell with both modalities.  That is not cosmetic.  MultiVI
aligns its two encoders with a penalty that is applied only where a cell has both
modalities (`mask12 = mask_expr & mask_acc`, module/_multivae.py), so with separate
rows nothing would ever align and z would simply separate by modality.  Paired
(10x multiome) cells are what makes MultiVI a multi-omic model rather than two
VAEs.  Where no barcodes are shared the rows are separate anyway, and the model
has no co-assayed cells to align on -- the report will show that in FOSCTTM.

That difference is also why FOSCTTM is computed two ways, both recorded in
metrics.json: with separate rows it is the cross-modal retrieval test of the deck;
with collapsed rows the cell's two modality-specific latents (z_expr, z_acc) are
retrieved against each other, which is the same question asked of the only two
representations the model exports per cell.

Held-out answer key (design doc section 4)
------------------------------------------
When an ATAC dataset shares barcodes with the anchor RNA dataset, a deterministic
fraction of those paired cells is frozen:

  * their real RNA is removed from the imputation reference;
  * for MultiVI the RNA modality of those cells is ALSO removed from training (a
    second fit), so the hidden values never enter the loss or the encoder;
  * steps 1 and 2 are scored against it (FOSCTTM, MSE, pFCS, pDES, rho,
    densification, downstream).

INPUTS -- exactly the ones config.yml already declares
------------------------------------------------------
  every dataset     pca_50.txt                          cells x N
  RNA datasets      counts.mtx + barcodes.txt + features.txt
  ATAC datasets     counts.mtx + barcodes.txt + features.txt   (multivi only;
                    the scSAGA arm still needs only the PCA)
  any dataset       labels.txt (optional; one cell-type label per PCA row)

Without labels.txt every label-based metric (ARI/NMI/silhouette/isolated-label F1/
graph connectivity, label transfer, pDES, downstream, regulon specificity) is
reported as unavailable, with the reason.  Without shared barcodes FOSCTTM and the
held-out imputation block are reported as unavailable.  Nothing is approximated
silently.

USAGE
-----
  python scripts/stepwise_eval.py --check
  python scripts/stepwise_eval.py -e paired --methods multivi,scsaga
  python scripts/stepwise_eval.py -e paired --only metrics
  python scripts/stepwise_eval.py --only grn --seeds 666,667,668

Outputs: results/stepwise/<experiment>/<method>/{integration,pipeline,heldout}
         results/stepwise/<experiment>/report.md  +  metrics.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run as base                                   # noqa: E402  (repo pipeline)
import stepwise_metrics as M                         # noqa: E402

STEPS = ["integrate", "impute", "grn", "metrics"]
METHODS = ["multivi", "scsaga"]
TARGET_SUM = 1e4                                     # run.py's CPM target
K_IMPUTE = 20                                        # run.py's K_NEIGHBOURS
NORM_TXT = ("log1p(x * 1e4 / cell_total) from raw counts.  MultiVI's px_scale is "
            "a within-cell proportion, so it is mapped through the same transform "
            "as log1p(1e4 * px_scale) and then per-gene affine-aligned to the "
            "measured reference rows (--no-scale-align disables the alignment).")


def log(msg=""):
    print(msg, flush=True)


def wjson(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)


# --------------------------------------------------------------------------- #
# data access -- always through run.py's Dataset, never re-derived
# --------------------------------------------------------------------------- #
def counts_orientation(ds):
    """-> 'ok' | 'transposed' | 'neither', plus the shapes, from the header only.

    The repo convention is CELLS x FEATURES.  The 10x / ArchR / Signac ATAC exports
    are usually FEATURES x CELLS, and nothing used to read ATAC counts (run.py
    exempts ATAC from its alignment check, and fix_counts.py skipped ATAC), so a
    transposed ATAC matrix was invisible until MultiVI.
    """
    rows, cols, _ = ds.mtx_shape()
    n_cells, n_feat = ds.n_cells_pca(), len(ds.load_features())
    if (rows, cols) == (n_cells, n_feat):
        return "ok", (rows, cols), (n_cells, n_feat)
    if (rows, cols) == (n_feat, n_cells):
        return "transposed", (rows, cols), (n_cells, n_feat)
    return "neither", (rows, cols), (n_cells, n_feat)


def counts_on_axis(ds, axis):
    """csr cells x len(axis) RAW counts, reindexed from the dataset's own axis.

    A counts.mtx stored features x cells is transposed in memory with a loud line,
    rather than failing a job that is already hours in.  `fix_counts.py` repairs the
    file itself; this keeps a first run from dying on it.
    """
    import scipy.io

    feats = ds.load_features()
    state, (rows, cols), (n_cells, n_feat) = counts_orientation(ds)
    if state == "neither":
        raise SystemExit(
            f"{ds.name}: counts.mtx is {rows} x {cols}, which is neither cells x "
            f"features ({n_cells} x {n_feat}) nor features x cells ({n_feat} x "
            f"{n_cells}).  barcodes.txt/features.txt do not describe this matrix.")
    X = scipy.io.mmread(ds.counts).tocsr().astype(np.float32)
    if state == "transposed":
        log(f"  [counts] {ds.name}: counts.mtx is {rows} x {cols} = features x cells; "
            f"transposing in memory.  Repair the file with: python fix_counts.py")
        X = X.T.tocsr()
    if X.shape != (n_cells, n_feat):
        raise SystemExit(f"{ds.name}: after orientation fix counts is {X.shape}, "
                         f"expected {(n_cells, n_feat)}")
    if list(feats) == list(axis):
        return X
    return base.reindex_counts_columns(X, feats, list(axis))


def lognorm_on_axis(ds, axis):
    """run.py's load_expression() transform, on an arbitrary gene axis.

    Totals are taken on the dataset's OWN matrix (as run.py does) and the log is
    applied before reindexing, so a dataset contributes exactly the numbers run.py
    would have produced for the genes it shares with the axis.
    """
    genes = ds.load_features()
    X = counts_on_axis(ds, genes)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    Y = X.multiply((TARGET_SUM / tot)[:, None]).tocsr()
    Y.data = np.log1p(Y.data)
    Y = Y.astype(np.float32)
    if list(genes) == list(axis):
        return Y
    return base.reindex_counts_columns(Y, genes, list(axis))


def load_labels(ds):
    """Optional labels.txt (or a config 'labels' path) -> list of labels or None."""
    path = getattr(ds, "labels", None)
    if path and os.path.exists(path):
        return base.read_lines(path)
    cand = os.path.join(os.path.dirname(ds.pca), "labels.txt")
    if os.path.exists(cand):
        return base.read_lines(cand)
    return None


def has_barcodes(ds):
    return bool(ds.barcodes) and os.path.exists(ds.barcodes)


def paired_map(rna_ds, atac_ds):
    """{atac_local_index: rna_local_index} for barcodes the two share.

    Returns {} when either side has no barcodes.txt: without barcode files there is
    no way to say which ATAC cell is which RNA cell, so there is no answer key.
    """
    if not (has_barcodes(rna_ds) and has_barcodes(atac_ds)):
        return {}
    rb = rna_ds.load_barcodes()
    ab = atac_ds.load_barcodes()
    pos = {b: i for i, b in enumerate(rb)}
    return {j: pos[b] for j, b in enumerate(ab) if b in pos}


# --------------------------------------------------------------------------- #
# experiment layout
# --------------------------------------------------------------------------- #
class Exp:
    """Dataset declaration order, anchor first, and the pipeline row space.

    `starts`/`sizes` describe the CELL space (one row per dataset cell) that the
    all-cells matrix and the held-out split live in.  Each method's embedding has
    its own row space, stored in integration_info.json as `row_index`.
    """

    root = ""

    def __init__(self, cfg, datasets, name):
        spec = cfg["experiments"][name]
        self.name = name
        self.datasets = datasets
        self.names = list(spec["datasets"])
        self.anchor = spec["anchor"]
        self.queries = list(spec.get("queries") or [])
        self.references = spec.get("references") or {}
        self.order = [self.anchor] + [n for n in self.names if n != self.anchor]
        self.ds = {n: datasets[n] for n in self.order}
        self.rna = [n for n in self.order if datasets[n].modality == "rna"]
        self.atac = [n for n in self.order if datasets[n].modality == "atac"]
        self.sizes, self.starts = {}, {}
        off = 0
        for n in self.order:
            d = self.ds[n]
            # never load the matrix here: --check must be able to report a missing
            # PCA instead of dying on it
            k = d.n_cells() if d.has_pca() else 0
            self.sizes[n], self.starts[n] = k, off
            off += k
        self.n_cells = off
        self.atac_queries = [q for q in self.queries if q in self.atac] or self.atac

    def gene_axis(self):
        """The anchor RNA dataset's own features -- run.py's convention."""
        return list(self.ds[self.anchor].load_features())

    def peak_axis(self, mode="union"):
        return base.shared_gene_axis([self.ds[a].load_features() for a in self.atac],
                                     mode)

    def out(self, method, *parts):
        return os.path.join(self.root, method, *parts)

    def n_real_rna(self):
        return sum(self.sizes[n] for n in self.rna)


def load_integ(exp, method):
    """-> (H, row_index, info, extras) for one method's embedding."""
    d = exp.out(method, "integration")
    info = json.load(open(f"{d}/integration_info.json"))
    H = np.load(f"{d}/joint_embedding_H.npy")
    rows = {n: np.asarray(v, dtype=int) for n, v in info["row_index"].items()}
    extras = {}
    for key, fn in (("expr", "joint_embedding_H_expr.npy"),
                    ("acc", "joint_embedding_H_acc.npy")):
        p = f"{d}/{fn}"
        if os.path.exists(p):
            extras[key] = np.load(p)
    return H, rows, info, extras


# --------------------------------------------------------------------------- #
# MultiVI
# --------------------------------------------------------------------------- #
def multivi_available():
    try:
        import scvi
        import muon                                        # noqa: F401
        from scvi.model import MULTIVI                      # noqa: F401
        return True, f"scvi-tools {scvi.__version__}"
    except Exception as exc:                                # noqa: BLE001
        return False, f"{type(exc).__name__}: {exc}"


def collapsed_index(exp):
    """One row per distinct barcode; {dataset: row ids per local cell}.

    A barcode shared by several datasets lands on one row, which is what makes a
    paired cell co-assayed for MultiVI.  Rows are numbered by first appearance in
    `exp.order`, so the anchor's cells come first.
    """
    key2row, rows, owner = {}, {}, {}
    for n in exp.order:
        rows[n] = np.empty(exp.sizes[n], dtype=int)
        for i, b in enumerate(exp.ds[n].load_barcodes()):
            k = b
            if k not in key2row:
                key2row[k] = len(key2row)
                owner[key2row[k]] = n
            rows[n][i] = key2row[k]
    return rows, len(key2row), owner


def _multivi_mdata(exp, gene_axis, peak_axis, drop_rna_cells=()):
    """Padded MuData + the row ids of every dataset cell.

    scvi >= 1.4 requires every cell to appear in every modality, so a modality a
    cell does not have is a ZERO ROW; MultiVI reads the per-cell modality mask from
    `rowsum > 0` (module/_multivae.py, _get_inference_input).  Rows of frozen cells
    whose RNA is dropped keep their ATAC and lose their RNA -- exactly an ATAC-only
    cell, which is what the held-out evaluation needs.
    """
    import anndata as ad
    import muon
    import scipy.sparse as sp

    rows, n, owner = collapsed_index(exp)
    drop = {(exp.anchor, int(i)) for i in drop_rna_cells}
    blocks_r, blocks_a = [], []
    batch = np.empty(n, dtype=object)
    batch[:] = ""
    for r, o in owner.items():
        batch[r] = o
    for name in exp.order:
        d = exp.ds[name]
        axis = gene_axis if d.modality == "rna" else peak_axis
        X = counts_on_axis(d, axis)
        keep = np.array([i for i in range(exp.sizes[name])
                         if (name, i) not in drop], dtype=int)
        if not len(keep):
            continue
        rid = rows[name][keep]
        coo = X[keep].tocoo()
        m = sp.coo_matrix((coo.data, (rid[coo.row], coo.col)),
                          shape=(n, len(axis))).tocsr()
        (blocks_r if d.modality == "rna" else blocks_a).append(m)
    Xrna = (sp.vstack(blocks_r).tocsr() if blocks_r
            else sp.csr_matrix((n, len(gene_axis)), np.float32))
    Xatac = (sp.vstack(blocks_a).tocsr() if blocks_a
             else sp.csr_matrix((n, len(peak_axis)), np.float32))
    Xrna = Xrna.astype(np.float32)
    Xatac = Xatac.astype(np.float32)
    Xrna.sum_duplicates()
    Xatac.sum_duplicates()

    rna = ad.AnnData(Xrna, obs={"batch": batch.copy()}, var={"gene": list(gene_axis)})
    atac = ad.AnnData(Xatac, obs={"batch": batch.copy()}, var={"peak": list(peak_axis)})
    rna.layers["counts"], atac.layers["counts"] = Xrna, Xatac
    return muon.MuData({"rna": rna, "atac": atac}), rows, owner


def _multivi_fit(exp, args, tag, drop_rna_cells=()):
    """Fit and save MultiVI.  tag='pipeline' keeps every cell's RNA."""
    import scvi
    from scvi.model import MULTIVI

    scvi.settings.verbosity = 40
    scvi.settings.seed = args.seed
    gene_axis, peak_axis = exp.gene_axis(), exp.peak_axis()
    mdata, rows, owner = _multivi_mdata(exp, gene_axis, peak_axis, drop_rna_cells)
    MULTIVI.setup_mudata(
        mdata, rna_layer="counts", atac_layer="counts", batch_key="batch",
        modalities={"rna_layer": "rna", "batch_key": "rna", "atac_layer": "atac"})
    model = MULTIVI(mdata, n_genes=len(gene_axis), n_regions=len(peak_axis),
                    n_hidden=args.n_hidden, n_latent=args.n_latent)
    t0 = time.time()
    model.train(max_epochs=args.epochs, early_stopping=True,
                enable_progress_bar=False)
    mdir = exp.out("multivi", "integration",
                   "multivi_model" if tag == "pipeline" else "multivi_model_heldout")
    model.save(mdir, overwrite=True, save_anndata=False)
    log(f"  [{tag}] MultiVI fit: {mdata.n_obs} cells "
        f"({len(rows[exp.anchor])} anchor + {sum(len(rows[a]) for a in exp.atac)} "
        f"ATAC barcodes), {len(gene_axis)} genes, {len(peak_axis)} peaks, "
        f"{args.epochs} max epochs, {time.time() - t0:.0f}s")
    return model, mdata, rows


def _multivi_load(exp, args, tag, drop_rna_cells=()):
    from scvi.model import MULTIVI

    mdir = exp.out("multivi", "integration",
                   "multivi_model" if tag == "pipeline" else "multivi_model_heldout")
    if not os.path.exists(mdir):
        return _multivi_fit(exp, args, tag, drop_rna_cells)
    mdata, rows, _ = _multivi_mdata(exp, exp.gene_axis(), exp.peak_axis(),
                                    drop_rna_cells)
    return MULTIVI.load(mdir, adata=mdata), mdata, rows


def _multivi_decode(model, mdata, rows, n_genes, tag):
    """Decoded expression for `rows` -> genes x cells on the repo's scale.

    px_scale is a within-cell proportion, so log1p(1e4 * px_scale) lands on the
    same footing as measured log1p(CPM/1e4) rows.  n_samples stays 1: scvi 1.4
    mishandles the n_samples>1 concatenation and returns a collapsed array.
    """
    rows = np.asarray(rows, dtype=int)
    out = model.get_normalized_expression(mdata, indices=rows, n_samples=1,
                                        return_mean=True, return_numpy=True)
    out = np.asarray(out, dtype=np.float32)
    if out.shape != (len(rows), n_genes):
        raise SystemExit(f"[{tag}] MultiVI decode gave {out.shape}, expected "
                         f"{(len(rows), n_genes)}")
    return np.log1p(TARGET_SUM * np.clip(out, 0, None)).T.astype(np.float32)


# --------------------------------------------------------------------------- #
# step 1: integrate
# --------------------------------------------------------------------------- #
def step_integrate(exp, method, args):
    odir = exp.out(method, "integration")
    os.makedirs(odir, exist_ok=True)

    if method == "scsaga":
        H, order, sizes = base.integrate(
            [exp.ds[n] for n in exp.order], exp.anchor, odir,
            (args.cfg.get("integration") or {}))
        rows = {n: np.arange(exp.starts[n], exp.starts[n] + exp.sizes[n])
                for n in exp.order}
        info = {"method": "scsaga", "row_space": "one row per dataset cell",
                "params": args.cfg.get("integration") or {}}
    else:
        ok, why = multivi_available()
        if not ok:
            raise SystemExit(
                f"MultiVI unavailable: {why}\n"
                "  build the env with: bash setup_stepwise.sh   (muon + scvi-tools)")
        model, mdata, rows = _multivi_fit(exp, args, "pipeline")
        z = np.asarray(model.get_latent_representation(), dtype=np.float32)
        np.save(f"{odir}/joint_embedding_H.npy", z)
        for key, modality in (("expr", "expression"), ("acc", "accessibility")):
            try:
                np.save(f"{odir}/joint_embedding_H_{key}.npy",
                        np.asarray(model.get_latent_representation(modality=modality),
                                   dtype=np.float32))
            except Exception as exc:                        # noqa: BLE001
                log(f"  [!] modality latent '{modality}' unavailable: {exc}")
        n = z.shape[0]
        for nme in exp.order:
            np.save(f"{odir}/aligned_{nme}.npy", z[rows[nme]])
        info = {"method": "multivi",
                "row_space": "one row per distinct barcode (paired cells collapsed)",
                "n_modality_specific_latents": 2,
                "params": {"epochs": args.epochs, "n_hidden": args.n_hidden,
                           "n_latent": args.n_latent, "seed": args.seed,
                           "n_genes": len(exp.gene_axis()),
                           "n_peaks": len(exp.peak_axis())}}
        if z.shape[0] != n:
            raise SystemExit(f"MultiVI z has {z.shape[0]} rows, expected {n}")

    H = np.load(f"{odir}/joint_embedding_H.npy")
    info.update({"anchor": exp.anchor, "order": exp.order, "sizes": exp.sizes,
                 "H_shape": list(H.shape),
                 "row_index": {k: [int(x) for x in v] for k, v in rows.items()},
                 "pairing_used_in_model": method == "multivi"})
    wjson(f"{odir}/integration_info.json", info)
    wjson(f"{odir}/blocks.json", {k: [int(x) for x in v] for k, v in rows.items()})
    log(f"  H {H.shape} -> {odir}")
    return H


# --------------------------------------------------------------------------- #
# step 2: imputation -- pipeline arm + frozen held-out arm
# --------------------------------------------------------------------------- #
def heldout_split(exp, args):
    """Deterministic frozen split of paired cells, in CELL space (dataset, local)."""
    pairs = {a: paired_map(exp.ds[exp.anchor], exp.ds[a]) for a in exp.atac_queries}
    pairs = {a: v for a, v in pairs.items() if v}
    if not pairs:
        return None
    rng = np.random.default_rng(args.split_seed)
    out = {"frac": args.holdout_frac, "seed": args.split_seed, "datasets": {}}
    for a, pm in pairs.items():
        keys = sorted(pm)
        k = int(round(args.holdout_frac * len(keys)))
        take = np.sort(rng.choice(len(keys), size=max(k, 0), replace=False))
        out["datasets"][a] = {
            "atac_local": [int(keys[i]) for i in take],
            "rna_local": [int(pm[keys[i]]) for i in take],
            "barcodes": [exp.ds[a].load_barcodes()[keys[i]] for i in take],
            "n_paired_with_anchor": int(len(keys))}
    out["n_heldout"] = sum(len(v["atac_local"]) for v in out["datasets"].values())
    out["n_paired"] = sum(v["n_paired_with_anchor"] for v in out["datasets"].values())
    out["rna_local"] = [r for v in out["datasets"].values() for r in v["rna_local"]]
    out["atac_local"] = [j for v in out["datasets"].values() for j in v["atac_local"]]
    return out


def step_impute(exp, method, args, split):
    odir = exp.out(method, "integration")
    sdir = exp.out(method, "pipeline")
    os.makedirs(sdir, exist_ok=True)
    H, rows, info, _ = load_integ(exp, method)
    genes = exp.gene_axis()

    rna_blocks = {n: lognorm_on_axis(exp.ds[n], genes) for n in exp.rna}
    ref_sparse = rna_blocks[exp.anchor].T.tocsr()             # genes x anchor cells
    ref_rows = rows[exp.anchor]
    q_rows = np.concatenate([rows[a] for a in exp.atac_queries])

    # ---- pipeline arm: every ATAC cell, full reference -------------------
    if method == "scsaga":
        imputed, mdist = M.knn_transfer(H[ref_rows], ref_sparse, H[q_rows],
                                        k=K_IMPUTE)
    else:
        model, mdata, _ = _multivi_load(exp, args, "pipeline")
        imputed = _multivi_decode(model, mdata, q_rows, len(genes), "pipeline")
        mdist = float("nan")
        if args.scale_align:
            raw = imputed.copy()
            imputed, sstats = M.align_scale(imputed, ref_sparse)
            log(f"  [pipeline] affine scale alignment: median |delta| "
                f"{np.median(np.abs(imputed - raw)):.4f}, "
                f"{sstats['frac_clipped'] * 100:.1f}% of entries clipped at 0")
    imputed = np.asarray(imputed, dtype=np.float32)
    log(f"  [pipeline] imputed {imputed.shape} genes x ATAC cells"
        + (f", mean nn distance {mdist:.3f}" if np.isfinite(mdist) else ""))
    np.save(f"{sdir}/imputed_expression_genes_x_atac.npy", imputed)
    np.save(f"{sdir}/imputed_expression_cells_x_genes.npy", imputed.T.copy())

    # ---- all-cells matrix: real RNA rows first, then imputed ATAC rows ----
    real = np.vstack([rna_blocks[n].toarray().astype(np.float32) for n in exp.rna])
    all_expr = np.vstack([real, imputed.T]).astype(np.float32)
    np.save(f"{sdir}/all_cells_gene_expression.npy", all_expr)
    mask = np.zeros(all_expr.shape, dtype=bool)
    mask[:real.shape[0], :] = True
    np.save(f"{sdir}/mask.npy", mask)
    with open(f"{sdir}/genes.txt", "w") as fh:
        fh.write("\n".join(genes) + "\n")
    with open(f"{sdir}/norm.txt", "w") as fh:
        fh.write(NORM_TXT + "\n")
    all_bc = [b for n in exp.rna for b in exp.ds[n].load_barcodes()]
    all_bc += [b for q in exp.atac_queries for b in exp.ds[q].load_barcodes()]
    with open(f"{sdir}/all_cells_barcodes.txt", "w") as fh:
        fh.write("\n".join(all_bc) + "\n")
    with open(f"{sdir}/n_cells.txt", "w") as fh:
        fh.write(str(all_expr.shape[0]))
    wjson(f"{sdir}/impute_info.json", {
        "method": method, "genes": len(genes),
        "rows_real_rna": int(real.shape[0]), "rows_imputed": int(imputed.shape[1]),
        "k": K_IMPUTE if method == "scsaga" else None,
        "mean_nn_distance": None if not np.isfinite(mdist) else float(mdist),
        "scale_align": bool(args.scale_align) and method == "multivi",
        "note": ("paired cells are imputed too, so their row is a reconstruction of "
                 "a measured cell, not a guess -- it duplicates the real RNA row "
                 "exactly as run.py's pipeline does")})
    log(f"  [pipeline] all-cells {all_expr.shape} "
        f"({real.shape[0]} real + {imputed.shape[1]} imputed) -> {sdir}")

    # ---- frozen held-out arm --------------------------------------------
    if split is None:
        log("  [heldout] no barcodes shared with the anchor -> step-1/2 "
            "answer-key metrics will be recorded as unavailable")
        return
    hdir = exp.out(method, "heldout")
    os.makedirs(hdir, exist_ok=True)
    wjson(f"{hdir}/split.json", {"heldout": split})
    r_local = np.asarray(split["rna_local"], dtype=int)
    q_h = np.concatenate([rows[a][np.asarray(split["datasets"][a]["atac_local"], int)]
                          for a in split["datasets"]])
    ref_keep = np.setdiff1d(np.arange(exp.sizes[exp.anchor]), r_local)

    if method == "scsaga":
        imp_h, _ = M.knn_transfer(H[ref_rows][ref_keep], ref_sparse[:, ref_keep],
                                  H[q_h], k=K_IMPUTE)
    else:
        model, mdata, _ = _multivi_load(exp, args, "heldout",
                                        drop_rna_cells=r_local)
        imp_h = _multivi_decode(model, mdata, q_h, len(genes), "heldout")
        if args.scale_align:
            imp_h, _ = M.align_scale(imp_h, ref_sparse[:, ref_keep])
    truth_h = np.asarray(ref_sparse[:, r_local].todense(), dtype=np.float32)
    np.save(f"{hdir}/imputed_heldout_genes_x_cells.npy",
            np.asarray(imp_h, dtype=np.float32))
    np.save(f"{hdir}/truth_heldout_genes_x_cells.npy", truth_h)
    log(f"  [heldout] {len(q_h)} frozen cells imputed from {len(ref_keep)} reference "
        f"cells" + ("; their RNA rows are excluded from the reference"
                    if method == "scsaga" else
                    "; their RNA rows are removed from the model's training data too"))


# --------------------------------------------------------------------------- #
# step 3: GRN
# --------------------------------------------------------------------------- #
def step_grn(exp, cfg, method, args):
    sdir = exp.out(method, "pipeline")
    gdir = f"{sdir}/grn"
    all_path = f"{sdir}/all_cells_gene_expression.npy"
    if not os.path.exists(all_path):
        raise SystemExit(f"no {all_path}; run --only impute first")
    genes = base.read_lines(f"{sdir}/genes.txt")
    g = cfg.get("grn") or {}
    reg = base.abspath(cfg["_root"], g.get("regulators") or "data/tf_only.txt")
    tgt = base.abspath(cfg["_root"], g.get("targets") or "data/trrust_tf.txt")

    seeds = [int(s) for s in str(args.seeds).split(",") if s.strip()]
    for i, seed in enumerate(seeds):
        out = gdir if i == 0 else f"{gdir}/seed{seed}"
        os.makedirs(out, exist_ok=True)
        if os.path.exists(f"{out}/network.tsv") and not args.force:
            log(f"  [grn] seed {seed}: network.tsv exists, skipping (--force to redo)")
            continue
        log(f"  [grn] GRNBoost2 seed {seed} -> {out}")
        base.run_grn(all_path, genes, reg, tgt, out,
                     n_workers=args.workers,
                     threads_per_worker=args.threads_per_worker,
                     seed=seed, scheduler=args.scheduler,
                     max_columns=args.max_columns)
        if i == 0:
            gt = {k: base.abspath(cfg["_root"], v)
                  for k, v in (cfg.get("ground_truth") or {}).items()}
            base.evaluate(f"{out}/network.tsv", gt, f"{out}/evaluation")


# --------------------------------------------------------------------------- #
# step 4: metrics
# --------------------------------------------------------------------------- #
def h_row_labels(exp, rows):
    """Cell-type labels for one method's H rows ('' where unknown)."""
    lab = np.array(["" for _ in range(max(rows[n].max() for n in rows) + 1)],
                   dtype=object)
    for n in exp.order:
        L = load_labels(exp.ds[n])
        if L is None:
            continue
        L = np.asarray(L, dtype=object)
        for i in range(exp.sizes[n]):
            r = rows[n][i]
            if lab[r] == "":
                lab[r] = L[i]
    return lab


def metrics_integration(exp, method, args, H, rows, extras, split, labels_row):
    out = {"method": method, "H_shape": list(H.shape),
           "row_space": ("collapsed barcodes" if method == "multivi"
                         else "one row per dataset cell")}

    # ---- FOSCTTM: grades the embedding, so it uses ALL paired cells ------
    if split is not None:
        a = next(iter(split["datasets"]))
        pm = paired_map(exp.ds[exp.anchor], exp.ds[a])
        keys = sorted(pm)
        r_rna = rows[exp.anchor][np.asarray([pm[j] for j in keys], int)]
        r_atac = rows[a][np.asarray(keys, int)]
        same = np.array_equal(r_rna, r_atac)
        variant, A, B = None, None, None
        if same and {"expr", "acc"} <= set(extras):
            variant = "z_expr vs z_acc of the same co-assayed cell"
            A, B = extras["expr"][r_rna], extras["acc"][r_rna]
        elif not same:
            variant = "RNA row vs ATAC row of the same cell"
            A, B = H[r_rna], H[r_atac]
        else:
            out["foscttm"] = {
                "status": "unavailable",
                "reason": "paired cells share one latent row and the model exported "
                          "no modality-specific latents, so no cross-modal pair exists"}
        if A is not None:
            order = np.arange(len(keys))
            val, f_a, f_b = M.foscttm(A, B, order, order)
            ctrl = M.foscttm_shuffled(A, B, order, order, seed=args.seed)
            out["foscttm"] = {
                "status": "ok", "foscttm": val, "shuffled_control": ctrl,
                "variant": variant, "n_pairs": int(len(keys)),
                "atac_dataset": a,
                "median_a_to_b": float(np.median(f_a)),
                "median_b_to_a": float(np.median(f_b)),
                "reference_range": "0.14-0.16 MMD-MA, 0.22 UnionCom, 0.5 random",
                "gate": "pass" if (val <= 0.22 and val < 0.75 * ctrl) else "fail"}
    else:
        out["foscttm"] = {"status": "unavailable",
                          "reason": "no ATAC dataset shares barcodes with the anchor"}

    # ---- label conservation + batch correction (all H rows) -------------
    keep = np.where(labels_row != "")[0] if labels_row is not None else np.array([], int)
    if len(keep) >= 30 and len(set(labels_row[keep])) >= 2:
        out["label_conservation"] = M.label_conservation(H[keep], labels_row[keep])
        rng = np.random.default_rng(args.seed)
        out["label_conservation"]["shuffled_control"] = M.label_conservation(
            H[keep], rng.permutation(labels_row[keep]))
        out["label_conservation"]["n_cells"] = int(len(keep))
    else:
        out["label_conservation"] = {
            "status": "unavailable", "n_cells_with_labels": int(len(keep)),
            "reason": ("needs labels.txt with >=2 cell types on >=30 cells; put "
                       "labels.txt next to pca_50.txt (one label per PCA row)")}

    # batch = the dataset a row came from (for a collapsed cell: the first one)
    batch = np.empty(H.shape[0], dtype=object)
    batch[:] = ""
    for n in exp.order:
        for r in rows[n]:
            if batch[r] == "":
                batch[r] = n
    out["batch_correction"] = M.batch_correction(H, batch)

    # ---- label transfer: learned on RNA -> ATAC half of the same cell ----
    if split is not None and len(split["atac_local"]):
        a = next(iter(split["datasets"]))
        lab_anchor = load_labels(exp.ds[exp.anchor])
        if lab_anchor is not None:
            r_local = np.asarray(split["datasets"][a]["rna_local"], int)
            q_local = np.asarray(split["datasets"][a]["atac_local"], int)
            ref_keep = np.setdiff1d(np.arange(exp.sizes[exp.anchor]), r_local)
            L = np.asarray(lab_anchor, dtype=object)
            out["label_transfer"] = M.label_transfer(
                H[rows[exp.anchor][ref_keep]], H[rows[a][q_local]],
                L[ref_keep], L[r_local])
            out["label_transfer"].update({
                "truth_from": "paired RNA cell's label",
                "scored_on": "frozen held-out cells only",
                "gate": "pass" if out["label_transfer"]["accuracy"] >= 0.85 else "fail"})
        else:
            out["label_transfer"] = {"status": "unavailable",
                                     "reason": "anchor RNA dataset has no labels.txt"}
    else:
        out["label_transfer"] = {"status": "unavailable", "reason": "needs paired cells"}
    return out


def metrics_imputation(exp, method, args, labels_row):
    hdir = exp.out(method, "heldout")
    pp, tp = (f"{hdir}/imputed_heldout_genes_x_cells.npy",
              f"{hdir}/truth_heldout_genes_x_cells.npy")
    if not (os.path.exists(pp) and os.path.exists(tp)):
        return {"status": "unavailable",
                "reason": "no held-out predictions: needs an ATAC dataset whose "
                          "barcodes match the anchor's, then --only impute"}
    pred, truth = np.load(pp), np.load(tp)
    hvg = M.hvg_mask(truth, n=1000)
    all_path = exp.out(method, "pipeline", "all_cells_gene_expression.npy")
    measured = None
    if os.path.exists(all_path):
        measured = np.asarray(np.load(all_path, mmap_mode="r")
                              [:exp.n_real_rna()]).T

    out = {"n_cells": int(truth.shape[1]), "n_genes": int(truth.shape[0]),
           "holdout_frac": json.load(open(f"{hdir}/split.json"))["heldout"]["frac"]}
    out.update(M.imputation_mse(pred, truth, hvg=hvg))
    out.update(M.gene_wise_rho(pred, truth, hvg=hvg))
    out.update(M.pfcs(pred, truth, n_hvg=100))
    out.update(M.densification(pred, truth, measured if measured is not None else truth))

    split = json.load(open(f"{hdir}/split.json"))["heldout"]
    lab = None
    if labels_row is not None:
        _, rows, _, _ = load_integ(exp, method)
        r_rows = rows[exp.anchor][np.asarray(split["rna_local"], int)]
        cand = np.asarray([labels_row[r] for r in r_rows], dtype=object)
        if not (cand == "").any() and len(set(cand)) >= 2:
            lab = cand
    if lab is not None:
        out.update(M.pdes(pred, truth, lab))
        out["downstream"] = M.downstream(pred, truth, lab)
    else:
        why = "needs >=2 cell types on the held-out cells (labels.txt)"
        out["pdes"] = {"status": "unavailable", "reason": why}
        out["downstream"] = {"status": "unavailable", "reason": why}

    out["shuffled_reference_control"] = _shuffled_control(exp, method, truth, hvg, args)
    ctrl = out.get("shuffled_reference_control") or {}
    gc = ctrl.get("gene_permuted_control") or {}
    # pFCS is gated against the GENE-permuted null, not the cell-permuted one: a
    # cell-permuted reference makes every query average a random subset, i.e. a
    # near-constant block, and a near-constant block can score high on pFCS from the
    # mean profile alone (the design doc's own warning about blurred matrices).
    out["gate"] = {
        "pfcs_above_gene_null": (
            "pass" if (np.isfinite(out.get("pfcs", np.nan))
                       and out["pfcs"] > gc.get("pfcs", 1.0)) else "fail"),
        "mse_below_control": _cmp(out, "shuffled_reference_control", "mse", "<"),
        "rho_in_published_range": _in_range(out.get("rho_mean"), 0.10, 0.40),
        "no_densification": (out.get("densification_ratio_vs_measured") is not None
                             and abs(out["densification_ratio_vs_measured"] - 1) < 0.5)}
    return out


def _cmp(d, ctrl_key, key, op):
    a = d.get(key)
    b = (d.get(ctrl_key) or {}).get(key)
    if a is None or b is None or not np.isfinite(a) or not np.isfinite(b):
        return "n/a"
    return "pass" if (a > b if op == ">" else a < b) else "fail"


def _in_range(x, lo, hi):
    if x is None or not np.isfinite(x):
        return "n/a"
    return "pass" if lo <= x <= hi else "fail"


def _shuffled_control(exp, method, truth, hvg, args):
    """Required negative control: same measurement with the reference permuted."""
    H, rows, _, _ = load_integ(exp, method)
    real = np.asarray(np.load(exp.out(method, "pipeline",
                                      "all_cells_gene_expression.npy"),
                              mmap_mode="r")[:exp.n_real_rna()])
    ref = np.ascontiguousarray(real[:exp.sizes[exp.anchor]].T)
    q_rows = np.concatenate([rows[a] for a in exp.atac_queries])
    rng = np.random.default_rng(args.seed + 1)
    perm = rng.permutation(exp.sizes[exp.anchor])
    imp, _ = M.knn_transfer(H[rows[exp.anchor]][perm], ref[:, perm], H[q_rows],
                            k=K_IMPUTE)
    split = json.load(open(exp.out(method, "heldout", "split.json")))["heldout"]
    q_h = np.concatenate([rows[a][np.asarray(split["datasets"][a]["atac_local"], int)]
                          for a in split["datasets"]])
    pos = {int(r): i for i, r in enumerate(q_rows)}
    sub = imp[:, np.asarray([pos[int(r)] for r in q_h], int)]
    out = {"note": "reference cells permuted with their expression; same k, same rows"}
    out.update(M.imputation_mse(sub, truth, hvg=hvg))
    out.update(M.gene_wise_rho(sub, truth, hvg=hvg))
    out.update(M.pfcs(sub, truth, n_hvg=100))

    # Second, sharper control for pFCS: permute the GENE rows of the reference.
    # Cell permutation alone produces a near-constant guess (every query averages a
    # random subset), and a near-constant block can score well on pFCS purely from
    # the mean profile.  Shuffling gene labels leaves every marginal intact and
    # destroys gene-gene structure, which is exactly what pFCS claims to measure.
    gp = rng.permutation(ref.shape[0])
    imp_g, _ = M.knn_transfer(H[rows[exp.anchor]], ref[gp], H[q_rows], k=K_IMPUTE)
    out_g = {"note": "reference GENE labels permuted; same k, same cells"}
    out_g.update(M.imputation_mse(imp_g[:, np.asarray([pos[int(r)] for r in q_h], int)],
                                  truth, hvg=hvg))
    out_g.update(M.gene_wise_rho(imp_g[:, np.asarray([pos[int(r)] for r in q_h], int)],
                                 truth, hvg=hvg))
    out_g.update(M.pfcs(imp_g[:, np.asarray([pos[int(r)] for r in q_h], int)],
                        truth, n_hvg=100))
    out["gene_permuted_control"] = out_g
    return out


def metrics_grn(exp, cfg, method, args, labels_row):
    import pandas as pd

    gdir = exp.out(method, "pipeline", "grn")
    main_path = f"{gdir}/network.tsv"
    if not os.path.exists(main_path):
        return {"status": "unavailable", "reason": f"no {main_path}"}
    seeds = [int(s) for s in str(args.seeds).split(",") if s.strip()]
    paths = [p for p in [main_path] + [f"{gdir}/seed{s}/network.tsv" for s in seeds[1:]]
             if os.path.exists(p)]
    nets = [pd.read_csv(p, sep="\t") for p in paths]
    main = nets[0]

    out = {"edges": int(len(main)), "seeds": seeds[:len(paths)],
           "stability": M.grn_stability(nets, top=args.stability_top),
           "vs_ground_truth": {}}
    truths = {}
    ev_path = f"{gdir}/evaluation/evaluation.json"
    ev = (json.load(open(ev_path)).get("ground_truths", {})
          if os.path.exists(ev_path) else {})
    for name, path in (cfg.get("ground_truth") or {}).items():
        full = base.abspath(cfg["_root"], path)
        if not os.path.exists(full):
            continue
        truths[name] = M.read_truth(full)
        r = M.grn_recovery(main, truths[name], cap=args.grn_cap)
        r["run_py_recovered"] = ev.get(name)
        out["vs_ground_truth"][name] = r
    if truths:
        tname = next(iter(truths))
        out["topology"] = M.grn_topology(main, truths[tname], n_hub=20,
                                        cap=args.grn_cap)
        out["topology"]["truth_used"] = tname
    else:
        out["topology"] = {"status": "unavailable", "reason": "no ground truth file"}
    out["perturbation"] = {
        "status": "not computed",
        "reason": "needs Perturb-seq / knockout data for these cells; none exists "
                  "for PBMC multiome. Recorded, not approximated."}

    # ---- regulon activity per cell type --------------------------------
    if labels_row is None or (labels_row == "").all():
        out["regulon_activity"] = {"status": "unavailable", "reason": "needs labels.txt"}
        return out
    sdir = exp.out(method, "pipeline")
    genes = base.read_lines(f"{sdir}/genes.txt")
    # all-cells rows: RNA cells (exp.rna order) then ATAC query cells
    labels = []
    for n in exp.rna:
        L = load_labels(exp.ds[n])
        labels += list(L) if L is not None else [""] * exp.sizes[n]
    lab_anchor = load_labels(exp.ds[exp.anchor])
    rna_pos = {b: i for n in exp.rna for i, b in enumerate(exp.ds[n].load_barcodes())}
    for q in exp.atac_queries:
        for b in exp.ds[q].load_barcodes():
            i = rna_pos.get(b)
            labels.append(lab_anchor[i] if (lab_anchor is not None and i is not None)
                          else "")
    y = np.asarray(labels, dtype=object)
    ok = y != ""
    if ok.sum() < 30 or len(set(y[ok])) < 2:
        out["regulon_activity"] = {"status": "unavailable",
                                   "reason": "needs labels on >=2 cell types"}
        return out
    X = np.asarray(np.load(f"{sdir}/all_cells_gene_expression.npy")[:len(y)][ok])
    out["regulon_activity"] = M.regulon_activity(
        X, genes, main, y[ok], top_edges=args.grn_cap, n_regulons=15)
    out["regulon_activity"]["n_cells_labelled"] = int(ok.sum())
    return out


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #
def _fmt(x, n=3):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, dict):
        return ", ".join(f"{k}:{v}" for k, v in list(x.items())[:6]) or "{}"
    if isinstance(x, str):
        return x
    if isinstance(x, (int, np.integer)):
        return f"{int(x)}"
    try:
        if not np.isfinite(x):
            return "n/a"
    except TypeError:
        return str(x)
    return f"{x:.{n}f}"


def write_report(exp, methods, payload):
    def get(m, keys):
        cur = payload.get(m, {})
        for k in keys:
            if not isinstance(cur, dict) or k not in cur:
                return None
            cur = cur[k]
        return cur

    lines = [f"# {exp.name}: step-wise measurement, " + " vs ".join(methods), "",
             f"datasets: {', '.join(exp.names)}  |  anchor: {exp.anchor}  |  "
             f"ATAC queries: {', '.join(exp.atac_queries)}", ""]
    variants = {m: (get(m, ("integration", "foscttm", "variant")) or "n/a")
                for m in methods}
    if len(set(variants.values())) > 1 or "n/a" not in variants.values():
        lines += ["FOSCTTM variant: " + "; ".join(f"{m} = {v}"
                                                  for m, v in variants.items()), ""]

    lines += ["| step | metric | " + " | ".join(methods) + " | what good looks like |",
              "|---|---|" + "---|" * len(methods) + "---|"]

    def row(step, label, keys, good, nd=3):
        lines.append(f"| {step} | {label} | "
                     + " | ".join(_fmt(get(m, keys), nd) for m in methods)
                     + f" | {good} |")

    for lbl, keys, good in [
            ("FOSCTTM", ("integration", "foscttm", "foscttm"), "0.14-0.22; 0.5 = random"),
            ("FOSCTTM, shuffled pairing", ("integration", "foscttm", "shuffled_control"), "≈0.5"),
            ("ARI (labels vs kmeans of H)", ("integration", "label_conservation", "ari"), "high, above control"),
            ("NMI", ("integration", "label_conservation", "nmi"), "high, above control"),
            ("isolated-label F1", ("integration", "label_conservation", "isolated_label_F1_macro"), "> 0.8"),
            ("graph connectivity", ("integration", "label_conservation", "graph_connectivity"), "> 0.8"),
            ("control ARI (labels shuffled)", ("integration", "label_conservation", "shuffled_control", "ari"), "≈0"),
            ("iLISI", ("integration", "batch_correction", "ilisi_mean"), "→ 1 (mixed)"),
            ("kBET-style acceptance", ("integration", "batch_correction", "kbet_acceptance"), "→ 0.95"),
            ("batch |ASW|", ("integration", "batch_correction", "batch_asw_abs"), "→ 0"),
            ("label transfer accuracy", ("integration", "label_transfer", "accuracy"), "> 0.85")]:
        row("1", lbl, keys, good)

    for lbl, keys, good in [
            ("MSE", ("imputation", "mse"), "low, below control"),
            ("MSE (1000 HVG)", ("imputation", "mse_hvg"), "low"),
            ("pFCS (top-100 HVG correlation)", ("imputation", "pfcs"), "high, above control"),
            ("control pFCS (shuffled reference)", ("imputation", "shuffled_reference_control", "pfcs"), "clearly below pFCS"),
            ("control MSE", ("imputation", "shuffled_reference_control", "mse"), "clearly above MSE"),
            ("control pFCS (gene labels shuffled)", ("imputation", "shuffled_reference_control", "gene_permuted_control", "pfcs"), "the null pFCS: must be below pFCS"),
            ("pDES (Welch t, not limma)", ("imputation", "pdes"), "high"),
            ("gene-wise rho, mean", ("imputation", "rho_mean"), "0.155-0.22; >0.4 = leak"),
            ("gene-wise rho, median", ("imputation", "rho_median"), "0.155-0.22"),
            ("control gene-wise rho", ("imputation", "shuffled_reference_control", "rho_mean"), "≈0"),
            ("densification vs measured rows", ("imputation", "densification_ratio_vs_measured"), "≈1, not inflated"),
            ("gene mean shift (median, sd units)", ("imputation", "gene_mean_shift_median_sd"), "small"),
            ("clf accuracy on imputed", ("imputation", "downstream", "imputed", "clf_accuracy"), "close to truth row"),
            ("clf macro-F1 on imputed", ("imputation", "downstream", "imputed", "clf_macro_f1"), "close to truth row"),
            ("clustering ARI on imputed", ("imputation", "downstream", "imputed", "clust_ari"), "close to truth row"),
            ("clf accuracy on truth", ("imputation", "downstream", "truth", "clf_accuracy"), "ceiling"),
            ("clustering ARI on truth", ("imputation", "downstream", "truth", "clust_ari"), "ceiling")]:
        row("2", lbl, keys, good)

    gt_names = []
    for m in methods:
        gt_names += list((get(m, ("grn", "vs_ground_truth")) or {}).keys())
    for name in dict.fromkeys(gt_names):
        row("3", f"AUPRC vs {name}", ("grn", "vs_ground_truth", name, "auprc"), "> prevalence")
        row("3", f"random baseline (prevalence) vs {name}", ("grn", "vs_ground_truth", name, "prevalence_random_auprc"), "reference")
        row("3", f"early precision vs {name}", ("grn", "vs_ground_truth", name, "early_precision"), "high")
        row("3", f"precision@10000 vs {name}", ("grn", "vs_ground_truth", name, "precision_at", "10000"), "high")
        row("3", f"recall@10000 vs {name}", ("grn", "vs_ground_truth", name, "recall_at", "10000"), "high")
        row("3", f"recovered vs {name} (run.py)", ("grn", "vs_ground_truth", name, "recovered"), "run.py number")

    for lbl, keys, good in [
            ("edge stability, mean Jaccard (top-10k, seeds)", ("grn", "stability", "jaccard_mean"), "differences must exceed this"),
            ("edges in every seed", ("grn", "stability", "edges_in_all_seeds"), "reported"),
            ("hub overlap vs truth", ("grn", "topology", "hub_overlap_at_k"), "high"),
            ("out-degree Spearman vs truth", ("grn", "topology", "spearman_out_degree_vs_truth"), "high"),
            ("symmetric edge fraction", ("grn", "topology", "symmetric_edge_fraction"), "low; high = smoothing"),
            ("regulon specificity, mean", ("grn", "regulon_activity", "specificity_mean"), "cell-type specific"),
            ("regulons peaking per cell type", ("grn", "regulon_activity", "regulons_peaking_in"), "matches lineages")]:
        row("3", lbl, keys, good)

    lines += ["", "The `truth` rows are the hidden measured RNA, which carries its own "
                  "dropout: a denoising imputer can beat them on a downstream task, so "
                  "read them as a reference point rather than a hard ceiling.", "",
              "Two negative controls are reported for step 2.  `shuffled reference` "
                  "permutes the reference cells (the design doc's control, and the one "
                  "MSE is gated on); it makes every guess a near-constant average. "
                  "`gene labels shuffled` permutes the reference's gene rows, which "
                  "leaves every marginal intact and destroys only the gene-gene "
                  "structure -- that is the null pFCS is gated on, because a "
                  "near-constant block can otherwise score well on pFCS.", "",
              "Gates (design doc sections 3-5):", "",
              "| gate | " + " | ".join(methods) + " |", "|---|" + "---|" * len(methods) + "|"]
    for gname, keys in (("step 1 FOSCTTM", ("integration", "foscttm", "gate")),
                        ("step 1 label transfer", ("integration", "label_transfer", "gate")),
                        ("step 2 pFCS > gene-permuted null", ("imputation", "gate", "pfcs_above_gene_null")),
                        ("step 2 MSE < control", ("imputation", "gate", "mse_below_control")),
                        ("step 2 rho in range", ("imputation", "gate", "rho_in_published_range")),
                        ("step 2 no densification", ("imputation", "gate", "no_densification"))):
        lines.append(f"| {gname} | " + " | ".join(_fmt(get(m, keys)) for m in methods) + " |")

    lines += ["", "`n/a` means the metric is unavailable for this method, or the input "
                  "it needs is absent.  metrics.json records the reason for every "
                  "one of them.  Nothing is approximated silently.", ""]

    for m in methods:
        p = payload.get(m, {})
        if "error" in p:
            lines += [f"## {m}", "", f"FAILED: {p['error']}", ""]
            continue
        imp, grn = p.get("imputation", {}), p.get("grn", {})
        lines += [f"## {m}", ""]
        if isinstance(imp, dict) and "densification_ratio_vs_truth" in imp:
            lines.append(f"- imputed row density vs hidden answer: "
                         f"{_fmt(imp['densification_ratio_vs_truth'])}x")
        if isinstance(imp, dict) and "frac_genes_mean_shift_gt_0.5sd" in imp:
            lines.append(f"- genes whose imputed mean is off by >0.5 sd: "
                         f"{100 * imp['frac_genes_mean_shift_gt_0.5sd']:.1f}%")
        st = (grn or {}).get("stability") or {}
        if st.get("pairwise"):
            lines.append(f"- seed-to-seed top-{st.get('top')} Jaccard: "
                         + ", ".join(_fmt(x) for x in st["pairwise"]))
        ra = (grn or {}).get("regulon_activity") or {}
        if ra.get("regulons"):
            lines += ["", "| regulon | targets | peak cell type | specificity |",
                      "|---|---|---|---|"]
            for r in ra["regulons"][:5]:
                lines.append(f"| {r['tf']} | {r['n_targets']} | "
                             f"{r['peak_cell_type']} | {_fmt(r['specificity'])} |")
        lines.append("")

    path = os.path.join(exp.root, "report.md")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"  report -> {path}")


# --------------------------------------------------------------------------- #
# check
# --------------------------------------------------------------------------- #
def check(exp, datasets, methods, args):
    ok = True
    log(f"experiment {exp.name}: {' + '.join(exp.names)}  anchor={exp.anchor}  "
        f"queries={exp.atac_queries}")
    bad = [n for n in exp.names if datasets[n].missing()]
    if bad:
        log(f"  MISSING run.py inputs: {bad}")
        ok = False
    pm = {a: len(paired_map(datasets[exp.anchor], datasets[a])) for a in exp.atac_queries}
    log(f"  barcodes shared with the anchor: {pm}")
    log("  -> " + ("held-out answer key available (FOSCTTM, MSE, pFCS, pDES, rho)"
                   if any(pm.values()) else
                   "no pairing: FOSCTTM and the held-out block will be 'unavailable'"))
    for n in exp.names:
        d = datasets[n]
        if not (d.has_expression() and d.has_pca()):
            continue
        state, got, want = counts_orientation(d)
        if state == "ok":
            log(f"  {n}: counts {got[0]} x {got[1]} = cells x features OK")
        elif state == "transposed":
            log(f"  {n}: counts {got[0]} x {got[1]} = FEATURES x CELLS (expected "
                f"{want[0]} x {want[1]}); stepwise_eval transposes it in memory, "
                f"`python fix_counts.py` repairs the file")
        else:
            log(f"  {n}: counts {got[0]} x {got[1]} matches neither {want[0]} x "
                f"{want[1]} nor transposed -- check barcodes/features")
            ok = False
    labs = {n: load_labels(datasets[n]) is not None for n in exp.names}
    log(f"  labels.txt present: {labs}")
    if not any(labs.values()):
        log("  -> label metrics will be 'unavailable' (put labels.txt next to pca_50.txt)")
    for m in methods:
        if m == "multivi":
            avail, why = multivi_available()
            log(f"  MultiVI backend: {why}")
            ok = ok and avail
            for a in exp.atac:
                if not datasets[a].has_expression():
                    log(f"  MultiVI needs ATAC peaks: {a} has no "
                        f"counts.mtx/barcodes.txt/features.txt")
                    ok = False
    for key in ("regulators", "targets"):
        p = base.abspath(args.cfg["_root"],
                         (args.cfg.get("grn") or {}).get(key) or f"data/{key}.txt")
        log(f"  GRN {key}: {p} {'OK' if os.path.exists(p) else 'MISSING'}")
        ok = ok and os.path.exists(p)
    return ok


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(ROOT, "config.yml"))
    ap.add_argument("-e", "--experiment", action="append", default=[])
    ap.add_argument("--methods", default="multivi,scsaga")
    ap.add_argument("--only", default=",".join(STEPS))
    ap.add_argument("--out", default=None,
                    help="results root (default <repo>/results/stepwise)")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--force", action="store_true", help="redo existing networks")
    ap.add_argument("--workers", type=int,
                    default=int(os.environ.get("GRN_WORKERS", "0")) or None)
    ap.add_argument("--threads-per-worker", type=int,
                    default=int(os.environ.get("GRN_THREADS_PER_WORKER", "1")))
    ap.add_argument("--scheduler", default=os.environ.get("GRN_SCHEDULER"))
    ap.add_argument("--max-columns", type=int,
                    default=int(os.environ.get("GRN_MAX_COLUMNS", "0")) or None)
    ap.add_argument("--seed", type=int, default=666)
    ap.add_argument("--seeds", default=os.environ.get("GRN_SEEDS", "666,667,668"),
                    help="comma list; the first seed owns network.tsv")
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--n-hidden", type=int, default=128)
    ap.add_argument("--n-latent", type=int, default=30,
                    help="latent width; 30 matches scSAGA's embedding width")
    ap.add_argument("--holdout-frac", type=float, default=0.2,
                    help="fraction of paired cells frozen as the answer key")
    ap.add_argument("--split-seed", type=int, default=0)
    ap.add_argument("--scale-align", dest="scale_align", action="store_true",
                    default=True,
                    help="per-gene affine alignment of MultiVI decodes to the "
                         "measured scale (default ON)")
    ap.add_argument("--no-scale-align", dest="scale_align", action="store_false")
    ap.add_argument("--grn-cap", type=int, default=50000,
                    help="score only the top-N edges (geneRNIB/BEELINE practice)")
    ap.add_argument("--stability-top", type=int, default=10000)
    ap.add_argument("--skip-metrics", nargs="*", default=[],
                    help="metric groups to skip: integration imputation grn")
    args = ap.parse_args()

    cfg_path = os.path.abspath(args.config)
    cfg = yaml.safe_load(open(cfg_path))
    root = os.path.abspath(os.environ.get("PIPELINE_ROOT") or os.path.dirname(cfg_path))
    cfg["_root"] = root
    args.cfg = cfg
    datasets = {n: base.Dataset(n, s or {}, root)
                for n, s in (cfg.get("datasets") or {}).items()}
    for n, s in (cfg.get("datasets") or {}).items():
        if (s or {}).get("labels"):
            datasets[n].labels = base.abspath(root, s["labels"])

    exps = cfg.get("experiments") or {}
    names = args.experiment or list(exps)
    unknown = [n for n in names if n not in exps]
    if unknown:
        raise SystemExit(f"unknown experiment(s) {unknown}; have {list(exps)}")
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    bad = [m for m in methods if m not in METHODS]
    if bad:
        raise SystemExit(f"unknown method(s) {bad}; choose from {METHODS}")
    steps = [s.strip() for s in args.only.split(",") if s.strip()]
    bad = [s for s in steps if s not in STEPS]
    if bad:
        raise SystemExit(f"unknown step(s) {bad}; choose from {STEPS}")

    out_root = os.path.abspath(args.out or os.path.join(root, "results", "stepwise"))
    if args.workers is None:
        import multiprocessing
        args.workers = int(os.environ.get("SLURM_CPUS_PER_TASK")
                           or min(multiprocessing.cpu_count(), 8))

    if args.check:
        ok = True
        for n in names:
            e = Exp(cfg, datasets, n)
            e.root = os.path.join(out_root, n)
            ok = check(e, datasets, methods, args) and ok
        return 0 if ok else 1

    log(f"config  {cfg_path}\nroot    {root}\nout     {out_root}\n"
        f"exps    {names}\nmethods {methods}\nsteps   {steps}\n"
        f"workers {args.workers}   seeds {args.seeds}")

    for name in names:
        exp = Exp(cfg, datasets, name)
        exp.root = os.path.join(out_root, name)
        os.makedirs(exp.root, exist_ok=True)
        log(f"\n{'=' * 74}\nEXPERIMENT {name}\n{'=' * 74}")
        split = heldout_split(exp, args)
        if split is not None:
            log(f"  answer key: {split['n_heldout']} of {split['n_paired']} paired "
                f"cells frozen ({split['frac']:.0%})")

        payload = {}
        for method in methods:
            log(f"\n--- {method}")
            t0 = time.time()
            try:
                if "integrate" in steps:
                    step_integrate(exp, method, args)
                if "impute" in steps:
                    step_impute(exp, method, args, split)
                if "grn" in steps:
                    step_grn(exp, cfg, method, args)
            except (Exception, SystemExit) as exc:                    # noqa: BLE001
                log(f"  [!] {method} FAILED: {exc}")
                payload[method] = {"error": str(exc)}
                continue
            log(f"  {method}: compute steps done in {time.time() - t0:.0f}s")

            if "metrics" in steps:
                H, rows, info, extras = load_integ(exp, method)
                labels_row = h_row_labels(exp, rows)
                pay = {"meta": {"method": method, "seeds": args.seeds,
                                "grn_cap": args.grn_cap,
                                "row_space": info.get("row_space"),
                                "labels": "labels.txt" if (labels_row != "").any()
                                          else "none"}}
                for group, fn in (("integration",
                                   lambda: metrics_integration(exp, method, args, H,
                                                               rows, extras, split,
                                                               labels_row)),
                                  ("imputation",
                                   lambda: metrics_imputation(exp, method, args,
                                                              labels_row)),
                                  ("grn",
                                   lambda: metrics_grn(exp, cfg, method, args,
                                                       labels_row))):
                    if group in args.skip_metrics:
                        continue
                    try:
                        pay[group] = fn()
                    except Exception as exc:                          # noqa: BLE001
                        pay[group] = {"error": f"{type(exc).__name__}: {exc}"}
                        log(f"  [!] {group} metrics failed: {exc}")
                wjson(exp.out(method, "metrics.json"), pay)
                payload[method] = pay
                log(f"  metrics -> {exp.out(method, 'metrics.json')}")

        # merge with whatever the other method(s) already wrote, so `--only metrics
        # --methods multivi` does not drop scsaga's column from the report
        for m in METHODS:                      # whatever exists on disk, in order
            p = exp.out(m, "metrics.json")
            if m not in payload and os.path.exists(p):
                try:
                    payload[m] = json.load(open(p))
                except Exception:                                     # noqa: BLE001
                    pass
        if payload:
            write_report(exp, [m for m in METHODS if m in payload], payload)
            wjson(os.path.join(exp.root, "metrics.json"), payload)
    log(f"\nDONE. results under {out_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
