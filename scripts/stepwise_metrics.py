#!/usr/bin/env python
"""Metrics for the step-wise multiome -> GRN experiment.

Every function here is pure (arrays in, numbers out) so it can be tested without
data.  Each metric names the workflow slide it implements:

  step 1 integration   -- slides 4 (workflow 1) and 11 (workflow 2)
      FOSCTTM, label conservation (ARI/NMI/silhouette/isolated-label F1/graph
      connectivity), batch correction (kBET-style acceptance, iLISI, batch ASW),
      label transfer.
  step 2 imputation    -- slide 6
      MSE, pFCS, pDES, gene-wise rho, densification, downstream
      (5 classification + 5 clustering metrics).
  step 3 GRN           -- slide 8
      AUPRC, early precision, stability, topology, perturbation (not
      computable without perturbation data), regulon activity.

Definitions follow the design document reports/pipeline-workflow-design.md
(sections 3-5), which in turn follows Liu et al. Nat Methods 2025 for MSE/pFCS/
pDES, BEELINE (Pratapa et al. 2020) for AUPRC/early precision, and scIB
(Luecken et al. 2022) for the label/batch metrics.

Run `python stepwise_metrics.py --self-test` for the unit checks (they build
synthetic data with a known answer and assert the metric returns it).
"""
from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #


def _as_float32(a):
    return np.ascontiguousarray(a, dtype=np.float32)


def knn_indices(X, k, exclude_self=True, query=None, chunk=4096):
    """k nearest neighbours in euclidean space.

    X      : reference points (n x d)
    query  : optional query points (m x d); neighbour indices then index X
    returns int array (n x k) or (m x k).  Self is dropped when query is None.
    """
    from sklearn.neighbors import NearestNeighbors

    X = _as_float32(X)
    kk = min(k + (1 if (exclude_self and query is None) else 0), X.shape[0])
    nn = NearestNeighbors(n_neighbors=kk, algorithm="auto").fit(X)
    Q = X if query is None else _as_float32(query)
    out = np.empty((Q.shape[0], kk), dtype=np.int32)
    for lo in range(0, Q.shape[0], chunk):
        _, idx = nn.kneighbors(Q[lo:lo + chunk], n_neighbors=kk)
        out[lo:lo + chunk] = idx
    if exclude_self and query is None:
        out = out[:, 1:]
    return out[:, :k] if out.shape[1] > k else out


def _drop_zero_variance(a, b):
    """Mask of genes usable for correlation.

    `a` and `b` are genes x observations, so every reduction is per gene (axis 1).
    """
    ok = np.isfinite(a).all(1) & np.isfinite(b).all(1)
    ok &= a.std(1) > 1e-9
    ok &= b.std(1) > 1e-9
    return ok


def hvg_mask(expr, n=100, min_cells=3):
    """Indices of the n most variable genes of `expr`, which is GENES x CELLS.

    The orientation is not guessed.  An earlier version transposed whenever
    genes > cells, which silently turned cell indices into gene indices on a
    900-gene x 90-cell held-out block -- so pFCS and gene-wise rho were computed
    over the first 90 genes of the axis.  Callers that hold cells x genes must
    transpose first.
    """
    e = np.asarray(expr, dtype=np.float32)
    if e.ndim != 2:
        raise ValueError(f"hvg_mask needs a 2-D genes x cells matrix, got {e.shape}")
    detected = (e > 0).sum(1)
    var = e.var(1, ddof=1)
    var[detected < min_cells] = -np.inf
    n = min(n, int(np.isfinite(var).sum()))
    return np.sort(np.argsort(-var)[:n])


# --------------------------------------------------------------------------- #
# step 1: integration
# --------------------------------------------------------------------------- #
def _frac_closer_chunked(Ha, Hb, partner, chunk=1024):
    """For each row of Ha: fraction of Hb rows closer than its partner row."""
    na, nb = Ha.shape[0], Hb.shape[0]
    frac = np.empty(na, dtype=np.float64)
    b_sq = (Hb.astype(np.float64) ** 2).sum(1)
    for lo in range(0, na, chunk):
        A = Ha[lo:lo + chunk].astype(np.float64)
        d2 = b_sq[None, :] - 2.0 * (A @ Hb.T.astype(np.float64)) + (A ** 2).sum(1)[:, None]
        d2 = np.maximum(d2, 0.0)
        rows = np.arange(A.shape[0])
        partner_d = d2[rows, partner[lo:lo + chunk]]
        frac[lo:lo + chunk] = (d2 < partner_d[:, None]).sum(1) / max(nb - 1, 1)
    return frac


def foscttm(H_rna, H_atac, partner_rna, partner_atac):
    """Fraction Of Samples Closer Than True Match, averaged over both directions.

    H_rna (nr x d), H_atac (na x d): per-modality embeddings.
    partner_rna[i] : row of H_atac that is cell i's true partner.
    partner_atac[j]: row of H_rna that is cell j's true partner.

    0 = every cell's nearest cross-modal neighbour is its true partner,
    0.5 = random.  Slide 4 reference values: 0.14-0.16 (MMD-MA), 0.22 (UnionCom).
    """
    a = _frac_closer_chunked(H_rna, H_atac, np.asarray(partner_rna))
    b = _frac_closer_chunked(H_atac, H_rna, np.asarray(partner_atac))
    return float(0.5 * (a.mean() + b.mean())), a, b


def foscttm_shuffled(H_rna, H_atac, partner_rna, partner_atac, seed=0):
    """Negative control: the same computation with the pairing permuted."""
    rng = np.random.default_rng(seed)
    pr = rng.permutation(partner_rna)
    pa = rng.permutation(partner_atac)
    val, _, _ = foscttm(H_rna, H_atac, pr, pa)
    return val


def label_conservation(H, labels, k=30, n_clusters=None, max_silhouette=20000,
                       seed=0):
    """Are cell types still separated in H?  (scIB label metrics.)

    ARI/NMI compare a clustering OF THE EMBEDDING with the given labels -- this
    is scIB's definition, and the clustering here is KMeans(k) rather than
    Leiden (leidenalg is not a dependency of this pipeline; recorded in the
    output as `cluster_algo`).

    silhouette        : mean silhouette of cells w.r.t. labels, in H.
    isolated_label_F1 : per label, precision/recall of a cell's kNN being the
                        same label; F1 per label, macro and size-weighted mean.
    graph_connectivity: for each label, size of the largest connected component
                        of the kNN graph induced on that label, over label size.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                                 silhouette_score)
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    H = _as_float32(H)
    labels = np.asarray(labels)
    uniq, y = np.unique(labels, return_inverse=True)
    n_lab = len(uniq)
    kk = min(k, max(2, H.shape[0] - 1))

    km = KMeans(n_clusters=n_clusters or n_lab, n_init=10,
                random_state=seed).fit_predict(H)
    idx = knn_indices(H, kk)
    same = (y[idx] == y[:, None])

    # isolated label F1
    f1, sizes = [], []
    for c in range(n_lab):
        cells = np.where(y == c)[0]
        if len(cells) < 2:
            continue
        nc = same[cells].sum(1).astype(np.float64)
        prec = nc / kk
        rec = nc / (len(cells) - 1)
        denom = prec + rec
        f1_c = np.where(denom > 0, 2 * prec * rec / np.maximum(denom, 1e-12), 0.0)
        f1.append(float(f1_c.mean()))
        sizes.append(len(cells))
    f1 = np.asarray(f1)
    sizes = np.asarray(sizes)

    # graph connectivity
    rows = np.repeat(np.arange(H.shape[0]), idx.shape[1])
    g = coo_matrix((np.ones(rows.size), (rows, idx.ravel())),
                   shape=(H.shape[0], H.shape[0])).tocsr()
    g = g.maximum(g.T)
    conn = []
    for c in range(n_lab):
        cells = np.where(y == c)[0]
        if len(cells) < 2:
            continue
        sub = g[cells][:, cells]
        ncomp, lab = connected_components(sub, directed=False)
        largest = np.bincount(lab).max()
        conn.append(largest / len(cells))

    sil = float("nan")
    if n_lab >= 2 and H.shape[0] > n_lab:
        n = H.shape[0]
        if n > max_silhouette:
            rs = np.random.default_rng(seed)
            take = rs.choice(n, max_silhouette, replace=False)
            sil = float(silhouette_score(H[take], y[take]))
        else:
            sil = float(silhouette_score(H, y))
    return {
        "ari": float(adjusted_rand_score(y, km)),
        "nmi": float(normalized_mutual_info_score(y, km)),
        "silhouette": sil,
        "isolated_label_F1_macro": float(f1.mean()) if f1.size else float("nan"),
        "isolated_label_F1_weighted":
            float((f1 * sizes).sum() / sizes.sum()) if f1.size else float("nan"),
        "graph_connectivity": float(np.mean(conn)) if conn else float("nan"),
        "n_labels": int(n_lab),
        "cluster_algo": "kmeans",
        "note": ("isolated_label_F1 uses scIB's definition (precision over the k "
                 "neighbours, recall over the label size), so it is capped by "
                 "k/(n_label-1) and must be read against the shuffled control"),
        "k": int(kk),
    }


def batch_correction(H, batch, k=50, max_silhouette=20000, seed=0):
    """Are datasets mixed in H?  (scIB batch metrics, no labels needed.)

    ilisi      : inverse Simpson index of the batch composition of each cell's
                 kNN, scaled to [0,1] (1 = perfectly mixed).  scIB reports the
                 median; this reports the mean over cells.
    kbet        : chi-square goodness of fit of each cell's kNN batch counts
                 against the global batch proportions; a cell is "rejected"
                 when p < 0.05.  acceptance = fraction accepted.  Under perfect
                 mixing this is ~0.95.  This is the kBET statistic computed
                 directly, not kBET's heuristic estimator, and is labelled as
                 such in the output.
    batch_asw   : |silhouette| with the batch as label (0 = indiscernible).
    cross_batch_nn : fraction of cells whose nearest neighbour is another batch.
    """
    from sklearn.metrics import silhouette_score
    from scipy.stats import chi2

    H = _as_float32(H)
    b = np.asarray(batch)
    keep = b != ""
    H = H[keep]
    b = b[keep]
    _, y = np.unique(b, return_inverse=True)
    B = int(y.max()) + 1
    if B < 2:
        return {"status": "unavailable",
                "reason": (f"only {B} batch in this row space, so there is nothing "
                           "to mix; a collapsed multiome row space has one batch per "
                           "dataset, not one per modality")}
    kk = min(k, max(2, H.shape[0] - 1))
    idx = knn_indices(H, kk)
    counts = np.zeros((H.shape[0], B), dtype=np.float64)
    np.add.at(counts, (np.repeat(np.arange(H.shape[0]), kk), y[idx].ravel()), 1)
    p = counts / kk
    glob = np.bincount(y, minlength=B) / len(y)

    ilisi = (1.0 / np.maximum((p ** 2).sum(1), 1e-12) - 1.0) / max(B - 1, 1)

    expected = kk * glob[None, :]
    stat = ((counts - expected) ** 2 / np.maximum(expected, 1e-12)).sum(1)
    accept = float((chi2.sf(stat, df=B - 1) > 0.05).mean())

    n = H.shape[0]
    if n > max_silhouette:
        rs = np.random.default_rng(seed)
        take = rs.choice(n, max_silhouette, replace=False)
        asw = float(abs(silhouette_score(H[take], y[take])))
    else:
        asw = float(abs(silhouette_score(H, y))) if B > 1 else 0.0

    return {
        "ilisi_mean": float(ilisi.mean()),
        "ilisi_median": float(np.median(ilisi)),
        "kbet_acceptance": accept,
        "kbet_statistic_mean": float(stat.mean()),
        "batch_asw_abs": asw,
        "cross_batch_nn_frac": float((y[idx[:, 0]] != y).mean()),
        "n_batches": B,
        "k": int(kk),
        "note": "kbet = chi2 acceptance on the kNN composition, not kBET's estimator",
    }


def label_transfer(H_ref, H_query, labels_ref, labels_query, k=30):
    """kNN label transfer from reference cells to query cells (slide 4, >0.85)."""
    from sklearn.metrics import accuracy_score, f1_score

    idx = knn_indices(_as_float32(H_ref), k, query=_as_float32(H_query))
    yref = np.asarray(labels_ref)
    uniq, yref_i = np.unique(yref, return_inverse=True)
    votes = yref_i[idx]
    pred = np.array([np.bincount(v, minlength=len(uniq)).argmax() for v in votes])
    yq = np.asarray(labels_query)
    keep = np.isin(yq, uniq)
    pred = uniq[pred[keep]]
    yq = yq[keep]
    return {
        "accuracy": float(accuracy_score(yq, pred)),
        "macro_f1": float(f1_score(yq, pred, average="macro")),
        "n_query": int(keep.sum()),
        "k": int(k),
    }


# --------------------------------------------------------------------------- #
# step 2: imputation  (all against held-out real RNA)
# --------------------------------------------------------------------------- #
def imputation_mse(pred, truth, hvg=None):
    """Mean squared error between guess and hidden RNA (lower better)."""
    d = np.asarray(pred, np.float64) - np.asarray(truth, np.float64)
    out = {"mse": float((d ** 2).mean()), "rmse": float(np.sqrt((d ** 2).mean()))}
    if hvg is not None and len(hvg):
        out["mse_hvg"] = float((d[hvg] ** 2).mean())
    return out


def gene_wise_rho(pred, truth, hvg=None):
    """Per-gene Pearson correlation against the hidden RNA (Polarbear/BABEL).

    Published ATAC->RNA ceiling: 0.155 (BABEL) to 0.22 (Polarbear) on SNARE-seq.
    """
    p = np.asarray(pred, np.float64)
    t = np.asarray(truth, np.float64)
    ok = _drop_zero_variance(p, t)
    if hvg is not None and len(hvg):
        sel = np.zeros(p.shape[0], bool)
        sel[np.asarray(hvg)] = True
        ok &= sel
    if ok.sum() < 2:
        return {"rho_mean": float("nan"), "rho_median": float("nan"),
                "n_genes_scored": 0}
    pp = p[ok] - p[ok].mean(1, keepdims=True)
    tt = t[ok] - t[ok].mean(1, keepdims=True)
    rho = (pp * tt).sum(1) / np.sqrt((pp ** 2).sum(1) * (tt ** 2).sum(1))
    return {"rho_mean": float(rho.mean()), "rho_median": float(np.median(rho)),
            "rho_p90": float(np.percentile(rho, 90)),
            "n_genes_scored": int(ok.sum())}


def pfcs(pred, truth, n_hvg=100):
    """pFCS: Pearson of the two gene-gene correlation matrices over the top-n HVGs.

    "Are the correlations the network model reads still the real ones?"  This is
    the metric to lead with (design doc section 4).  HVGs are selected on the
    HIDDEN TRUTH and the same genes are used for both matrices, so the two
    correlation matrices are over one gene set by construction.
    """
    hvg = hvg_mask(np.asarray(truth, np.float32), n=n_hvg)
    if len(hvg) < 3:
        return {"pfcs": float("nan"), "n_hvg": int(len(hvg))}
    cp = np.corrcoef(np.asarray(pred, np.float64)[hvg])
    ct = np.corrcoef(np.asarray(truth, np.float64)[hvg])
    iu = np.triu_indices(len(hvg), 1)
    a, b = cp[iu], ct[iu]
    good = np.isfinite(a) & np.isfinite(b)
    if good.sum() < 3:
        return {"pfcs": float("nan"), "n_hvg": int(len(hvg))}
    return {"pfcs": float(np.corrcoef(a[good], b[good])[0, 1]),
            "n_hvg": int(len(hvg))}


def pdes(pred, truth, labels, min_cells=10):
    """pDES: are differential-expression statistics preserved?

    Liu et al. recompute limma per cell type.  limma is R-only, so this uses the
    Welch t statistic of each gene (group vs the rest), which is the same
    quantity limma's moderated t estimates; the metric is the Pearson
    correlation between the truth's and the imputed matrix's (gene x group)
    statistic tables.  Recorded as `pdes_tstat` so it is never confused with a
    limma run.
    """
    p = np.asarray(pred, np.float64)
    t = np.asarray(truth, np.float64)
    y = np.asarray(labels)
    groups = [g for g in np.unique(y) if (y == g).sum() >= min_cells]
    if len(groups) < 2:
        return {"pdes": float("nan"), "n_groups": len(groups)}

    def stat_table(m):
        cols = []
        for g in groups:
            a = m[:, y == g]
            b = m[:, y != g]
            va = a.var(1, ddof=1) / a.shape[1]
            vb = b.var(1, ddof=1) / b.shape[1]
            se = np.sqrt(np.maximum(va + vb, 1e-12))
            cols.append((a.mean(1) - b.mean(1)) / se)
        return np.vstack(cols).T                      # genes x groups

    sp, st = stat_table(p), stat_table(t)
    ok = _drop_zero_variance(sp, st)
    if ok.sum() < 3:
        return {"pdes": float("nan"), "n_groups": len(groups)}
    a, b = sp[ok].ravel(), st[ok].ravel()
    return {"pdes": float(np.corrcoef(a, b)[0, 1]), "n_groups": len(groups),
            "n_genes_scored": int(ok.sum()),
            "variant": "welch-t (limma substitute)"}


def densification(pred, truth, measured):
    """Are guessed rows denser than measured rows?  (no inflation is the gate.)

    `measured` is the real-RNA block (genes x cells).  Also reports the per-gene
    scale check from the design doc: mean and variance shift between the two
    blocks, in units of the measured block's own sd.
    """
    p = np.asarray(pred)
    t = np.asarray(truth)
    m = np.asarray(measured)
    nz_p = float((p > 0).mean())
    nz_t = float((t > 0).mean())
    nz_m = float((m > 0).mean())
    sd = m.std(1, ddof=1) + 1e-9
    dmean = np.abs(p.mean(1) - m.mean(1)) / sd
    dvar = np.abs(p.var(1, ddof=1) - m.var(1, ddof=1)) / (sd ** 2)
    return {
        "nonzero_frac_imputed": nz_p,
        "nonzero_frac_hidden_truth": nz_t,
        "nonzero_frac_measured_reference": nz_m,
        "densification_ratio_vs_measured": nz_p / max(nz_m, 1e-12),
        "densification_ratio_vs_truth": nz_p / max(nz_t, 1e-12),
        "gene_mean_shift_median_sd": float(np.median(dmean)),
        "gene_var_shift_median_sd": float(np.median(dvar)),
        "frac_genes_mean_shift_gt_0.5sd": float((dmean > 0.5).mean()),
        "frac_genes_var_shift_gt_0.5sd": float((dvar > 0.5).mean()),
    }


def downstream(pred, truth, labels, n_pca=50, seed=0):
    """5 classification + 5 clustering metrics on the imputed matrix.

    Liu et al.'s downstream block.  Computed on PCA(n_pca) of each matrix, so the
    row is "does the imputed data still classify / cluster like the real data".
    Both the imputed and the truth matrix are scored; the truth row is the
    ceiling, not a competing method.
    """
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (accuracy_score, adjusted_rand_score,
                                 completeness_score, f1_score,
                                 homogeneity_score, normalized_mutual_info_score,
                                 precision_score, recall_score, silhouette_score)
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    y = np.asarray(labels)
    cells_x_genes = np.asarray(pred, np.float32).T
    truth_x_genes = np.asarray(truth, np.float32).T

    def block(m):
        k = min(n_pca, min(m.shape) - 1)
        z = PCA(n_components=k, random_state=seed).fit_transform(m)
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(max_iter=2000))
        cv = StratifiedKFold(5, shuffle=True, random_state=seed)
        yp = cross_val_predict(clf, z, y, cv=cv)
        _, yk = np.unique(y, return_inverse=True)
        km = KMeans(n_clusters=len(np.unique(y)), n_init=10,
                    random_state=seed).fit_predict(z)
        return {
            "clf_accuracy": float(accuracy_score(y, yp)),
            "clf_macro_f1": float(f1_score(y, yp, average="macro")),
            "clf_weighted_f1": float(f1_score(y, yp, average="weighted")),
            "clf_precision_macro": float(precision_score(y, yp, average="macro",
                                                         zero_division=0)),
            "clf_recall_macro": float(recall_score(y, yp, average="macro",
                                                   zero_division=0)),
            "clust_ari": float(adjusted_rand_score(yk, km)),
            "clust_nmi": float(normalized_mutual_info_score(yk, km)),
            "clust_silhouette": float(silhouette_score(z, km)),
            "clust_homogeneity": float(homogeneity_score(yk, km)),
            "clust_completeness": float(completeness_score(yk, km)),
        }

    return {"imputed": block(cells_x_genes), "truth": block(truth_x_genes)}


# --------------------------------------------------------------------------- #
# step 3: GRN
# --------------------------------------------------------------------------- #
def read_truth(path):
    """TF,TARGET csv -> deduplicated set of (TF, TARGET) pairs."""
    import csv

    pairs = set()
    with open(path) as fh:
        for row in csv.DictReader(fh):
            tf = (row.get("TF") or "").strip()
            tg = (row.get("TARGET") or row.get("target") or "").strip()
            if tf and tg:
                pairs.add((tf, tg))
    return pairs


def _ranked(net_df, cap):
    import pandas as pd

    tfcol = "TF" if "TF" in net_df.columns else net_df.columns[0]
    tgcol = "target" if "target" in net_df.columns else net_df.columns[1]
    net = net_df.drop_duplicates(subset=[tfcol, tgcol], keep="first")
    net = net.sort_values("importance", ascending=False).head(cap)
    return net, tfcol, tgcol


def grn_recovery(net_df, truth, ks=(100, 1000, 10000), cap=50000):
    """The repo's existing evaluation (run.py evaluate()) plus AUPRC.

    run.py reports `recovered` (how many truth edges appear anywhere in the
    network), precision@K and recall@K.  The deck asks for AUPRC and early
    precision, which are the ranked versions of the same thing; `prevalence` is
    what a random edge list scores, i.e. the negative control for AUPRC.
    """
    from sklearn.metrics import average_precision_score

    net, tfcol, tgcol = _ranked(net_df, cap)
    y = np.fromiter(((tf, tg) in truth for tf, tg in
                     zip(net[tfcol].values, net[tgcol].values)), bool, len(net))
    out = {
        "edges_total": int(len(net_df)),
        "edges_scored": int(len(net)),
        "truth_edges": int(len(truth)),
        "recovered": int((y).sum()),
        "auprc": float(average_precision_score(y, net["importance"].values))
        if y.any() else float("nan"),
        "prevalence_random_auprc": float(y.mean()) if len(y) else float("nan"),
        "early_precision": float(y[:max(int(y.sum()), 1)].mean()) if y.any()
        else float("nan"),
        "precision_at": {}, "recall_at": {},
    }
    for K in ks:
        top = y[:min(K, len(y))]
        out["precision_at"][str(K)] = float(top.mean()) if len(top) else float("nan")
        out["recall_at"][str(K)] = float(top.sum() / len(truth)) if truth else float("nan")
    return out


def grn_stability(net_dfs, top=10000):
    """Jaccard of the top-N edges across seeds (any reported difference must
    exceed this).  Also returns an edge-set size check."""
    sets, sizes = [], []
    for df in net_dfs:
        net, tfcol, tgcol = _ranked(df, top)
        sets.append(set(zip(net[tfcol], net[tgcol])))
        sizes.append(len(net))
    pairs = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            u = sets[i] | sets[j]
            pairs.append(len(sets[i] & sets[j]) / len(u) if u else float("nan"))
    out = {"n_seeds": len(sets), "top": top,
           "jaccard_mean": float(np.nanmean(pairs)) if pairs else float("nan"),
           "jaccard_min": float(np.nanmin(pairs)) if pairs else float("nan"),
           "pairwise": [float(x) for x in pairs],
           "edges_per_seed": [int(s) for s in sizes]}
    if len(sets) > 1:
        common = set.intersection(*sets)
        out["edges_in_all_seeds"] = int(len(common))
    return out


def grn_topology(net_df, truth, n_hub=20, cap=50000):
    """Hubs and degree structure, referenced to the truth edge list.

    The deck's version ("comparable to the planted network") needs a simulation.
    With a curation truth the comparable statement is: does the inferred
    out-degree ranking of regulators agree with the truth's, and do the same TFs
    come out as hubs?  Both degree vectors are computed on the SAME target
    universe (the targets the network scored) so the comparison is not an
    artefact of truth covering more genes.  Also reports the symmetric-edge
    fraction, which is the deck's signature of smoothing-inflated correlation.
    """
    from scipy.stats import spearmanr

    net, tfcol, tgcol = _ranked(net_df, cap)
    edges = list(zip(net[tfcol], net[tgcol]))
    targets = set(net[tgcol])
    tfs = sorted(set(net[tfcol]))
    inf_out = {tf: 0 for tf in tfs}
    for tf, tg in edges:
        inf_out[tf] += 1
    tru = {(a, b) for a, b in truth if b in targets}
    tru_out = {tf: 0 for tf in tfs}
    for tf, _ in tru:
        if tf in tru_out:
            tru_out[tf] += 1
    a = np.array([inf_out[t] for t in tfs], float)
    b = np.array([tru_out[t] for t in tfs], float)
    n_hub = min(n_hub, len(tfs))
    rho = float(spearmanr(a, b).statistic) if a.std() > 0 and b.std() > 0 else float("nan")
    hub_inf = set(np.array(tfs)[np.argsort(-a)[:n_hub]])
    hub_tru = set(np.array(tfs)[np.argsort(-b)[:n_hub]])
    eset = set(edges)
    sym = sum((tg, tf) in eset for tf, tg in eset) / max(len(eset), 1)
    return {
        "spearman_out_degree_vs_truth": rho,
        "hub_overlap_at_k": len(hub_inf & hub_tru) / max(n_hub, 1),
        "hub_k": int(n_hub),
        "n_tfs": int(len(tfs)),
        "out_degree_mean": float(a.mean()), "out_degree_median": float(np.median(a)),
        "out_degree_max": float(a.max()) if a.size else float("nan"),
        "frac_targets_with_in_degree_gt_0":
            float(len({tg for _, tg in edges}) / max(len(targets), 1)),
        "symmetric_edge_fraction": float(sym),
        "top_hubs_inferred": sorted(hub_inf),
        "n_targets_scored": int(len(targets)),
    }


def regulon_activity(expr, genes, net_df, labels, top_edges=50000,
                     n_regulons=15, min_targets=10):
    """AUCell activity of the top regulons, per cell type (slide 8).

    A regulon here is a TF plus its `min_targets`..N highest-ranked targets in
    the network.  AUCell score = rank-based AUC of the regulon's gene set in
    each cell (1 - normalised mean rank, the standard area-under-the-recovery-
    curve definition).  Specificity = (max cell-type mean - second max) / max,
    so a regulon active everywhere scores ~0 and a lineage-specific one ~1.
    """
    import pandas as pd

    net, tfcol, tgcol = _ranked(net_df, top_edges)
    y = np.asarray(labels)
    uniq = list(dict.fromkeys(y.tolist()))
    if len(uniq) < 2:
        return {"error": "need >=2 cell types"}

    gpos = {g: i for i, g in enumerate(genes)}
    X = np.asarray(expr, np.float32)                     # cells x genes
    E = X.T                                              # genes x cells
    n_genes, n_cells = E.shape

    # rank each gene in each cell from the TOP: rank 0 = highest expressed gene,
    # which is the direction AUCell's recovery curve integrates over.
    order = np.argsort(-E, axis=0, kind="stable")
    rank = np.empty_like(order)
    np.put_along_axis(rank, order, np.arange(n_genes)[:, None].repeat(n_cells, 1), axis=0)

    rows = []
    for tf, grp in net.groupby(tfcol, sort=False):
        tg = [g for g in grp[tgcol].tolist() if g in gpos]
        if len(tg) < min_targets:
            continue
        idx = np.array([gpos[g] for g in tg])
        m = len(idx)
        # AUC = 1 - (mean top-rank of the set)/(n_genes): 1 = the set is always at
        # the top of the cell's ranking, 0 = always at the bottom.
        auc = 1.0 - rank[idx].mean(0) / n_genes
        means = {str(c): float(auc[y == c].mean()) for c in uniq}
        vals = np.array(list(means.values()))
        srt = np.sort(vals)[::-1]
        spec = float((srt[0] - srt[1]) / srt[0]) if len(srt) > 1 and srt[0] > 0 else float("nan")
        rows.append({"tf": tf, "n_targets": m,
                     "peak_cell_type": max(means, key=means.get),
                     "specificity": spec,
                     "mean_auc_by_cell_type": means})
        if len(rows) >= n_regulons:
            break
    if not rows:
        return {"error": f"no TF has >= {min_targets} targets in the network"}
    specs = np.array([r["specificity"] for r in rows], float)
    peaking = {}
    for r in rows:
        peaking[r["peak_cell_type"]] = peaking.get(r["peak_cell_type"], 0) + 1
    return {"n_regulons": len(rows),
            "specificity_mean": float(np.nanmean(specs)),
            "specificity_median": float(np.nanmedian(specs)),
            "regulons_peaking_in": peaking,
            "regulons": rows}


# --------------------------------------------------------------------------- #
# shared transfer used by both methods (same math as run.py impute())
# --------------------------------------------------------------------------- #
def knn_transfer(H_ref, ref_expr, H_query, k=20):
    """exp(-distance) weighted average of the k nearest reference cells.

    Identical arithmetic to run.py's reverse-imputeKNN (K_NEIGHBOURS=20,
    weights exp(-d) row-normalised, imputed = ref_expr @ W.T); `--self-test`
    verifies that against run.impute on the repo's own layout.

    ref_expr : genes x n_ref.  Returns genes x n_query.
    """
    from scipy.sparse import coo_matrix

    kk = min(k, H_ref.shape[0])
    idx = knn_indices(H_ref, kk, query=H_query)
    d = np.linalg.norm(_as_float32(H_query)[:, None, :]
                       - _as_float32(H_ref)[idx], axis=2)
    e = np.exp(-d)
    W = e / e.sum(1, keepdims=True)
    rows = np.repeat(np.arange(H_query.shape[0]), kk)
    Ws = coo_matrix((W.ravel(), (rows, idx.ravel())),
                    shape=(H_query.shape[0], H_ref.shape[0])).tocsr()
    out = ref_expr @ Ws.T
    if hasattr(out, "toarray"):          # sparse reference in -> dense genes x cells
        out = out.toarray()
    return np.asarray(out, dtype=np.float32), float(d.mean())


def _gene_moments(x):
    """Per-gene mean and sd of a genes x cells matrix, sparse or dense.

    Sparse input is the normal case here (the measured reference is a CSR count
    matrix); densifying 36k genes x 10k cells just to take two moments would cost
    most of a gigabyte for nothing.
    """
    if hasattr(x, "tocsr"):
        a = x.tocsr()
        n = a.shape[1]
        s1 = np.asarray(a.sum(axis=1)).ravel().astype(np.float64)
        s2 = np.asarray(a.multiply(a).sum(axis=1)).ravel().astype(np.float64)
        mu = s1 / max(n, 1)
        var = np.maximum(s2 / max(n, 1) - mu ** 2, 0.0)
        # ddof=1, to match the dense branch exactly
        sd = np.sqrt(var * n / max(n - 1, 1))
    else:
        a = np.asarray(x, np.float64)
        mu = a.mean(1)
        sd = a.std(1, ddof=1)
    return mu[:, None], sd[:, None]


def align_scale(pred, reference):
    """Per-gene affine map of `pred` onto the reference block's moments.

    MultiVI decodes on its own normalised scale, so its ATAC rows are not on the
    same footing as measured log1p-CPM rows; a model that can tell the blocks
    apart learns the block, not the biology (design doc, step 2).  Fitted on
    MEASURED reference cells only -- never on the held-out answer key -- and it
    is affine, so it leaves Pearson correlations (pFCS, gene-wise rho) untouched
    and only repairs mean/variance.  Fit is a no-op for a block already on scale.
    """
    p = np.asarray(pred, np.float64)
    p_mu = p.mean(1, keepdims=True)
    p_sd = p.std(1, ddof=1, keepdims=True)
    r_mu, r_sd = _gene_moments(reference)
    good = (p_sd.ravel() > 1e-8) & (r_sd.ravel() > 1e-8) & np.isfinite(p_mu.ravel())
    out = p.copy()
    if good.any():
        out[good] = (p[good] - p_mu[good]) / p_sd[good] * r_sd[good] + r_mu[good]
    out[~good] = np.clip(out[~good], 0, None)
    frac_clipped = float((out < 0).mean())
    out = np.clip(out, 0, None).astype(np.float32)
    stats = {"frac_clipped": frac_clipped,
             "n_genes_rescaled": int(good.sum()),
             "n_genes_total": int(len(good))}
    return out, stats


# --------------------------------------------------------------------------- #
# self-test
# --------------------------------------------------------------------------- #
def _self_test():
    import tempfile, os, sys
    ok = True

    def chk(name, cond, detail=""):
        nonlocal ok
        print(f"  [{'ok' if cond else 'FAIL'}] {name} {detail}")
        ok = ok and bool(cond)

    rng = np.random.default_rng(0)

    # ---- FOSCTTM: aligned pairs vs random pairs
    n, d = 300, 8
    base = rng.normal(size=(n, d))
    Ha = base + 0.01 * rng.normal(size=(n, d))
    Hb = base + 0.01 * rng.normal(size=(n, d))
    idx = np.arange(n)
    v, _, _ = foscttm(Ha, Hb, idx, idx)
    ctrl = foscttm_shuffled(Ha, Hb, idx, idx, seed=1)
    chk("foscttm aligned ~0", v < 0.01, f"{v:.4f}")
    chk("foscttm shuffled ~0.5", 0.4 < ctrl < 0.6, f"{ctrl:.4f}")
    rnd_a, rnd_b = rng.normal(size=(n, d)), rng.normal(size=(n, d))
    chk("foscttm independent ~0.5", 0.4 < foscttm(rnd_a, rnd_b, idx, idx)[0] < 0.6)

    # ---- label metrics: two well separated clusters
    lab = np.repeat(["A", "B", "C"], 100)
    Hc = np.vstack([rng.normal(loc, 0.3, size=(100, d)) for loc in
                    [[0, 0, 0, 0, 0, 0, 0, 0], [10, 0, 0, 0, 0, 0, 0, 0],
                     [0, 10, 0, 0, 0, 0, 0, 0]]])
    lc = label_conservation(Hc, lab)
    chk("label ARI=1 on separated clusters", lc["ari"] > 0.99, f"ari={lc['ari']:.3f}")
    chk("graph connectivity ~1", lc["graph_connectivity"] > 0.95,
        f"{lc['graph_connectivity']:.3f}")
    lcs = label_conservation(Hc, rng.permutation(lab))
    chk("isolated-label F1 > 2x shuffled control",
        lc["isolated_label_F1_macro"] > 2 * max(lcs["isolated_label_F1_macro"], 1e-3),
        f"{lc['isolated_label_F1_macro']:.3f} vs {lcs['isolated_label_F1_macro']:.3f}")
    chk("label metrics collapse on shuffled labels", lcs["ari"] < 0.1,
        f"shuffled ari={lcs['ari']:.3f}")

    # ---- batch metrics: one mixed cloud vs two separated clouds
    b = np.repeat([0, 1], 200)
    mixed = rng.normal(size=(400, 5))
    sep = np.vstack([rng.normal(0, 1, (200, 5)), rng.normal(8, 1, (200, 5))])
    bm, bs = batch_correction(mixed, b), batch_correction(sep, b)
    chk("ilisi high when mixed", bm["ilisi_mean"] > 0.9, f"{bm['ilisi_mean']:.3f}")
    chk("ilisi low when separated", bs["ilisi_mean"] < 0.2, f"{bs['ilisi_mean']:.3f}")
    chk("kbet accepts mixing, rejects separation",
        bm["kbet_acceptance"] > 0.8 and bs["kbet_acceptance"] < 0.2,
        f"{bm['kbet_acceptance']:.2f} vs {bs['kbet_acceptance']:.2f}")

    # ---- label transfer
    lt = label_transfer(Hc, Hc + 0.01 * rng.normal(size=Hc.shape), lab, lab)
    chk("label transfer accurate", lt["accuracy"] > 0.99, f"{lt['accuracy']:.3f}")

    # ---- imputation metrics
    # genes > cells must still index genes: 5 high-variance genes in a 300x20 block
    big = rng.normal(0, 0.01, size=(300, 20)).astype(np.float32)
    big[[7, 42, 101, 200, 299]] += rng.normal(0, 5.0, size=(5, 20))  # scale, not shift
    h = hvg_mask(big, n=5, min_cells=2)
    chk("hvg_mask indexes genes when genes > cells",
        sorted(h.tolist()) == [7, 42, 101, 200, 299], str(h.tolist()))

    truth = np.abs(rng.normal(size=(60, 200))).astype(np.float32)
    truth[truth < 0.3] = 0
    pred_exact = truth.copy()
    pred_noise = np.abs(rng.normal(size=truth.shape)).astype(np.float32)
    m_ex, m_no, p_ex, p_no = (imputation_mse(pred_exact, truth),
                              imputation_mse(pred_noise, truth),
                              pfcs(pred_exact, truth), pfcs(pred_noise, truth))
    chk("mse 0 for exact copy", m_ex["mse"] == 0.0)
    chk("pfcs 1 for exact copy", abs(p_ex["pfcs"] - 1) < 1e-6, f"{p_ex['pfcs']:.3f}")
    chk("pfcs ~0 for noise", abs(p_no["pfcs"]) < 0.25, f"{p_no['pfcs']:.3f}")
    rho = gene_wise_rho(pred_exact, truth)
    chk("rho 1 for exact copy", abs(rho["rho_mean"] - 1) < 1e-6, f"{rho['rho_mean']:.3f}")
    dsl = densification(pred_exact, truth, truth)
    chk("densification 1 for exact copy",
        abs(dsl["densification_ratio_vs_measured"] - 1) < 1e-6)
    ylab = np.repeat(["A", "B"], 100)
    pdes_ex = pdes(pred_exact, truth, ylab)
    chk("pdes 1 for exact copy", abs(pdes_ex["pdes"] - 1) < 1e-6, f"{pdes_ex['pdes']:.3f}")
    ds = downstream(pred_exact, truth, ylab)
    chk("downstream ceiling: imputed == truth",
        abs(ds["imputed"]["clf_accuracy"] - ds["truth"]["clf_accuracy"]) < 0.15,
        f"{ds['imputed']['clf_accuracy']:.2f} vs {ds['truth']['clf_accuracy']:.2f}")

    # ---- scale alignment is correlation preserving
    scaled, sstats = align_scale(pred_noise, truth)
    chk("align_scale keeps gene correlation",
        abs(gene_wise_rho(scaled, truth)["rho_mean"]
            - gene_wise_rho(pred_noise, truth)["rho_mean"]) < 0.02,
        f"frac_clipped={sstats['frac_clipped']:.3f}")
    chk("align_scale fixes the block's scale",
        abs(float(scaled.mean()) - float(truth.mean()))
        < 0.25 * float(truth.std()), f"{scaled.mean():.3f} vs {truth.mean():.3f}")

    # ---- GRN metrics
    import pandas as pd
    truth_set = {(f"TF{i}", f"G{j}") for i in range(5) for j in range(20)}
    rows = []
    for i in range(5):
        for j in range(20):
            rows.append((f"TF{i}", f"G{j}", 10.0 - 0.01 * j))     # true edges first
        for j in range(20, 60):
            rows.append((f"TF{i}", f"G{j}", 0.5))                # decoys after
    perfect = pd.DataFrame(rows, columns=["TF", "target", "importance"])
    g = grn_recovery(perfect, truth_set, ks=(100,))
    chk("auprc 1 for perfect ranking", abs(g["auprc"] - 1) < 1e-9, f"{g['auprc']:.3f}")
    chk("recovered counts truth", g["recovered"] == 100, str(g["recovered"]))
    chk("prevalence = random baseline", abs(g["prevalence_random_auprc"] - 100 / 300) < 1e-9,
        f"{g['prevalence_random_auprc']:.3f}")
    rev = perfect.copy()
    rev["importance"] = rng.random(len(rev))          # same edges, random ranking
    chk("auprc ~prevalence for a random ranking",
        abs(grn_recovery(rev, truth_set)["auprc"]
            - grn_recovery(rev, truth_set)["prevalence_random_auprc"]) < 0.08)
    st = grn_stability([perfect, perfect, perfect], top=300)
    chk("stability 1 for identical networks", abs(st["jaccard_mean"] - 1) < 1e-12)
    topo = grn_topology(perfect, truth_set, n_hub=3)
    chk("topology hub overlap=1", abs(topo["hub_overlap_at_k"] - 1) < 1e-12,
        f"{topo['hub_overlap_at_k']:.2f} (k={topo['hub_k']})")

    # ---- regulon activity: a regulon built from cluster-A markers
    genes = [f"G{i}" for i in range(200)]
    X = np.zeros((200, 200), np.float32)
    X[:100, :20] = rng.normal(5, 0.1, size=(100, 20))     # cluster A markers
    X[100:, 20:40] = rng.normal(5, 0.1, size=(100, 20))   # cluster B markers
    X += np.abs(rng.normal(0, 0.01, size=X.shape)).astype(np.float32)
    net = pd.DataFrame([("TF0", genes[j], 10.0 - 0.01 * j) for j in range(20)],
                       columns=["TF", "target", "importance"])
    ra = regulon_activity(X, genes, net, np.repeat(["A", "B"], 100))
    chk("aucell regulon peaks in the right cell type",
        ra["regulons"][0]["peak_cell_type"] == "A",
        str(ra["regulons"][0]["peak_cell_type"]))

    # ---- knn_transfer must reproduce run.py's impute() arithmetic
    try:
        sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
        import run as base
        rng2 = np.random.default_rng(3)
        Href, Hq = rng2.normal(size=(80, 6)), rng2.normal(size=(25, 6))
        expr = rng2.normal(size=(30, 80)).astype(np.float32)   # genes x cells
        got, _ = knn_transfer(Href, expr, Hq, k=20)
        H = np.vstack([Href, Hq]).astype(np.float32)
        sizes = {"ref": Href.shape[0], "q": Hq.shape[0]}
        with tempfile.TemporaryDirectory() as td:
            base.K_NEIGHBOURS = 20
            ref, _ = base.impute(H, ["ref", "q"], sizes, expr, Href, ["q"], td,
                                 {"q": ["c"] * Hq.shape[0]})
        chk("knn_transfer == run.impute",
            ref.shape == got.shape and np.allclose(ref, got, atol=1e-6),
            f"max|d|={np.abs(ref - got).max():.2e}" if ref.shape == got.shape else "shape")
    except Exception as exc:                                  # noqa: BLE001
        chk("knn_transfer vs run.impute", False, f"skipped: {exc}")

    print("\nSELF-TEST", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        sys.exit(_self_test())
    print(__doc__)
