#!/usr/bin/env python
"""Flatten results/stepwise/*/*/metrics.json into one tidy table (CSV or XLSX).

One row per (experiment, method). A second sheet holds one row per (experiment,
method, seed) so the stability claim can be checked, and a third states what every
column means, which direction is better, and which deck slide it comes from.

    python scripts/metrics_to_table.py                      # -> results/stepwise/
    python scripts/metrics_to_table.py --xlsx               # one workbook, 3 sheets
    python scripts/metrics_to_table.py --out ~/table.csv

Long format on purpose: adding a third method adds rows, not columns, and every
column stays filterable/pivotable.  A metric is never written without its control
next to it -- a pFCS or MSE on its own says nothing.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def get(d, *path, default=None):
    cur = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


# (column, path into metrics.json, decimals)   -- order is the order in the sheet
RUNS = [
    # ---- identity -----------------------------------------------------------
    ("experiment",              ("__exp__",), None),
    ("method",                  ("meta", "method"), None),
    ("row_space",               ("meta", "row_space"), None),
    ("labels",                  ("meta", "labels"), None),
    ("seeds",                   ("meta", "seeds"), None),
    ("grn_cap",                 ("meta", "grn_cap"), None),
    ("n_cells_anchor",          ("meta", "n_cells", "ANCHOR"), None),
    ("n_cells_atac",            ("meta", "n_cells", "ATAC"), None),
    ("n_genes_axis",            ("imputation", "n_genes"), None),
    ("n_heldout_cells",         ("imputation", "n_cells"), None),
    ("holdout_frac",            ("imputation", "holdout_frac"), 2),
    ("grn_columns",             ("__grn_run_info__", "columns"), None),
    ("grn_regulators",          ("__grn_run_info__", "regulators"), None),
    ("grn_targets",             ("__grn_run_info__", "targets"), None),
    ("grn_edges",               ("grn", "edges"), None),
    ("grn_edges_scored",        ("grn", "vs_ground_truth", "PBMC-TRRUST", "edges_scored"), None),

    # ---- step 1: integration ------------------------------------------------
    ("foscttm",                 ("integration", "foscttm", "foscttm"), 3),
    ("foscttm_control",         ("integration", "foscttm", "shuffled_control"), 3),
    ("foscttm_variant",         ("integration", "foscttm", "variant"), None),
    ("n_pairs",                 ("integration", "foscttm", "n_pairs"), None),
    ("ari",                     ("integration", "label_conservation", "ari"), 3),
    ("ari_control",             ("integration", "label_conservation", "shuffled_control", "ari"), 3),
    ("nmi",                     ("integration", "label_conservation", "nmi"), 3),
    ("silhouette",              ("integration", "label_conservation", "silhouette"), 3),
    ("isolated_label_f1",       ("integration", "label_conservation", "isolated_label_F1_macro"), 3),
    ("graph_connectivity",      ("integration", "label_conservation", "graph_connectivity"), 3),
    ("ilisi",                   ("integration", "batch_correction", "ilisi_mean"), 3),
    ("kbet_acceptance",         ("integration", "batch_correction", "kbet_acceptance"), 3),
    ("batch_asw_abs",           ("integration", "batch_correction", "batch_asw_abs"), 3),
    ("label_transfer_acc",      ("integration", "label_transfer", "accuracy"), 3),
    ("label_transfer_f1",       ("integration", "label_transfer", "macro_f1"), 3),

    # ---- step 2: imputation -------------------------------------------------
    ("mse",                     ("imputation", "mse"), 4),
    ("mse_hvg",                 ("imputation", "mse_hvg"), 4),
    ("mse_control_cells",       ("imputation", "shuffled_reference_control", "mse"), 4),
    ("pfcs",                    ("imputation", "pfcs"), 3),
    ("pfcs_gene_null",          ("imputation", "shuffled_reference_control", "gene_permuted_control", "pfcs"), 3),
    ("pfcs_control_cells",      ("imputation", "shuffled_reference_control", "pfcs"), 3),
    ("pdes",                    ("imputation", "pdes"), 3),
    ("rho_mean",                ("imputation", "rho_mean"), 3),
    ("rho_median",              ("imputation", "rho_median"), 3),
    ("rho_control_cells",       ("imputation", "shuffled_reference_control", "rho_mean"), 3),
    ("densification_ratio",     ("imputation", "densification_ratio_vs_measured"), 3),
    ("genes_mean_shift_gt_0.5sd", ("imputation", "frac_genes_mean_shift_gt_0.5sd"), 3),
    ("clf_accuracy_imputed",    ("imputation", "downstream", "imputed", "clf_accuracy"), 3),
    ("clf_accuracy_truth",      ("imputation", "downstream", "truth", "clf_accuracy"), 3),
    ("clust_ari_imputed",       ("imputation", "downstream", "imputed", "clust_ari"), 3),
    ("clust_ari_truth",         ("imputation", "downstream", "truth", "clust_ari"), 3),

    # ---- step 3: network, vs the repo's two truth sets ----------------------
    ("auprc_trrust",            ("grn", "vs_ground_truth", "PBMC-TRRUST", "auprc"), 4),
    ("prevalence_trrust",       ("grn", "vs_ground_truth", "PBMC-TRRUST", "prevalence_random_auprc"), 4),
    ("early_precision_trrust",  ("grn", "vs_ground_truth", "PBMC-TRRUST", "early_precision"), 4),
    ("recovered_trrust",        ("grn", "vs_ground_truth", "PBMC-TRRUST", "recovered"), None),
    ("precision10k_trrust",     ("grn", "vs_ground_truth", "PBMC-TRRUST", "precision_at", "10000"), 4),
    ("recall10k_trrust",        ("grn", "vs_ground_truth", "PBMC-TRRUST", "recall_at", "10000"), 4),
    ("auprc_blood",             ("grn", "vs_ground_truth", "PBMC-Blood", "auprc"), 4),
    ("prevalence_blood",        ("grn", "vs_ground_truth", "PBMC-Blood", "prevalence_random_auprc"), 4),
    ("early_precision_blood",   ("grn", "vs_ground_truth", "PBMC-Blood", "early_precision"), 4),
    ("recovered_blood",         ("grn", "vs_ground_truth", "PBMC-Blood", "recovered"), None),
    ("jaccard_3seed",           ("grn", "stability", "jaccard_mean"), 3),
    ("edges_in_all_seeds",      ("grn", "stability", "edges_in_all_seeds"), None),
    ("hub_overlap",             ("grn", "topology", "hub_overlap_at_k"), 3),
    ("hub_k",                   ("grn", "topology", "hub_k"), None),
    ("spearman_out_degree",     ("grn", "topology", "spearman_out_degree_vs_truth"), 3),
    ("symmetric_edge_frac",     ("grn", "topology", "symmetric_edge_fraction"), 3),
    ("regulon_specificity",     ("grn", "regulon_activity", "specificity_mean"), 3),

    # ---- gates --------------------------------------------------------------
    ("gate_foscttm",            ("integration", "foscttm", "gate"), None),
    ("gate_label_transfer",     ("integration", "label_transfer", "gate"), None),
    ("gate_pfcs",               ("imputation", "gate", "pfcs_above_gene_null"), None),
    ("gate_mse",                ("imputation", "gate", "mse_below_control"), None),
    ("gate_rho",                ("imputation", "gate", "rho_in_published_range"), None),
    ("gate_densification",      ("imputation", "gate", "no_densification"), None),

    # ---- cost (a speed claim is meaningless without the hardware) ----------
    ("sec_integrate",           ("meta", "seconds", "integrate"), 1),
    ("sec_impute",              ("meta", "seconds", "impute"), 1),
    ("sec_grn",                 ("meta", "seconds", "grn"), 1),
    ("sec_metrics",             ("meta", "seconds", "metrics"), 1),
    ("grn_workers",             ("meta", "workers"), None),
]


def rows_for(root, out_rows, seed_rows):
    for exp in sorted(os.listdir(root)):
        edir = os.path.join(root, exp)
        if not os.path.isdir(edir):
            continue
        for method in sorted(os.listdir(edir)):
            mpath = os.path.join(edir, method, "metrics.json")
            if not os.path.exists(mpath):
                continue
            d = json.load(open(mpath))
            # cell counts come from the integration step's own record: it is written
            # by every version of the pipeline, while meta.n_cells is newer
            info_path = os.path.join(edir, method, "integration", "integration_info.json")
            info = json.load(open(info_path)) if os.path.exists(info_path) else {}
            sizes = get(info, "sizes", default={}) or get(d, "meta", "n_cells",
                                                          default={}) or {}
            anc = info.get("anchor") or next(iter(sizes), None)
            anchor = sizes.get(anc) if anc else None
            atac = sum(v for k, v in sizes.items() if k != anc)
            row = {"experiment": exp}
            # GRN column counts, from run.py's own run_info.txt
            ri_path = os.path.join(edir, method, "pipeline", "grn", "run_info.txt")
            ri = {}
            if os.path.exists(ri_path):
                for line in open(ri_path):
                    if ":" in line:
                        k, v = line.split(":", 1)
                        ri[k.strip()] = v.strip()
            row["__grn_run_info__"] = ri
            for col, path, nd in RUNS:
                if path == ("__exp__",):
                    v = exp
                elif len(path) >= 3 and path[0] == "meta" and path[1] == "n_cells" \
                        and path[2] == "ANCHOR":
                    v = anchor
                elif len(path) >= 3 and path[0] == "meta" and path[1] == "n_cells" \
                        and path[2] == "ATAC":
                    v = atac
                elif path[0] == "__grn_run_info__":
                    v = ri.get(path[1])
                else:
                    v = get(d, *path)
                if isinstance(v, float) and nd is not None:
                    v = round(v, nd)
                row[col] = "" if v is None else v
            out_rows.append(row)

            # per-seed step-3 numbers, from the seed directories
            gdir = os.path.join(edir, method, "pipeline", "grn")
            seeds = str(get(d, "meta", "seeds") or "").split(",")
            for i, s in enumerate(seeds):
                net = os.path.join(gdir, "network.tsv" if i == 0 else f"seed{s}/network.tsv")
                if not os.path.exists(net):
                    continue
                ap_val = get(d, "grn", "vs_ground_truth", "PBMC-TRRUST", "auprc")
                seed_rows.append({
                    "experiment": exp, "method": method, "seed": s,
                    "network": os.path.relpath(net, root),
                    "edges": sum(1 for _ in open(net)) - 1,
                    "auprc_trrust": "" if i else (round(ap_val, 4)
                                                  if isinstance(ap_val, float) else ""),
                })


DEFINITIONS = [
    # column, what it asks, better, ground truth, deck slide
    ("foscttm", "for each cell, the fraction of cross-modal cells closer than its true partner", "lower", "paired cells", "4 / 11"),
    ("foscttm_control", "the same with the pairing permuted", "≈0.5", "paired cells", "4 / 11"),
    ("ari", "agreement between a clustering of the embedding and the cell-type labels", "higher", "labels.txt", "4"),
    ("ari_control", "the same with labels shuffled", "≈0", "labels.txt", "4"),
    ("nmi", "mutual information between embedding clusters and labels", "higher", "labels.txt", "4"),
    ("silhouette", "label separation in the embedding", "higher", "labels.txt", "4"),
    ("isolated_label_f1", "per-label kNN precision/recall; capped by k/n_label", "higher", "labels.txt", "4"),
    ("graph_connectivity", "largest connected component within each label's kNN subgraph", "higher", "labels.txt", "4"),
    ("ilisi", "batch mixing in each cell's neighbourhood, 1 = perfectly mixed", "higher", "dataset ids", "4"),
    ("kbet_acceptance", "chi-square acceptance of local batch composition (not kBET's estimator)", "→0.95", "dataset ids", "4"),
    ("batch_asw_abs", "|silhouette| with the batch as label, 0 = indiscernible", "lower", "dataset ids", "4"),
    ("label_transfer_acc", "kNN label from RNA, scored on the ATAC half of frozen paired cells", "higher (>0.85)", "paired cells + labels", "4"),
    ("mse", "squared error against the hidden RNA, on log1p CPM/1e4", "lower (but see control)", "held-out paired cells", "6"),
    ("mse_control_cells", "the same with reference cells permuted: a blurred guess wins this", "reference", "held-out paired cells", "6"),
    ("pfcs", "Pearson of the two gene-gene correlation matrices over the top-100 HVG", "higher", "held-out paired cells", "6"),
    ("pfcs_gene_null", "the same with reference gene labels permuted; the null pFCS is gated on", "≈0", "held-out paired cells", "6"),
    ("pdes", "correlation of per-gene differential-expression statistics (Welch t, not limma)", "higher", "held-out cells + labels", "6"),
    ("rho_mean", "mean per-gene Pearson against the hidden RNA; 0.155-0.22 is the published ATAC->RNA ceiling", "0.155-0.22", "held-out paired cells", "6"),
    ("densification_ratio", "nonzero fraction of imputed rows over measured rows", "≈1", "reference block", "6"),
    ("clf_accuracy_imputed", "5-fold logistic regression on PCA of the imputed held-out rows", "near truth", "labels", "6"),
    ("clust_ari_imputed", "KMeans on PCA of the imputed held-out rows vs labels", "near truth", "labels", "6"),
    ("auprc_trrust", "average precision of the ranked edge list against the truth set", "> prevalence", "PBMC-TRRUST", "8"),
    ("prevalence_trrust", "AUPRC a random ranking would score (= true-edge density in the list)", "reference", "PBMC-TRRUST", "8"),
    ("early_precision_trrust", "precision among the top-|truth| edges", "higher", "PBMC-TRRUST", "8"),
    ("recovered_trrust", "truth edges present anywhere in the network (run.py's own number)", "higher", "PBMC-TRRUST", "run.py"),
    ("jaccard_3seed", "mean Jaccard of the top-10k edges across seeds; differences must exceed it", "reported", "none", "8"),
    ("hub_overlap", "overlap of the top-k regulators by out-degree with the truth's", "higher", "PBMC-TRRUST", "8"),
    ("spearman_out_degree", "rank correlation of out-degree with the truth's, same target universe", "higher", "PBMC-TRRUST", "8"),
    ("symmetric_edge_frac", "fraction of edges present in both directions; high = smoothing", "low", "none", "8"),
    ("regulon_specificity", "(max - second) / max of a regulon's mean AUCell across cell types", "higher", "labels", "8"),
    ("sec_integrate", "wall time of the integration step (scSAGA OT+svds, or MultiVI train)", "reported", "none", "cost"),
    ("sec_impute", "wall time of imputation (kNN transfer, or MultiVI decode + the frozen refit)", "reported", "none", "cost"),
    ("sec_grn", "wall time of all GRN seeds together", "reported", "none", "cost"),
    ("grn_workers", "dask workers used for the GRN; memory scales with this, not with CPUs", "reported", "none", "cost"),
]


def write_csv(path, rows, columns=None):
    import csv
    cols = columns or list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {path}  ({len(rows)} rows x {len(cols)} cols)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(ROOT, "results", "stepwise"))
    ap.add_argument("--out", default=None, help="CSV path (default <results>/stepwise_runs.csv)")
    ap.add_argument("--xlsx", action="store_true", help="one workbook with 3 sheets")
    a = ap.parse_args()
    if not os.path.isdir(a.results):
        raise SystemExit(f"no results at {a.results}")
    runs, seeds = [], []
    rows_for(a.results, runs, seeds)
    if not runs:
        raise SystemExit(f"no metrics.json under {a.results}/*/*/")
    cols = [c for c, _, _ in RUNS]          # RUNS already starts with experiment, method
    defs = [{"column": c, "what it asks": q, "better": b,
             "ground truth": g, "deck slide": s} for c, q, b, g, s in DEFINITIONS]

    base = a.out or os.path.join(a.results, "stepwise_runs.csv")
    if not a.xlsx:
        write_csv(base, runs, cols)
        write_csv(os.path.join(os.path.dirname(base), "stepwise_per_seed.csv"), seeds,
                  ["experiment", "method", "seed", "network", "edges", "auprc_trrust"])
        write_csv(os.path.join(os.path.dirname(base), "stepwise_definitions.csv"), defs,
                  ["column", "what it asks", "better", "ground truth", "deck slide"])
        print("\nopen the CSVs in Excel, or re-run with --xlsx for one workbook")
        return 0

    try:
        import pandas as pd
    except ImportError:
        raise SystemExit("--xlsx needs pandas: pip install pandas openpyxl")
    path = base.replace(".csv", ".xlsx")
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        pd.DataFrame(runs)[cols].to_excel(xl, sheet_name="runs", index=False)
        pd.DataFrame(seeds).to_excel(xl, sheet_name="per_seed", index=False)
        pd.DataFrame(defs).to_excel(xl, sheet_name="definitions", index=False)
    print(f"wrote {path}  (sheets: runs {len(runs)}, per_seed {len(seeds)}, "
          f"definitions {len(defs)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
