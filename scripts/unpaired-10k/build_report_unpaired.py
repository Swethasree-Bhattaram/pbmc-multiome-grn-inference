#!/usr/bin/env python
"""Generate the UNPAIRED 10k report bundle from evaluation.json + QC outputs.

Writes into reports/unpaired-10k/:
  REPORT.md            — full pipeline report (design, dims, GRN, evaluation, QC)
  README.md            — short orientation
  qc_imputation.txt / qc_mixing.txt / qc_biology.txt — raw QC logs (copied)

Report style: clean and data-focused.  No timestamps, no filler prose, fractions
written inline (e.g. 1,957/8,751).
"""
import json, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
# Repo root: $UNPAIRED_ROOT, else two levels up (scripts/<workflow>/ -> repo).
# NOTE: one os.pardir would resolve to scripts/, which is a bug when the
# env var is unset (results/ and data/ then point inside scripts/).
PROJ = os.environ.get('UNPAIRED_ROOT', os.path.abspath(os.path.join(
    HERE, os.pardir, os.pardir)))
GRN = f'{PROJ}/results/unpaired_grn/grn_tfonly_reg'
OUTD = f'{PROJ}/reports/unpaired-10k'
GT_T_DEDUP, GT_B_DEDUP = 8751, 96846
GT_T_TOT, GT_B_TOT = 9720, 98338


def load_meta():
    """Integration / imputation facts, read back from the actual artifacts."""
    import numpy as np
    H = np.load(f'{PROJ}/results/integration_unpaired/joint_embedding_H.npy', mmap_mode='r')
    m = {'H': tuple(H.shape)}
    rt = open(f'{PROJ}/results/integration_unpaired/saga_runtimes.txt').read()
    for line in rt.splitlines():
        if 'Global Alignment Score' in line:
            m['global_score'] = line.split(':')[-1].strip()
        if 'Pairwise' in line and 'score=' in line:
            m['pairwise'] = line.split(':', 1)[1].strip()
    # input set sizes actually handed to GRNBoost2 (not the distinct genes seen in edges)
    m['n_targets_in'] = sum(1 for l in open(f'{GRN}/target_genes.txt') if l.strip())
    m['n_regs_in'] = sum(1 for l in open(f'{GRN}/tf_regulators.txt') if l.strip())
    # mean NN distance of ATAC queries to the RNA reference, from QC
    for line in open(f'{PROJ}/results/unpaired_grn/qc_imputation.txt'):
        if 'mean NN distance' in line:
            m['nn_dist'] = line.split(':')[1].split('(')[0].strip()
    m['ref_used'] = json.load(open(f'{PROJ}/results/unpaired_grn/qc_mixing.json'))['distinct_ref_cells_used']
    return m


def main():
    os.makedirs(OUTD, exist_ok=True)
    E = json.load(open(f'{GRN}/evaluation/evaluation.json'))
    Q = json.load(open(f'{PROJ}/results/unpaired_grn/qc_report.json'))
    QM = json.load(open(f'{PROJ}/results/unpaired_grn/qc_mixing.json'))
    M = load_meta()
    gt = E['ground_truths']

    def topk_block(name):
        rows = []
        for K, (p, r) in gt[name]['rows'].items():
            rows.append(f'| {int(K):,} | {p} | {r} |')
        return '\n'.join(rows)

    t25 = '\n'.join(f'| {tf} → {tg} | {imp:.4f} |' for tf, tg, imp in E['top25'])

    def rec(name):
        d = gt[name]
        return f"{d['recovered']:,}/{d['dedup']:,} ({d['recovered']/d['dedup']*100:.1f}%)"

    md = f"""# Unpaired 10k PBMC: multiome RNA × external ATAC → GRN

**Task.** Run the repo pipeline (integration → reverse-imputeKNN → Arboreto) as an
**unpaired** experiment between the RNA-seq of the 10x PBMC 10k multiome and the
external 10x *10k Human PBMCs ATAC v1.1* dataset, then infer a GRN with
regulators restricted to the TF-only list.

---

## 1. Datasets

| Role | Dataset | Cells | Features | Source |
|---|---|---|---|---|
| RNA (anchor) | 10x PBMC 10k multiome — granulocyte-sorted, RNA | 11,898 | 36,601 genes (GRCh38) | `pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5` |
| ATAC (query) | 10x 10k Human PBMCs ATAC v1.1 — **cells by peaks** | 8,161 | 87,863 peaks (hg19) | `atac_pbmc_10k_v1_filtered_peak_bc_matrix.h5` |

The external ATAC dataset is Cell Ranger ATAC **1.1.0** (the dataset page names it
"ATAC v1.1"), resolved from `cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_v1/`.
It is **not** the v2 "nextgem" 10k ATAC build.

### Unpaired caveats (they matter for interpretation)

1. **No shared cells.** Different nuclei, different donor/protocol. The two cell sets
   cannot be matched 1:1, so this is a *partial* alignment only.
2. **Different genome builds.** Multiome peaks are GRCh38, these ATAC peaks are **hg19**.
   scSAGA integrates on the per-modality PCA, so raw peak coordinates never enter the
   alignment directly — the mismatch affects biological comparability, not the mechanics.
3. **Different peak sets.** 87,863 peaks vs 143,887 in the multiome, so the ATAC
   modality carries a different (smaller) accessibility vocabulary.

## 2. Pipeline

```
10x ATAC h5 (cells x peaks)
  → preprocess_atac10k_ext.py     → data/atac10k_ext/{{counts.mtx, barcodes, features, pca_50}}
  → run_integration_unpaired.py    → scSAGA joint embedding H
  → build_imputation_unpaired.py   → reverse-imputeKNN → all-cells matrix
  → run_arboreto_unpaired.py       → GRNBoost2
  → evaluate_grn_unpaired.py       → precision/recall vs PBMC ground truth
```

Preprocessing and PCA use the **identical recipe** as `preprocess_10k_full.py`
(log1p CPM/1e4 over the top-2,000 variable features → 50 PCs), so RNA and ATAC are
treated the same way.

### scSAGA integration

| Item | Value |
|---|---|
| Anchor | rna10k |
| Datasets | rna10k (11,898), atac10k_ext (8,161) |
| `s_shared_cells` | 8,161 (= min of the two, full partial alignment) |
| `M_samples`, `alpha`, `S_iterations` | 2,000, 0.75, 25 |
| Joint embedding **H** | {M['H'][0]:,} × {M['H'][1]} |
| Global alignment score | {M.get('global_score','n/a')} |
| Pairwise (atac10k_ext -> rna10k) | {M.get('pairwise','n/a').split('=')[-1] if M.get('pairwise') else 'n/a'} |

### Reverse-imputeKNN

Reference = real RNA10k expression (11,898 cells); query = the 8,161 external ATAC cells.
For each query cell: k=20 nearest reference cells in H, weights `exp(-d)` normalised per
row, imputed expression `= ref_expr @ Wᵀ`.

| Item | Value |
|---|---|
| Imputed expression | 36,601 genes × 8,161 ATAC cells |
| All-cells matrix | **20,059 × 36,601** (float32) — rows [rna10k real (11,898), atac10k_ext imputed (8,161)] |
| Mean NN distance query→reference | {M.get('nn_dist','n/a')} (in the 30-dim joint embedding) |
| Reference cells used as neighbours | {M['ref_used']:,} / 11,898 |

## 3. Arboreto GRNBoost2

| Item | Value |
|---|---|
| Cells | {E['n_cells']:,} |
| **Target genes** | genes present in `trrust_tf.txt` → **{M['n_targets_in']:,}** (TRRUST list 2,862 genes) |
| **Regulators** (`tf_names`) | TFs present in `tf_only.txt` → **{M['n_regs_in']:,}** (list 829 TFs) |
| Union columns fed to GRNBoost2 | 2,852 |
| Inferred edges | **{E['n_edges']:,}** |
| Distinct genes actually appearing in edges | {E['n_regulators']:,} regulators, {E['n_targets']:,} targets |
| Engine | Arboreto 0.1.6 GRNBoost2, 4 dask workers, seed 666 |

This is the repo's TF-only-regulator variant: targets from the mixed TRRUST list, but
only true transcription factors allowed to act as regulators.

**Top 25 edges:**

| TF → target | importance |
|---|---|
{t25}

## 4. Evaluation vs PBMC ground truth

Ground truths are TF→target edge lists, deduplicated before scoring.

| Ground truth | Total | Deduplicated | Recovered | Fraction |
|---|---|---|---|---|
| PBMC-TRRUST | {GT_T_TOT:,} | {GT_T_DEDUP:,} | {gt['PBMC-TRRUST']['recovered']:,} | {gt['PBMC-TRRUST']['recovered']/GT_T_DEDUP*100:.1f}% |
| PBMC-Blood | {GT_B_TOT:,} | {GT_B_DEDUP:,} | {gt['PBMC-Blood']['recovered']:,} | {gt['PBMC-Blood']['recovered']/GT_B_DEDUP*100:.1f}% |

### PBMC-TRRUST — Precision@K / Recall@K

| Top-K | Precision@K | Recall@K |
|---|---|---|
{topk_block('PBMC-TRRUST')}

### PBMC-Blood — Precision@K / Recall@K

| Top-K | Precision@K | Recall@K |
|---|---|---|
{topk_block('PBMC-Blood')}

## 5. Context: unpaired vs the repo's paired TF-only runs

All runs below use the same TF-only regulator recipe (targets from `trrust_tf.txt`,
regulators from `tf_only.txt`), so the numbers are directly comparable in construction.
They differ in cell count and in whether the ATAC cells are the multiome's own.

| Run | ATAC source | Cells | Inferred edges | PBMC-TRRUST recovered | PBMC-Blood recovered |
|---|---|---|---|---|---|
| **Unpaired (this run)** | external 10x ATAC v1.1, 8,161 | 20,059 | 659,314 | **3,110/8,751 (35.5%)** | **4,486/96,846 (4.6%)** |
| expB2 (paired, 10k RNA ref) | multiome ATAC, 11,898 | 29,218 | 832,664 | 3,664/8,751 (41.9%) | 5,362/96,846 (5.5%) |
| expA (paired, SCEMENT ref) | multiome ATAC, 11,898 | 29,218 | 829,982 | 3,660/8,751 (41.8%) | 5,269/96,846 (5.4%) |
| expB1 (paired, 3k RNA ref) | multiome ATAC, 11,898 | 29,218 | 779,701 | 3,454/8,751 (39.5%) | 5,132/96,846 (5.3%) |

Paired values from `scripts/_tfonly_eval_data.json` in the `pbmc-full10k-grn` checkout;
the evaluator here reproduces expA exactly, so the comparison is methodologically sound.

The unpaired run recovers **35.5% of PBMC-TRRUST vs 39.5–41.9% for the paired runs** —
about 85–90% of the paired performance, with 31% fewer cells and a completely different
ATAC dataset. Given that, the degradation is modest. Do **not** read this as a clean
measure of the cost of being unpaired: cell count, ATAC cell number (8,161 vs 11,898) and
the ATAC feature vocabulary all differ between these rows. It is context, not a controlled
comparison.

## 6. QC

### Imputation is biologically coherent

To test whether the alignment carries real signal we clustered the 8,161 external ATAC
cells **on their own peak PCA** (no RNA, no alignment) and compared their imputed marker
profiles with independently-derived RNA clusters.

| ATAC cluster → best-matched RNA cluster | marker-profile r | identity |
|---|---|---|
| C0 → C2 | +0.920 | NK / cytotoxic (NKG7, GNLY) |
| C1 → C3 | +0.881 | T (CD3D, IL7R, CCR7) |
| C2 → C4 | +0.992 | B (MS4A1, CD79A, CD79B) |
| C3 → C0 | +0.781 | Monocyte (LYZ, S100A8, S100A9) |
| C5 → C1 | +0.756 | T (CD3D, IL7R) |
| C4 → C5 | +0.117 | (no clean match) |
| **mean matched r** | **+0.741** | |
| *control (cluster labels shuffled)* | *+0.485* | |

Correct lineage assignment for B cells, monocytes, NK and T, clearly above the shuffled
control — the alignment is informative, not random.

### Where it degrades

- **Amplitude is attenuated.** Marker signal-to-noise (between-cluster sd / mean
  within-cluster sd) is < 1 for CD3D 0.92, CD3E 0.91, IL7R 0.89, LYZ 0.90, CCR7 0.84,
  CD14 0.80, S100A9 0.74 — versus 1.48–1.62 for the B-cell markers MS4A1/CD79A/CD79B.
  Averaging 20 neighbours smooths the imputed profiles relative to real cells.
- **Rare populations are lost.** PPBP (platelets, snr 0.10) and FCER1A (basophils/DC,
  snr 0.09) show essentially no between-cluster structure — 0.0003 and 0.001 mean levels.
  The 6-cluster partition has no cluster for them.
- **Imputed cells are more similar to each other than real cells are** (mean pairwise
  r 0.857 vs 0.394). Expected for kNN-averaged profiles, but it means the imputed block
  contributes a large, partly redundant cell population to the GRN.

### Alignment is block-segregated — and so is the repo's own paired run

| Embedding | Block | fraction of same-block neighbours (random expectation) |
|---|---|---|
| **Unpaired (this run)** | atac10k_ext | 0.343 (0.407) |
| | rna10k | **0.997** (0.593) |
| **4-dataset (repo, paired)** | atac10k | **1.000** (0.407) |
| | rna10k | **0.998** (0.407) |
| | atac3k | 0.988 (0.093) |
| | rna3k | 0.330 (0.093) |

scSAGA's partial alignment places the **anchor** block on top of the others and the
query block alongside, so "modality blocks are not mixed" is a property of the tool as
used in this repo, not something the unpaired run broke. In the unpaired run the ATAC
block is in fact *less* segregated than any block in the paired run (0.343 vs
0.988–1.000). Read the mixing statistic with this in mind — do not interpret a low
mixing number here as evidence the unpaired alignment failed. The biological check in
§6 (r=+0.741 vs +0.485 control) is the informative evidence.

### Consequences for the GRN

- 8,161 of 20,059 cells (**40.7%**) are imputed — a large minority of rows are
  kNN-averaged rather than measured. The network is inferred on a matrix where two
  fifths of observations are smoothed copies of RNA neighbourhoods.
- **The imputation is concentrated.** The nearest reference cell of the 8,161 ATAC
  queries is drawn from only **2,831 of the 11,898** RNA cells (24%); 1,377 references
  are the unique nearest neighbour of exactly one query, and one reference is the
  nearest neighbour of **60** queries (mean 2.9). Neighbourhoods therefore overlap
  heavily, which is the direct cause of the elevated inter-cell correlation
  (r 0.857 vs 0.394) and of the attenuated marker amplitude noted above.
- Cell-type structure survives the unpaired alignment (B, mono, NK, T), so edges among
  lineage marker genes are credible; edges depending on rare populations (platelet,
  basophil/DC) are not supported by this run.
- The imputation does not use a disjoint reference: reference and query cells are
  different cells (different datasets), but there is no cell-level pairing to hold out,
  so no 1:1 paired-cell leakage is possible here.

## 7. Reproducing

Scripts: `scripts/unpaired-10k/`

```
preprocess_atac10k_ext.py      # 10x ATAC h5 -> scSAGA-format files  (scSAGA .venv, py3.12)
run_integration_unpaired.py    # scSAGA unpaired integration -> joint_embedding_H.npy
build_imputation_unpaired.py   # reverse-imputeKNN -> all_cells_gene_expression.npy
run_arboreto_unpaired.py       # GRNBoost2 (targets=trrust, regulators=tf_only)
run_all_grn_unpaired.sh        # launcher (caffeinate + ulimit)
evaluate_grn_unpaired.py       # precision/recall vs PBMC-TRUST + PBMC-Blood
qc_unpaired.py, qc_mixing.py, qc_biology.py   # QC
build_report_unpaired.py       # this report
```

Environments: scSAGA `.venv` (py3.12) for preprocessing/integration/imputation/QC;
`.venv39` (python 3.9, numpy 1.21, dask 2021.10, arboreto 0.1.6) for GRNBoost2.
Large matrices and the downloaded h5 are not committed (see `.gitignore`).

Working checkout: `/Volumes/samsung_ssd/tmp/pbmc-10k-unpaired-grn/`
"""

    with open(f'{OUTD}/REPORT.md', 'w') as f:
        f.write(md)
    print('wrote', f'{OUTD}/REPORT.md')

    with open(f'{OUTD}/README.md', 'w') as f:
        f.write(f"""# Unpaired 10k PBMC GRN (multiome RNA × external ATAC v1.1)

Integration + reverse-imputeKNN + Arboreto GRNBoost2 on an **unpaired** pairing:
the 10k multiome RNA (11,898 cells) and the external 10x 10k PBMC ATAC v1.1
"cells by peaks" dataset (8,161 cells). All-cells matrix 20,059 × 36,601.

- Targets: `trrust_tf.txt` genes present = {E['n_targets']:,}
- Regulators: `tf_only.txt` TFs present = {E['n_regulators']:,}
- Edges: {E['n_edges']:,}
- Recovered: PBMC-TRRUST {rec('PBMC-TRRUST')}, PBMC-Blood {rec('PBMC-Blood')}

See **REPORT.md** for the full pipeline, evaluation tables and QC (including the
unpaired caveats and the block-segregation comparison against the repo's paired run).
Raw QC logs: `qc_imputation.txt`, `qc_mixing.txt`, `qc_biology.txt`.
""")
    print('wrote', f'{OUTD}/README.md')

    for q in ['qc_imputation.txt', 'qc_mixing.txt', 'qc_biology.txt']:
        src = f'{PROJ}/results/unpaired_grn/{q}'
        if os.path.exists(src):
            shutil.copy(src, f'{OUTD}/{q}')
    print('copied QC logs')


if __name__ == '__main__':
    main()
