# Step-wise experiment: MultiVI vs scSAGA, measured at every step

Two of the workflows in `workflow.pptx` are run over the same datasets, and each
step is scored with the metrics the deck lists for that step — before the next
step is allowed to run.

```
workflow 2   MultiVI          raw RNA counts + raw ATAC peak counts + modality mask
                              -> joint latent -> decoded RNA for the ATAC cells
              GRN             -> GRNBoost2 -> network

workflow 1   scSAGA          per-dataset PCA -> joint embedding H
              reverse-imputeKNN   -> k=20, exp(-distance) average of reference RNA
              GRN             -> GRNBoost2 -> network
```

Everything from imputation onwards is the same code for both, so a difference in
the network is attributable to the embedding and the imputation, not to the GRN
step: same gene axis, same TF list, same target list, same seeds, same scorer.

New files only. `run.py`, `config.yml`, `scement.py` and the original pipeline are
untouched.

| file | what it is |
|---|---|
| `scripts/stepwise_eval.py` | the experiment (integrate -> impute -> grn -> metrics) |
| `scripts/stepwise_metrics.py` | every metric, as pure functions; `--self-test` checks them |
| `scripts/make_smoke_data.py` | a tiny paired dataset with a planted network, for a safe first run |
| `setup_stepwise.sh` | adds `muon` + `scvi-tools` to the `setup.sh` env |
| `run_stepwise.slurm` | PACE Phoenix submission |
| `README_STEPWISE.md` | this file |

## What is measured, per step

The metrics are the ones in the deck (slides 4, 6, 8, 11). Each one is computed
per method, on that method's own artifacts.

**Step 1 — integration** (slides 4 and 11)

| metric | question |
|---|---|
| FOSCTTM (+ shuffled-pairing control) | is each cell nearest its true cross-modal partner? 0 perfect, 0.5 random, 0.14–0.22 published |
| ARI, NMI, silhouette, isolated-label F1, graph connectivity (+ shuffled-label control) | are cell types still separated? |
| iLISI, kBET-style acceptance, batch \|ASW\| | are datasets mixed? |
| label transfer accuracy / macro-F1 | can a label learned on RNA be assigned to the ATAC half of the same cell? (> 0.85) |

**Step 2 — imputation** (slide 6)

| metric | question |
|---|---|
| MSE, MSE over 1000 HVG | how far is the guess from the hidden value? |
| pFCS | is the gene–gene correlation structure preserved? (the metric to lead with) |
| pDES | are differential-expression statistics preserved? |
| gene-wise ρ (mean, median) | published ATAC→RNA ceiling 0.155–0.22; > 0.4 suggests leakage |
| densification | nonzero fraction of guessed vs measured rows; no inflation |
| downstream: 5 classification + 5 clustering metrics | does the imputed matrix classify/cluster like the real one? |

Two negative controls run with it: reference cells permuted (the design doc's
control, what MSE is gated on) and reference **gene labels** permuted (the null
pFCS is gated on). The second exists because permuting cells makes every guess a
near-constant average, and a near-constant block can score well on pFCS from its
mean profile alone.

**Step 3 — gene network** (slide 8)

| metric | question |
|---|---|
| recovered / precision@K / recall@K vs PBMC-TRRUST and PBMC-Blood | the number `run.py` already reports (the current way) |
| AUPRC, early precision, prevalence baseline | the ranked version of the same comparison (BEELINE) |
| stability: Jaccard of the top-10k edges across 3 seeds | any difference between methods must exceed it |
| topology: hub overlap, out-degree Spearman vs truth, symmetric-edge fraction | are hubs and degree structure recovered? |
| regulon activity: AUCell per cell type, and specificity | is the top regulon cell-type specific? |
| perturbation | **not computed** — needs Perturb-seq/knockout data that does not exist for PBMC multiome. The key is written to `metrics.json` with that reason, not faked. |

Slide 11's "run-to-run agreement" is the seed-to-seed Jaccard above (three GRNBoost2
seeds on identical input); it is not a second MultiVI fit.

Every metric that cannot be computed is written as
`{"status": "unavailable", "reason": "..."}`. Nothing is silently approximated.

## Inputs: the same four files, plus two optional ones

`config.yml` is unchanged and is still the only file you edit. The script reads
the datasets it declares.

| file | needed by |
|---|---|
| `pca_50.txt` (cells × N) | every dataset, both methods |
| `counts.mtx` + `barcodes.txt` + `features.txt` | RNA datasets (as before) |
| the same three files for **ATAC** datasets | **MultiVI only** — it needs peaks. The scSAGA arm still needs only the PCA. |
| `labels.txt` (one cell-type label per PCA row) | optional; without it every label-based metric is reported unavailable |

`python scripts/stepwise_eval.py --check` tells you, per experiment, which of
these are present, how many barcodes the anchor shares with each ATAC dataset
(that is what enables FOSCTTM and the held-out block), and whether the MultiVI
backend is importable.

## The held-out answer key

When an ATAC dataset shares barcodes with the anchor RNA dataset, 20% of those
paired cells are frozen deterministically (`--holdout-frac`, `--split-seed`):

- their real RNA is removed from the imputation reference;
- for MultiVI their RNA modality is also removed from training — a second fit —
  so the hidden values never enter the loss or the encoder;
- FOSCTTM, MSE, pFCS, pDES, ρ, densification and downstream are scored against them.

Without shared barcodes those metrics are reported unavailable, not approximated.

## Run on PACE Phoenix

```bash
# once, on a LOGIN node (compute nodes have no network): home is quota'd at 10 GB,
# so put the three things that default to home on scratch first
export TMPDIR=/storage/scratch1/<n>/<user>/tmp
export PIP_CACHE_DIR=$TMPDIR/pip
export CONDA_PKGS_DIRS=$TMPDIR/pkgs

cd /storage/.../pbmc-multiome-grn-inference
bash setup_stepwise.sh                      # setup.sh env + muon + scvi-tools
python scripts/stepwise_eval.py --check     # paths, pairing, labels, backend
```

```bash
# the run (CPU node; 16 workers, 64 GB)
sbatch --account=<acct> --qos=inferno run_stepwise.slurm

# useful variants
SEEDS="666"            sbatch ... run_stepwise.slurm    # one network per method first
EXPERIMENTS="paired"   sbatch ... run_stepwise.slurm
METHODS="scsaga"       sbatch ... run_stepwise.slurm
STEPS="metrics"        sbatch ... run_stepwise.slurm    # rescore without recomputing
GRN_WORKERS=8 sbatch --cpus-per-task=8 --mem=32G ... run_stepwise.slurm
```

Interactively, the same thing without SLURM:

```bash
python scripts/stepwise_eval.py -e paired --methods multivi,scsaga --workers 8
```

Steps can be re-run independently: `--only integrate`, `--only impute`,
`--only grn`, `--only metrics`. Each step reuses what the previous one wrote,
including the saved MultiVI model, so iterating on the network never refits the
model.

### Cost (estimates, not measurements on your data)

- scSAGA: minutes. MultiVI: 2 fits (full + held-out) of a few hundred epochs —
  tens of minutes on CPU for ~10k cells.
- GRNBoost2 dominates: runtime scales with the number of **columns**, not targets,
  and it runs once per seed and twice (once per method). With 16 workers expect a
  few hours per seed at the repo's scale; the 3-seed loop is what fills a 24 h job.
  Start with `SEEDS="666"`, then add the other seeds with `STEPS=grn SEEDS="667,668"`.
- Memory: GRN workers × ~0.35 GB + ~2 GB (8 → 32G, 16 → 64G, 20 → 96G), plus the
  count matrices MultiVI holds.

### Smoke test (five minutes, no real data needed)

```bash
python scripts/stepwise_metrics.py --self-test     # metrics against known answers
python scripts/make_smoke_data.py --out results/smoke/data
python scripts/stepwise_eval.py \
    --config results/smoke/config_smoke.yml --epochs 100 --seeds 666 \
    --workers 4 --out results/smoke/results
```

This builds two datasets of the same 320 cells from one shared latent state, with
labels, a planted TF→target network and a ground-truth edge list, and runs the
whole experiment on them. It proves the plumbing. **Its numbers are not a
benchmark of either method** — 320 cells and 240 genes are far below the regime
the methods were designed for.

## Three things to know before reading the numbers

**1. MultiVI collapses paired barcodes into one cell, scSAGA does not.** MultiVI
aligns its two encoders with a penalty applied only where a cell has *both*
modalities (`mask12 = mask_expr & mask_acc`). With one row per modality nothing
would ever align and the latent would simply separate by modality — measured:
FOSCTTM 0.43 with separate rows vs 0.00 with co-assayed rows on the same data. So
paired multiome cells are what make MultiVI a multi-omic model rather than two
VAEs. scSAGA keeps one block of rows per dataset, as `run.py` has always done.

A consequence worth knowing when you read the imputed block: for a *paired* ATAC
cell, MultiVI's decoded value is a reconstruction (the model saw that cell's RNA),
while scSAGA's is a kNN average from other cells. Both are written to the matrix
and both are flagged by `mask.npy`; `pipeline/impute_info.json` records which is
which. Step 2's metrics never touch these rows — they use the frozen arm.

**2. FOSCTTM therefore has two variants**, and `metrics.json` records which one
was used: `RNA row vs ATAC row of the same cell` (scSAGA, and MultiVI on unpaired
data) or `z_expr vs z_acc of the same co-assayed cell` (MultiVI on paired data,
where the two modality-specific latents of one cell are retrieved against each
other — the same question asked of the only two representations the model exports
per cell).

**3. Batch metrics are `n/a` for MultiVI in a fully paired experiment**, because a
collapsed row space has one batch per dataset rather than one per modality: there
is no cross-modality batch effect left to measure for cells that are literally the
same cells. They come back as soon as the experiment mixes a multiome dataset with
a single-modality one.

Two more substitutions worth knowing: pDES uses a Welch t statistic per gene
because limma is R-only, and `kbet_acceptance` is the chi-square acceptance on the
kNN composition rather than kBET's own estimator. Both are labelled as such in
`metrics.json`.

## Output

```
results/stepwise/<experiment>/
├── multivi/  and  scsaga/
│   ├── integration/
│   │   ├── joint_embedding_H.npy        cells x 30, float32
│   │   ├── joint_embedding_H_expr.npy   MultiVI only: per-modality latents
│   │   ├── joint_embedding_H_acc.npy
│   │   ├── aligned_<dataset>.npy        that dataset's rows of H
│   │   ├── blocks.json / integration_info.json   dataset -> row ids, params
│   │   └── multivi_model/               so --only impute/grn never refits
│   ├── pipeline/
│   │   ├── all_cells_gene_expression.npy  real RNA rows + imputed ATAC rows
│   │   ├── imputed_expression_*.npy, mask.npy, genes.txt, norm.txt
│   │   └── grn/network.tsv + seed<k>/network.tsv + evaluation/
│   ├── heldout/ imputed vs hidden answer key, split.json
│   └── metrics.json
├── metrics.json     both methods side by side
└── report.md        the table you actually read
```
