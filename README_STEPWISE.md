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

In one table, per step, what it produces and what grades it:

| step | what it produces | validated by (deck slide) |
|---|---|---|
| 1 integration | `H`: cells × 30 joint embedding (scSAGA), or the joint latent `z` plus per-modality `z_expr`/`z_acc` (MultiVI) | FOSCTTM, label conservation, batch correction, label transfer (4 / 11) |
| 2 imputation | a gene × ATAC-cell matrix, then the all-cells matrix `X` + `mask` (measured / guessed) | MSE, pFCS, pDES, gene-wise rho, densification, downstream (6 / 11) |
| 3 GRN | `network.tsv`: TF → target → importance, one per seed | AUPRC, early precision, stability, topology, perturbation, regulon activity (8), plus `run.py`'s recovered / precision@K / recall@K |

They run in that order because each one protects the next: step 1 answers *is `H` a
space where "nearest neighbour" means the same cell* (if not, step 2's average is
meaningless), step 2 answers *does the matrix the booster reads still carry the real
gene–gene correlation structure, on one scale* (if not, the network is confidently
wrong), step 3 answers *is the network better than random, and better than its own
seed noise* (if not, there is no result to report, whatever the AUPRC says).

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

## Run on PACE Phoenix — step by step

Eight steps, in order. Steps 1–3 run on a **login node**; steps 4–6 are the job.

If you have done it before, this is the whole thing:

```bash
export TMPDIR=/storage/scratch1/<n>/<user>/tmp PIP_CACHE_DIR=$TMPDIR/pip CONDA_PKGS_DIRS=$TMPDIR/pkgs
cd /storage/scratch1/<n>/<user>/pbmc-multiome-grn-inference
bash setup_stepwise.sh                          # once
conda activate scgrn && python scripts/stepwise_eval.py --check
SEEDS="666" sbatch --account=<acct> --qos=inferno run_stepwise.slurm
```

Then `results/stepwise/<experiment>/report.md`.

### 1. Put the repo on scratch, not in home

Home is quota'd at 10 GB and torch + scvi-tools alone will not fit. Everything —
repo, env, TMPDIR, pip cache — goes on scratch.

```bash
ssh <user>@phoenix.pace.gatech.edu

export TMPDIR=/storage/scratch1/<n>/<user>/tmp      # <n> = your scratch number
mkdir -p "$TMPDIR"

cd /storage/scratch1/<n>/<user>
git clone https://github.com/Swethasree-Bhattaram/pbmc-multiome-grn-inference.git
cd pbmc-multiome-grn-inference
git log --oneline -1          # should read: Step-wise experiment: MultiVI vs scSAGA ...
```

### 2. Build the environment (login node, once, ~10 min)

```bash
export TMPDIR=/storage/scratch1/<n>/<user>/tmp
export PIP_CACHE_DIR=$TMPDIR/pip
export CONDA_PKGS_DIRS=$TMPDIR/pkgs      # without these three you get Errno 122

bash setup_stepwise.sh
```

It runs the repo's own `setup.sh` (numpy/scipy/sklearn/pandas/scanpy/anndata,
arboreto+dask, scSAGA+torch) and then adds `muon` + `scvi-tools`. The last lines
must read:

```
  MULTIVI      importable  (...)
  scSAGA       importable
  arboreto     importable (+ modern-dask patch applied)
  all imports OK
```

`Errno 122` during the install means one of the three exports above did not take.
`muon`/`scvi-tools` failing only breaks the MultiVI arm, not scSAGA.

### 3. Check the inputs (login node, seconds, no compute)

The env is a conda env named `scgrn` (or `ENV_NAME` if you changed it). Activate it
before running anything by hand — the slurm script activates it itself.

```bash
conda activate scgrn          # if "conda: command not found": module load anaconda3
python scripts/stepwise_eval.py --check
```

One block per experiment. It must print, for the experiment you care about:

```
experiment paired: rna10k + atac10k  anchor=rna10k  queries=['atac10k']
  barcodes shared with the anchor: {'atac10k': 5000}     <- non-zero = answer key works
  -> held-out answer key available (FOSCTTM, MSE, pFCS, pDES, rho)
  labels.txt present: {'rna10k': True, 'atac10k': True}
  MultiVI backend: scvi-tools 1.4.2
  GRN regulators: .../data/tf_only.txt OK
  GRN targets: .../data/trrust_tf.txt OK
```

Exit code is 0 only when every check passes. Read the failure table below before
submitting if anything is flagged.

| `--check` says | what it means | what to do |
|---|---|---|
| `MISSING run.py inputs: [...]` | `pca_50.txt`, or an RNA dataset's `counts.mtx`/`barcodes.txt`/`features.txt`, is not there | fix the path in `config.yml` |
| `MultiVI needs ATAC peaks: <ds> has no counts.mtx/...` | MultiVI reads raw peaks; the ATAC dir currently has only the PCA | add the three files, or run `METHODS="scsaga"` |
| `counts N x M = FEATURES x CELLS (expected M x N)` | the matrix is transposed relative to the repo convention (cells x features). Usual for an ATAC export from 10x/ArchR/Signac | nothing to do: the script transposes it in memory. `python fix_counts.py` repairs the file itself (it now covers ATAC datasets too) |
| `counts N x M matches neither ...` | `barcodes.txt`/`features.txt` do not describe this matrix | check the three files belong together |
| `barcodes shared with the anchor: {..., 0}` | no ATAC cell matches an RNA cell by barcode | then FOSCTTM and the held-out block are reported unavailable — everything else still runs |
| `labels.txt present: {..., False}` | no cell-type labels | add `labels.txt` next to `pca_50.txt` (one label per PCA row, same order) or accept that the label metrics read `unavailable` |

The four-file layout, per dataset, is exactly what `run.py` already expects —
MultiVI simply also wants it for the ATAC side. `counts.mtx` must be cells x features;
a features x cells file is accepted and transposed in memory, but `python fix_counts.py`
is the thing to run once on the data:

```
data/3k_rna/    pca_50.txt  counts.mtx  barcodes.txt  features.txt  labels.txt(optional)
data/3k_atac/   pca_50.txt  counts.mtx  barcodes.txt  features.txt   <- last three new, for MultiVI
```

`labels.txt` is one line per PCA row. If you already have an annotation in an
AnnData or a CSV, write it out in that row order:

```bash
python - <<'PY'
import anndata as ad
a = ad.read_h5ad("your.h5ad")
open("data/3k_rna/labels.txt", "w").write("\n".join(a.obs["cell_type"].astype(str)) + "\n")
PY
```

### 4. Submit

```bash
sbatch --account=<acct> --qos=inferno run_stepwise.slurm
# e.g.
sbatch --account=gts-saluru8-coda20 --qos=inferno run_stepwise.slurm
```

`--account` and `--qos` are deliberately not baked in. Defaults: 1 node,
16 CPUs, 64 GB, 24 h, both methods, all experiments in `config.yml`, three seeds.

**First submission: one seed.** The 3-seed loop runs once per method and is what
fills the wall clock. Prove the input is right with

```bash
SEEDS="666" sbatch --account=<acct> --qos=inferno run_stepwise.slurm
```

then add the other two seeds without recomputing anything:

```bash
STEPS="grn" SEEDS="667,668" sbatch --account=<acct> --qos=inferno run_stepwise.slurm
STEPS="metrics" sbatch --account=<acct> --qos=inferno run_stepwise.slurm
```

Other overrides, all on the `sbatch` line:

```bash
EXPERIMENTS="paired"   ...   # one experiment (default: all in config.yml)
METHODS="scsaga"       ...   # one workflow (default: multivi,scsaga)
METHODS="multivi"      ...   # the other
EPOCHS=100             ...   # MultiVI max epochs (default 200, early stopping on)
FORCE=1                ...   # overwrite existing network.tsv
GRN_WORKERS=8 sbatch --cpus-per-task=8 --mem=32G ...   # small-memory variant
```

### 5. Watch it

```bash
squeue -u $USER
tail -f slurm-stepwise-<jobid>.out
```

The log goes, per experiment and per method, through these lines — if the last
one you see is much older than the elapsed time, look at the next line down for a
`[!]`:

```
EXPERIMENT paired
--- multivi
  [pipeline] MultiVI fit: 5422 cells (2711 anchor + 2711 ATAC barcodes), ... 
  [pipeline] affine scale alignment: median |delta| ...
  [pipeline] imputed (36601, 2711) genes x ATAC cells
  [pipeline] all-cells (5422, 36601) (2711 real + 2711 imputed)
  [heldout] 542 frozen cells imputed from 2169 reference cells; their RNA rows are removed from the model's training data too
  [grn] GRNBoost2 seed 666 -> .../pipeline/grn
  GRNBoost2: <N> edges in <t> min
  metrics -> .../metrics.json
--- scsaga
  ...
  report -> .../report.md
```

A failure inside one method is caught and printed as
`[!] multivi FAILED: <reason>`; the other method still runs and the report still
gets written. A missing `network.tsv` at the end is the one thing that fails the
job (`set -e` plus the final `find`).

### 6. Read the result

```bash
less results/stepwise/<experiment>/report.md     # the table you actually read
less results/stepwise/<experiment>/metrics.json  # every number, plus why anything is unavailable
```

`report.md` is one table: metric by step, one column per method, and a column
saying what a good value looks like. `n/a` always has a reason in `metrics.json`.
The slurm job also prints the first 60 lines of every report when it finishes.

### 7. Cost, so you can size the job

- scSAGA: minutes. MultiVI: two fits (full + frozen) — tens of minutes on CPU for
  ~10k cells at 200 epochs, dominated by the ATAC peak count.
- GRNBoost2 dominates, and it runs once per seed and twice (once per method).
  Runtime scales with **columns**, not targets. With 16 workers expect a few
  hours per seed at the repo's scale, so 3 seeds × 2 methods is what fills 24 h.
- Memory: GRN workers × ~0.35 GB + ~2 GB (8 → 32G, 16 → 64G, 20 → 96G), plus the
  count matrices MultiVI holds. If the job is killed for memory, drop
  `--cpus-per-task` and `--mem` together, in that ratio.

### 8. Smoke test first, if you want a five-minute rehearsal

On any machine, with no real data:

```bash
python scripts/stepwise_metrics.py --self-test     # metrics against known answers
python scripts/make_smoke_data.py --out results/smoke/data
python scripts/stepwise_eval.py \
    --config results/smoke/config_smoke.yml --epochs 100 --seeds 666 \
    --workers 4 --out results/smoke/results
```

This builds two datasets of the same 320 cells from one shared latent state, with
labels, a planted TF→target network and a ground-truth edge list. It proves the
plumbing end to end in about a minute. **Its numbers are not a benchmark of either
method** — 320 cells and 240 genes are far below the regime the methods were
designed for.

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
