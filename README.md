# Multimodal GRN inference

```
scSAGA integration -> reference RNA -> reverse-imputeKNN -> GRNBoost2 -> evaluation
   (joint embedding)   (one or all RNA)  (RNA -> ATAC)      (network)     (vs truth)
```

You give it datasets and a config. It gives you, for each experiment, an inferred
gene regulatory network and its recovery against ground truth — plus a comparison
across reference choices when there is more than one RNA dataset.

Nothing in the code is specific to PBMC, to an organism, to a modality or to a
number of datasets. **`config.yml` is the only file you edit.**

---

## Quick start

```bash
bash setup.sh                  # once: creates one environment, installs everything
python run.py --check          # verify your paths (no compute)
python run.py                  # run every experiment in config.yml
```

On PACE Phoenix:

```bash
bash setup.sh                                    # once, on a login node
cd <repo>                                        # run_grn.slurm resolves the
sbatch --account=<acct> --qos=inferno run_grn.slurm   # repo from $SLURM_SUBMIT_DIR
```

That is the whole interface.

Two things that are easy to get wrong the first time:

- **`EXPERIMENTS` defaults to every experiment in config.yml**, including the
  bundled PBMC demo entries. If you only add your own experiment, the untouched
  `ab` / `paired` / `unpaired` entries will fail `--check` on missing demo data.
  Either delete the demo experiments, or submit with `EXPERIMENTS="my_exp"`.
- **`setup.sh` must run on a login node.** Compute nodes usually have no network
  and it pip-installs scSAGA from GitHub.

---

## What you provide

Each dataset needs a PCA matrix. **RNA datasets additionally need raw counts**,
because the matrix handed to GRNBoost2 stacks the real RNA cells on top of the
imputed ATAC cells — those real rows have to come from somewhere.

| file | who needs it | format |
|---|---|---|
| `pca_50.txt` | **every** dataset | cells x N, whitespace separated. Any N works. |
| `counts.mtx` | **RNA only** | **cells x features**, MatrixMarket. Raw counts. |
| `barcodes.txt` | **RNA only** | one per line; row *i* of counts is cell *i* |
| `features.txt` | **RNA only** | one per line; column *j* of counts is feature *j* |

ATAC datasets need **only** the PCA. Their expression is imputed, never read.

`python run.py --check` verifies the files exist **and** that they agree:
counts rows == PCA rows == barcodes, and counts columns == features. All four
files must describe the same cells in the same order — a mismatch there
silently misaligns the integration, so it is checked up front.

Two ways to declare a dataset:

```yaml
datasets:
  # all four files in one directory, under the standard names
  my_rna:  {modality: rna,  dir: /scratch/me/my_rna}

  # or give every path explicitly (e.g. the PCA lives elsewhere)
  my_atac:
    modality: atac
    pca: /scratch/me/pcas/my_atac_50pc.txt

  # mixed: `dir` for the counts, `pca` for an external PCA
  other_rna:
    modality: rna
    dir: /scratch/me/other_rna
    pca: /scratch/me/pcas/other_50pc.txt
```

To use a different dataset, change the paths and run. No code changes.

---

## Experiments: declare, don't code

An experiment is one integration plus one or more **reference strategies**. A
strategy answers: *which real RNA expression do we propagate onto the ATAC cells?*

```yaml
experiments:
  ab:
    datasets:  [rna3k, atac3k, rna10k, atac10k]   # enters the integration
    anchor:    rna3k                              # everything aligns TO this
    queries:   [atac3k, atac10k]                  # ATAC cells to impute
    references:
      combined: {combine: all_rna}   # combine every RNA dataset into one
      rna3k:    {dataset: rna3k}     # use 3k RNA alone
      rna10k:   {dataset: rna10k}    # use 10k RNA alone
```

**N RNA datasets give N+1 strategies** — one combined, plus one per individual
dataset. Add a third RNA dataset and a third per-dataset strategy appears with
no code change.

The all-cells matrix is *always* the same real RNA cells plus the imputed ATAC
cells, so the strategies differ **only** in which reference fed the imputation —
which is exactly what isolates the effect of the reference choice.

---

## Output

```
results/<experiment>/
├── integration/
│   ├── joint_embedding_H.npy        anchor block first, then declaration order
│   ├── aligned_<name>.npy           per-dataset block of H
│   └── integration_info.json        anchor, order, sizes, params
├── <strategy>/
│   ├── all_cells_gene_expression.npy   real RNA rows + imputed ATAC rows
│   ├── genes.txt / all_cells_barcodes.txt / n_cells.txt
│   ├── imputed_expression_genes_x_atac.npy
│   └── grn/
│       ├── network.tsv              TF, target, importance
│       ├── run_info.txt             cells, columns, workers, seed, wall time
│       └── evaluation/
│           ├── evaluation_summary.txt   recovered edges vs each ground truth
│           ├── evaluation.json
│           └── top_edges.csv
├── comparison.md                    (experiment with >1 reference)
└── comparison.json
```

`comparison.md` reports network size per strategy, pairwise top-100 overlap
(shared + Jaccard), and recovery vs each ground truth with the best marked.

---

## Running part of the pipeline

Each step reuses what the earlier ones wrote, so you can iterate on the GRN
without redoing scSAGA:

```bash
python run.py -e ab --only integrate              # scSAGA only
python run.py -e ab --only impute                 # reference + imputation
python run.py -e ab --only grn,evaluate           # iterate here freely
python run.py -e ab --only grn --max-columns 3    # smoke test
```

Steps: `integrate`, `impute`, `grn`, `evaluate`.

---

## GRN: how much work this is

`grnboost2()` accepts no target-gene list. It forwards to `create_graph()`, whose
`target_genes` defaults to `'all'`, so **one regression is fitted per column** of
the matrix — including the regulator columns. Runtime scales with **columns**,
not targets:

```
~2,827 targets + ~816 regulators -> union ~2,852 columns -> ~2,852 regressions
```

That is why the GRN step dominates the run, and why the only way to make a test
cheap is to shrink the column set (`--max-columns`).

Parallelism is a runtime knob (Arboreto takes a pre-built dask client):

```bash
GRN_WORKERS=16 python run.py                    # or --workers 16
```

Measured on a real matrix (2,000 cells x 400 genes, seed 666, identical 21,998
edges at every setting):

| workers | wall | speedup |
|---|---|---|
| 1 | 92.0 s | 1.00× |
| 2 | 47.2 s | 1.95× |
| 4 | 25.2 s | 3.65× |
| **8** | **18.6 s** | **4.95×** |

Gains flatten past ~8–16 (the slowest single gene bounds the wall time).
**Memory scales with worker count**, because the regulator matrix is broadcast to
every worker: ~0.35 GB per worker + ~2 GB.

| workers | `--mem` |
|---|---|
| 8 | 32G |
| **16** | **64G** |
| 20 | 96G |

---

## On PACE Phoenix

`run_grn.slurm` carries no hardcoded paths — it runs the `config.yml` next to it,
wherever you cloned the repo. Submit from the repo directory or point at it:

```bash
sbatch --account=<acct> --qos=inferno run_grn.slurm
sbatch --account=<acct> --qos=inferno --chdir=/storage/.../my-repo run_grn.slurm

EXPERIMENTS="ab" sbatch ... run_grn.slurm       # one experiment
ONLY=grn,evaluate sbatch ... run_grn.slurm      # reuse earlier steps
GRN_WORKERS=8 sbatch --cpus-per-task=8 --mem=32G ... run_grn.slurm
```

`--account` and `--qos` are deliberately not baked in. Keep the repo on scratch,
not home (home is 20 GB). Use `--qos=inferno` for long runs; `embers` is free but
preemptible.

Build the environment **once on a login node** — compute nodes usually have no
network, so `setup.sh` (which pip-installs scSAGA from GitHub) cannot run inside
a job.

---

## Environment

One environment for everything.

```bash
bash setup.sh                  # conda env named scgrn
ENV_NAME=mygrn bash setup.sh   # different name
USE_VENV=1 bash setup.sh       # plain ./.venv, no conda
```

scSAGA already requires python>=3.10 and depends on anndata/scanpy, so the
integration, the reference combination (ComBat) and the GRN step can share one
interpreter. Python 3.11 is used by default.

### Why arboreto works on a modern stack

arboreto's final release (0.1.6) predates modern dask. In
`arboreto/core.py:create_graph` it calls `from_delayed(delayed_meta_dfs, ...)`
*unconditionally*, but only fills `delayed_meta_dfs` when `include_meta=True`.
With the default `include_meta=False` that is `from_delayed([])`, which old dask
tolerated and current dask rejects with `TypeError: Must supply at least one
delayed object`.

`grn_compat.py` guards that single call; the value is discarded when
`include_meta=False`. Nothing else in arboreto touches a changed dask API.

Verified equivalence — fixed input (300 cells x 40 genes, seed 666), **identical
390-edge network** on both stacks:

| | python | numpy | pandas | dask | sklearn |
|---|---|---|---|---|---|
| legacy | 3.9.6 | 1.21.5 | 1.4.4 | 2021.10.0 | 1.1.3 |
| modern | 3.11.16 | 2.4.6 | 3.0.6 | 2026.8.0 | 1.9.1 |

Edge sets match exactly (390/390). Edge *importances* differ in the low digits
because sklearn's `GradientBoostingRegressor` changed between 1.1 and 1.9 — that
is inherent to the sklearn version, not to the patch.

### Reference combination

The combined reference uses the vendored SCEMENT `sct_sparse` ComBat in
`scement.py`, which reproduces the published numbers. This is deliberate:
**scanpy's `pp.combat` is a different implementation** — on real PBMC counts
(5,422 cells x 36,601 genes, 2 batches) it differs from `sct_sparse` by
mean |Δ| 9.5e-4, with 10.6% of nonzeros off by more than 1e-3. Swapping in
scanpy would silently change published results, so the vendored version stays.

---

## Layout

```
config.yml        the only file you edit
run.py            the whole pipeline (integration, reference, impute, grn, evaluate)
grn_compat.py     makes arboreto work on modern dask
scement.py        vendored SCEMENT ComBat (sct_sparse)
setup.sh          one environment
requirements.txt  dependencies
run_grn.slurm     PACE Phoenix submission
data/             gene lists + ground truth (+ the bundled PBMC demo files)
results/          outputs
reports/          earlier write-ups
```

---

## Notes

- **No patched scSAGA.** Upstream's CLI computes the joint embedding and discards
  it; this pipeline calls `Saga.run_multi()` directly and saves `H` itself, so any
  upstream clone works.
- **Imputation** is `exp(-distance)` over the k=20 nearest reference cells,
  row-normalised. The impute step prints the mean neighbour distance and how many
  distinct reference cells were used — a low count means a collapsed alignment.
- **Ground truth is optional.** Any missing file is skipped with a note.
- Integration params (`integration:` in config.yml) are optional; the defaults
  are the values used for the published runs. `s_shared_cells` defaults to the
  smallest dataset; for paired multiome data you usually want the full cell count.
