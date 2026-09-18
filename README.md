# Multimodal GRN inference

Integration → imputation → gene regulatory network, in one config-driven
pipeline. Runs on PACE Phoenix or on a laptop, unchanged.

```
scSAGA integration  ->  reverse-imputeKNN  ->  Arboreto GRNBoost2  ->  evaluation
   (joint embedding H)     (RNA -> ATAC)          (network)            (vs ground truth)
```

The only file you edit is **`config.yml`**. Nothing is specific to PBMC, to a
particular modality, or to a particular number of datasets.

---

## Quick start

```bash
# 1. environments (once) -- creates scmint, scement, grn39
module load anaconda3 && bash envs/setup_envs.sh

# 2. scSAGA checkout (any clone of upstream works -- no patched fork needed)
git clone https://github.com/AluruLab/scSAGA.git tools/scSAGA

# 3. see what the config defines, and validate your inputs
python pipeline/run.py --list
python pipeline/run.py --check

# 4. run
python pipeline/run.py --all                 # every experiment in config.yml
python pipeline/run.py --experiment ab       # just one
```

On PACE submit the included SLURM file instead of step 4:

```bash
sbatch --account=<acct> --qos=inferno slurm_grn.slurm
```

---

## Datasets: you provide them

You supply the data; this repo does not download or format anything. Each dataset
needs a directory containing:

| file | needed for | format |
|---|---|---|
| `pca_50.txt` | **integration** (scSAGA) | cells × 50 PCs, whitespace-separated |
| `counts.mtx` | **RNA reference** only | raw counts, cells × features, MatrixMarket |
| `barcodes.txt` | RNA reference only | one barcode per line, order matches `counts` rows |
| `features.txt` | RNA reference only | one feature per line, order matches `counts` columns |

**Already have PCs?** Point straight at them — nothing else is required for
integration, because scSAGA consumes only the PCA file:

```yaml
datasets:
  my_rna:
    modality: rna
    dir: /scratch/me/expr          # counts.mtx barcodes.txt features.txt
    pca: /scratch/me/my_50pc.txt   # your own PCs; any cells x N file
```

If you have *only* PCs and no expression, that dataset can be integrated but
cannot serve as an RNA reference (there is no real expression to propagate) —
`--check` tells you which datasets are usable for what.

### Building the files from 10x h5 (optional helper)

`scripts/prepare_data.sh` builds `data/<name>/{pca_50.txt,counts.mtx,barcodes.txt,
features.txt}` from 10x h5 files. The PCA is 50 PCs on log1p(CPM/1e4) of the
top-2000 variable features. Skip it entirely if you already have the files.

---

## Experiments: declare, don't code

An experiment is one integration plus one or more **reference strategies**. A
strategy answers "which real RNA expression do we propagate onto the ATAC cells?"

```yaml
experiments:
  ab:
    datasets:  [rna3k, atac3k, rna10k, atac10k]   # what enters the integration
    anchor:    rna3k                              # everything aligns TO this
    queries:   [atac3k, atac10k]                  # ATAC cells to impute
    references:
      A:  {combine: all_rna}      # SCEMENT-combine every RNA dataset
      B1: {dataset: rna3k}        # use 3k RNA alone
      B2: {dataset: rna10k}       # use 10k RNA alone
```

**More RNA datasets require no code change.** With N RNA datasets in `datasets`
you get N+1 strategies — one SCEMENT-combined, plus one per individual dataset.
Add a fourth RNA dataset and a fourth strategy appears automatically.

The all-cells matrix is always the same real RNA cells plus the imputed ATAC
cells, so experiments differ **only** in which reference fed the imputation —
which is exactly what isolates the effect of the reference choice.

---

## Output layout

```
results/<experiment>/
├── integration/
│   ├── joint_embedding_H.npy        anchor block first, then datasets order
│   ├── aligned_<name>.npy            per-dataset block of H
│   └── integration_info.json         anchor, order, sizes, params
└── <strategy>/
    ├── all_cells_gene_expression.npy   real RNA rows + imputed ATAC rows
    ├── imputed_expression_genes_x_atac.npy
    ├── genes.npy / genes.txt / all_cells_barcodes.txt / n_cells.txt
    └── grn_tfonly_reg/
        ├── network.tsv                 TF, target, importance
        ├── tf_regulators.txt / target_genes.txt / all_columns.txt
        ├── run_info.txt                cells, columns, workers, seed, wall time
        └── evaluation/
            ├── evaluation_summary.txt  recovered edges vs each ground truth
            ├── evaluation.json
            └── top_edges.csv
```

`reports/` and `results/` from earlier runs are kept as-is.

---

## GRN specification

Set in `config.yml`:

```yaml
grn:
  regulators: data/tf_only.txt      # TFs -> GRNBoost2 tf_names
  targets:    data/trrust_tf.txt    # genes -> regression targets
```

- **Regulators** = TFs from `regulators` present in the matrix (~816)
- **Targets** = genes from `targets` present in the matrix (~2,827)
- The matrix handed to GRNBoost2 is the **union** of both (~2,852 columns),
  because every regulator must also be a column

This is identical on PACE and on a laptop; there is no platform-specific branch.

---

## Parallelism

Arboreto accepts a pre-built dask Client, so the worker count is a runtime knob —
see the [Arboreto user guide](https://arboreto.readthedocs.io/en/latest/userguide.html#running-with-a-custom-dask-client).

```bash
GRN_WORKERS=16 python pipeline/run.py --experiment ab --stage grn
GRN_SCHEDULER=tcp://host:8786 python pipeline/run.py ... --stage grn   # external
```

Measured on a real matrix (2,000 cells × 400 genes, seed 666, **identical
21,998 edges** at every setting):

| workers | wall | speedup |
|---|---|---|
| 1 | 92.0 s | 1.00× |
| 2 | 47.2 s | 1.95× |
| 4 | 25.2 s | 3.65× |
| **8** | **18.6 s** | **4.95×** |

Gains flatten past ~8–16 because work splits per target gene and the last gene
bounds the wall time.

**Memory** scales with worker *count* (not CPUs) because the TF matrix is
broadcast to every worker — `n_TFs × n_cells × 4 bytes`:

| workers | `--mem` |
|---|---|
| 8 | 32G |
| **16** | **64G** |
| 20 | 96G |

---

## Stages

Running a later stage reuses what earlier stages already wrote, so you can
iterate on the GRN without redoing scSAGA.

```bash
--stage integrate    scSAGA only            -> joint_embedding_H.npy
--stage impute       references + imputation
--stage grn          GRNBoost2
--stage evaluate     scoring
```

```bash
python pipeline/run.py --experiment ab --stage integrate    # once
python pipeline/run.py --experiment ab --stage grn         # iterate freely
```

---

## Environments

`envs/setup_envs.sh` builds all three. They cannot be merged:

| env | python | used for | why separate |
|---|---|---|---|
| `scmint` | 3.12 | integration, imputation, evaluation | torch/pot/geosketch |
| `scement` | 3.11 | SCEMENT combined reference | needs anndata/scanpy |
| `grn39` | 3.9 | GRNBoost2 | arboreto 0.1.6 needs numpy 1.21 + dask 2021.10 |

`grn39` pins **`click==8.0.4`**: `distributed 2021.10` imports
`click._unicodefun`, which click 8.1 removed, so with click ≥ 8.1 the
`dask-worker` / `dask-scheduler` CLIs die with an ImportError.

`pipeline/stages.sh` activates the right env per stage automatically.

---

## Environment variables

| var | meaning |
|---|---|
| `PIPELINE_ROOT` | repo root for relative paths in config.yml (default: config's dir) |
| `SCAGA_REPO` | scSAGA checkout (default `<root>/tools/scSAGA`) |
| `SCEMENT_PYTHON` | interpreter with anndata/scanpy, for the SCEMENT step |
| `GRN_WORKERS`, `GRN_THREADS_PER_WORKER` | dask workers; keep product = `--cpus-per-task` |
| `GRN_SCHEDULER` | attach to an external dask scheduler |
| `GRN_MAX_TARGETS` | debug: cap target genes (default off = all 2,827) |

---

## Layout

```
config.yml                 <- the only file you edit
pipeline/
  run.py                   entry point (--list, --check, stages)
  engine.py                integrate / build_reference / impute / run_grn / evaluate
  scement_combine.py       SCEMENT batch integration (runs in the scement env)
  scement_sparse.py        vendored numpy2-safe SCEMENT sct_sparse
  stages.sh                env-per-stage driver
envs/setup_envs.sh         conda envs
scripts/                   optional 10x h5 -> pca/counts prep, raw data download
slurm_grn.slurm            PACE Phoenix submission
data/                      your datasets + gene lists + ground truth
results/, reports/         outputs (earlier runs preserved)
```

---

## Notes

- **No patched scSAGA needed.** Upstream's CLI computes the joint embedding and
  discards it. This pipeline calls `Saga.run_multi()` directly and saves `H`
  itself, so any upstream clone works — and `aligned_*.npy` blocks are indexed by
  the integration order, so downstream never depends on the CLI's output.
- **Imputation weight**: `exp(-distance)` over the k=20 nearest reference cells,
  row-normalised. `impute` reports the mean neighbour distance and how many
  distinct reference cells were used — a low count means a collapsed alignment.
- **Ground truth is optional per file.** Any path missing from `ground_truth` is
  skipped with a printed note.
