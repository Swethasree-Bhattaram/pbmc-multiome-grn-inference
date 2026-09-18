# PACE Phoenix run bundle — scSAGA → reverse-imputeKNN → Arboreto GRNBoost2

Everything needed to run the PBMC multiome GRN pipeline on the PACE **Phoenix**
cluster from a clean checkout: environments, data download, per-step SLURM
scripts and a submit-the-whole-chain script.

```
pipeline/
├── pace.env                    # ALL site settings -- source this first
├── README_PACE.md              # this file
├── envs/
│   ├── 02_setup_envs.sh        # conda envs: scmint (py3.12), grn39 (py3.9)
│   └── 03_patch_scsaga.py      # clone AluruLab/scSAGA + apply the SAVE-H patch
├── scripts/
│   ├── 01_download_data.sh     # 4 raw 10x datasets, MD5-verified
│   ├── 00_split_10x_h5.py      # h5 -> scSAGA-format counts/barcodes/features/pca
│   ├── 00b_prepare_4k_rna.py   # 4k v2 tar.gz -> data/4k_rna (RNA baseline)
│   ├── 04_run_integration.py   # scSAGA integration -> joint_embedding_H.npy
│   ├── 05_build_imputation.py  # reverse-imputeKNN -> all-cells matrix
│   ├── 06_run_grn.py           # Arboreto GRNBoost2 (multi-worker dask client)
│   ├── 07_evaluate_grn.py      # Precision@K / Recall@K vs PBMC ground truth
│   ├── 99_preflight.sh         # verify setup before submitting
│   └── _bench_arboreto_workers.py, _smoke_arboreto_multiclient.py
├── pace/
│   ├── submit_all.sh           # submit the chain with afterok dependencies
│   └── slurm/
│       ├── _common.sh          # shared boilerplate (workspace + conda + banner)
│       ├── 00_download.slurm
│       ├── 01_preprocess.slurm
│       ├── 02_integration.slurm
│       ├── 03_imputation.slurm
│       ├── 04_grn.slurm
│       └── 05_evaluate.slurm
└── requirements.txt            # exact package versions, per env
```

---

## 1. Quick start

```bash
# ---- on the PACE login node -------------------------------------------------
git clone <this-repo> ~/scratch/pbmc-grn-pipeline && cd ~/scratch/pbmc-grn-pipeline
cd pipeline

# 1. point the bundle at your workspace and charge account
export PACE_ROOT=$HOME/scratch/pbmc-grn
export PACE_ACCOUNT=$(pace-quota | ...)       # your MAM account, e.g. GT-xxxx
source pace.env                                # or edit pace.env once and forget

# 2. build the two conda envs  (login node, ~10-15 min)
module load anaconda3
bash envs/02_setup_envs.sh

# 3. clone + patch scSAGA (REQUIRED -- see §4)
python envs/03_patch_scsaga.py
export SCAGA_REPO=$PACE_ROOT/tools/scSAGA

# 4. fetch the raw data (316 MB; fine on a login node)
bash scripts/01_download_data.sh

# 5. preflight -- catches configuration mistakes without burning queue time
bash scripts/99_preflight.sh

# 6. submit everything
bash pace/submit_all.sh unpaired
```

`submit_all.sh` prints the job IDs, chains each step with `afterok`, and tells
you the exact `tail -f` command for the running job.

Prefer manual control? Each step is its own script:

```bash
sbatch --account=$PACE_ACCOUNT --qos=embers pace/slurm/01_preprocess.slurm
sbatch --account=$PACE_ACCOUNT --qos=embers pace/slurm/02_integration.slurm
sbatch --account=$PACE_ACCOUNT --qos=embers pace/slurm/03_imputation.slurm
sbatch --account=$PACE_ACCOUNT --qos=embers pace/slurm/04_grn.slurm
sbatch --account=$PACE_ACCOUNT --qos=embers pace/slurm/05_evaluate.slurm
```

---

## 2. Datasets downloaded (step 0)

All four are the standard 10x Genomics public downloads, MD5-verified against the
copies used for the published runs in `pbmc-multiome-grn-inference`.

| File | Dataset | Size | MD5 |
|---|---|---|---|
| `3k_multiome.h5` | PBMC granulocyte-sorted **3k multiome** (ARC 2.0.0), 2,711 cells | 38,844,318 | `e326066b51ec8975197c29a7f911a4fd` |
| `10k_multiome.h5` | PBMC granulocyte-sorted **10k multiome** (ARC 2.0.0), 11,898 cells | 192,125,528 | `df86844b99161b9487090d91e644745e` |
| `atac10k_v1.1.h5` | **10k Human PBMCs ATAC v1.1** "cells by peaks" (hg19), 8,161 cells | 67,041,999 | `5ee75b4d7b5d70945a3c33a78c63582d` |
| `pbmc4k_…tar.gz` | **PBMC 4k** (10x v2, Cell Ranger 2.1.0), 4,340 cells | 18,423,814 | `f61f4deca423ef0fa82d63fdfa0497f7` |

URLs (in `scripts/01_download_data.sh`):

```
https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5
https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_v1/atac_pbmc_10k_v1_filtered_peak_bc_matrix.h5
https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz
```

**Naming trap:** the external ATAC dataset page says *v1.1* but the 10x sample
slug is `atac_pbmc_10k_v1` under pipeline version directory `1.1.0`. It is hg19,
while the multiome is GRCh38 — that mismatch is expected and does not affect
scSAGA, which integrates on each dataset's own PCA.

`01_download_data.sh` is idempotent (re-running verifies and skips) and resumes
partial downloads with `curl -C -`. `FORCE=1` re-downloads.

### What step 1 produces

`00_split_10x_h5.py` turns each h5 into `counts.mtx` (cells × features),
`barcodes.txt`, `features.txt` and `pca_50.txt` (50 PCs on log1p(CPM/1e4) of the
top-2000 variable features).

| Output | Contents | Needed by |
|---|---|---|
| `data/3k_rna`, `data/3k_atac` | 2,711 cells × 36,601 / 98,319 features | 3k and 3k+6k workflows |
| `data/10k_rna`, `data/10k_atac` | 11,898 cells × 36,601 / 143,887 features | **unpaired** and paired |
| `data/6k_rna`, `data/6k_atac` | deterministic 2,711-cell subsample of 10k (seed 0, same cells both modalities) | 3k+6k workflow only — skip with `SKIP_6K=1` |
| `data/atac10k_ext` | 8,161 cells × 87,863 peaks | **unpaired** |
| `data/4k_rna` | 4,340 cells × 33,694 genes | RNA-only baseline — skip with `SKIP_4K=1` |

The unpaired workflow needs only `10k_rna`, `atac10k_ext` plus the repo's
`data/{trrust_tf,tf_only}.txt` and `data/ground_truth/*.csv`, so
`SKIP_6K=1 SKIP_4K=1` makes a minimal run.

**Verified equivalence.** Each generated directory was hashed against the
originals used for the published results — all four files byte-identical for
`3k_rna`, `10k_rna`, `10k_atac`, `atac10k_ext`; `4k_rna` identical after a
`barcodes.tsv` → `barcodes.txt` rename. The PCA step in the splitter is also
numerically identical to the repo's two different PCA recipes (identical feature
selection; max |PC difference| 4.3e-12), so no downstream number shifts.

---

## 3. SLURM specification — which resources and why

| Step | Script | CPUs | Mem | Time | Binding constraint |
|---|---|---|---|---|---|
| download | `00_download.slurm` | 2 | 8G | 1:00:00 | network only; also fine on a login node |
| preprocess | `01_preprocess.slurm` | 8 | 64G | 2:00:00 | 10k h5 = 180,488 × 11,898, 127.5M nonzeros; sparse-variance PCA keeps peak RAM low |
| integration | `02_integration.slurm` | 8 | 64G | 2:00:00 | scSAGA is CPU-only; published unpaired run = 157 s of solver time |
| imputation | `03_imputation.slurm` | 8 | 96G | 4:00:00 | **memory**, not CPU: all-cells matrix is 20,059 × 36,601 float32 = 2.70 GB, plus a 1.07 GB reference |
| **GRN** | `04_grn.slurm` | **8** | 64G | 6:00:00 | wall time; scales with dask workers (see §5). Input subset is only 20,059 × 2,852 (0.23 GB) |
| evaluate | `05_evaluate.slurm` | 4 | 16G | 1:00:00 | trivial, reads the network .tsv |

All scripts use `--nodes=1 --ntasks-per-node=1`. Notes:

- **Queue.** Default `--qos=embers` (free, preemptible after a guaranteed first
  hour). For the multi-hour GRN run use `PACE_QOS=inferno` (spends credits) or
  `GRN_JOBQUEUE=1` to spread workers across jobs. `submit_all.sh` passes
  `--account` / `--qos` / `--partition` to every job from `pace.env`.
- **`--output`/`--error` are relative**, so the bundle writes them into
  `logs/slurm_<jobname>-%j.{out,err}` in the submit directory. That is why the
  SLURM scripts deliberately do **not** set `--chdir` — `_common.sh` does the
  `cd` after resolving the workspace, which keeps `logs/` next to the bundle.
- **Storage.** Put `PACE_ROOT` on `~/scratch` (15 TB) or project storage. Home is
  20 GB; the all-cells matrix alone is 2.7 GB. Scratch is purged after 60 idle
  days — the pipeline is fully regenerable from `raw/`, so that is safe.
- `conda activate` happens inside `_common.sh`; batch jobs never read `~/.bashrc`.

---

## 4. The scSAGA SAVE-H patch (required)

Upstream `AluruLab/scSAGA` computes the 30-dim joint embedding inside
`Saga.run_multi()` and then **discards it** — the released code never writes
`joint_embedding_H.npy`, and neither the `AluruLab` repo nor its fork contains
the modified file. The reverse-imputeKNN step is built entirely on that file
(and on nothing else from scSAGA), so the published runs used a locally modified
`scmint/scsaga.py`.

`envs/03_patch_scsaga.py` reproduces exactly that modification against a fresh
clone: right after the existing `save_joint_embedding_plot(...)` call it writes

```
aligned_<name>.npy / aligned_<name>_barcodes.txt     for every dataset
joint_embedding_H.npy                                (anchor block first)
joint_embedding_barcodes.txt
```

It is idempotent, verifies the result compiles, and installs the patched module
under both `scmint.scsaga` and `scmint.scsaga_saveH` so either import works.

```bash
python envs/03_patch_scsaga.py           # clone + patch into $SCAGA_REPO
python envs/03_patch_scsaga.py --check   # exit 0 only if the patch is present
```

`02_integration.slurm` runs `--check` first and `04_run_integration.py` aborts if
`joint_embedding_H.npy` is missing, so an unpatched checkout fails loudly instead
of silently poisoning downstream steps.

Verified on a fresh clone of upstream HEAD `b0d1438`: a small integration run
produced `joint_embedding_H.npy` of the correct shape.

---

## 5. Running Arboreto in parallel — the answer

**Yes.** Arboreto takes a pre-built `distributed.Client` through
`client_or_address`, so worker count is a runtime choice, not a code change
(see the [Arboreto user guide, "Running with a custom Dask
Client"](https://arboreto.readthedocs.io/en/latest/userguide.html#running-with-a-custom-dask-client)).
`06_run_grn.py` implements three modes:

```bash
# (a) single node, N worker processes  -- the default
GRN_WORKERS=8 GRN_THREADS_PER_WORKER=1 python scripts/06_run_grn.py --experiment unpaired

# (b) attach to a scheduler you started yourself
dask-scheduler --port 8786 &
dask-worker tcp://<node>:8786 --nthreads 8 &
GRN_SCHEDULER=tcp://<node>:8786 python scripts/06_run_grn.py --experiment unpaired

# (c) multi-node: the job submits its own SLURM worker jobs
GRN_JOBQUEUE=1 GRN_NODE_WORKERS=2 python scripts/06_run_grn.py --experiment unpaired \
    --slurm-account $PACE_ACCOUNT --slurm-qos embers
```

Keep `GRN_WORKERS × GRN_THREADS_PER_WORKER == --cpus-per-task`. Measured on the
real matrix shape (2,000 cells × 400 genes from the 10k multiome RNA, seed 666 —
**identical 21,998 edges at every setting**, so this is pure speedup):

| dask workers | cores | wall | speedup |
|---|---|---|---|
| 1 | 1 | 92.0 s | 1.00× |
| 2 | 2 | 47.2 s | 1.95× |
| 4 | 4 | 25.2 s | 3.65× |
| **8** | **8** | **18.6 s** | **4.95×** |

Diminishing returns past ~8: GRNBoost2 splits work per **target gene**, so the
last gene to finish bounds the wall time. That is why `04_grn.slurm` defaults to
8 CPUs — going to 16/32 buys little for a single node, and multi-node (mode c) is
the better lever for a large speedup.

Reference point: the published unpaired run (20,059 × 2,852, **4 workers**) took
2 h 10 m.

`run_info.txt` in the output directory records the mode, worker count, cores,
seed and wall time of every run, so the parallelism actually used is auditable.

---

## 6. Outputs

```
$PACE_ROOT/
├── data/                          prepared matrices (§2)
├── results/
│   ├── integration_unpaired/
│   │   ├── joint_embedding_H.npy          20,059 × 30
│   │   ├── joint_embedding_barcodes.txt
│   │   ├── aligned_{rna10k,atac10k_ext}.npy
│   │   ├── T_atac10k_ext_to_rna10k.npy    transport plan (11898 × 8161)
│   │   ├── saga_runtimes.txt              alignment score + timings
│   │   └── scsaga_config.yml              the resolved config actually used
│   └── unpaired_grn/
│       ├── all_cells_gene_expression.npy  20,059 × 36,601 float32 (2.7 GB)
│       ├── imputed_expression_genes_x_atac.npy  36,601 × 8,161
│       ├── all_cells_barcodes.txt, genes.txt, genes.npy, n_cells.txt
│       └── grn_tfonly_reg/
│           ├── grnboost2_network.tsv       inferred edges (TF, target, importance)
│           ├── tf_regulators.txt (816), target_genes.txt (2,827), all_columns.txt (2,852)
│           ├── run_info.txt                seed / workers / wall time
│           └── evaluation/{evaluation_summary.txt,top_edges.csv,evaluation.json}
└── logs/slurm_*.{out,err}
```

Sanity values from a real full-size local run of this bundle (unpaired):

| Quantity | Value |
|---|---|
| joint embedding H | 20,059 × 30 |
| scSAGA global alignment score | 0.564 (published 0.558; run-to-run variation) |
| all-cells matrix | 20,059 × 36,601 |
| regulators / targets / union columns | 816 / 2,827 / 2,852 |
| GRN input subset | 20,059 × 2,852 (0.23 GB) |

---

## 7. Paired vs unpaired

`PACE_EXPERIMENT` selects the workflow; everything else is identical.

| | unpaired (default) | paired |
|---|---|---|
| RNA | 10k multiome RNA, 11,898 cells | same |
| ATAC | **external** 10k ATAC v1.1, 8,161 cells (hg19, different nuclei) | the multiome's own ATAC, 11,898 cells |
| `s_shared_cells` | 8,161 (`min(n_rna, n_atac)` → partial alignment) | 11,898 |
| H | 20,059 × 30 | 23,796 × 30 |
| all-cells | 20,059 × 36,601 | 23,796 × 36,601 |

```bash
PACE_EXPERIMENT=paired bash pace/submit_all.sh paired
```

Note the published unpaired run used `s_shared_cells = 8,161` (the full smaller
dataset, i.e. a deliberately *partial* alignment, since the two modalities share
no nuclei). For a tol-converged run, raise the GW iterations — set
`S_iterations` in `scripts/04_run_integration.py` (`build_config`) from 25 to 200
and tighten `gw_epsilon` to 1e-7.

---

## 8. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `ImportError: cannot import name '_unicodefun' from 'click'` | `click ≥ 8.1` with `distributed 2021.10`. `02_setup_envs.sh` pins `click=8.0.4`; otherwise `conda install -n grn39 click=8.0.4` |
| `scSAGA finished but joint_embedding_H.npy is missing` | SAVE-H patch not applied → `python envs/03_patch_scsaga.py` |
| `H has N rows but expected …` in step 3 | integration ran with a different dataset set than the imputation expects; re-run step 2 with the same `PACE_EXPERIMENT` |
| dask worker warnings / `Event loop stopped before Future completed` at exit | harmless teardown noise from distributed 2021.10 after results are written; the job's exit code is what matters |
| `Sinkhorn … numerical errors at iteration N` | expected for a partial (unpaired) alignment; the run converges and writes its score |
| job killed at the time limit | GRN step: raise `--time`, use `PACE_QOS=inferno`, or `GRN_JOBQUEUE=1` to parallelise across nodes |
| OOM in step 3 | lower `GENE_CHUNK` (default 4,000) in `scripts/05_build_imputation.py`, or raise `--mem` |
| `SLURMCluster … unexpected keyword argument 'account'` | dask-jobqueue 0.7.2 uses `job_extra`, not `account=`; `06_run_grn.py` already does this |

## 9. Reproducing the upstream repo's own scripts

This bundle is additive — it does not modify `pbmc-multiome-grn-inference`. To
drive that repo directly instead, its scripts read `PBSC4K_ROOT` /
`UNPAIRED_ROOT` / `SCAGA_REPO`, and this bundle's data layout matches what they
expect (`data/{3k,6k,10k}_{rna,atac}`, `data/atac10k_ext`, `data/4k_rna`).
Portability fixes applied to that repo alongside this bundle are listed in
`PORTABILITY_FIXES.md`.
