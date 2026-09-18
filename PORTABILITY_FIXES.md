# Portability fixes applied to this repo (for PACE Phoenix)

The scripts in this repository were written against a single macOS workstation
and contained hardcoded absolute paths plus one broken relative path. These
fixes make the repo runnable from a clean `git clone` on any machine (PACE
Phoenix in particular) without changing any numerical behaviour.

Every item below is a **path-resolution change only** — no algorithm, parameter,
seed, threshold or output format was touched. All edited files were re-compiled
after editing, and the pipeline was re-run end to end to confirm identical
shapes and identical regulator/target/edge counts.

## Fixes

| File | Was | Now |
|---|---|---|
| `scripts/unpaired-10k/*.py` (9 files) | `PROJ = os.environ.get('UNPAIRED_ROOT', os.path.join(HERE, os.pardir))` — **off-by-one**: with `UNPAIRED_ROOT` unset this resolves to `scripts/`, so `data/` and `results/` point inside `scripts/`. Only ever worked because the env var was always exported | two `os.pardir` levels, i.e. `scripts/<workflow>/` → repo root |
| `scripts/unpaired-10k/run_arboreto_unpaired.py` | `TR_F = f'{PROJ}/trrust_tf.txt'` — **broken path**: the file lives at `data/trrust_tf.txt`, and there is no `trrust_tf.txt` at the repo root, so the script could only ever work in a workspace that happened to have a loose copy | `TR_F = os.environ.get('TRRUST_FILE', f'{PROJ}/data/trrust_tf.txt')` |
| `scripts/unpaired-10k/run_arboreto_unpaired.py` | `TF_F` fixed to `data/tf_only.txt` | overridable via `TF_ONLY_FILE` |
| `scripts/unpaired-10k/qc_mixing.py` | `OLD = '/Volumes/samsung_ssd/tmp/pbmc-full10k-grn/...'` | `os.environ.get('PAIRED_H', f'{PROJ}/results/integration_10k/joint_embedding_H.npy')`; the comparison is skipped with a printed note when absent |
| `scripts/unpaired-10k/run_integration_unpaired.py` | `SCAGA_REPO` default `/Volumes/samsung_ssd/tmp/scSAGA` | `$SCAGA_REPO`, else `<repo parent>/scSAGA` |
| `scripts/multi-dataset-4x/run_integration_10k.py` | same hardcoded scSAGA default | `$SCAGA_REPO`, else `<repo parent>/scSAGA` |
| `scripts/multi-dataset-4x/preprocess_10k_full.py` | `PROJ = os.path.dirname(HERE)` — **off-by-one**: resolved to `scripts/`, so `PROJ/raw/10k.h5` and `PROJ/data` pointed inside `scripts/` | `os.environ.get('PBSC4K_ROOT', <two levels up>)` |
| `scripts/multi-dataset-4x/preprocess_10k.py` | same off-by-one | same fix |
| `scripts/common/extract_data.py` | `RAW`/`OUT` hardcoded under `/Users/sbhattaram/pbmc3k_analysis` | `$PBMC3K_ROOT` (default `~/pbmc3k_analysis`), `$PBMC3K_H5`, `$PBMC3K_INPUT` |
| `scripts/common/run_arboreto_single.py` | `DOWN`/`OUT`/`TF_F` hardcoded absolute paths | `$PBMC3K_ROOT` + `$SINGLE_DOWN_DIR` / `$SINGLE_GRN_DIR` / `$TRRUST_FILE` |
| `scripts/common/evaluate_grn_single.py` | `BASE` hardcoded | `$PBMC3K_ROOT` (default `~/pbmc3k_analysis`) |
| `scripts/common/run_scsaga.py` | `sys.path.insert(0, "/Volumes/samsung_ssd/tmp/scSAGA")`, `CFG` hardcoded | `$SCAGA_REPO`, `$SCSAGA_CONFIG` |
| `scripts/single-dataset-3k/reverse_impute.py` | `RES`/`IN`/`DOWN` hardcoded | `$PBMC3K_ROOT` + `$SCSAGA_OUT_DIR` / `$SCSAGA_IN_DIR` / `$SINGLE_DOWN_DIR` |
| `scripts/rna-baseline/write_report_tfonly.py` | `PROJ` hardcoded to `/Users/sbhattaram/pbmc-multiome-grn-inference`; `ORIG` hardcoded to a macOS scratch dir | `$PBSC4K_ROOT` (default: two levels up); `ORIG` via `$ORIG_GRN_DIR` |
| `scripts/rna-baseline-3k10k/write_report_tfonly_3k10k.py` | `PROJ` hardcoded | `$PBSC4K_ROOT` (default: two levels up) |

## Known issue that is NOT a path bug

`scmint/scsaga_saveH.py` (the scSAGA module that saves the joint embedding `H`)
was never committed to either `AluruLab/scSAGA` or its fork. Upstream
`scmint/scsaga.py` computes the 30-dim aligned embeddings inside
`Saga.run_multi()` and then discards them, so a fresh clone cannot produce
`joint_embedding_H.npy`, which the reverse-imputeKNN step depends on.

`pace-phoenix/envs/03_patch_scsaga.py` reproduces the required modification
against a fresh upstream clone and is idempotent + self-verifying. See
`pace-phoenix/README_PACE.md` §4.

## Pre-existing minor issues (left as-is, noted for completeness)

- `scripts/unpaired-10k/scsaga_unpaired_config.yml` and
  `scripts/multi-dataset-4x/scsaga_4dataset_config.yml` are *records* of the
  configs used for published runs and contain absolute macOS paths. They are
  overwritten by the runner scripts on the next run (the runners write a
  resolved config), so they are documentation rather than live configuration.
  `pace-phoenix` regenerates them per run via `--experiment`.
- `scripts/unpaired-10k/run_all_grn_unpaired.sh` is a macOS wrapper (uses
  `caffeinate`, absolute `.venv39` path). `pace-phoenix/scripts/06_run_grn.py`
  replaces its role portably; the original is kept for reproduction on macOS.
