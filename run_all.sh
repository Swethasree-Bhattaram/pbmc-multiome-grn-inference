#!/bin/bash
# =============================================================================
#  RUN EVERYTHING  --  ONE command, start to finish
# =============================================================================
#
#      bash run_all.sh
#
#  That is the whole interface.  It does, in order:
#
#      1. create the three conda envs (skips any that already exist)
#      2. install scSAGA from GitHub (pip, no manual clone)
#      3. validate config.yml and every dataset
#      4. integrate -> impute -> GRN -> evaluate, for every experiment
#
#  You only ever edit config.yml.  Point it at your datasets, where each one
#  declares its modality and its file paths:
#
#      datasets:
#        my_rna:
#          modality: rna
#          pca:      /path/pca_50.txt        # REQUIRED  (cells x 50)
#          counts:   /path/counts.mtx        # RNa datasets only
#          barcodes: /path/barcodes.txt      # RNA datasets only
#          features: /path/features.txt      # RNA datasets only
#        my_atac:
#          modality: atac
#          pca:      /path/pca_50.txt        # REQUIRED; nothing else needed
#
#  WHY RNA NEEDS counts/barcodes/features: the ATAC cells' expression is
#  imputed, but the matrix that goes into GRNBoost2 also stacks the REAL RNA
#  cells on top.  Those come from counts.mtx.  ATAC datasets never need it.
#
#  On PACE, submit slurm_grn.slurm instead of running this directly; the SLURM
#  file calls this script.  On a laptop run it directly.
# =============================================================================
set -euo pipefail

ROOT="${PIPELINE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CONFIG="${CONFIG:-$ROOT/config.yml}"
cd "$ROOT"

SC_ENV="${SC_ENV:-scmint}"      # integration / imputation / evaluation
SCE_ENV="${SCE_ENV:-scement}"   # SCEMENT combined reference
GRN_ENV="${GRN_ENV:-grn39}"     # Arboreto GRNBoost2

GRN_WORKERS="${GRN_WORKERS:-${SLURM_CPUS_PER_TASK:-8}}"
GRN_THREADS_PER_WORKER="${GRN_THREADS_PER_WORKER:-1}"

EXPERIMENTS="${EXPERIMENTS:-}"          # blank = every experiment in config.yml
SKIP_ENVS="${SKIP_ENVS:-0}"             # 1 = assume envs already exist
SKIP_DATA="${SKIP_DATA:-0}"             # 1 = don't build data/ from raw/

echo "============================================================"
echo " multimodal GRN pipeline"
echo " root        : $ROOT"
echo " config      : $CONFIG"
echo " experiments : ${EXPERIMENTS:-<all in config>}"
echo " GRN workers : $GRN_WORKERS x $GRN_THREADS_PER_WORKER threads"
echo " host        : $(hostname)   date: $(date -Is)"
echo "============================================================"

# --- 0. conda ---------------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
  echo "--- loading anaconda3 module"
  module load anaconda3 2>/dev/null || true
fi
command -v conda >/dev/null 2>&1 || {
  echo "ERROR: conda not found.  On PACE: module load anaconda3" >&2; exit 1; }
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

# --- 1. environments --------------------------------------------------------
if [ "$SKIP_ENVS" != "1" ]; then
  echo; echo "### [1/4] environments"
  bash envs/setup_envs.sh
else
  echo; echo "### [1/4] environments -- skipped (SKIP_ENVS=1)"
fi

# --- 2. scSAGA --------------------------------------------------------------
echo; echo "### [2/4] scSAGA"
if conda run -n "$SC_ENV" python -c "import scmint.scsaga" 2>/dev/null; then
  echo "    already installed"
else
  echo "    installing from GitHub (no manual clone needed)"
  conda run --no-capture-output -n "$SC_ENV" \
    pip -q install "scsaga @ git+https://github.com/AluruLab/scSAGA.git"
  conda run -n "$SC_ENV" python -c \
    "import scmint.scsaga; print('    scmint.scsaga OK')"
fi

# --- 3. validate ------------------------------------------------------------
EXP_ARGS=(--all)
if [ -n "$EXPERIMENTS" ]; then
  EXP_ARGS=(); for e in $EXPERIMENTS; do EXP_ARGS+=(--experiment "$e"); done
fi

echo; echo "### [3/4] validate config + datasets"
conda run --no-capture-output -n "$SC_ENV" \
  python pipeline/run.py --config "$CONFIG" --check "${EXP_ARGS[@]}"

# --- 4. run -----------------------------------------------------------------
echo; echo "### [4/4] integrate -> impute -> GRN -> evaluate"
SC_ENV="$SC_ENV" SCE_ENV="$SCE_ENV" GRN_ENV="$GRN_ENV" \
PIPELINE_ROOT="$ROOT" \
GRN_WORKERS="$GRN_WORKERS" GRN_THREADS_PER_WORKER="$GRN_THREADS_PER_WORKER" \
  bash pipeline/stages.sh "$CONFIG" "${EXP_ARGS[@]}"

echo; echo "### done $(date -Is)"
find "$ROOT/results" -name network.tsv 2>/dev/null | sort | sed 's/^/  network:   /'
find "$ROOT/results" -name evaluation_summary.txt 2>/dev/null | sort | sed 's/^/  evaluated: /'
