#!/bin/bash
# Drive the pipeline stages that need DIFFERENT python environments.
#
#   pipeline/stages.sh [config.yml] [--experiment NAME | --all]
#
# Why this exists: the stages cannot share one interpreter.
#   integrate, impute, evaluate   sklearn/scipy/h5py            -> $SC_ENV  (py3.12)
#   SCEMENT combined reference    anndata/scanpy               -> $SCE_ENV (py3.11)
#   GRNBoost2                     arboreto 0.1.6 + dask 2021.10 -> $GRN_ENV (py3.9)
#
# The SCEMENT step is launched as a subprocess by engine.build_reference using
# $SCEMENT_PYTHON, so the main env never needs anndata.
set -euo pipefail

CONFIG="${1:-config.yml}"
shift || true
EXP_ARGS=("$@")
[ ${#EXP_ARGS[@]} -eq 0 ] && EXP_ARGS=(--all)

SC_ENV="${SC_ENV:-scmint}"
SCE_ENV="${SCE_ENV:-scement}"
GRN_ENV="${GRN_ENV:-grn39}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "$HERE/.." && pwd)}"

if ! declare -f module >/dev/null 2>&1; then module() { return 0; }; fi
module load anaconda3 2>/dev/null || true
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda unavailable -- run 'module load anaconda3' first" >&2
  exit 1
fi
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

# Resolve the SCEMENT interpreter once and hand it to the pipeline.
export SCEMENT_PYTHON="$(conda run -n "$SCE_ENV" python -c 'import sys;print(sys.executable)' 2>/dev/null || true)"
if [ -z "$SCEMENT_PYTHON" ]; then
  echo "WARNING: env '$SCE_ENV' not found; any SCEMENT ('combine: all_rna')" >&2
  echo "         reference strategy will fail.  Run envs/02_setup_envs.sh" >&2
else
  echo "SCEMENT_PYTHON=$SCEMENT_PYTHON"
fi

run() {   # run <env> <label> <args...>
  local env="$1"; shift
  local label="$1"; shift
  echo
  echo "--- [$env] $label"
  local t0; t0=$(date +%s)
  if conda run --no-capture-output -n "$env" "$@"; then
    echo "    OK $label ($(( $(date +%s) - t0 ))s)"
  else
    echo "    FAIL $label" >&2
    exit 1
  fi
}

# integrate already ran in the SLURM script before calling this, but including it
# here makes stages.sh usable standalone; run.py skips work that is already done.
run "$SC_ENV" "integrate (scSAGA joint embedding)" \
    python "$HERE/run.py" --config "$CONFIG" "${EXP_ARGS[@]}" --stage integrate

run "$SC_ENV" "impute (reverse-imputeKNN -> all-cells matrices)" \
    python "$HERE/run.py" --config "$CONFIG" "${EXP_ARGS[@]}" --stage impute

run "$GRN_ENV" "grn (Arboreto GRNBoost2)" \
    python "$HERE/run.py" --config "$CONFIG" "${EXP_ARGS[@]}" --stage grn

run "$SC_ENV" "evaluate (vs ground truth)" \
    python "$HERE/run.py" --config "$CONFIG" "${EXP_ARGS[@]}" --stage evaluate
