#!/bin/bash
# Drive the pipeline stages that need DIFFERENT python environments.
#
#   pipeline/stages.sh [config.yml] [--all | --experiment NAME] [--stage S ...]
#
# Why this exists: the stages cannot share one interpreter.
#   integrate, impute, evaluate   sklearn/scipy/torch/h5py      -> $SC_ENV  (py3.12)
#   SCEMENT combined reference    anndata/scanpy                -> $SCE_ENV (py3.11)
#   GRNBoost2                     arboreto 0.1.6 + dask 2021.10 -> $GRN_ENV (py3.9)
#
# Each of SC_ENV / SCE_ENV / GRN_ENV may be either a conda env NAME or a full
# path to a python interpreter.  Paths work without conda, which makes this
# usable on a laptop that has plain venvs.
#
# The SCEMENT step is launched as a subprocess by engine.build_reference using
# $SCEMENT_PYTHON, so the main env never needs anndata.
#
# Stage selection: pass --stage one or more times; the default is all four.
set -euo pipefail

CONFIG="${1:-config.yml}"
shift || true
ARGS=("$@")
[ ${#ARGS[@]} -eq 0 ] && ARGS=(--all)

# Split --stage flags out so we know which stages to drive.  Anything else is
# forwarded to run.py verbatim.
STAGES=()
PASS_ARGS=()
i=0
while [ $i -lt ${#ARGS[@]} ]; do
  a="${ARGS[$i]}"
  if [ "$a" = "--stage" ]; then
    i=$((i+1)); STAGES+=("${ARGS[$i]}")
  else
    PASS_ARGS+=("$a")
  fi
  i=$((i+1))
done
[ ${#STAGES[@]} -eq 0 ] && STAGES=(integrate impute grn evaluate)
has_stage() { for s in "${STAGES[@]}"; do [ "$s" = "$1" ] && return 0; done; return 1; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "$HERE/.." && pwd)}"

# --- resolve the three interpreters -----------------------------------------
# Accept either an env name (resolved via conda) or a direct interpreter path.
resolve_py() {   # resolve_py <value> -> prints an interpreter path
  local v="$1"
  if [ -x "$v" ] && [ ! -d "$v" ]; then echo "$v"; return 0; fi
  if [ -x "$v/bin/python" ]; then echo "$v/bin/python"; return 0; fi
  if command -v conda >/dev/null 2>&1; then
    conda run -n "$v" python -c 'import sys;print(sys.executable)' 2>/dev/null && return 0
  fi
  return 1
}

if ! command -v conda >/dev/null 2>&1; then
  if declare -f module >/dev/null 2>&1; then module load anaconda3 2>/dev/null || true; fi
  if command -v conda >/dev/null 2>&1; then :; fi
fi

# Plain variables, not associative arrays: macOS ships bash 3.2, which has no
# `declare -A`.  Keeping this portable means the same script runs on the macmini
# and on PACE.
PY_SC=""; PY_SCE=""; PY_GRN=""
resolve_into() {  # resolve_into <VARNAME> <value>
  if p=$(resolve_py "$2"); then
    eval "$1=\$p"; echo "  $1 python: $p"
  else
    echo "ERROR: cannot resolve interpreter for $1 ('$2')." >&2
    echo "       Give a conda env name or a path to a python." >&2
    exit 1
  fi
}
resolve_into PY_SC  "$SC_ENV"
resolve_into PY_SCE "$SCE_ENV"
resolve_into PY_GRN "$GRN_ENV"
export SCEMENT_PYTHON="${SCEMENT_PYTHON:-$PY_SCE}"

# --- run one stage ----------------------------------------------------------
run() {   # run <interpreter> <label> <stage> <extra args...>
  local py="$1"; shift
  local label="$1"; shift
  local stage="$1"; shift
  echo; echo "--- [$label]"
  local t0; t0=$(date +%s)
  if "$py" "$HERE/run.py" --config "$CONFIG" "${PASS_ARGS[@]}" \
        --stage "$stage" "$@"; then
    echo "    OK $label ($(( $(date +%s) - t0 ))s)"
  else
    echo "    FAIL $label" >&2
    exit 1
  fi
}

has_stage integrate && run "$PY_SC"  "integrate (scSAGA joint embedding)"      integrate
has_stage impute    && run "$PY_SC"  "impute (reverse-imputeKNN -> all-cells)" impute
has_stage grn       && run "$PY_GRN" "grn (Arboreto GRNBoost2)"                grn \
                          --max-targets "${GRN_MAX_TARGETS:-0}" \
                          --max-regulators "${GRN_MAX_REGULATORS:-0}"
has_stage evaluate  && run "$PY_SC"  "evaluate (vs ground truth)"              evaluate
exit 0
