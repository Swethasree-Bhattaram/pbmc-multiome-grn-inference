#!/bin/bash
# Shared SLURM boilerplate for every step in this pipeline.
#
# Not submitted directly.  Each step script in pace/slurm/ sources it after its
# own #SBATCH block:
#
#     ROOT_OVERRIDE=/path/to/workspace   # optional, else $PACE_ROOT/$SCRATCH
#     source "$(dirname "$0")/_common.sh"
#
# It resolves the workspace, loads conda, activates the requested env and prints
# a short banner so every .out file records which code and env ran.

set -euo pipefail

# ---- workspace -------------------------------------------------------------
# Resolution order: explicit SBATCH --chdir / ROOT_OVERRIDE, then PACE_ROOT,
# then $SCRATCH, then a workspace next to this script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUNDLE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"        # .../pipeline

if [ -n "${ROOT_OVERRIDE:-}" ]; then
  export PACE_ROOT="$ROOT_OVERRIDE"
elif [ -z "${PACE_ROOT:-}" ]; then
  if [ -n "${SCRATCH:-}" ]; then
    export PACE_ROOT="$SCRATCH/pbmc-grn"
  else
    export PACE_ROOT="$BUNDLE_DIR/workspace"
  fi
fi

# pull in the site config if the workspace / bundle carries one
for cand in "$PACE_ROOT/pace.env" "$BUNDLE_DIR/pace.env" "$PACE_ROOT/pipeline/pace.env"; do
  if [ -f "$cand" ]; then
    # shellcheck disable=SC1090
    source "$cand"
    break
  fi
done

export PACE_DATA="${PACE_DATA:-$PACE_ROOT/data}"
export PACE_RAW="${PACE_RAW:-$PACE_ROOT/raw}"
export SCAGA_REPO="${SCAGA_REPO:-$PACE_ROOT/tools/scSAGA}"

# ---- conda -----------------------------------------------------------------
# 'module' is a shell function on PACE; provide a no-op fallback so this file can
# be linted/run off-cluster.
if ! declare -f module >/dev/null 2>&1; then
  module() { return 0; }
fi
module load anaconda3

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

ENV_NAME="${ENV_NAME:?ENV_NAME must be set before sourcing _common.sh}"
conda activate "$ENV_NAME"

# ---- banner ----------------------------------------------------------------
echo "==================== $(basename "$0") ===================="
echo "host        : $(hostname)"
echo "slurm job   : ${SLURM_JOB_ID:-<none>}  account=${SLURM_JOB_ACCOUNT:-<none>} qos=${SLURM_JOB_QOS:-<none>}"
echo "date        : $(date -Is)"
echo "PACE_ROOT   : $PACE_ROOT"
echo "PACE_DATA   : $PACE_DATA"
echo "conda env   : $ENV_NAME  ($(command -v python))"
echo "python      : $(python -V 2>&1)"
echo "======================================================"
echo

mkdir -p "$PACE_ROOT/logs" "$PACE_ROOT/results"
