#!/bin/bash
# Submit the whole pipeline as a SLURM dependency chain, so you can walk away.
#
#   bash pace/submit_all.sh                          # unpaired, default resources
#   bash pace/submit_all.sh paired                   # the paired variant
#   PACE_ACCOUNT=GT-xxxx bash pace/submit_all.sh     # charge account
#   SKIP_DOWNLOAD=1 bash pace/submit_all.sh          # data already in $PACE_RAW
#   SKIP_6K=1 SKIP_4K=1 bash pace/submit_all.sh      # only what unpaired needs
#
# Chain:  download -> preprocess -> integration -> imputation -> grn -> evaluate
# Each step starts only if the previous one succeeded (afterok).
#
# Everything is configurable through the environment; see pace.env for the full
# list of knobs.  Prints the job IDs and the log-file paths when done.

set -euo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BUNDLE_DIR"

# ---- site config -----------------------------------------------------------
if [ -f ./pace.env ]; then
  # shellcheck disable=SC1091
  source ./pace.env
fi

EXPERIMENT="${1:-${PACE_EXPERIMENT:-unpaired}}"
export PACE_EXPERIMENT="$EXPERIMENT"
export PACE_ROOT="${PACE_ROOT:-$BUNDLE_DIR/workspace}"
export PACE_DATA="${PACE_DATA:-$PACE_ROOT/data}"
export PACE_RAW="${PACE_RAW:-$PACE_ROOT/raw}"

mkdir -p "$PACE_ROOT/logs" "$PACE_ROOT/results" logs

# ---- optional account / qos overrides -------------------------------------
SB_OPTS=()
[ -n "${PACE_ACCOUNT:-}" ]   && SB_OPTS+=("--account=${PACE_ACCOUNT}")
[ -n "${PACE_QOS:-}" ]       && SB_OPTS+=("--qos=${PACE_QOS}")
[ -n "${PACE_PARTITION:-}" ] && SB_OPTS+=("--partition=${PACE_PARTITION}")
SB_OPTS+=("--export=ALL")

submit() {   # submit <script> [dependency]
  local script="$1" dep="${2:-}"
  local args=("sbatch" "${SB_OPTS[@]}")
  [ -n "$dep" ] && args+=("--dependency=afterok:$dep")
  # shellcheck disable=SC2068
  local out
  out="$(${args[@]} "$script")"
  echo "$out" >&2
  echo "$out" | grep -oE '[0-9]+$'
}

echo "=== project     : $PACE_ROOT"
echo "=== experiment  : $EXPERIMENT"
echo "=== account     : ${PACE_ACCOUNT:-<none: set PACE_ACCOUNT>}"
echo "=== qos         : ${PACE_QOS:-<none>}"
echo "=== workers     : GRN_WORKERS=${GRN_WORKERS:-auto} x threads=${GRN_THREADS_PER_WORKER:-1}"
echo

# ---- build the chain -------------------------------------------------------
PREV=""

if [ "${SKIP_DOWNLOAD:-0}" != "1" ]; then
  J_DL=$(submit pace/slurm/00_download.slurm)
  echo "download     : $J_DL"
  PREV="$J_DL"
else
  echo "download     : skipped (SKIP_DOWNLOAD=1)"
fi

J_PRE=$(submit pace/slurm/01_preprocess.slurm "$PREV")
echo "preprocess   : $J_PRE"
J_INT=$(submit pace/slurm/02_integration.slurm "$J_PRE")
echo "integration  : $J_INT"
J_IMP=$(submit pace/slurm/03_imputation.slurm "$J_INT")
echo "imputation   : $J_IMP"
J_GRN=$(submit pace/slurm/04_grn.slurm "$J_IMP")
echo "grnboost2    : $J_GRN"
J_EVL=$(submit pace/slurm/05_evaluate.slurm "$J_GRN")
echo "evaluate     : $J_EVL"

cat <<EOF

=== submitted. chain: ${J_DL:-skip} -> $J_PRE -> $J_INT -> $J_IMP -> $J_GRN -> $J_EVL

Watch it:
  squeue -u \$USER
  tail -f $PACE_ROOT/logs/slurm_grn-grnboost2-$J_GRN.out

Cancel everything:
  scancel $J_EVL $J_GRN $J_IMP $J_INT $J_PRE ${J_DL:-}

NOTE: 06_run_grn.py is the long step.  The embers queue is free but preemptible
after a guaranteed first hour -- for a multi-hour GRN run either use inferno
(PACE_QOS=inferno) or set GRN_JOBQUEUE=1 to spread workers over several jobs.
EOF
